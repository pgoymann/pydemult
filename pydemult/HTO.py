import sys
import argparse
import logging
import re
import regex as mrab
import time
import subprocess
import heapq
from functools import partial
from mputil import lazy_map
from collections import defaultdict
from mutationhash import mutationhash
from worker import entryfunc
from buffered_reader import buffered_blob
import csv
from scipy import sparse
import scipy.io
import json
import numpy as np
import os

def _find_bc(chunk, mutationhash, regex,  keep_empty = False):
    start = time.time()
    barcodes_result_dict  = defaultdict(set)
    count = 0
    line = 0
    for entry in entryfunc(chunk):
        line += 1
        # Ignore empty sequences
        if entry[1].decode('utf-8') == '' and not keep_empty:
            continue

        is_unmatched = False
        match = regex.match(entry[1].decode('utf-8'))

        if match is not None:
            try:
                bc_match = match.group()
                origin = mutationhash[bc_match]
                if len(origin) > 1:
                    is_unmatched = True
                else:
                    barcode = list(origin)[0]
                    count += 1

                    barcodes_result_dict[barcode].add(line)
            except KeyError:
                    is_unmatched = True
    parsing = time.time()

    return((count, parsing-start,barcodes_result_dict,line))
    

def _find_hto(chunk, mutationhash, regex, keep_empty = False):

    start = time.time()
    hto_result_dict  = {}
    count = 0
    line = 0
     
    for entry in entryfunc(chunk):
        line += 1
        # Ignore empty sequences
        if entry[1].decode('utf-8') == '' and not keep_empty:
            continue

        is_unmatched = False
        match = regex.match(entry[1].decode('utf-8'))

        if match is not None:
            try:
                bc_match = match.group()
                origin = mutationhash[bc_match]

                if len(origin) > 1:
                    is_unmatched = True
                else:
                    hto_hash = list(origin)[0]
                    
                    hto_result_dict[line] = hto_hash
                    count += 1
            except KeyError:
                is_unmatched = True
    parsing = time.time()
    
    return((count, parsing-start,hto_result_dict, line))

def count():
    parser = argparse.ArgumentParser(description='Demultiplex samples based fastq files from hash tag oligo data')
    parser.add_argument('--out', '-o', help='Out put folder')
    parser.add_argument('--reference', '-r', help='Tab-separated reference file containing hash tag sequences and names')
    parser.add_argument('--whitelist', '-w', help='Cell barcode whitelist of allowed / known cell barcodes')
    parser.add_argument('--barcode-regex', '-b', help = 'Regular expression to parse cell barcode (CB) from barcode sequences', default = '[ATGCN]{16}', type = str)
    parser.add_argument('--hashtag-regex', '-c', help = 'Regular expression to parse hash tag sequences (HTO) from hash tag sequences', default = '^.{10}\K([ATGCN]{15})', type = str)
    parser.add_argument('--barcode-edit-distance', help='Maximum allowed edit distance for barcodes', metavar = '1', type=int, default = 1)
    parser.add_argument('--hashtag-edit-distance', help='Maximum allowed edit distance for hash tag oligos', metavar = '2', type=int, default = 1)
    parser.add_argument('--edit-alphabet', help='The alphabet that is used to created edited barcodes / hash tag sequences', choices=['N', 'ACGT', 'ACGTN'], default = "ACGTN", type = str, metavar = "ACGTN")
    parser.add_argument('--buffer-size', help="Buffer size for the FASTQ reader (in Bytes). Must be large enough to contain the largest entry.", type = int, default = 40000000, metavar = '40000000')
    parser.add_argument('--threads', '-t', help='Number of threads to use for multiprocessing.', type=int, metavar='1', default=1)
    parser.add_argument('--barcode-sequences', '-i', help='FASTQ file containing barcode sequences', metavar='input_BC.fastq.gz', type=str)
    parser.add_argument('--hashtag-sequences', '-j', help='FASTQ file containing hash tag sequences', metavar='input_HT.fastq.gz', type=str)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    logger = logging.getLogger(__name__)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    sh.setLevel(logging.DEBUG)
    logger.addHandler(sh)

    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARNING)

    logger.debug('Working on {} using {} threads'.format(args.barcode_sequences, args.threads))

    #
    # Create regular expression for barcode / hash tag oligo parsing from sequences
    #
    c_barcode_regex = re.compile(args.barcode_regex)
    c_hashtag_regex = mrab.compile(args.hashtag_regex)

    #
    # Create mutationhash for barcodes from whitelist
    #
    logger.info('Creating mutationhash for barcodes from whitelist')
    with open(args.whitelist, 'r') as file:
        barcodes = [line.rstrip('\n') for line in file]
    barcodes_dict = {}
    for i in barcodes:
        barcodes_dict[i] = []
    logger.debug('Whitelist contains the following barcodes: ' + ','.join(barcodes))

    logger.debug('Creating mutation hash with edit distance {} for {} barcodes.'.format(args.barcode_edit_distance, len(barcodes)))
    barcode_mutationhash = mutationhash(strings = barcodes, nedit = args.barcode_edit_distance, alphabet = list(args.edit_alphabet), log = logger)

    #
    # Construct a list of barcode indices from barcode reads
    #
    logger.info('Construct a list of barcode indices from barcode reads')
    bufsize = args.buffer_size
    
    zcat = subprocess.Popen(['zcat', args.barcode_sequences], stdout=subprocess.PIPE, bufsize = bufsize)
    with zcat.stdout as fh:
        blob_generator = buffered_blob(fh, bufsize)
        find_bc = partial(_find_bc, mutationhash = barcode_mutationhash, regex = c_barcode_regex)
        barcode_data = lazy_map(find_bc, blob_generator, n_cpus = args.threads)
    for d in barcode_data:
        logger.info('{}\t{} seconds'.format(d[0], d[1]))
    print(barcode_data[0])

    #
    # Create mutationhash for barcodes from whitelist
    #
    logger.info('Create mutationhash for barcodes from whitelist')
    hashes = []
    hashes_names = []
    with open(args.reference, 'r') as file:
        for line in file:
            hashes.append(line.rstrip('\n').split('\t')[1])
            hashes_names.append(line.rstrip('\n').split('\t')[0])
    hashes_names.append('unmaped')
    logger.debug('Whitelist contains the following barcodes: ' + ','.join(hashes))

    logger.debug('Creating mutation hash with edit distance {} for {} barcodes.'.format(args.hashtag_edit_distance, len(hashes)))
    hashtag_mutationhash = mutationhash(strings = hashes, nedit = args.hashtag_edit_distance, alphabet = list(args.edit_alphabet), log = logger)


    #
    # Construct a list of hto indices from cDNA reads
    #
    logger.info('Construct a list of hto indices from cDNA reads')

    bufsize = args.buffer_size
    zcat = subprocess.Popen(['zcat', args.hashtag_sequences], stdout=subprocess.PIPE, bufsize = bufsize)
    with zcat.stdout as fh:
        blob_generator = buffered_blob(fh, bufsize)
        find_hto = partial(_find_hto, mutationhash = hashtag_mutationhash, regex = c_hashtag_regex)
        hto_data = lazy_map(find_hto, blob_generator, n_cpus = args.threads)

    for d in hto_data:
        print('{}\t{} seconds'.format(d[0], d[1]))

    hash_count = 0
    hash_result = {}
    for i in hashes:
        hash_result[i] = hash_count
        hash_count += 1

    barcode_result = []
    hash_result_for_barcode = []
    cell_number = []
    mtx_matrix = []
    dict_all = {}
    print(hto_data[0][2])
    for k in barcode_data[0][2].keys():
        barcode_result.append(k)
        count_result = [0,0,0,0,0]
        for a in barcode_data[0][2][k]:
            if a in hto_data[0][2]:
                count_result[hash_result[hto_data[0][2][a]]] += 1
            else:
                count_result[4] += 1
        dict_all[k] =  count_result
        cell_number.append(len(barcode_data[0][2][k]))
        if hashes_names[count_result.index(max(count_result))] != 'unmaped':
            hash_result_for_barcode.append(hashes_names[count_result.index(max(count_result))])
        else:
            hash_result_for_barcode.append(hashes_names[heapq.nlargest(2, range(len(count_result)), key=count_result.__getitem__)[1]] + '_unmaped')

        mtx_matrix.append(count_result)
    

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    
    with open(args.out + 'info.txt', 'w') as file:
        file.write(json.dumps(dict_all)) # use `json.loads` to do the reverse

    with open(args.out + 'barcodes.csv','w') as f:
        for barcode,hashes_name,cell_number in zip(barcode_result,hash_result_for_barcode, cell_number):
            f.write(barcode + '\t')
            f.write(str(cell_number) + '\t')
            f.write("%s\n" % hashes_name)
    hashes.append('unmaped')
    with open(args.out + 'features.csv','w') as f:
        for item in hashes:
            f.write("%s\n" % item)

    mat_coo = sparse.coo_matrix((np.array(mtx_matrix).astype(np.float)))
    scipy.io.mmwrite(args.out + 'matrix.mtx', mat_coo)

def main():
    count()

if __name__ == '__main__':
    main()