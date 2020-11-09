import sys
import argparse
import logging
import re
import regex as mrab
import time
import subprocess
import heapq
from collections import Counter
import multiprocessing
from functools import reduce
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
import Levenshtein
from solo import hashsolo
import anndata

def barcode_merger(accumulator, element):
    for key, value in element.items():
        accumulator[key] = accumulator.get(key, set()).union(value)
    return accumulator

#assigend hashes to cells 

def hash_solo(adata):
    cell_hashing_data = anndata.read(adata)
    hashsolo.hashsolo(cell_hashing_data)
    return cell_hashing_data.obs

#demultiplexing functions

def find_best_match_shift(TAG_seq, tags, maximum_distance):
    """
    Find the best match from the list of tags with sliding window.
    Compares the Levenshtein distance between tags and the trimmed sequences.
    The tag and the sequence must have the same length.
    If no matches found returns 'unmapped'.
    We add 1
    Args:
        TAG_seq (string): Sequence from R1 already start trimmed
        tags (dict): A dictionary with the TAGs as keys and TAG Names as values.
        maximum_distance (int): Maximum distance given by the user.
    Returns:
        best_match (string): The TAG name that will be used for counting.
    """
    best_match = 'unmapped'
    best_score = maximum_distance
    shifts = range(0,len(TAG_seq) - len(max(tags,key=len)))

    for shift in shifts:
        for tag, name in tags.items():
            score = Levenshtein.hamming(tag, TAG_seq[shift:len(tag)+shift])
            if score == 0:
                #Best possible match
                return(name)
            elif score <= best_score:
                best_score = score
                best_match = name
                return(best_match)
    return(best_match)

def _find_bc(queue, chunk, mutationhash, regex,  keep_empty = False):
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
    queue.put(barcodes_result_dict)
    
def _find_hto(queue,chunk, mutationhash, regex, keep_empty = False):
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
    queue.put(hto_result_dict)
    

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
    parser.add_argument('--buffer-size', help="Buffer size for the FASTQ reader (in Bytes). Must be large enough to contain the largest entry.", type = int, default = 100000018, metavar = '100000018')
    parser.add_argument('--threads', '-t', help='Number of threads to use for multiprocessing.', type=int, metavar='1', default=1)
    parser.add_argument('--barcode-sequences', '-i', help='FASTQ file containing barcode sequences', metavar='input_BC.fastq.gz', type=str)
    parser.add_argument('--hashtag-sequences', '-j', help='FASTQ file containing hash tag sequences', metavar='input_HT.fastq.gz', type=str)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--silent', '-s', action='store_true')

    args = parser.parse_args()
    # Set Logger
    logger = logging.getLogger(__name__)
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    sh.setLevel(logging.DEBUG)
    logger.addHandler(sh)

    #set logger level
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.silent:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    logger.info('Working on Barcodefile {} and Hashfile {} using {} threads and a buffsize of {}'.format(args.barcode_sequences, 
        args.hashtag_sequences, args.threads, args.buffer_size))

    bufsize = args.buffer_size

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    elif os.path.isfile(outdir ):
        sys.exit('Their is already a folder ' + outdir)

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
    logger.info('Creating mutation hash with edit distance {} for {} barcodes.'.format(args.barcode_edit_distance, len(barcodes)))
    

    barcode_mutationhash = mutationhash(strings = barcodes, nedit = args.barcode_edit_distance, 
                                            alphabet = list(args.edit_alphabet), log = logger)

    #
    # Construct a list of barcode indices from barcode reads
    #

    logger.info('Construct a list of barcode indices from barcode reads')
    logger.info('Start reading barcodes file.')
    
    zcat = subprocess.Popen(['zcat', args.barcode_sequences], stdout=subprocess.PIPE, bufsize = bufsize)
    queue = multiprocessing.Queue()#Que vor result from multiprocesses
    barcode_results = []#storting results from que
    with zcat.stdout as fh:
        blob_generator = buffered_blob(fh, bufsize)#genert blob by bufsize
        procs_counter = 0
        procs = []
        for index in blob_generator:
            procs_counter +=1
            proc = multiprocessing.Process(target=_find_bc, args=(queue,index, barcode_mutationhash, c_barcode_regex))
            procs.append(proc)
            proc.start()
            if procs_counter == args.threads:
                for i in procs:
                    barcode_results.append( queue.get())
                for i in procs:
                    i.join()
                procs = []
                procs_counter = 0

    for i in procs:
        barcode_results.append( queue.get())

    for i in procs:
        i.join()#wait until Process have finished
    logger.info('Finish reading barcodes')

    #
    # Create mutationhash for barcodes from whitelist
    #
    logger.info('Reading hastags reference file')
    hashes = []
    hashes_names = []
    with open(args.reference, 'r') as file:
        for line in file:
            hashes.append(line.rstrip('\n').split(',')[0])
            hashes_names.append(line.rstrip('\n').split(',')[1])
    
    hashes_names.append('unmaped')#count all reads which didn't contain a hash sequence

    logger.debug('Whitelist contains the following : ' + ','.join(hashes))
    logger.info('Creating mutation hash with edit distance {} for {} barcodes.'.format(args.hashtag_edit_distance, len(hashes)))
    hashtag_mutationhash = mutationhash(strings = hashes, nedit = args.hashtag_edit_distance,
                                             alphabet = list(args.edit_alphabet), log = logger)


    #
    # Construct a list of hto indices from cDNA reads
    #
    logger.info('Construct a list of hto indices from cDNA reads')
    logger.info('Start reading hashtag file')
    queue = multiprocessing.Queue()#Que vor result from multiprocesses    
    hto_results = []
    zcat = subprocess.Popen(['zcat', args.hashtag_sequences], stdout=subprocess.PIPE, bufsize = bufsize)
    with zcat.stdout as fh:
        blob_generator = buffered_blob(fh, bufsize)
        procs_counter = 0
        procs = []
        for index in blob_generator:
            procs_counter +=1
            proc = multiprocessing.Process(target=_find_hto, args=(queue,index, hashtag_mutationhash, c_hashtag_regex))
            procs.append(proc)
            proc.start()
            if procs_counter == args.threads:
                for i in procs:
                    hto_results.append( queue.get())
                for i in procs:
                    i.join()
                procs = []
                procs_counter = 0

    for i in procs:
        hto_results.append( queue.get())

    for i in procs:
        i.join()#wait until Process have finished
    logger.info('Finished reading hashtag file')

    logger.info('Strat merging results')

    #merging Hash tag data
    hto_data = {}
    for i in hto_results:
        hto_data = {**hto_data, **i}

    #merging barcodes data
    barcode_data = reduce(barcode_merger, barcode_results, {})
    
    logger.info('Finished merging results')

    logger.info('Creating results Matrix')
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

    for barcode in barcode_data:
        barcode_result.append(barcode)
        count_result = [0] * len(hashes_names)
        for a in barcode_data[barcode]:
            if a in hto_data:
                count_result[hash_result[hto_data[a]]] += 1
            else:
                count_result[-1] += 1
        dict_all[k] =  count_result
        cell_number.append(len(barcode_data[barcode]))

        if hashes_names[count_result.index(max(count_result))] != 'unmaped':
            hash_result_for_barcode.append(hashes_names[count_result.index(max(count_result))])
        else:
            hash_result_for_barcode.append(hashes_names[heapq.nlargest(2, range(len(count_result)), key=count_result.__getitem__)[1]] + '_unmaped')

        mtx_matrix.append(count_result)

    
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

    logging.info('Assign hahes to Cells')

def main():
    count()

if __name__ == '__main__':
    main()