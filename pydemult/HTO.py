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
import csv
import pandas as pd
from scipy import sparse
import scipy.io
import json
import numpy as np
import os
import Levenshtein
from solo import hashsolo
import anndata as ad

#self imports 
from buffered_reader import buffered_blob

from mutationhash import mutationhash
from worker import entryfunc


def barcode_merger(accumulator, element):
    for key, value in element.items():
        accumulator[key] = accumulator.get(key, set()).union(value)
    return accumulator

#assigend hashes to cells 

def hash_solo(adata):
    print(adata)
    cell_hashing_data = ad.read('/mnt/workspace/pgoyman/SC_analysis/test_pydemult/test_data/test_data/out/adata_sub_hash.h5ad')
    #hashsolo.hashsolo(cell_hashing_data)
    #return cell_hashing_data.obs

#demultiplexing functions

def find_best_match_shift(TAG_seq, tags, maximum_distance):
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
                return (best_match)
    return(best_match)

def _find_bc(chunk,barcodes_result_dict, mutationhash, regex,  keep_empty = False):

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

                    t = barcodes_result_dict[barcode]
                    t.append(line)
                    barcodes_result_dict[barcode] = t
            except KeyError:
                    is_unmatched = True

    
def _find_hto(queue,chunk, mutationhash, regex, dict_hashes, inv_dict_hash, sliding_window_hemming_distance, sliding_window = False,keep_empty = False):
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
                if sliding_window:
                    
                    bc_match = find_best_match_shift(entry[1].decode('utf-8'), inv_dict_hash, sliding_window_hemming_distance)
                    if  bc_match != 'unmapped':
                        origin = mutationhash[dict_hashes[bc_match]]
                        if len(origin) > 1:
                            is_unmatched = True
                        else:
                            hto_hash = list(origin)[0]  
                            hto_result_dict[line] = hto_hash
                            count += 1
                else:
                    is_unmatched = True
    queue.put(hto_result_dict)
    

def count():
    parser = argparse.ArgumentParser(description='Demultiplex samples based fastq files from hash tag oligo data')
    parser.add_argument('--out', '-o', help='Out put folder', default = 'out/')
    parser.add_argument('--reference', '-r', help='Tab-separated reference file containing hash tag sequences and names')
    parser.add_argument('--whitelist', '-w', help='Cell barcode whitelist of allowed / known cell barcodes')
    parser.add_argument('--barcode-regex', '-b', help = 'Regular expression to parse cell barcode (CB) from barcode sequences', default = '[ATGCN]{16}', type = str)
    parser.add_argument('--hashtag-regex', '-c', help = 'Regular expression to parse hash tag sequences (HTO) from hash tag sequences', default = '^.{10}\K([ATGCN]{15})', type = str)
    parser.add_argument('--barcode-edit-distance', help='Maximum allowed edit distance for barcodes', metavar = '1', type=int, default = 1)
    parser.add_argument('--hashtag-edit-distance', help='Maximum allowed edit distance for hash tag oligos', metavar = '2', type=int, default = 1)
    parser.add_argument('--edit-alphabet', help='The alphabet that is used to created edited barcodes / hash tag sequences', choices=['N', 'ACGT', 'ACGTN'], default = "ACGTN", type = str, metavar = "ACGTN")
    parser.add_argument('--buffer-size', help="Buffer size for the FASTQ reader (in Bytes). Must be large enough to contain the largest entry.", type = int, default = 100000018, metavar = '100000018')#100000018
    parser.add_argument('--threads', '-t', help='Number of threads to use for multiprocessing.', type=int, metavar='1', default=1)
    parser.add_argument('--barcode-sequences', '-i', help='FASTQ file containing barcode sequences', metavar='input_BC.fastq.gz', type=str)
    parser.add_argument('--hashtag-sequences', '-j', help='FASTQ file containing hash tag sequences', metavar='input_HT.fastq.gz', type=str)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--sliding-window-hemming-distance', help='Maximum allowed edit distance for hash tag oligos', metavar = '2', type=int, default = 2)
    parser.add_argument('--sliding-window', action='store_true')
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
    outdir = args.out
    if outdir[-1] != '/':
        outdir += '/'

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    elif os.path.isfile(outdir):
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
    manager = multiprocessing.Manager()
    barcodes_dict = manager.dict()
    
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

    barcode_results = []#storting results from que
    with zcat.stdout as fh:
        blob_generator = buffered_blob(fh, bufsize)#genert blob by bufsize
        procs_counter = 0
        procs = []
        for index in blob_generator:
            procs_counter +=1
            proc = multiprocessing.Process(target=_find_bc, args=(index,barcodes_dict, barcode_mutationhash, c_barcode_regex))
            procs.append(proc)
            proc.start()
            if procs_counter == args.threads:

                for i in procs:
                    i.join()
                procs = []
                procs_counter = 0

    for i in procs:
        i.join()#wait until Process have finished
    
    logger.info('Finish reading barcodes')

    #
    # Create mutationhash for barcodes from whitelist
    #
    logger.info('Reading hastags reference file')
    hashes = []
    hashes_names = []
    #read_hatag file
    hash_frame = pd.read_csv(args.reference)
    hashes_names = hash_frame['name'].to_list()
    hashes = hash_frame['sequence'].to_list()
    hash_dict = {key:val for key,val in zip(hashes_names,hashes) }
    inv_dict_hash = {v: k for k, v in hash_dict.items()}

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
            
            #_find_hto(queue,index, hashtag_mutationhash, c_hashtag_regex)
            procs_counter +=1
            proc = multiprocessing.Process(target=_find_hto, args=(queue,index, hashtag_mutationhash, 
                    c_hashtag_regex, hash_dict,inv_dict_hash, args.sliding_window_hemming_distance, args.sliding_window))
            procs.append(proc)
            proc.start()
            if procs_counter == args.threads:
                for i in procs:
                    hto_results.append( queue.get())
                for i in procs:
                    i.join()
                procs = []
                procs_counter = 0
    #get last values from que
    for i in procs:
        hto_results.append( queue.get())
    #wait for porcesses
    for i in procs:
        i.join()#wait until Process have finished"""
    logger.info('Finished reading hashtag file')

    logger.info('Strat merging results')

    #merging Hash tag data
    hto_data = {}
    for i in hto_results:
        hto_data = {**hto_data, **i}

    #merging barcodes data
    barcode_data = barcodes_dict#reduce(barcode_merger, barcode_results, {})
    
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
        dict_all[barcode] =  count_result
        cell_number.append(len(barcode_data[barcode]))

        if hashes_names[count_result.index(max(count_result))] != 'unmaped':
            hash_result_for_barcode.append(hashes_names[count_result.index(max(count_result))])
        else:
            hash_result_for_barcode.append(hashes_names[heapq.nlargest(2, range(len(count_result)), key=count_result.__getitem__)[1]] + '_unmaped')
        mtx_matrix.append(count_result)


    print(hashes)
    print(hashes_names)
    logging.info('Assign hahes to Cells')
    mat_coo = sparse.coo_matrix((np.array(mtx_matrix).astype(np.float)))
    scipy.io.mmwrite(args.out + 'matrix.mtx', mat_coo)

    obs_data_farme = pd.DataFrame(np.column_stack([hash_result_for_barcode  ,cell_number]),index=barcode_result,columns=['Hash names by max count' , 'Number found barcodes'])

    #obs_data_farme = pd.DataFrame({'Barcodes':barcode_result,'Hash names' : hash_result_for_barcode , 'Number found barcodes' : cell_number})
    var_data_farme = pd.DataFrame(index=hashes_names)
    obs_data_farme.to_csv(args.out + 'barcode.csv')
    var_data_farme.to_csv(args.out + 'feature.csv')
    print(hashes_names)
    adata = ad.AnnData(sparse.csr_matrix(mat_coo),  var=var_data_farme , obs=obs_data_farme, dtype='int32')
    adata.write(filename=args.out + 'adata.h5ad', compression='gzip')
    hashs = adata.var.index[adata.var.index.str.contains('HashTag*')].tolist()
    adata_subset = adata[:, hashs]
    adata_subset.write(filename=args.out + 'adata_sub_hash.h5ad',as_dense='X')
    print(args.out + 'adata_sub_hash.h5ad')
    hash_solo(args.out + 'adata_sub_hash.h5ad')
 
def main():
    count()

if __name__ == '__main__':
    main()