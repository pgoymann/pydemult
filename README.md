# Streamed and parallel demultiplexing of fastq files

## Quickstart

```
pydemult --fastq input.fastq.gz --barcodes barcodes.txt --threads 4 --writer-threads 16
```

## Requirements and usage

`pydemult` allows you to demultiplex fastq files in a streamed and parallel way. It expects that a **sample barcode** can be matched by a regular expression from the first line of each fastq entry and that sample barcodes are known in advance.

Suppose we have a file containing **sample barcodes** like this:

```
Sample  Barcode
sample1 CTTCAA
sample2 CAACAA
sample3 GTACGG
```

and a typical entry in the fastq file looks like this:

```
@HWI-ST808:140:H0L10ADXX:1:1101:8463:2:NNNNNN:CTTCCA
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATGATGCTGTGAGTTCC
+
@CCDDDDFHHHHHJIJFDDDDDDDDDBDDDDDBB0@B#####################
```

Since the sample barcode is six bases long, we have to set the corresponding `--barcode-regex` option to `(.*):(?P<CB>[ATGCN]{6}` in the call

```
pydemult --fastq input.fastq.gz --barcodes barcodes.txt --barcode-regex "(.*):(?P<CB>[ATGCN]{6}"
```

### Barcode and UMI regular expressions

By default, `pydemult` parses the read name for the cell barcode with regular expressions. Cell barcodes are indicated by a capturing group called `CB`, while (optional) UMIs are indicated by a capturing group called `UMI`. Some examples include:

- `(.*):(?P<CB>[ATGCN]{11})`, for a cell barcode of length 11 that is present after the last colon of the read name.
- `(.*):CELL_(?P<CB>[ATGCN]{10}):UMI_(?P<UMI>[ATGCN]{8})`, for a cell barcode of length 10, followed by a UMI sequence of length 8. For DropSeq data preprocessed by the [umis](https://github.com/vals/umis) tool, a regex like this is advisable.

### Output

`pydemult` will create a compressed fastq file for each **sample barcode**, with the filename taken from the corresponding *Sample* column entry of `barcodes.txt`.  

### A note on multithreading

`pydemult` divides its work into a demultiplexing and output part. The main thread streams the input and lazily distributes data blobs (of size `--buffer-size`) across `n` different demultiplexing threads (set with `--threads`), where the actual work happens. Demultiplexed input is then sent over to `m` threads for writing into individual output files (set with `--writer-threads`). Reading and demultiplexing are fast and CPU-bound operations, while output speed is determined by how fast data can be written to the underlying file system. In our experience, output is much slower than demultiplexing itself and requires proportionally more cores to speed up the runtime. We obtained best results when distributing output to three threads for each demultiplexing thread (`1:3` ratio of `--threads` to `--writer-threads`).  

# HTO detection

HTOCount is able to demultiplex hash SC experiments and counts hashes per cell.

# Installation

```
conda create --name pydemult  python=3.8
conda activate pydemult
python setup.py install
```
## Quick start

```HTOcount --reference Hash_tag.csv --whitelist barcodes.tsv --threads {cores} --barcode-sequences BC.fastq.gz --hashtag-sequences cDNA.fastq.gz ```

## Functions

- ```--reference```  Reference file with hashtag sequences
- ```--whitelist``` Whitelist which contains cell barcodes

- ```--barcode-sequences```  Fastqfile with barcodes
- ```--barcode-regex```  regex where HTOCount searches for barcode default first 16 bp 
- ```--barcode-edit-distance``` Number of allowed Mutations in barcode

- ```--hashtag-sequences``` Fastqfile with hashes
- ```--hashtag-regex```  regex where HTOCount searches for hash default between bp 10 and 15
- ```--hashtag-edit-distance``` Number of allowed Mutations in hash
- ```--sliding-window``` Boolean searches for hash in hole read for best match
- ```--sliding-window-hemming-distance``` Number of allowed Mutations in hash

- ```--buffer-size```
- ```--threads```

## License

The project is licensed under the MIT license. See the `LICENSE` file for details.
