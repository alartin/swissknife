#!/usr/bin/env python

"""

FASTQ SAMPLER 

Author: Xin Wu

GrandOmics Copyright 2019

Usage: 
    fastq_sampler.py  [-b bases] [-k logk]  IN_FASTQ  OUT_FASTQ_GZ

FASTQ Sampler is a script for sampling reads in FASTQ format.

Arguments:
    IN_FASTQ        input file in FASTQ format
    OUT_FASTQ_GZ    output file in FASTQ and GZ format

Options:
    -h --help               Print help
    -b --bases=bases        Total number of base pairs sampled from original FASTQ
    -k --logk=logk          Print log info every k reads [default: 100] 
    -v --version            Version

"""

from docopt import docopt
from Bio import SeqIO
from random import randint
import gzip, logging

def main(in_fastq: str, out_fastq: str, bases: int, logk: int):
    fq = SeqIO.index(in_fastq, 'fastq')
    out_fq = gzip.open(out_fastq, 'wt')
    
    keys = list(fq.keys())
    sampled_ids = set()
    sampled_reads = []
    accumulated_bases = 0
    while(accumulated_bases < bases):
        read_id = keys[randint(0, len(keys) - 1)]
        if read_id not in sampled_ids:
            read = fq.get(read_id)
            sampled_reads.append(read)
            accumulated_bases += len(read) 
            if len(sampled_reads) % logk == 0: 
                logger.info('[Accumulated Bases]  {acc}'.format(acc=accumulated_bases))
        else:
            logger.info('[Duplicated sampled to skip]  {id}'.format(id=read_id))
            continue
    

    SeqIO.write(sequences=sampled_reads, handle=out_fq, format='fastq')
    out_fq.close()
    


if __name__ == '__main__':
    args = docopt(__doc__, version='fastq_sampler version 0.1')
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('sampler')
    logger.info(args)
   
    in_fastq = args['IN_FASTQ']
    out_fastq = args['OUT_FASTQ_GZ']
    bases = int(args['--bases'])
    logk = int(args['--logk'])
    main(in_fastq, out_fastq, bases, logk)
    
