#!/usr/bin/env python
"""

Gene Coordinates Fetcher 

Author: Xin Wu

GrandOmics Copyright 2019

Usage: 
    gene_coordinates_fetcher.py  [-e email] [-a assembly]  GENE_FILE  OUTPUT_FILE

Gene Coordinates Fetcher is a script to get HUMAN chr positions by gene name.
Note: gene name with '-' is not allowd in NCBI and it just keep the first part before '-'.
Such as NREP-EPB41L4A will be replaced by NREP

Arguments:
    GENE_FILE                 Plain file with one gene name per line
    OUTPUT_FILE               Output file with one gene per line

Options:
    -h --help                 Print help
    -e --email=email          User email [default: usr@nextomics.com] 
    -a --assembly=assembly    Assembly of hg38 or hg19 [default: hg19] 
    -v --version              Version

"""

from docopt import docopt
from Bio import Entrez
import logging
import time
import collections
import json

# qname is the query name, name is the standard name
Gene = collections.namedtuple('Gene', 'qname name id assembly chr start stop') 

def main(infile: str, assembly: str, outfile: str):
    #gene = query('OBFC1', 'hg19') 
    #print(gene)
    names = readNames(infile)
    genes = list()
    for name in names:
        gene = query(name, assembly)
        genes.append(gene)
    write(genes, outfile)
    logger.info('[Entrez] query done.')
    

def readNames(infile: str):
    names = set()
    with open(infile) as f:
        for line in f:
            name = line.strip()    
            # keep only the first part if there is a '-' in gene name since NCBI can not find it
            # and we found MBIP and NREP are good for MBIP-NKX2-1 and NREP-EPB41L4A 
            name = name.split('-')[0]
            names.add(name)
        return names

def write(genes: list, outfile: str):
     with open(outfile, 'w') as f:
         header = ','.join(['query_name','name','id','assembly','chr','start','stop'])
         f.write(header)
         for gene in genes:
           row = ','.join([gene.qname, gene.name, gene.id, gene.assembly, gene.chr, str(gene.start), str(gene.stop)])
           f.write(row) 
           f.write('\n')
    

def query(qname: str, assembly: str):
    # sleep 2 seconds to for throttling of NCBI Entrez 
    time.sleep(1)
    # info = Entrez.esearch(db = "gene",term = "OBFC1[sym] AND human[ORGN]", retmode="json")
    term = qname + "[sym] AND human[ORGN]"
    logger.info('[Entrez] search the first gene id for {qname}'.format(qname=qname))
    res = Entrez.esearch(db = "gene",term = term, retmode="json" )
    res = json.load(res)
    # just use 1st id
    time.sleep(1)
    ids = res['esearchresult']['idlist']
    if len(ids) == 0:
        logger.info('[Entrez] can not find any id for {qname}, skip ...'.format(qname=qname))
        return None
    else:
        id = ids[0]
        logger.info('[Entrez] get gene summary for {id}'.format(id=id))
        res =  Entrez.esummary(db="gene",id=id, retmode="json")
        res = json.load(res)
        # refseq assembly for hg19 patch 13
        acc = 'GCF_000001405.25'
        if assembly == 'hg38':
            # refseq assembly for hg38 patch 13
            acc = 'GCF_000001405.39'
        name =  res['result'][id]['name']
        loc_hist = res['result'][id]['locationhist']
        chr = res['result'][id]['chromosome']
        loc = next(item for item in loc_hist if item['assemblyaccver'] == acc)
        start, stop = loc['chrstart'], loc['chrstop'] 
        gene = Gene(qname=qname, name=name, id=id, assembly=assembly, chr=chr, start=start, stop=stop) 
        logger.info(str(gene))
        return gene
     
    


if __name__ == '__main__':
    args = docopt(__doc__, version='version 0.1')
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('fetcher')
    logger.info(args)
   
    genefile = args['GENE_FILE']
    outfile = args['OUTPUT_FILE']
    email = args['--email']
    assembly = args['--assembly']
    Entrez.email = email
    main(genefile, assembly, outfile)
    
