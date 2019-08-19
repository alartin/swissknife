#!/usr/bin/env python

"""

ClinVar Searcher 

Author: Xin Wu

GrandOmics Copyright 2019

Usage: 
    clinvar_searcher.py [-m max] GENE 

ClinVar Searcher is a script for listing all variants given a gene. 

Arguments:
    Gene        Gene name. e.g SLC26A4

Options:
    -m --max=max            Max count of result for a query [default: 500]
    -h --help               Print help
    -v --version            Version

"""

from docopt import docopt
from Bio import Entrez
import time
import logging
import collections
import json


Variant = collections.namedtuple('Variant', 'accession name vtype assembly band start stop') 

def main(gene: str):
    
    logger.info('[Entrez] search ClinVar for {id}'.format(id=gene))
    term = gene + '[Gene]'
    query = Entrez.esearch(db='clinvar', term=term, retmax=500)
    result = Entrez.read(query)
    count = result['Count']
    ids = result['IdList']
    logger.info('[Count] {count}'.format(count=count))
    #logger.info('[Ids] {ids}'.format(ids=ids))
    for id in ids:
        var = getVariant(id)
        logger.info(var) 



def getVariant(id: str):

    # sleep 2 seconds to for throttling of NCBI Entrez 
    time.sleep(2)
    logger.info('[Entrez] get ClinVar variant for {id}'.format(id=id))
    result = Entrez.esummary(db='clinvar', id=id, retmode='json')
    res = json.load(result)
    vid = res['result'][id] 
    acc, vtype, name = vid['accession'], vid['obj_type'], vid['title']
    # all locs of the 1st var
    vlocs = vid['variation_set'][0]['variation_loc']
    # current loc
    cvloc = next(l for l in vlocs if l['status'] == 'current') 
    assembly, band, start, stop = cvloc['assembly_name'], cvloc['band'], cvloc['start'], cvloc['stop']
    variant = Variant(accession=acc, name=name, vtype=vtype, assembly=assembly, band=band, start=start, stop=stop)
    return variant
    
     
    
    


if __name__ == '__main__':
    args = docopt(__doc__, version='clinvar_searcher version 0.1')
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('clinvar')
    logger.info(args)
   
    gene = args['GENE']
    Entrez.email = 'robot@grandomics.com'
    main(gene)
    
