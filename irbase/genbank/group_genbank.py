#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
import logging
import time
import os

from Bio import SeqIO

from irbase.utils import open_compressed

def calculate_year(journal):
    if journal == 'Unpublished':
        return 0
    elif journal.startswith('Submitted ('):
        date = journal[len('Submitted ('):journal.index(')')]
        day, month, year = date.split('-')
        return int(year)
    else:
        if journal.find('(') != -1:
            year = journal[journal.rindex('(') + 1:journal.rindex(')')]
            return int(year)
        else:
            return 0

def calculate_key_year(record):
    pubmed_id = None
    pubmed_year = None
    pubmed_title = None
    for reference in record.annotations['references']:
        if reference.pubmed_id != '':
            pubmed_id = int(reference.pubmed_id)
            pubmed_title = reference.title
            try:
                pubmed_year = calculate_year(reference.journal)
            except:
                #print('===>', file=sys.stderr)
                #print(reference.journal, file=sys.stderr)
                #print('<===', file=sys.stderr)
                raise

    if pubmed_id is not None:
        return pubmed_id, pubmed_year, pubmed_title
    else:
        last_reference = record.annotations['references'][-1]
        key = (last_reference.authors, last_reference.title, last_reference.journal)
        year = calculate_year(last_reference.journal)
        return key, year, last_reference.title

def main():
    parser = argparse.ArgumentParser(description='extract immune receptor sequences from Genbank records',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # directory to store the batches in
    parser.add_argument('batch_dirname', metavar='dir', help='name for the batch directory')
    # input files
    parser.add_argument('genbank_filenames', metavar='genbank-file', nargs='+', help='the file with the Genbank records')
    # options
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])
    start_time = time.time()

    # check if the batch directory already exists
    if os.path.exists(args.batch_dirname):
        logging.error('batch direcotry %s already exists', args.batch_dirname)
        return 10
    else:
        os.mkdir(args.batch_dirname)

    data = {}
    key_title_map = {}
    record_counts = 0
    
    # read in all the Genbank records
    for filename in args.genbank_filenames:
        with open_compressed(filename, 'rt') as genbank_handle:
            for record in SeqIO.parse(genbank_handle, 'genbank'):
                record_counts += 1

                # store the data by year and key
                key, year, title = calculate_key_year(record)
                if year not in data:
                    data[year] = {}
                if key not in data[year]:
                    data[year][key] = []
                data[year][key].append(record) 
                key_title_map[key] = title

    # make a directory for each year
    for year in data:
        os.mkdir('%s/%04d' % (args.batch_dirname, year))

        # for each key, make a file with the Genbank records
        for key in data[year]:
            title = key_title_map[key]
            if type(key) is int:
                filename = '%s/%04d/pmid_%d.genbank' % (args.batch_dirname, year, key)
            else:
                filename = '%s/%04d/hash_%d.genbank' % (args.batch_dirname, year, abs(hash(key)))
            with open(filename, 'wt') as handle:
                SeqIO.write(data[year][key], handle, 'genbank')
                for record in data[year][key]:
                    if type(key) is int:
                        print(year, 'pmid_%d' % key, record.id, title, record.description, sep='\t')
                    else:
                        print(year, 'hash_%d' % abs(hash(key)), record.id, title, record.description, sep='\t')
                        
    elapsed_time = time.time() - start_time
    logging.info('elapsed time %s', time.strftime('%H hours, %M minutes, %S seconds', time.gmtime(elapsed_time)))
    
if __name__ == '__main__':
    sys.exit(main())
