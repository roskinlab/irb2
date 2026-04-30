#!/usr/bin/env python

import sys
import argparse
import logging
import time

from Bio import SeqIO
from irbase.utils import open_compressed


def genbank_filter(it, organism=None, max_size=None):
    for record in it:
        if (organism is None) or (record.annotations['organism'] == organism):
            sequence_length = len(record.seq)
            if (max_size is None) or (sequence_length <= max_size):
                yield record

def main():
    parser = argparse.ArgumentParser(description='subset Genbank files bases on organism and length',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # files
    parser.add_argument('genbank_filenames', metavar='genbank_file', nargs='+', help='the file with the Genbank records')
    parser.add_argument('--organism', '-o', metavar='O', type=str, default='Homo sapiens', help='only process records with the given organism')
    parser.add_argument('--max-size', '-m', metavar='M', type=int, default=50000, help='ignore sequences longer than this')

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    start_time = time.time()

    processed_record_count = 0

    for genbank_filename in args.genbank_filenames:
        logging.info('processing %s', genbank_filename)
        with open_compressed(genbank_filename, 'rt') as genbank_file:
            records = SeqIO.parse(genbank_file, 'genbank')
            filtered_records = genbank_filter(records,
                                              organism=args.organism,
                                              max_size=args.max_size)
            SeqIO.write(filtered_records, sys.stdout, 'genbank')

    elapsed_time = time.time() - start_time
    logging.info('elapsed time %s', time.strftime('%H hours, %M minutes, %S seconds', time.gmtime(elapsed_time)))
    
if __name__ == '__main__':
    sys.exit(main())
