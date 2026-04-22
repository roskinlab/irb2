#!/usr/bin/env python

import sys
import argparse
import logging
import time


import fastavro

from irbase.utils import open_compressed

def main():
    parser = argparse.ArgumentParser(description='list the parses and lineages present in the given seqrec file')
    # input files
    parser.add_argument('seqrec_filename', metavar='seqrec.avro', help='the Avro file with the sequence records')
    parser.add_argument('inspect_count', metavar='N', type=int, default=1, nargs='?', help='number of records to inspect')
    # options
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])

    with open_compressed(args.seqrec_filename, 'rb') as seq_record_handle:
        seq_record_reader = fastavro.reader(seq_record_handle)

        parses = set()
        lineages = set()
        for record, c in zip(seq_record_reader, range(args.inspect_count)):
            parses.update(record['parses'].keys())
            lineages.update(record['lineages'].keys())
        
    if len(parses) > 0:
        print('parses:')
        for p in parses:
            print('\t', p)
    if len(lineages) > 0:
        print('lineages:')
        for l in lineages:
            print('\t', l)

if __name__ == '__main__':
    sys.exit(main())