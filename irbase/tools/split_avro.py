#!/usr/bin/env python

import logging
import time
import sys
import argparse
from itertools import batched, count

import fastavro
from fastavro import writer
from fastavro.read import BLOCK_READERS

from irbase.utils import open_compressed

def main():
    parser = argparse.ArgumentParser(description='split an Avro file into batches', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # input files
    parser.add_argument('avro_filename', metavar='file.avro', help='the Avro files to batch')
    parser.add_argument('batch_filename_prefix', metavar='prefix', help='filename prefix of each batch')
    parser.add_argument('batch_filename_suffix', metavar='suffix', help='filename suffix of each batch')
    # options
    parser.add_argument('--batch-size', '-b', metavar='B', type=int, default=50000,
            help='the number of records per batch')
    parser.add_argument('--count-length', '-n', metavar='N', type=int, default=3,
            help='field width for the batch number')
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')
    # options for avro
    avro_group = parser.add_argument_group('Avro options')
    avro_group.add_argument('--codec', choices=BLOCK_READERS.keys(), default='snappy', help='compression codec to use')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])
    start_time = time.time()

    with open_compressed(args.avro_filename, 'rb') as read_handle:
        record_reader = fastavro.reader(read_handle)

        for b, n in zip(batched(record_reader, args.batch_size), count()):
            with open(args.batch_filename_prefix + f'{n:0{args.count_length}}' + args.batch_filename_suffix, 'wb') as write_handle:
                writer(write_handle, record_reader.writer_schema, b, codec=args.codec)

if __name__ == '__main__':
    sys.exit(main())