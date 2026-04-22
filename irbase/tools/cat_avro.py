#!/usr/bin/env python

import logging
import sys
import argparse

import fastavro
from fastavro import writer
from fastavro.read import BLOCK_READERS

from irbase.utils import open_compressed

def chain_avro(filenames):
    for filename in filenames:
        with open_compressed(filename, 'rb') as handle:
            for record in fastavro.reader(handle):
                yield record

def main():
    parser = argparse.ArgumentParser(description='split an Avro file into batches', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # input files
    parser.add_argument('avro_filenames', metavar='file.avro', nargs='+', help='the Avro files to concatinate')
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')
    # options for avro
    avro_group = parser.add_argument_group('Avro options')
    avro_group.add_argument('--codec', choices=BLOCK_READERS.keys(), default='snappy', help='compression codec to use')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])

    # get the schema from the first file
    with open_compressed(args.avro_filenames[0], 'rb') as read_handle:
        record_reader = fastavro.reader(read_handle)
        output_schema = record_reader.writer_schema

    records = chain_avro(args.avro_filenames)
    writer(sys.stdout.buffer, output_schema, records, codec=args.codec)

if __name__ == '__main__':
    sys.exit(main())