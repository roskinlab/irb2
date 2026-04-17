#!/usr/bin/env python

import sys
import argparse
import logging
import time
import itertools
import io

import fastavro

from irbase.utils import open_compressed

def main():
    parser = argparse.ArgumentParser(description='convert a sequence records in Avro file format to FASTA or FASTQ')
    # input files
    parser.add_argument('seqrec_filenames', metavar='seqrec.avro', nargs='+', help='the Avro file with the sequence records')
    # options
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')
    # output options
    output_format = parser.add_argument_group('output options').add_mutually_exclusive_group()
    output_format.add_argument('--fasta', '-a', default=True, action='store_true', help='output a FASTA file')
    output_format.add_argument('--fastq', '-q', action='store_false', dest='fasta', help='output a FASTQ file')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])
    start_time = time.time()

    for input_filename in args.seqrec_filenames:
        with open_compressed(input_filename, 'rb') as seq_record_handle:
            seq_record_reader = fastavro.reader(seq_record_handle)

            for record in seq_record_reader:
                if args.fasta:
                    print('>%s\n%s' % (record['name'], record['sequence']['sequence']))
                else:
                    print('@%s\n%s\n+\n%s' % (record['name'], record['sequence']['sequence'], record['sequence']['qual']))

    elapsed_time = time.time() - start_time
    logging.info('elapsed time %s', time.strftime('%H hours, %M minutes, %S seconds', time.gmtime(elapsed_time)))

if __name__ == '__main__':
    sys.exit(main())