#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
import logging
import time

from Bio import SeqIO
from fastavro import parse_schema
from fastavro import writer

from irbase.utils import open_compressed
from irbase.schemata.avro import SEQUENCE_RECORD

def seqrec_from_genbank(genbank_record):
    sequence_annotations = {'organism':    genbank_record.annotations['organism'],
                            'description': genbank_record.description}
    sequence = {'sequence': str(genbank_record.seq),
                'qual': '',
                'annotations': sequence_annotations}
    return {'name': 'genbank:' + genbank_record.id,
            'source': 'Genbank',
            'subject': None,
            'sample': None,
            'sequence': sequence}

def genbank_filter_chain(filenames, organism=None, max_length=None):
    processed_read_count = 0
    for genbank_filename in filenames:
        logging.info('processing %s', genbank_filename)
        with open_compressed(genbank_filename, 'rt') as genbank_file:
            for record in SeqIO.parse(genbank_file, 'genbank'):
                if (organism is None) or (organism == record.annotations['organism']):
                    sequence_length = len(record.seq)
                    if (max_length is None) or (sequence_length <= max_length):
                        processed_read_count += 1
                        yield seqrec_from_genbank(record)
                    else:
                        logging.info('record %s (%s) is too big, %d > %s, ignoring',
                                record.id, record.description, sequence_length, max_length)

    logging.info('processed %s reads', processed_read_count)


def main():
    parser = argparse.ArgumentParser(description='loads Genbank sequences',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # input files
    parser.add_argument('genbank_filenames', metavar='genbank_file', nargs='+', help='the file with the Genbank records')
    # options
    parser.add_argument('--organism', '-o', metavar='O', help='only process records with the given organism')
    parser.add_argument('--max-length', '-m', metavar='L', type=int, default=50000, help='ignore sequences longer than this')
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')
    # options for avro
    avro_group = parser.add_argument_group('Avro options')
    avro_group.add_argument('--codec', choices=BLOCK_READERS.keys(), default='snappy', help='compression codec to use')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])
    start_time = time.time()

    # setup schema
    avro_schema = parse_schema(SEQUENCE_RECORD)

    genbank_records = genbank_filter_chain(args.genbank_filenames, args.organism, args.max_length)
    writer(sys.stdout.buffer, avro_schema, genbank_records, codec=args.codec)

    elapsed_time = time.time() - start_time
    logging.info('elapsed time %s', time.strftime('%H hours, %M minutes, %S seconds', time.gmtime(elapsed_time)))
    
if __name__ == '__main__':
    sys.exit(main())
