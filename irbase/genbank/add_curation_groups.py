#!/usr/bin/env python

import sys
import argparse
import logging

import fastavro
from fastavro.read import BLOCK_READERS
from collections import defaultdict

from irbase.utils import open_compressed
from irbase.schemata.avro import SEQUENCE_RECORD

def curation_group_annotator(curation_data, seq_record_iter):
    for seq_record in seq_record_iter:
        genbank, name = seq_record['name'].split(':')
        assert genbank == 'genbank'

        year, group = curation_data[name]
        seq_record['sequence']['annotations']['curation_year'] = year
        seq_record['sequence']['annotations']['curation_group'] = group

        yield seq_record

def main():
    parser = argparse.ArgumentParser(description='',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('seqrec_filename', metavar='seqrec.avro', help='the Avro file to add annotations to')
    parser.add_argument('curation_groups_filename', metavar='groups_table.tsv', help='TSV file with the curation groups')
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')
    # column options
    column_group = parser.add_argument_group('Curation column options')
    column_group.add_argument('--curation-year',      '-y', metavar='N', default=0, help='which column index has the curation year')
    column_group.add_argument('--curation-group',     '-g', metavar='N', default=1, help='which column index has the curation group')
    column_group.add_argument('--curation-accession', '-a', metavar='N', default=2, help='which column index has the accesss id')
    # options for avro
    avro_group = parser.add_argument_group('Avro options')
    avro_group.add_argument('--codec', choices=BLOCK_READERS.keys(), default='snappy', help='compression codec to use')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])

    output_schema = fastavro.parse_schema(SEQUENCE_RECORD)

    # load in curation group data
    curation_data = {}
    with open_compressed(args.curation_groups_filename, 'rt') as read_handle:
        for row in read_handle:
            row = row.split('\t')
            curation_year      = int(row[args.curation_year])
            curation_group     = row[args.curation_group]
            curation_accession = row[args.curation_accession]
            
            curation_data[curation_accession] = (curation_year, curation_group)

    with open_compressed(args.seqrec_filename, 'rb') as seq_record_handle:
        seq_record_reader = fastavro.reader(seq_record_handle)

        annotator = curation_group_annotator(curation_data, seq_record_reader)
        fastavro.writer(sys.stdout.buffer, output_schema, annotator, codec=args.codec)

if __name__ == '__main__':
    sys.exit(main())
