#!/usr/bin/env python

import sys
import argparse
import logging
import time
import csv
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from fastavro import parse_schema, writer
from fastavro.read import BLOCK_READERS

from irbase.utils import open_compressed
from irbase.schemata.avro import SEQUENCE_RECORD

def make_seq_records(sequences, subject=None, sample=None, source=None):
    for name, sequence, qual in sequences:
        sequence = {'sequence': sequence,
                    'qual':     qual}

        record = {'name': name,
                  'subject': subject,
                  'sample': sample,
                  'source': source,
                  'sequence': sequence
                 }
        
        yield record

def fasta_chain(fasta_filenames):
    for filename in fasta_filenames:
        with open_compressed(filename, 'rt') as handle:
            logging.info('processing FASTA %s', filename)
            for id, seq in SimpleFastaParser(handle):
                yield id, seq, None

def fastq_chain(fasta_filenames):
    for filename in fasta_filenames:
        with open_compressed(filename, 'rt') as handle:
            logging.info('processing FASTQ %s', filename)
            yield from FastqGeneralIterator(handle)

def main():
    parser = argparse.ArgumentParser(description='convert FASTA/FASTQ file into SeqRec', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # paths
    parser.add_argument('fastx_filenames', metavar='sequences.fastx', type=Path, nargs='+', help='FASTA/FASTQ to convert')
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')
    # input options
    input_format = parser.add_argument_group('input options').add_mutually_exclusive_group()
    input_format.add_argument('--fasta', '-a', default=True, action='store_true', help='output a FASTA file')
    input_format.add_argument('--fastq', '-q', action='store_false', dest='fasta', help='output a FASTQ file')
    # values for metadata
    meta_group = parser.add_argument_group('metadata options')
    meta_group.add_argument('--subject', metavar='subject', default=None, help='the subject to label this data')
    meta_group.add_argument('--sample',  metavar='sample',  default=None, help='the sample to label this data')
    meta_group.add_argument('--source',  metavar='source',  default=None, help='the source to label this data')
    # options for avro
    avro_group = parser.add_argument_group('Avro options')
    avro_group.add_argument('--codec', choices=BLOCK_READERS.keys(), default='snappy', help='compression codec to use')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])
    start_time = time.time()
    
    output_schema = parse_schema(SEQUENCE_RECORD)

    if args.fasta:
        input_sequences = fasta_chain(args.fastx_filenames) 
    else:
        input_sequences = fastq_chain(args.fastx_filenames)

    seq_records = make_seq_records(input_sequences, subject=args.subject, sample=args.sample, source=args.source)

    writer(sys.stdout.buffer, output_schema, seq_records, codec=args.codec)

    elapsed_time = time.time() - start_time
    logging.info('elapsed time %s', time.strftime('%H hours, %M minutes, %S seconds', time.gmtime(elapsed_time)))
    
if __name__ == '__main__':
    sys.exit(main())
