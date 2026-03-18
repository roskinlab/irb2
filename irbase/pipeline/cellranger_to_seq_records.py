#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
import logging
import time
import csv
from pathlib import Path

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from fastavro import parse_schema, writer

from irbase.utils import open_compressed
from irbase.schemata.avro import SEQUENCE_RECORD

def make_seq_records(cell_umis, cell_contigs, sequences, annotations, subject=None, sample=None, source=None, drop_contig=False):
    for (cell, chain) in cell_umis:
        name = cell_contigs[(cell, chain)]
        if drop_contig:
            name = name.split('_')[0]

        contig_seq, contig_qual = sequences[cell_contigs[(cell, chain)]]
        # form the sequence objects with annotations
        sequence = {'sequence': contig_seq,
                    'qual':     contig_qual,
                    'annotations': annotations[(cell, chain)], 
                    'ranges': {}}

        record = {'name': name,
                  'subject': subject,
                  'sample': sample,
                  'source': source,
                  'sequence': sequence,
                  'parses': {},
                  'lineages': {}
                 }
        
        yield record

def main():
    parser = argparse.ArgumentParser(description='select chain from Cell Ranger VDJ outout with the most UMIs and convert to IRB format', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # paths
    parser.add_argument('cell_ranger_vdj_path', metavar='cellranger_outs_dir', default='.', type=Path, help='directory with the output of Cell Ranger VDJ pipeline')
    # options to processes data
    parser.add_argument('--drop-contig', '-d', action='store_true', default=False, help='drop the contig name from the read name')
    parser.add_argument('--group', '-g', metavar='G', default='all', help='which set of VDJ contigs or annotations to load')
    # subset based on chain type
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--heavy-only',  action='store_true', help='only parse heavy-chains')
    group.add_argument('--light-only',  action='store_true', help='only parse light-chains')
    group.add_argument('--both-chains', action='store_true', help='parse both chains')
    # values for metadata
    meta_group = parser.add_argument_group('metadata options')
    meta_group.add_argument('--subject', metavar='subject', default=None, help='the subject to label this data')
    meta_group.add_argument('--sample',  metavar='sample',  default=None, help='the sample to label this data')
    meta_group.add_argument('--source',  metavar='source',  default=None, help='the source to label this data')
    

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    start_time = time.time()

    logging.info(args)

    if args.light_only:
        chains_to_parse = {'IGK', 'IGL'}
    elif args.heavy_only:
        chains_to_parse = {'IGH'}
    else:
        chains_to_parse = {'IGH', 'IGK', 'IGL'}

    logging.info(','.join(chains_to_parse))

    output_schema = parse_schema(SEQUENCE_RECORD)

    with open_compressed(args.cell_ranger_vdj_path / (args.group + '_contig_annotations.csv'), 'rt') as filtered_contig, \
         open_compressed(args.cell_ranger_vdj_path / (args.group + '_contig.fastq'), 'rt') as filtered_sequences:
        
        cell_umi_count = {} # max UMI count for this cell, chain
        annotations    = {} # seq. annotations for cell, chaim with max UMI
        cell_contig    = {} # contig name with the max UMI count

        for row in csv.DictReader(filtered_contig):
            chain = row['chain']
            # some basic filteres
            if row['is_cell'] == 'true' and row['high_confidence'] == 'true' and \
               row['productive'] == 'true' and (chain in chains_to_parse):
                # extract columns
                cell = row['barcode']
                umis = int(row['umis'])
                # for each cell, record the contig with the most UMIs
                if (cell, chain) in cell_umi_count:
                    if umis >= cell_umi_count[(cell, chain)]:
                        cell_umi_count[(cell, chain)] = umis
                        cell_contig[(cell, chain)] = row['contig_id']
                        annotations[(cell, chain)] = {
                            'cr_umis': row['umis'],
                            'cr_reads': row['reads'],
                            'cr_chain': row['chain'],
                            'cr_clonotype_id': row['raw_clonotype_id']}
                else:
                    cell_umi_count[(cell, chain)] = umis
                    cell_contig[(cell, chain)] = row['contig_id']
                    annotations[(cell, chain)] = {
                            'cr_umis': row['umis'],
                            'cr_reads': row['reads'],
                            'cr_chain': row['chain'],
                            'cr_clonotype_id': row['raw_clonotype_id']}

        # load the sequences for each contig
        sequences = {}
        for contig_id, contig_seq, contig_qual in FastqGeneralIterator(filtered_sequences):
            sequences[contig_id] = (contig_seq, contig_qual)

        seq_records = make_seq_records(cell_umi_count, cell_contig, sequences, annotations,
                subject=args.subject, sample=args.sample, source=args.source, drop_contig=args.drop_contig)

        writer(sys.stdout.buffer, output_schema, seq_records, codec='bzip2')

    elapsed_time = time.time() - start_time
    logging.info('elapsed time %s', time.strftime('%H hours, %M minutes, %S seconds', time.gmtime(elapsed_time)))
    
if __name__ == '__main__':
    sys.exit(main())
