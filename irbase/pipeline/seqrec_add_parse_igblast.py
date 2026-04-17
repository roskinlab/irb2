#!/usr/bin/env python

import sys
import argparse
import logging
import time
from collections import defaultdict

import fastavro
from fastavro.read import BLOCK_READERS
from Bio import SeqIO

from irbase.utils import open_compressed, make_range, make_named_range
from irbase.parsers.igblast import IgBLASTParser
from irbase.schemata.avro import SEQUENCE_RECORD

def get_padding(seq):
    length = len(seq)
    seq = seq.lstrip('-')
    start_padding = length - len(seq)
    length = len(seq)
    seq = seq.rstrip('-')
    stop_padding = length - len(seq)
    return start_padding, seq, stop_padding

def igblast_annotator(germline_lengths, seq_record_iter, igblast_iter, parse_label, min_v_score, min_j_score):
    for record, parse in zip(seq_record_iter, igblast_iter):
        # make sure the sequence record name matches the IgBLAST record name
        assert record['name'] == parse.query_name, '%s != %s' % (record['name'], parse.query_name)
        assert parse_label not in record['parses']

        if not parse:
            record['parses'][parse_label] = None
            yield record
        else:
            # form the basic parse record
            parse_record = {'modified_sequence': None,
                            'chain': parse.chain_type,
                            'has_stop_codon': parse.stop_codon,
                            'v_j_in_frame': parse.v_j_in_frame,
                            'positive_strand': parse.strand == '+',
                            'v_frame_shift': parse.v_frame_shift,
                            'alignments': [],
                            'annotations': {}}

            # add the query sequence to the list of alignments
            query_alignment = parse.alignment_lines[0]
            assert query_alignment.segment_type == 'Q'
            parse_record['alignments'].append(
                    {'type': 'Q',
                    'name': '',
                    'length': parse.query_length,
                    'score': float('nan'),
                    'e_value': float('nan'),
                    'range': make_range(query_alignment.start, query_alignment.end),
                    'padding': make_range(0, 0),
                    'alignment': query_alignment.line
                    })
            
            # name, start, and stop for the top scoring segment type
            top_vdj_ranges = {}

            # process the alignments and keep best scores for each segment
            best_scores = defaultdict(float)
            for align_line in parse.alignment_lines[1:]:
                segment_type = align_line.segment_type
                align_score = parse.significant_alignments[align_line.name]
                start_padding, trimmed_line, stop_padding = get_padding(align_line.line)
                parse_record['alignments'].append(
                        {'type': segment_type,
                        'name': align_line.name,
                        'length': germline_lengths[align_line.name],
                        'score': align_score.bit_score,
                        'e_value': align_score.e_value,
                        'range': make_range(align_line.start, align_line.end),
                        'padding': make_range(start_padding, stop_padding),
                        'alignment': trimmed_line
                        })
                # save the best score for each segment type
                best_scores[segment_type] = max(best_scores[segment_type], align_score.bit_score)
                # store the name and range of the first observed segment type
                if segment_type not in top_vdj_ranges:
                    top_vdj_ranges[segment_type] = make_named_range(align_line.name, start_padding, start_padding + len(trimmed_line))

            # if the scores aren't good enough, return a null parse
            if best_scores['V'] < min_v_score or best_scores['J'] < min_j_score:
                record['parses'][parse_label] = None
                yield record
            else:
                # add the ranges for the V-, D-, J-segment
                for segment_type, named_range in top_vdj_ranges.items():
                    assert segment_type not in parse_record['annotations']
                    parse_record['annotations'][segment_type] = named_range

                # add in the ranges for the regions
                for region_name, region_range in parse.alignment_regions.items():
                    assert region_name not in parse_record['annotations']
                    parse_record['annotations'][region_name] = make_range(region_range.start, region_range.stop)

                # store the parses
                record['parses'][parse_label] = parse_record

                yield record


def igblast_chain(igblast_filenames):
    for filename in igblast_filenames:
        with open_compressed(filename, 'rt') as igblast_handle:
            logging.info('processing parsed in %s', filename)
            igblast_parse_reader = IgBLASTParser(igblast_handle)
            for record in igblast_parse_reader:
                yield record


def main():
    parser = argparse.ArgumentParser(description='load IgBLAST annotations into an Avro sequence record',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # input files
    parser.add_argument('parse_label', metavar='label', help='the parse label to use for the parse')
    parser.add_argument('repertoire_filenames', metavar='repertoire-file', nargs=4, help='the V(D)JC repertoire file used in IgBLAST')
    parser.add_argument('seqrec_filename', metavar='seq_record.avro', help='the Avro file with the sequence records')
    parser.add_argument('igblast_output_filenames', metavar='parse.igblast', nargs='+', help='the output of IgBLAST to parse and attach to the sequence record')
    # options
    parser.add_argument('--min-v-score', metavar='S', type=float, default=70.0, help='the minimum score for the V-segment')
    parser.add_argument('--min-j-score', metavar='S', type=float, default=26.0, help='the minimum score for the J-segment')
    parser.add_argument('--log-level', '-l', choices=[i.lower() for i in reversed(logging.getLevelNamesMapping())],
        default='error', help='set logging level')
    # options for avro
    avro_group = parser.add_argument_group('Avro options')
    avro_group.add_argument('--codec', choices=BLOCK_READERS.keys(), default='snappy', help='compression codec to use')

    args = parser.parse_args()
    logging.basicConfig(level=logging.getLevelNamesMapping()[args.log_level.upper()])
    start_time = time.time()

    output_schema = fastavro.parse_schema(SEQUENCE_RECORD)

    logging.info('calculating V(D)J repertoire lengths')
    germline_lengths = {}
    for rep_filename in args.repertoire_filenames:
        with open(rep_filename, 'rt') as rep_handle:
            for record in SeqIO.parse(rep_handle, 'fasta'):
                germline_lengths[record.id] = len(record)

    logging.info('adding parses to sequence records')

    with open_compressed(args.seqrec_filename, 'rb') as seq_record_handle:
        seq_record_reader = fastavro.reader(seq_record_handle)
        igblast_parse_reader = igblast_chain(args.igblast_output_filenames)

        annotator = igblast_annotator(germline_lengths, seq_record_reader, igblast_parse_reader, args.parse_label,
                                    args.min_v_score, args.min_j_score)

        fastavro.writer(sys.stdout.buffer, output_schema, annotator, codec=args.codec)

    elapsed_time = time.time() - start_time
    logging.info('elapsed time %s', time.strftime('%H hours, %M minutes, %S seconds', time.gmtime(elapsed_time)))
    
if __name__ == '__main__':
    sys.exit(main())
