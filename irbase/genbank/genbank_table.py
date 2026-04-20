#!/usr/bin/env python

import sys
import argparse
import csv
from collections import defaultdict

import fastavro
from Bio.Seq import Seq
from test.test_dataclasses.dataclass_module_1_str import annotations

from irbase.utils import open_compressed

def best_vdj_score(parse):
    '''extracts the best segment name and score given a parse record'''
    best_v_name  = None
    best_v_score = None
    best_d_name  = None
    best_d_score = None
    best_j_name  = None
    best_j_score = None

    if parse is not None:
        for align in parse['alignments']:
            if align['type'] == 'V':
                if best_v_score is None:
                    best_v_name = align['name']
                    best_v_score = align['score']
            elif align['type'] == 'D':
                if best_d_score is None:
                    best_d_name = align['name']
                    best_d_score = align['score']
            elif align['type'] == 'J':
                if best_j_score is None:
                    best_j_name = align['name']
                    best_j_score = align['score']

    return best_v_name, best_v_score, \
           best_d_name, best_d_score, \
           best_j_name, best_j_score
      
def get_inframe_aa_seq(parse, anchor_region) -> None | str:
    if parse is None or anchor_region not in parse['annotations']:
        return None

    assert parse['alignments'][0]['type'] == 'Q'
    sequence = parse['alignments'][0]['alignment']

    # get the start position of given anchor region, which should be in frame
    region_start = parse['annotations'][anchor_region]['start']

    # get the sequence before and after the anchor point
    prefix = sequence[:region_start].replace('-', '')
    postfix = sequence[region_start:].replace('-', '')

    # pad sequences to be in frame
    prefix = 'N' * (-len(prefix) % 3) + prefix
    assert len(prefix) % 3 == 0
    postfix = postfix + 'N' * (-len(postfix) % 3)

    # translate the sequence
    assert len(postfix) % 3 == 0
    sequence = Seq(prefix + postfix)
    translated_sequence = str(sequence.translate())

    return translated_sequence

def get_cdr3_seq(parse) -> None | str:
    if parse is None or 'CDR3' not in parse['annotations']:
        return None, None

    assert parse['alignments'][0]['type'] == 'Q'
    sequence = parse['alignments'][0]['alignment']

    # get the start and stop position of the CDR3
    cdr3_start = parse['annotations']['CDR3']['start']
    cdr3_stop = parse['annotations']['CDR3']['stop']

    cdr3_nt = sequence[cdr3_start: cdr3_stop].replace('-', '')
    cdr3_nt += (-len(cdr3_nt) % 3) * 'N'
    assert len(cdr3_nt) % 3 == 0, cdr3_nt

    return cdr3_nt, str(Seq(cdr3_nt).translate())

def main():
    parser = argparse.ArgumentParser(description='',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('parse_label', metavar='label', help='the parse label to use for the parse')
    parser.add_argument('filenames', metavar='file', nargs='+', help='the Avro file to read')
    parser.add_argument('-a', '--anchor-region', metavar='region', default= 'CDR3', help= 'anchor region for CDR parsing and framework adjustment')
    args = parser.parse_args()

    writer = csv.DictWriter(sys.stdout, fieldnames=[
        'accession', 'description',
        #'curation_group', 'curation_year',
        'v_name', 'v_score', 'd_name', 'd_score', 'j_name', 'j_score',
        'sequence_nt', 'sequence_aa', 'cdr3_nt', 'cdr3_aa'
        ])
    writer.writeheader()

    for filename in args.filenames:
        with open_compressed(filename, 'rb') as read_handle:
            reader = fastavro.reader(read_handle)

            for record in reader:
                parse = record['parses'][args.parse_label]
                if parse is not None:

                    description = None
                    if 'description' in record['sequence']['annotations']:
                        description = record['sequence']['annotations']['description']

                    v_name, v_score, d_name, d_score, j_name, j_score = best_vdj_score(parse)

                    assert parse['alignments'][0]['type'] == 'Q'
                    sequence_nt = parse['alignments'][0]['alignment']

                    try:
                        sequence_aa = get_inframe_aa_seq(parse, args.anchor_region)
                        cdr3_nt, cdr3_aa = get_cdr3_seq(parse)
                    except:
                        print('error with ' + record['name'], file=sys.stderr)
                        raise 

                    row = {'accession': record['name'],
                           'description': description,
                           #'curation_group': record['sequence']['annotations']['curation_group'],
                           #'curation_year': record['sequence']['annotations']['curation_year'],
                           'sequence_nt': sequence_nt,
                           'sequence_aa': sequence_aa,
                           'v_name': v_name,
                           'v_score': v_score,
                           'd_name': d_name,
                           'd_score': d_score,
                           'j_name': j_name,
                           'j_score': j_score,
                           'cdr3_nt': cdr3_nt,
                           'cdr3_aa': cdr3_aa
                        }

                    writer.writerow(row)


if __name__ == '__main__':
    sys.exit(main())
