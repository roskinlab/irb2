#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse

import fastavro
from collections import defaultdict

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


def main():
    parser = argparse.ArgumentParser(description='',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('parse_label', metavar='label', help='the parse label to use for the parse')
    parser.add_argument('filenames', metavar='file', nargs='+', help='the Avro file to read')
    args = parser.parse_args()

    print('accession', 'description', 'v_name', 'v_score', 'd_name', 'd_score', 'j_name', 'j_score', sep='\t')

    for filename in args.filenames:
        with open_compressed(filename, 'rb') as read_handle:
            reader = fastavro.reader(read_handle)

            for record in reader:
                name = record['name']
                assert name.startswith('genbank:')
                accession = name.split(':')[1]

                parse = record['parses'][args.parse_label]
                v_name, v_score, d_name, d_score, j_name, j_score = best_vdj_score(parse)
                description = None
                if 'description' in record['sequence']['annotations']:
                    description = record['sequence']['annotations']['description']

                print(accession, description, v_name, v_score, d_name, d_score, j_name, j_score, sep='\t')


if __name__ == '__main__':
    sys.exit(main())
