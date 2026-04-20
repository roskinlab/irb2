#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
import logging
import time
import os
import io
import csv
from itertools import zip_longest
from collections import defaultdict

from Bio import SeqIO
import xlsxwriter

from irbase.utils import open_compressed


feature_whitelist = set([
    ('CDS', 'experiment'),
    ('CDS', 'function'),
    ('CDS', 'note'),
    ('CDS', 'product'),
    ('C_region', 'gene'),
    ('C_region', 'note'),
    ('C_region', 'product'),
    ('mat_peptide', 'product'),
    ('misc_feature', 'note'),
    ('mRNA', 'product'),
    ('source', 'bio_material'),
    ('source', 'cell_line'),
    ('source', 'cell_type'),
    ('source', 'clone'),
    ('source', 'clone_lib'),
    ('source', 'dev_stage'),
    ('source', 'isolate'),
    ('source', 'isolation_source'),
    ('source', 'lab_host'),
    ('source', 'mol_type'),
    ('source', 'note'),
    ('source', 'serotype'),
    ('source', 'sex'),
    ('source', 'specimen_voucher'),
    ('source', 'sub_clone'),
    ('source', 'tissue_lib'),
    ('source', 'tissue_type')
])

def get_features(record):
    source_features = []
    other_features  = []
    for feature in record.features:
        for key, value in feature.qualifiers.items():
            if (feature.type, key) in feature_whitelist:
                if feature.type == 'source':
                    source_features.append('%s=%s' % (key, '|'.join(value)))
                else:
                    other_features.append('%s:%s=%s' % (feature.type, key, '|'.join(value)))
            if feature.type == 'CDS' and key == 'db_xref':
                for v in value:
                    if v.startswith('HSSP:') or v.startswith('PDB:'):
                        other_features.append('%s:%s=%s' % (feature.type, key, v))

    return source_features, other_features

def equal_references(refs1, refs2):
    for r1, r2 in zip_longest(refs1, refs2):
        if r1.authors != r2.authors:
            return False
        if r1.title != r2.title:
            return False
        if r1.journal != r2.journal:
            # allow 2 mismatches for the journal entry
            mismatch_count = 0
            for c1, c2, in zip_longest(r1.journal, r2.journal):
                if c1 != c2:
                    mismatch_count += 1
                    if mismatch_count > 2:
                        return False
    return True

def get_master_references(records):
    all_references = set()

    direct_journal_dates = defaultdict(set)
    direct_author = {}
    direct_title = {}

    for record in records:
        for reference in record.annotations['references']:
            if reference.title == 'Direct Submission':
                # trim and extract the date from the "journal"
                journal =  reference.journal
                assert journal.startswith('Submitted (')
                date = journal[journal.index('(') + 1:journal.index(')')]
                journal = journal[journal.index(')') + 1:]

                direct_journal_dates[journal].add(date)
                direct_author[journal] = reference.authors
                direct_title[journal]  = reference.title
            else:
                all_references.add((reference.authors, reference.title, reference.journal, reference.pubmed_id))

    for journal in direct_journal_dates:
        all_references.add((direct_author[journal], direct_title[journal],
                            'Submitted (' + ','.join(direct_journal_dates[journal]) + ')' + journal, ''))

    return all_references

def write_references(workbook, worksheet, references, current_row=0):
    #format_dark       = workbook.add_format({'bg_color': '#CCCCCC'})
    #format_dark_bold  = workbook.add_format({'bold': True, 'bg_color': '#CCCCCC'})
    #format_light      = workbook.add_format({'bg_color': '#EEEEEE'})
    #format_light_bold = workbook.add_format({'bold': True, 'bg_color': '#EEEEEE'})

    #formats = [format_light, format_dark]
    #formats_bold = [format_light_bold, format_dark_bold]

    worksheet.write_string(current_row, 0, 'Authors')
    worksheet.write_string(current_row, 1, 'Title')
    worksheet.write_string(current_row, 2, 'Journal')
    worksheet.write_string(current_row, 3, 'PMID')

    for authors, title, journal, pubmed_id in references:
        current_row += 1

        worksheet.write_string(current_row, 0, authors)
        worksheet.write_string(current_row, 1, title)
        worksheet.write_string(current_row, 2, journal)
        worksheet.write_url(current_row, 3, 'https://www.ncbi.nlm.nih.gov/pubmed/' + pubmed_id, string=pubmed_id)

    current_row += 2
    worksheet.write_string(current_row, 0, 'Notes')

def write_curation_header(workbook, worksheet, header_row=0):
    # set the header formatting
    #format_merge_dark = workbook.add_format({'bold': True, 'center_across': True, 'bg_color': '#CCCCCC'})
    #format_merge_light = workbook.add_format({'bold': True, 'center_across': True, 'bg_color': '#EEEEEE'})
    #worksheet.merge_range(header_row,  0, header_row,  3, 'Genbank Features',  format_merge_dark)
    #worksheet.merge_range(header_row,  4, header_row,  7, 'IgBLAST Features',  format_merge_light)
    #worksheet.merge_range(header_row,  8, header_row,  9, 'Ann. Notes',        format_merge_dark)
    #worksheet.merge_range(header_row, 10, header_row, 14, 'Subject Features',  format_merge_light)
    #worksheet.merge_range(header_row, 15, header_row, 19, 'Sample Features',   format_merge_dark)
    #worksheet.merge_range(header_row, 20, header_row, 24, 'Cell Origin',       format_merge_light)
    #worksheet.merge_range(header_row, 25, header_row, 27, 'Sequence Method',   format_merge_dark)
    #worksheet.merge_range(header_row, 28, header_row, 34, 'Antibody Features', format_merge_light)
    #worksheet.write(      header_row, 35,                 'Notes',             format_merge_dark)
    #worksheet.set_row(header_row + 1, 30)

    #worksheet.set_column( 4,  4,  12.0)
    #worksheet.set_column( 5,  5,   5.5)
    #worksheet.set_column( 6,  6,  12.0)
    #worksheet.set_column( 7,  7,   5.5)

    #header_row += 1

    format_bold_dark = workbook.add_format({'bold': True, 'bg_color': '#CCCCCC'})
    format_bold_light = workbook.add_format({'bold': True, 'bg_color': '#EEEEEE'})


    # Genbank Features
    worksheet.write_string(header_row,  0, 'Accession', format_bold_dark)
    worksheet.write_string(header_row,  1, 'Description', format_bold_dark)
    worksheet.write_string(header_row,  2, 'Source Features', format_bold_dark)
    worksheet.write_string(header_row,  3, 'Other Features', format_bold_dark)
    # IgBLAST Features
    worksheet.write_string(header_row,  4, 'V-segment', format_bold_light)
    worksheet.write_string(header_row,  5, 'V-score', format_bold_light)
    worksheet.write_string(header_row,  6, 'J-segment', format_bold_light)
    worksheet.write_string(header_row,  7, 'J-score', format_bold_light)
    # Annotation fields
    worksheet.write_string(header_row,  8, 'Seq. Cat.', format_bold_dark)
    worksheet.write_string(header_row,  9, 'Ann. Level', format_bold_dark)
    # Subject Features
    worksheet.write_string(header_row, 10, 'Subject Label', format_bold_light)
    worksheet.write_string(header_row, 11, 'Sex', format_bold_light)
    worksheet.write_string(header_row, 12, 'Race', format_bold_light)
    worksheet.write_string(header_row, 13, 'Ethnicity', format_bold_light)
    worksheet.write_string(header_row, 14, 'Genotype', format_bold_light)
    # Sample Features
    worksheet.write_string(header_row, 15, 'Sample Label', format_bold_dark)
    worksheet.write_string(header_row, 16, 'Age', format_bold_dark)
    worksheet.write_string(header_row, 17, 'Disease', format_bold_dark)
    worksheet.write_string(header_row, 18, 'Time Point', format_bold_dark)
    worksheet.write_string(header_row, 19, 'Intervention', format_bold_dark)
    # Cell Origin
    worksheet.write_string(header_row, 20, 'Tissue', format_bold_light)
    worksheet.write_string(header_row, 21, 'Isolation\nMethod', format_bold_light)
    worksheet.write_string(header_row, 22, 'Cell Type', format_bold_light)
    worksheet.write_string(header_row, 23, 'Cell Type\nAssayed', format_bold_light)
    worksheet.write_string(header_row, 24, 'Immortalized', format_bold_light)
    # Sequence Method
    worksheet.write_string(header_row, 25, 'Template', format_bold_dark)
    worksheet.write_string(header_row, 26, 'Amp.\nMethod', format_bold_dark)
    worksheet.write_string(header_row, 27, 'Seq.\n Method', format_bold_dark)
    # Antibody Features
    worksheet.write_string(header_row, 28, 'Antibody Label', format_bold_light)
    worksheet.write_string(header_row, 29, 'Isotype', format_bold_light)
    worksheet.write_string(header_row, 30, 'Isotype\nAssayed', format_bold_light)
    worksheet.write_string(header_row, 31, 'Specificity', format_bold_light)
    worksheet.write_string(header_row, 32, 'Specificity\nAssayed', format_bold_light)
    worksheet.write_string(header_row, 33, 'Binding\nAffinity', format_bold_light)
    worksheet.write_string(header_row, 34, 'Autoantibody', format_bold_light)
    # Notes
    worksheet.write_string(header_row, 35, 'Notes', format_bold_dark)

    header_row += 1

    worksheet.freeze_panes(header_row, 1)

    return header_row

def write_curation_row(workbook, worksheet, records, igblast_annotations, current_row):
    format_dark = workbook.add_format({'bg_color': '#CCCCCC'})
    format_dark_wrap = workbook.add_format({'bg_color': '#CCCCCC', 'text_wrap': True})
    format_light = workbook.add_format({'bg_color': '#EEEEEE'})
    format_red = workbook.add_format({'bg_color': '#EE0000'})

    for record in records:
        # output the name with URL
        name = record.id.split('.')[0]
        if worksheet.hlink_count < 65530:
            worksheet.write_url(current_row, 0, 'https://www.ncbi.nlm.nih.gov/nuccore/' + name, string=name, cell_format=format_dark)
        else:
            worksheet.write(current_row, 0, name, format_dark)

        # clean up and write description
        description = record.description
        if description.startswith('Homo sapiens '):
            description = description[len('Homo sapiens '):]
        worksheet.write_string(current_row, 1, description, format_dark_wrap)

        # output features about the source and other features
        source_features, other_features = get_features(record)
        worksheet.write_string(current_row, 2, '\n'.join(source_features), format_dark)
        worksheet.write_string(current_row, 3, '\n'.join(other_features), format_dark)

        # outout IgBLAST annotations
        igblast_annotations
        v_name, v_score, d_name, d_score, j_name, j_score = igblast_annotations[name]

        if v_name is None:
            if worksheet.hlink_count < 65530:
                worksheet.write_url(current_row,  4, 'https://www.ncbi.nlm.nih.gov/igblast/igblast.cgi?germline_db_V=IG_DB%2Fimgt.Homo_sapiens.V.f.orf.p&germline_db_D=IG_DB%2Fimgt.Homo_sapiens.D.f.orf&germline_db_J=IG_DB%2Fimgt.Homo_sapiens.J.f.orf&germline_db_C=IG_DB%2Fimgt.Homo_sapiens.C&num_alignments_V=3&num_alignments_D=3&num_alignments_J=3&num_alignments_C=2&translation=true&domain=imgt&num_clonotype=0&CMD=request&queryseq=' \
                        + name, string='-', cell_format=format_red)
            else:
                worksheet.write(current_row,  4, '-', format_red)
            worksheet.write(current_row,  5, v_score, format_red)
            worksheet.write(current_row,  6, '-', format_red)
            worksheet.write(current_row,  7, j_score, format_red)
        else:
            if worksheet.hlink_count < 65530:
                worksheet.write_url(current_row,  4, 'https://www.ncbi.nlm.nih.gov/igblast/igblast.cgi?germline_db_V=IG_DB%2Fimgt.Homo_sapiens.V.f.orf.p&germline_db_D=IG_DB%2Fimgt.Homo_sapiens.D.f.orf&germline_db_J=IG_DB%2Fimgt.Homo_sapiens.J.f.orf&CMD=request&queryseq=' \
                        + name, string=v_name, cell_format=format_light)
            else:
                worksheet.write(current_row,  4, v_name, format_light)
            worksheet.write(current_row,  5, round(v_score, 1), format_light)
            worksheet.write(current_row,  6, j_name, format_light)
            worksheet.write(current_row,  7, round(j_score, 1), format_light)

        # color code the groups of columns
        worksheet.write(current_row,  8, '', format_dark)
        worksheet.write(current_row,  9, '', format_dark)

        worksheet.write(current_row, 10, '', format_light)
        worksheet.write(current_row, 11, '', format_light)
        worksheet.write(current_row, 12, '', format_light)
        worksheet.write(current_row, 13, '', format_light)
        worksheet.write(current_row, 14, '', format_light)

        worksheet.write(current_row, 15, '', format_dark)
        worksheet.write(current_row, 16, '', format_dark)
        worksheet.write(current_row, 17, '', format_dark)
        worksheet.write(current_row, 18, '', format_dark)
        worksheet.write(current_row, 19, '', format_dark)

        worksheet.write(current_row, 20, '', format_light)
        worksheet.write(current_row, 21, '', format_light)
        worksheet.write(current_row, 22, '', format_light)
        worksheet.write(current_row, 23, '', format_light)
        worksheet.write(current_row, 24, '', format_light)

        worksheet.write(current_row, 25, '', format_dark)
        worksheet.write(current_row, 26, '', format_dark)
        worksheet.write(current_row, 27, '', format_dark)

        worksheet.write(current_row, 28, '', format_light)
        worksheet.write(current_row, 29, '', format_light)
        worksheet.write(current_row, 30, '', format_light)
        worksheet.write(current_row, 31, '', format_light)
        worksheet.write(current_row, 32, '', format_light)
        worksheet.write(current_row, 33, '', format_light)
        worksheet.write(current_row, 34, '', format_light)

        worksheet.write(current_row, 35, '', format_dark)

        worksheet.set_row(current_row, 60)

        current_row += 1

    #worksheet.set_column( 8, 12,  9.5)
    #worksheet.set_column(13, 17, 10.5)
    #worksheet.set_column(18, 22, 12.5)
    #worksheet.set_column(23, 25,  9.5)
    #worksheet.set_column(26, 32, 10.5)

    return current_row

def write_genbank_records(workbook, worksheet, records, current_row=0):
    for record in records:
        genbank_string = io.StringIO()

        SeqIO.write(record, genbank_string, 'genbank')

        worksheet.write_string(current_row, 0, genbank_string.getvalue())
        worksheet.set_row(current_row, 500)

        current_row += 1

    format_small_wrap = workbook.add_format({'font_size': 8, 'text_wrap': True})
    worksheet.set_column(0, 0, 65, format_small_wrap)

    return current_row

def get_igblast_annotations(igblast_table_filename, names, v_cutoff, j_cutoff):
    max_v_score = 0.0
    max_j_score = 0.0
    annotations = {}
    good_count = 0
    with open_compressed(igblast_table_filename, 'rt') as igblast_table_handle:
        for record in csv.DictReader(igblast_table_handle, delimiter='\t'):
            name = record['accession'].split('.')[0]
            if name in names:
                v_name = None if record['v_name'] == 'None' else record['v_name']
                d_name = None if record['d_name'] == 'None' else record['d_name']
                j_name = None if record['j_name'] == 'None' else record['j_name']

                v_score = 0.0 if record['v_score'] == 'None' else float(record['v_score'])
                d_score = 0.0 if record['d_score'] == 'None' else float(record['d_score'])
                j_score = 0.0 if record['j_score'] == 'None' else float(record['j_score'])

                max_v_score = max(max_v_score, v_score)
                max_j_score = max(max_j_score, j_score)

                if v_score >= v_cutoff and j_score >= j_cutoff:
                    good_count += 1

                annotations[name] = (v_name, v_score, d_name, d_score, j_name, j_score)
    return max_v_score, max_j_score, good_count, annotations

def main():
    parser = argparse.ArgumentParser(description='generate a sheet for Genbank immune receptor annotation',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # input file
    parser.add_argument('genbank_filename', metavar='genbank-file', help='the file with the Genbank records')
    parser.add_argument('igblast_table', metavar='igblast-tabel.tsv', help='the file with the IgBLAST annotations')
    # options
    parser.add_argument('--min-v-score', metavar='S', type=float, default=70.0, help='the minimum score for the V-segment')
    parser.add_argument('--min-j-score', metavar='S', type=float, default=26.0, help='the minimum score for the J-segment')

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    start_time = time.time()

    record_counts = 0

    excel_filename = args.genbank_filename
    if '.' in excel_filename:
        excel_filename = excel_filename[:excel_filename.rindex('.')]
    excel_filename += '.xlsx'

    stat_filename = args.genbank_filename
    if '.' in stat_filename:
        stat_filename = stat_filename[:stat_filename.rindex('.')]
    stat_filename += '.stat'

    current_row = 0
    
    # read in the Genbank records
    with open_compressed(args.genbank_filename, 'rt') as genbank_handle:
        # load all the records
        records = list(SeqIO.parse(genbank_handle, 'genbank'))
        names = set(r.id.split('.')[0] for r in records)
        max_v_score, max_j_score, good_count, igblast_annotations = get_igblast_annotations(args.igblast_table, names, v_cutoff=args.min_v_score, j_cutoff=args.min_j_score)
        if max_v_score < args.min_v_score or max_j_score < args.min_j_score:
            logging.info('group does not appear to contain IgH sequences (V-score: %f, J-score: %f', max_v_score, max_j_score)
        else:
            logging.info('group does appear to contain IgH sequence (V-score: %f, J-score: %f', max_v_score, max_j_score)
            # get a unique references list for all records
            references = get_master_references(records)

            # create the workbook
            workbook = xlsxwriter.Workbook(excel_filename)
            curation_worksheet = workbook.add_worksheet('Curation')
            references_worksheet = workbook.add_worksheet('References')
            
            # write the header
            current_row = write_curation_header(workbook, curation_worksheet)

            # write curation annotation
            current_row = write_curation_row(workbook, curation_worksheet, records, igblast_annotations, current_row)

            # write the references to the sheet
            write_references(workbook, references_worksheet, references)

            workbook.close()

            with open(stat_filename, 'wt') as stat_handle:
                all_title = ';'.join(r[1] for r in references)
                print(good_count, all_title, sep='\t', file=stat_handle)

    logging.info('wrote %s rows', current_row)
    elapsed_time = time.time() - start_time
    logging.info('elapsed time %s', time.strftime('%H hours, %M minutes, %S seconds', time.gmtime(elapsed_time)))
    
if __name__ == '__main__':
    sys.exit(main())
