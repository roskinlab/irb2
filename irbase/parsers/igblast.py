import re
import sys

_re_space_sep     = re.compile(' {1,}')
_re_double_space  = re.compile(' {2,}')
_re_region_ranges = re.compile('<-*(?P<label>[^-]*)-?[^-]*-*>')

region_labels = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']

def upper_mismatch(t):
    q, v, d, j = t
    if v != ' ':
        if v != '.' and v != '-': return q.upper()
    elif d != ' ':
        if d != '.' and d != '-': return q.upper()
    elif j != ' ':
        if j != '.' and j != '-': return q.upper()
    return q.lower()

class IgBLASTRecord:
    class RegionInfo:
        def __init__(self, start, end, length, match_count=None, mismatch_count=None, gap_count=None, percent_identity=None):
            self.start            = start
            self.end              = end
            self.length           = length
            self.match_count      = match_count
            self.mismatch_count   = mismatch_count
            self.gap_count        = gap_count
            self.percent_identity = percent_identity
        def __repr__(self):
            return '[%d, %d] %f%%' % (self.start, self.end, 100 * self.percent_identity)
    class AlignmentScore:
        def __init__(self, bit_score, e_value):
            self.bit_score = bit_score
            self.e_value = e_value
        def __repr__(self):
            return '%s: bit score %f, E value %g' % (self.target_name, self.bit_score, self.e_value)
    class AlignmentLine:
        def __init__(self, segment_type, name, start, end, line):
            self.segment_type = segment_type
            self.name = name
            self.start = start
            self.end   = end
            self.line  = line
        def __repr__(self):
            return '%s %s %d %s %d' % (self.segment_type, self.name, self.start, self.line, self.end)

    def __init__(self, block):
        self.regions = {}
        self.regions_order = []
        self.alignment_regions = {}
        self.alignment_regions_definition = None
        self.is_vdj = None
        self.is_vj  = None
        self.strand = '+'   # the default
        self.v_frame_shift = None
        self.alignment_lines = []

        subblocks = block.split('\n\n\n')
        if subblocks[1] == '***** No hits found *****':
            self._parse_query_length_no_hits(subblocks[0])
        else:
            query_length_sig_align, domain_classification, rearrangement_junction_align_summary, alignments, statistics \
                    = subblocks
            self._parse_query_length_sig_align(query_length_sig_align)
            try:
                self._parse_domain_classification(domain_classification)
                self._parse_rearrangement_junction_align_summary(rearrangement_junction_align_summary)
                self._parse_alignments(alignments)
            except:
                print('error parsing', self.query_name, file=sys.stderr)
                raise
    def __nonzero__(self):
        return len(self.significant_alignments) > 0
    def _parse_query_length_sig_align(self, block):
        # break block into subblocks
        query_name, length_and_header, sig_alignments = block.split('\n\n')

        # extract the query name
        if not query_name.startswith('Query= '):
            raise ValueError(query_name)
        self.query_name = ' '.join(query_name[7:].split('\n')) # join-split long query names

        # get the length and validate the header
        length, header, title = length_and_header.split('\n')

        if not length.startswith('Length='):
            raise ValueError(length)
        self.query_length = int(length[7:])

        if header.split() != ['Score', 'E'] or title.split() != ['Sequences', 'producing', 'significant', 'alignments:', '(Bits)', 'Value']:
            raise ValueError(header + title)

        # parse the list of alignments and their score
        self._parse_sig_alignments(sig_alignments)
    def _parse_query_length_no_hits(self, block):
        # break block into subblocks
        query_name, length = block.split('\n\n')

        # extract the query name
        if not query_name.startswith('Query= '):
            raise ValueError(query_name)
        self.query_name = ''.join(query_name[7:].split('\n')) # join-split long query names

        if not length.startswith('Length='):
            raise ValueError(length)
        self.query_length = int(length[7:])

        self.significant_alignments = {}
    def _parse_sig_alignments(self, block):
        self.significant_alignments = {}
        for alignment in block.split('\n'):
            target_name, target_bit_score, target_e_value = alignment.split()
            self.significant_alignments[target_name] = self.AlignmentScore(float(target_bit_score), float(target_e_value))
    def _parse_domain_classification(self, block):
        if not block.startswith('Domain classification requested: '):
            raise ValueError(block)
        self.domain_classification = block[33:]
    def _parse_rearrangement_junction_align_summary(self, blocks):
        # process each summary table
        for block in blocks.split('\n\n'):
            # split block into header and data
            if block.count('\n') == 0:
                header, data = block, None
            else:
                header, data = block.split('\n', 1)

            # process each header type
            if header == 'Note that your query represents the minus strand of a V gene and has been converted to the plus strand. The sequence positions refer to the converted sequence. ':
                self.strand = '-'
                # no data
            elif header == 'V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Top C gene match, Chain type, stop codon, V-J frame, Productive, Strand, V Frame shift).  Multiple equivalent top matches, if present, are separated by a comma.':
                if self.is_vdj is not None or self.is_vj is not None:   # make sure rearrangement type not already defined
                    raise ValueError
                self.is_vdj = True
                self.is_vj  = False
                self._parse_vdjc_rearrangement_summary(header, data, has_v_frame_shift=True)
            elif header == 'V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Top C gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches, if present, are separated by a comma.':
                if self.is_vdj is not None or self.is_vj is not None:   # make sure rearrangement type not already defined
                    raise ValueError
                self.is_vdj = True
                self.is_vj  = False
                self._parse_vdjc_rearrangement_summary(header, data, has_v_frame_shift=False)
            elif header == 'V-(D)-J rearrangement summary for query sequence (Top V gene match, Top J gene match, Top C gene match, Chain type, stop codon, V-J frame, Productive, Strand, V Frame shift).  Multiple equivalent top matches, if present, are separated by a comma.':
                if self.is_vdj is not None or self.is_vj is not None:   # make sure rearrangement type not already defined
                    raise ValueError
                self.is_vdj = False
                self.is_vj  = True
                self._parse_vjc_rearrangement_summary(header, data, has_v_frame_shift=True)
            elif header == 'V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches, if present, are separated by a comma.':
                if self.is_vdj is not None or self.is_vj is not None:   # make sure rearrangement type not already defined
                    raise ValueError
                self.is_vdj = True
                self.is_vj  = False
                self._parse_vdj_rearrangement_summary(header, data, has_v_frame_shift=False)
            elif header == 'V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand, V Frame shift).  Multiple equivalent top matches, if present, are separated by a comma.':
                if self.is_vdj is not None or self.is_vj is not None:   # make sure rearrangement type not already defined
                    raise ValueError
                self.is_vdj = True
                self.is_vj  = False
                self._parse_vdj_rearrangement_summary(header, data, has_v_frame_shift=True)
            elif header == 'V-(D)-J rearrangement summary for query sequence (Top V gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches, if present, are separated by a comma.':
                if self.is_vdj is not None or self.is_vj is not None:   # make sure rearrangement type not already defined
                    raise ValueError
                self.is_vdj = False
                self.is_vj  = True
                self._parse_vj_rearrangement_summary(header, data)
            elif header == 'V-(D)-J rearrangement summary for query sequence (Top V gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand, V Frame shift).  Multiple equivalent top matches, if present, are separated by a comma.':
                if self.is_vdj is not None or self.is_vj is not None:   # make sure rearrangement type not already defined
                    raise ValueError
                self.is_vdj = False
                self.is_vj  = True
                self._parse_vj_rearrangement_summary(header, data, has_v_frame_shift=True)
            elif header == 'V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself':
                if not self.is_vdj and self.is_vj:      # check for mismatch of rearrangement type
                    raise ValueError
                self._parse_vdj_junction_details(header, data)
            elif header == 'V-(D)-J junction details based on top germline gene matches (V end, V-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself':
                if self.is_vdj and not self.is_vj:      # check for mismatch of rearrangement type
                    raise ValueError
                self._parse_vj_junction_details(header, data)
            elif header == 'Sub-region sequence details (nucleotide sequence, translation, start, end)':
                self._parse_subregion_details(header, data)
            elif header == 'Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)':
                self._parse_alignment_summary(header, data)
            else:
                # unknown block type
                raise ValueError(header)

    def _parse_vdjc_rearrangement_summary(self, header, data, has_v_frame_shift):
        if data.count('\n') != 0:
            raise ValueError(data)

        if has_v_frame_shift:
            top_v_segment_matches, top_d_segment_matches, top_j_segment_matches, top_c_segment_matches, chain_type, stop_codon, v_j_frame, productive, stand, v_frame_shift = data.split('\t')
        else:
            top_v_segment_matches, top_d_segment_matches, top_j_segment_matches, top_c_segment_matches, chain_type, stop_codon, v_j_frame, productive, stand = data.split('\t')
        
        # turn the top hits into lists
        if top_v_segment_matches == 'N/A':
            self.top_v_segment_matches = [None]
        else:
            self.top_v_segment_matches = top_v_segment_matches.split(',')
        if top_d_segment_matches == 'N/A':
            self.top_d_segment_matches = [None]
        else:
            self.top_d_segment_matches = top_d_segment_matches.split(',')
        if top_j_segment_matches == 'N/A':
            self.top_j_segment_matches = [None]
        else:
            self.top_j_segment_matches = top_j_segment_matches.split(',')

        if top_c_segment_matches == 'N/A':
            self.top_c_segment_matches = [None]
        else:
            self.top_c_segment_matches = top_c_segment_matches.split(',')

        if chain_type == 'N/A':
            self.chain_type = None
        else:
            self.chain_type = chain_type

        self._parse_stop_codon(stop_codon)
        self._parse_v_j_frame(v_j_frame)
        self._parse_productive(productive)
        if has_v_frame_shift:
            self._parse_v_frame_shift(v_frame_shift)

        assert self.strand == stand
    
    def _parse_vjc_rearrangement_summary(self, header, data, has_v_frame_shift):
        if data.count('\n') != 0:
            raise ValueError(data)

        if has_v_frame_shift:
            top_v_segment_matches, top_j_segment_matches, top_c_segment_matches, chain_type, stop_codon, v_j_frame, productive, stand, v_frame_shift = data.split('\t')
        else:
            top_v_segment_matches, top_j_segment_matches, top_c_segment_matches, chain_type, stop_codon, v_j_frame, productive, stand = data.split('\t')
        
        # turn the top hits into lists
        if top_v_segment_matches == 'N/A':
            self.top_v_segment_matches = [None]
        else:
            self.top_v_segment_matches = top_v_segment_matches.split(',')

        self.top_d_segment_matches = []

        if top_j_segment_matches == 'N/A':
            self.top_j_segment_matches = [None]
        else:
            self.top_j_segment_matches = top_j_segment_matches.split(',')

        if top_c_segment_matches == 'N/A':
            self.top_c_segment_matches = [None]
        else:
            self.top_c_segment_matches = top_c_segment_matches.split(',')

        if chain_type == 'N/A':
            self.chain_type = None
        else:
            self.chain_type = chain_type

        self._parse_stop_codon(stop_codon)
        self._parse_v_j_frame(v_j_frame)
        self._parse_productive(productive)
        if has_v_frame_shift:
            self._parse_v_frame_shift(v_frame_shift)

        assert self.strand == stand

    def _parse_vdj_rearrangement_summary(self, header, data, has_v_frame_shift):
        if data.count('\n') != 0:
            raise ValueError(data)

        if has_v_frame_shift:
            top_v_segment_matches, top_d_segment_matches, top_j_segment_matches, chain_type, stop_codon, v_j_frame, productive, stand, v_frame_shift = data.split('\t')
        else:
            top_v_segment_matches, top_d_segment_matches, top_j_segment_matches, chain_type, stop_codon, v_j_frame, productive, stand = data.split('\t')
        
        # turn the top hits into lists
        if top_v_segment_matches == 'N/A':
            self.top_v_segment_matches = [None]
        else:
            self.top_v_segment_matches = top_v_segment_matches.split(',')
        if top_d_segment_matches == 'N/A':
            self.top_d_segment_matches = [None]
        else:
            self.top_d_segment_matches = top_d_segment_matches.split(',')
        if top_j_segment_matches == 'N/A':
            self.top_j_segment_matches = [None]
        else:
            self.top_j_segment_matches = top_j_segment_matches.split(',')

        if chain_type == 'N/A':
            self.chain_type = None
        else:
            self.chain_type = chain_type

        self._parse_stop_codon(stop_codon)
        self._parse_v_j_frame(v_j_frame)
        self._parse_productive(productive)
        if has_v_frame_shift:
            self._parse_v_frame_shift(v_frame_shift)

        assert self.strand == stand

    def _parse_vj_rearrangement_summary(self, header, data, has_v_frame_shift):
        if data.count('\n') != 0:
            raise ValueError(data)
        if has_v_frame_shift:
            top_v_segment_matches, top_j_segment_matches, chain_type, stop_codon, v_j_frame, productive, stand, v_frame_shift = data.split('\t')
        else:
            top_v_segment_matches, top_j_segment_matches, chain_type, stop_codon, v_j_frame, productive, stand = data.split('\t')

        
        # turn the top hits into lists
        if top_v_segment_matches == 'N/A':
            self.top_v_segment_matches = [None]
        else:
            self.top_v_segment_matches = top_v_segment_matches.split(',')

        self.top_d_segment_matches = []

        if top_j_segment_matches == 'N/A':
            self.top_j_segment_matches = [None]
        else:
            self.top_j_segment_matches = top_j_segment_matches.split(',')

        if chain_type == 'N/A':
            self.chain_type = None
        else:
            self.chain_type = chain_type

        self._parse_stop_codon(stop_codon)
        self._parse_v_j_frame(v_j_frame)
        self._parse_productive(productive)

        if has_v_frame_shift:
            self._parse_v_frame_shift(v_frame_shift)

        assert self.strand == stand

    def _parse_stop_codon(self, stop_codon):
        if stop_codon == 'Yes':
            self.stop_codon = True
        elif stop_codon == 'No':
            self.stop_codon = False
        elif stop_codon == 'N/A':
            self.stop_codon = None
        else:
            raise ValueError(stop_codon)

    def _parse_v_frame_shift(self, v_frame_shift):
        if v_frame_shift == 'Yes':
            self.v_frame_shift = True
        elif v_frame_shift == 'No':
            self.v_frame_shift = False
        elif v_frame_shift == 'N/A':
            self.v_frame_shift = None
        else:
            raise ValueError(v_frame_shift)

    def _parse_v_j_frame(self, v_j_frame):
        if v_j_frame == 'In-frame':
            self.v_j_in_frame = True
        elif v_j_frame == 'Out-of-frame':
            self.v_j_in_frame = False
        elif v_j_frame == 'N/A':
            self.v_j_in_frame = None
        else:
            raise ValueError(v_j_frame)
    def _parse_productive(self, productive):
        if productive == 'Yes':
            self.productive = True
        elif productive == 'No':
            self.productive = False
        elif productive == 'N/A':
            self.productive = None
        else:
            raise ValueError(productive)

    def _parse_vdj_junction_details(self, header, data):
        if data.count('\n') != 0:
            raise ValueError(data)
        v_end_seq, v_d_junction_seq, d_segment_seq, d_j_junction_seq, j_start_seq, _ = data.split('\t')

        assert _ == ''
        if v_end_seq   == 'N/A': v_end_seq   = None
        if j_start_seq == 'N/A': j_start_seq = None
        
        self.v_end_seq = v_end_seq

        if v_d_junction_seq == 'N/A':
            self.v_d_junction_seq, self.v_d_junction_overlap = None, None
        else:
            self.v_d_junction_seq, self.v_d_junction_overlap = self._parse_overlap_seq(v_d_junction_seq)

        self.d_segment_seq = d_segment_seq if d_segment_seq != 'N/A' else None

        if d_j_junction_seq == 'N/A':
            self.d_j_junction_seq, self.d_j_junction_overlap = None, None
        else:
            self.d_j_junction_seq, self.d_j_junction_overlap = self._parse_overlap_seq(d_j_junction_seq)

        self.j_start_seq = j_start_seq

    def _parse_vj_junction_details(self, header, data):
        if data.count('\n') != 0:
            raise ValueError(data)
        v_end_seq, v_j_junction_seq, j_start_seq, _ = data.split('\t')

        assert _ == ''
        if v_end_seq   == 'N/A': v_end_seq   = None
        if j_start_seq == 'N/A': j_start_seq = None
        
        self.v_end_seq = v_end_seq

        if v_j_junction_seq == 'N/A':
            self.v_j_junction_seq, self.v_j_junction_overlap = None, None
        else:
            self.v_j_junction_seq, self.v_j_junction_overlap = self._parse_overlap_seq(v_j_junction_seq)

        self.j_start_seq = j_start_seq

    def _parse_overlap_seq(self, s):
        if s[0] == '(' and s[-1] == ')':
            return s[1:-1], True
        else:
            return s, False
    def _parse_subregion_details(self, header, data):
        if data.count('\n') != 0:
            raise ValueError(data)
        region_label, region_seq_nt, region_seq_aa, region_start, region_end, _ = data.split('\t')

        assert _ == ''
        if region_label != 'CDR3':
            raise ValueError('only CDR3 special sub-regions are currently suppported. Found ' + region_label)

        region_start = int(region_start) - 1    # change positions to 0-based,
        region_end   = int(region_end)          # half open ranges

        self.cdr3_seq_nt = region_seq_nt
        self.cdr3_seq_aa = region_seq_aa
        self.cdr3_start  = region_start
        self.cdr3_end    = region_end

        self.regions['CDR3'] = self.RegionInfo(region_start, region_end, region_end - region_start)
    def _parse_alignment_summary(self, header, data):
        domain_class = self.domain_classification.upper()

        for region_def in data.split('\n'):
            if region_def.count('\t') == 7:
                region_label, region_start, region_end, region_length, match_count, mismatch_count, gap_count, percent_identity = region_def.split('\t')
            else:
                region_label, region_start, region_end, region_length, match_count, mismatch_count, gap_count, percent_identity, _, _ = region_def.split('\t')

            region_start     = int(region_start) - 1 if region_start != 'N/A' else None # change positions to 0-based,
            region_end       = int(region_end) if region_end != 'N/A' else None         # half open ranges
            region_length    = int(region_length) if region_length != 'N/A' else None
            match_count      = int(match_count) if match_count != 'N/A' else None
            mismatch_count   = int(mismatch_count) if mismatch_count != 'N/A' else None
            gap_count        = int(gap_count) if gap_count != 'N/A' else None
            percent_identity = float(percent_identity) if percent_identity != 'N/A' else None

            # check that region names match the given domain classification
            if region_label == 'Total': # except for Total
                assert region_start is None and region_end is None
                self.regions['total'] = self.RegionInfo(region_start, region_end, region_length, match_count, mismatch_count, gap_count, percent_identity)
            elif region_label.startswith('CDR3-') and region_label.endswith(' (germline)'):
                if region_label[5:-11] != domain_class:
                    raise ValueError('regions label does not match domain classification, ' + region_label)
                self.regions['CDR3-germline'] = self.RegionInfo(region_start, region_end, region_length, match_count, mismatch_count, gap_count, percent_identity)
            else:
                region_label, region_domain_class = region_label.split(' ')[0].split('-')
                if region_domain_class != domain_class:
                    raise ValueError('regions label does not match domain classification, ' + region_domain_class)
                if region_label in self.regions:
                    raise ValueError('region ' + region_label + ' defined twice')
                self.regions[region_label] = self.RegionInfo(region_start, region_end, region_length, match_count, mismatch_count, gap_count, percent_identity)
                self.regions_order.append(region_label)

        # CDR3 is added out of order, check for it here to added it to the end
        if 'CDR3' in  self.regions:
            self.regions_order.append('CDR3')

    def _parse_alignments(self, block):
        header, alignments = block.split('\n\n')
        if header != 'Alignments':
            raise ValueError(header)
        
        alignments = alignments.split('\n')
        # if region line is missing
        first_nonblank = alignments[0].lstrip() # first non-space is the query
        if first_nonblank.startswith('Query_') or first_nonblank.startswith('lcl|Query_'):
            regions_def_line = None
            query_line, alignment_lines = alignments[0], alignments[1:]
        elif len(first_nonblank) == 0:
            regions_def_line = None
            query_line, alignment_lines = alignments[1], alignments[2:]
        else:
            regions_def_line, query_line, alignment_lines = alignments[0], alignments[1], alignments[2:]
        # parse the query line to get the basic structure
        begin_query_sep, query_start_sep, start_seq_sep, seq_end_sep = list(_re_space_sep.finditer(query_line))
        pre_slice    = slice(None,                  begin_query_sep.end())
        target_slice = slice(begin_query_sep.end(), query_start_sep.end())
        start_slice  = slice(query_start_sep.end(), start_seq_sep.end())
        align_slice  = slice(start_seq_sep.end(),   seq_end_sep.end())
        end_slice    = slice(seq_end_sep.end(),     None)
        
        # parse the first line to get the rest of the structure
        segment_ident_sep, ident_target_sep = list(_re_double_space.finditer(alignment_lines[0][pre_slice]))
        segment_slice = slice(None,                    segment_ident_sep.end())
        ident_slice   = slice(segment_ident_sep.end(), ident_target_sep.end())

        query_number_label = query_line[target_slice].strip()
        if query_number_label.startswith('Query_'):
            if self.strand != '+':
                raise ValueError('mismatch between strand definitions')
        elif query_number_label.startswith('lcl|Query_'):
            if not query_number_label.endswith('_reversed'):
                raise ValueError('expected _reversedsuffix', query_number_label)
            if self.strand != '-':
                raise ValueError('mismatch between strand definitions')
        else:
            raise ValueError(query_number_label)
        query_start_pos = int(query_line[start_slice]) - 1  # convert to zero-based ranges
        query_end_pos   = int(query_line[end_slice])
        self.alignment_lines.append(self.AlignmentLine('Q', query_number_label, query_start_pos, query_end_pos, query_line[align_slice].strip()))

        # now that we have processed the query line, we can look at the regions line
        if regions_def_line is not None:
            self._process_regions_definition_line(regions_def_line[align_slice])
        else:
            self._process_regions_definition_line(None)
        
        # parse each line
        for a in alignment_lines:
            segment_type   = a[segment_slice].strip()
            ident_column   = a[ident_slice]
            target_name    = a[target_slice].strip()
            start_position = int(a[start_slice]) - 1    # int already ignores extra whitespace, convert to zero-based ranges
            alignment      = a[align_slice].strip()
            end_position   = int(a[end_slice])

            # break the precent identity line into its parts
            ident_percent, ident_numerator, ident_denominator = self._parse_ident(ident_column)

            assert segment_type in ['V', 'D', 'J', 'C'], 'unknown segment segment type: ' + segment_type
            assert not self.alignment_regions_definition or len(alignment) == len(self.alignment_regions_definition), 'aligment length does not match legnth of regions string'
            self.alignment_lines.append(self.AlignmentLine(segment_type, target_name, start_position, end_position, alignment))
    def _parse_ident(self, ident):
        percent, numerator_denominator = ident.strip().split()
        if not percent.endswith('%'):
            raise ValueError
        percent = float(percent[:-1])

        if numerator_denominator[0] != '(' or numerator_denominator[-1] != ')':
            raise ValueError(numerator_denominator)
        numerator, denominator = numerator_denominator[1:-1].split('/')
        numerator   = int(numerator)
        denominator = int(denominator)
        
        return percent, numerator_denominator, denominator

    def _guess_region_name(self, regions, index):
        current_region_label = regions[index].group('label')
        if index + 1 < len(regions):
            # get the next label
            next_region = regions[index + 1]
            next_region_label = next_region.group('label')
            # assume that the current region is the next region's previous region
            prev_region_label = region_labels[region_labels.index(next_region_label) - 1]
            assert prev_region_label.startswith(current_region_label)
            return prev_region_label
        elif index > 0:
            # get the previous label
            prev_region = regions[index - 1]
            prev_region_label = prev_region.group('label')
            # assume that the current region is the previous region's next region
            next_region_label = region_labels[region_labels.index(prev_region_label) + 1]
            if next_region_label.startswith(current_region_label):
                return next_region_label
            else:
                next_region_label = region_labels[region_labels.index(prev_region_label) + 2]
                assert next_region_label.startswith(current_region_label), 'next: ' + next_region_label + ' next region: ' + current_region_label
                return next_region_label
        elif len(self.regions) == 2:    # if there is only one region (not countint total) in the alignment summary table, use full label from there
            for r in self.regions:
                if r.startswith(current_region_label):
                    return r
                else:
                    assert r == 'total'
            return None
        else:
            return None
    def _process_regions_definition_line(self, regions_def_line):
        self.alignment_regions_definition = regions_def_line
        if self.alignment_regions_definition is not None:
            # get list of regions from regions def. line
            regions = list(_re_region_ranges.finditer(self.alignment_regions_definition))
            for r_idx in range(len(regions)):
                r = regions[r_idx]
                label = r.group('label')
                if label in ['F', 'FR', 'C', 'CD', 'CDR']:
                    label = self._guess_region_name(regions, r_idx)
                assert label not in ['F', 'FR', 'C', 'CD', 'CDR'], '%s: truncated label %s' % (self.query_name, label) # check for truncated labels
                region_slice = slice(r.start(), r.end())
                self.alignment_regions[label] = region_slice
                region_string = self.alignment_regions_definition[region_slice]
                assert region_string.startswith('<') and region_string.endswith('>')
    def _get_mutation_capped_query(self, v_segment=None, d_segment=None, j_segment=None):
        capped_alignment = self.alignment_query.line

        if v_segment is None:
            v_segment = self.get_v_segment()
        if d_segment is None:
            d_segment = self.get_d_segment()
        if j_segment is None:
            j_segment = self.get_j_segment()

        if v_segment is None:
            v_alignment = ' ' * len(capped_alignment)
        else:
            v_alignment = self.alignment_v_segments[v_segment].line

        if d_segment is None:
            d_alignment = ' ' * len(capped_alignment)
        else:
            d_alignment = self.alignment_d_segments[d_segment].line

        if j_segment is None:
            j_alignment = ' ' * len(capped_alignment)
        else:
            j_alignment = self.alignment_j_segments[j_segment].line

        assert len(capped_alignment) == len(v_alignment)
        assert len(capped_alignment) == len(d_alignment)
        assert len(capped_alignment) == len(j_alignment)

        capped_alignment = ''.join(map(upper_mismatch,
                                   zip(capped_alignment, v_alignment, d_alignment, j_alignment)))

        return capped_alignment
    def __bool__(self):
        return len(self.significant_alignments) > 0
    def get_v_segment(self):
        return self.top_v_segment_matches[0]
    def get_d_segment(self):
        if self.is_vdj:
            return self.top_d_segment_matches[0]
        else:
            assert len(self.top_d_segment_matches) == 0
            return None
    def get_j_segment(self):
        return self.top_j_segment_matches[0]
    def get_region_nt_sequence(self, region):
        if region in self.alignment_regions:
            return self.alignment_query.line[self.alignment_regions[region]]
        else:
            return None
    def get_cdr3_nt(self):
        if 'CDR3' in self.alignment_regions:
            cdr3 = self.alignment_query.line[self.alignment_regions['CDR3']]
            assert cdr3.replace('-', '') == self.cdr3_seq_nt
            return self.cdr3_seq_nt
        else:
            return None

class IgBLASTParser:
    def __init__(self, handle):
        self.handle = handle

        self.header = None
        self.footer = None

        self._read_header()

    def _read_header(self):
        handle_readline = self.handle.readline

        # gather all the header lines until we hit 3 blank new lines
        header = ''
        line = handle_readline()
        new_line_count = 0
        while line and not line.startswith('Query= ') and not line.startswith('Total queries ='):
            header += line
            if line == '\n':
                new_line_count += 1
            else:
                new_line_count = 0
            if new_line_count == 3:
                break
            line = handle_readline()
        self.header = header    # store the header
        if line.startswith('Total queries ='):  # empty file
            self._read_footer(line)
    def _read_footer(self, first_line):
        handle_readline = self.handle.readline

        footer = first_line
        line = handle_readline()    
        while line:
            footer += line
            line = handle_readline()
        self.footer = footer
    def __iter__(self):
        return self
    def __next__(self):
        return self.next()
    def next(self):
        handle_readline = self.handle.readline

        # if we've already read the footer, the file is empty and we should stop
        if self.footer is not None:
            raise StopIteration
        
        # skip blank lines
        line = handle_readline()    
        while line == '\n':
            line = handle_readline()    

        # if this is the footer
        if line.startswith('Total queries = '):
            self._read_footer(line)
            raise StopIteration

        # start storing the record
        block = line
        line = handle_readline()
        # look for last line
        while not line.startswith('Effective search space used: '):
            block += line
            line = handle_readline()
        block += line
        return IgBLASTRecord(block)
