#cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False

from __future__ import division, print_function, absolute_import
from pywfa cimport WFA_wrap as wfa
from dataclasses import dataclass


__all__ = ["WavefrontAligner", "clip_cigartuples", "cigartuples_to_str", "elide_mismatches_from_cigar"]


cdef int[89] codes
codes[:] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,0,0,0,0,9,0,2,0,0,0,5,1,
 0,0,0,0,3,0,6,0,0,4,0,0,0,0,8]


@dataclass
class AlignmentResult:
    """Holds the result of an alignment
    """
    pattern_length: int
    text_length: int
    pattern_start: int
    pattern_end: ind
    text_start: int
    text_end: int
    cigartuples: object
    score: int
    pattern: str
    text: str
    status: int
    def __init__(self, pl, tl, ps, pe, ts, te, ct, s, p, t, status):
        self.pattern_length = pl
        self.text_length = tl
        self.pattern_start = ps
        self.pattern_end = pe
        self.text_start = ts
        self.text_end = te
        self.cigartuples = ct
        self.score = s
        self.pattern = p
        self.text = t
        self.status = status

    def __repr__(self):
        return str(self.__dict__)

    @property
    def aligned_pattern(self):
        """Returns the pattern sequence aligned by the cigar
        Returns
        -------
        str
            Aligned pattern sequence
        Notes
        -----
        This will return `None` if `suppress_sequences` was True when this
        object was created
        """
        if self.pattern:
            return self._get_aligned_sequence(self.pattern,
                                              self.cigartuples,
                                              self.pattern_start, self.pattern_end,
                                              "D")

    @property
    def aligned_text(self):
        """Returns the text sequence aligned by the cigar
        Returns
        -------
        str
            Aligned text sequence
        Notes
        -----
        This will return `None` if `suppress_sequences` was True when this
        object was created
        """
        if self.text:
            return self._get_aligned_sequence(self.text,
                                              self.cigartuples,
                                              self.text_start,
                                              self.text_end,
                                              "I")

    def _get_aligned_sequence(self, sequence, tuple_cigar, begin, end,
                              gap_type):
        aligned_sequence = []
        seq = sequence[begin:end]
        index = 0
        for length, mid in tuple_cigar:
            if mid == gap_type:
                aligned_sequence += ['-' * length]
            else:
                aligned_sequence += [seq[index:index + length]]
                index += length
        aligned_sequence += [seq[index:end - begin]]
        return "".join(aligned_sequence)


cpdef clip_cigartuples(object align_result, int min_aligned_bases_left=5, int min_aligned_bases_right=5):
    """Returns cigartuples with blocks of aligned bases < threshold removed from each end
    Returns
    -------
    list
        cigartuples
    """
    ct = align_result.cigartuples
    if not ct:
        return align_result

    cdef int i, j, start_l, end_l, left_clip, right_clip
    i = 0
    text_start = 0
    pattern_start = 0
    for i in range(len(ct)):
        if ct[i][0] == 0:
            if ct[i][1] >= min_aligned_bases_left:
                break
            else:
                text_start += ct[i][1]
                pattern_start += ct[i][1]
        elif ct[i][0] == 2:  # deletion
            pattern_start += ct[i][1]
        elif ct[i][0] == 8 :  # mismatch
            text_start += ct[i][1]
            pattern_start += ct[i][1]
        elif ct[i][0] == 1: # insertion
            text_start += ct[i][1]

    text_end = align_result.text_length
    pattern_end = align_result.pattern_length
    j = len(ct) - 1
    for j in range(len(ct) -1, -1, -1):
        if ct[j][0] == 0 :
            if ct[j][1] >= min_aligned_bases_right:
                break
            else:
                text_end -= ct[j][1]
                pattern_end -= ct[j][1]
        elif ct[j][0] == 2:
            pattern_end -= ct[j][1]
        elif ct[j][0] == 8:
            pattern_end -= ct[j][1]
            text_end -= ct[j][1]
        elif ct[j][0] == 1:
            text_end -= ct[j][1]

    modified = []
    if align_result.text_start + text_start > 0:
        modified.append((4, align_result.text_start + text_start))
    modified += ct[i:j+1]
    if align_result.text_length - text_end > 0:
        modified.append((4, align_result.text_length - text_end))
    align_result.cigartuples = modified

    align_result.text_start = text_start
    align_result.text_end = text_end

    align_result.pattern_start = pattern_start
    align_result.pattern_end = pattern_end

    return align_result


cpdef elide_mismatches_from_cigar(cigartuples):
    """Returns cigartuples with mismatched 'X' merged into aligned blocks 'M'
    Returns
    -------
    list
        cigartuples
    """
    if not cigartuples:
        return []
    modified = []
    cdef int l
    block = 0
    for opp, l in cigartuples:
        if opp != 8 and opp != 0:
            if block:
                modified.append((0, block))
                block = 0
            modified.append((opp, l))
        else:
            block += l
    if block:
        modified.append((0, block))
    return modified


cpdef cigartuples_to_str(cigartuples):
    """Returns string format of cigartuples
    Returns
    -------
    str
        cigar in string format
    """
    if not cigartuples:
        return ""
    cdef int l
    str_codes = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "B"]
    cigarstring = ""
    for opp, l in cigartuples:
        cigarstring += f"{l}{str_codes[opp]}"
    return cigarstring


cdef class WavefrontAligner:
    """Wrapper class for WFA2-lib. If a pattern is supplied, it will be cached for re-use
    """
    def __init__(self,
                 pattern=None,
                 distance="affine",
                 int match=0,
                 int mismatch=4,
                 int gap_opening=6,
                 int gap_extension=2,
                 int gap_opening2=24,
                 int gap_extension2=1,
                 scope="full",
                 span="ends-free",
                 int pattern_begin_free=0,
                 int pattern_end_free=0,
                 int text_begin_free=0,
                 int text_end_free=0,
                 heuristic=None,
                 # min_k=-10,
                 # max_k=10,
                 int min_wavefront_length=10,
                 int max_distance_threshold=50,
                 int steps_between_cutoffs=1,
                 int xdrop=20,
                 ):

        self.pattern_len = 0
        self.text_len = 0
        if pattern:
            self._pattern = pattern

        # could get a malloc version working
        # self.attributes = <wfa.wavefront_aligner_attr_t* > malloc(sizeof(wfa.wavefront_aligner_attr_default))
        self.attributes = &wfa.wavefront_aligner_attr_default

        if distance == "affine":
            self.attributes.distance_metric = wfa.gap_affine
            self.attributes.affine_penalties.match = match
            self.match_score = match
            self.attributes.affine_penalties.mismatch = mismatch
            self.attributes.affine_penalties.gap_opening = gap_opening
            self.attributes.affine_penalties.gap_extension = gap_extension
        elif distance == "affine2p":
            self.attributes.distance_metric = wfa.gap_affine_2p
            self.attributes.affine2p_penalties.match = match
            self.match_score = match
            self.attributes.affine2p_penalties.mismatch = mismatch
            self.attributes.affine2p_penalties.gap_opening1 = gap_opening
            self.attributes.affine2p_penalties.gap_extension1 = gap_extension
            self.attributes.affine2p_penalties.gap_opening2 = gap_opening2
            self.attributes.affine2p_penalties.gap_extension2 = gap_extension2
        else:
            print(NotImplementedError(f'{distance} distance not implemented'))
            # raise NotImplementedError(f'{distance} distance not implemented')
        if scope == "full":
            self.attributes.alignment_scope = wfa.compute_alignment
        elif scope == "score":
            self.attributes.alignment_scope = wfa.compute_score
            self.score_only = True
        else:
            print(ValueError(f'{scope} scope not understood'))
            # raise ValueError(f'{scope} scope not understood')

        self.attributes.alignment_form.pattern_begin_free = pattern_begin_free
        self.attributes.alignment_form.pattern_end_free = pattern_end_free
        self.attributes.alignment_form.text_begin_free = text_begin_free
        self.attributes.alignment_form.text_end_free = text_end_free
        if span == "ends-free":
            self.attributes.alignment_form.span = wfa.alignment_endsfree

        elif span == "end-to-end":
            self.attributes.alignment_form.span = wfa.alignment_end2end
        else:
            print(NotImplementedError(f'{span} span not implemented'))
            # raise NotImplementedError(f'{span} span not implemented')

        if heuristic is None:
            self.attributes.heuristic.strategy = wfa.wf_heuristic_none
        elif heuristic == "adaptive":
            self.attributes.heuristic.strategy = wfa.wf_heuristic_wfadaptive
            self.attributes.heuristic.min_wavefront_length = min_wavefront_length
            self.attributes.heuristic.max_distance_threshold = max_distance_threshold
            self.attributes.heuristic.steps_between_cutoffs = steps_between_cutoffs
        elif heuristic == "X-drop":
            self.attributes.heuristic.strategy = wfa.wf_heuristic_xdrop
            self.attributes.heuristic.xdrop = xdrop
            self.attributes.heuristic.steps_between_cutoffs = steps_between_cutoffs
        else:
            print(NotImplementedError(f'{heuristic} heuristic not implemented'))
            # raise NotImplementedError(f'{heuristic} heuristic not implemented')

        self.wf_aligner = wfa.wavefront_aligner_new(self.attributes)

    def wavefront_align(self, text, pattern=None):
        """The main alignment function. Returns alignment score
        Returns
        -------
        int
            alignment score
        """
        cdef bytes p
        if pattern is not None:
            p = pattern.encode('ascii')
            self._pattern = pattern
        else:
            p = self._pattern.encode('ascii')
        cdef bytes t = text.encode('ascii')
        self.text_len = len(t)
        self.pattern_len = len(p)
        wfa.wavefront_align(self.wf_aligner, p, <size_t>len(p), t, <size_t>len(text))
        return self.wf_aligner.cigar.score

    @property
    def status(self):
        return self.wf_aligner.align_status.status

    @property
    def score(self):
        return self.wf_aligner.cigar.score

    @property
    def cigarstring(self):
        cdef wfa.cigar_t* cigar
        cdef char last_op
        cdef int last_op_length, i, length

        cigar = &self.wf_aligner.cigar

        # Check null CIGAR
        if cigar.begin_offset >= cigar.end_offset:
            return ""
        # Print operations
        last_op = cigar.operations[cigar.begin_offset]
        last_op_length = 1
        length = 1
        result = ""
        for i in range(cigar.begin_offset+1, cigar.end_offset):
            if cigar.operations[i] == last_op:
                last_op_length += 1
                length += 1
            else:
                result += f"{length}{chr(last_op)}"
                length = 1
            last_op = cigar.operations[i]
            last_op_length = 1
        result += f"{length}{chr(last_op)}"
        return result

    @property
    def cigartuples(self):
        cdef wfa.cigar_t* cigar
        cdef char last_op
        cdef int last_op_length, i, length

        cigar = &self.wf_aligner.cigar

        # Check null CIGAR
        if cigar.begin_offset >= cigar.end_offset:
            return []
        # Print operations
        last_op = cigar.operations[cigar.begin_offset]
        last_op_length = 1
        length = 1
        # codes = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8, "B": 9}
        result = []
        for i in range(cigar.begin_offset+1, cigar.end_offset):
            if cigar.operations[i] == last_op:
                last_op_length += 1
                length += 1
            else:
                result.append((codes[<int>last_op], length))
                length = 1
            last_op = cigar.operations[i]
            last_op_length = 1
        result.append((codes[<int>last_op], length))
        return result

    @property
    def locations(self):
        if self.score_only:
            return [0, 0, 0, 0]
        cigartuples = self.cigartuples
        if not cigartuples or self.text_len == 0 or self.pattern_len == 0:
            return [0, 0, 0, 0]
        cdef int pattern_start, pattern_end, text_start, text_end, i, j

        ct = cigartuples
        text_start = 0
        pattern_start = 0
        for i in range(len(cigartuples)):
            if ct[i][0] == 0:
                if ct[i][1] >= 1:
                    break
                else:
                    text_start += ct[i][1]
                    pattern_start += ct[i][1]
            elif ct[i][0] == 2:  # deletion
                pattern_start += ct[i][1]
            elif ct[i][0] == 8 :  # mismatch
                text_start += ct[i][1]
                pattern_start += ct[i][1]
            elif ct[i][0] == 1: # insertion
                text_start += ct[i][1]

        text_end = self.text_len
        pattern_end = self.pattern_len
        j = len(ct) - 1
        for j in range(len(ct) -1, -1, -1):
            if ct[j][0] == 0 :
                if ct[j][1] >= 1:
                    break
                else:
                    text_end -= ct[j][1]
                    pattern_end -= ct[j][1]
            elif ct[j][0] == 2:
                pattern_end -= ct[j][1]
            elif ct[j][0] == 8:
                pattern_end -= ct[j][1]
                text_end -= ct[j][1]
            elif ct[j][0] == 1:
                text_end -= ct[j][1]

        return pattern_start, pattern_end, text_start, text_end

    def __call__(self, text, pattern=None, clip_cigar=False, min_aligned_bases_left=1, min_aligned_bases_right=1, elide_mismatches=False,
                 supress_sequences=False):

        if pattern is None:
            p = self._pattern
            if not p:
                raise ValueError("pattern is None")
            lp = len(self._pattern)
            score = self.wavefront_align(text)
        else:
            lp = len(pattern)
            p = pattern
            score = self.wavefront_align(text, pattern)

        ct = self.cigartuples
        locs = self.locations
        status = self.status
        if supress_sequences:
            res = AlignmentResult(lp, len(text), locs[0], locs[1], locs[2], locs[3], ct, score, "", "", status)
        else:
            res = AlignmentResult(lp, len(text), locs[0], locs[1], locs[2], locs[3], ct, score, p, text, status)
        if not self.score_only:
            if clip_cigar:
                res = clip_cigartuples(res, min_aligned_bases_left, min_aligned_bases_right)
            if elide_mismatches:
                res.cigartuples = elide_mismatches_from_cigar(res.cigartuples)
        return res

    def __dealloc__(self):
        wfa.wavefront_aligner_delete(self.wf_aligner)
