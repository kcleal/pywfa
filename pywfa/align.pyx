#cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False

from pywfa cimport WFA_wrap as wfa
from dataclasses import dataclass
from libc.stdio cimport stdout, FILE, fopen, fclose, fputs
from libc.limits cimport INT_MAX

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
        data = ['score',
                'pattern_start',
                'pattern_end',
                'text_start',
                'text_end',
                'cigartuples',
                'pattern',
                'text']
        d = self.__dict__
        s = f""
        for k in data:
            s += f"    {k}: {d[k]}\n"
        return s

    def __str__(self):
        score = "Score: %d" % self.score
        if self.pattern and self.cigartuples:
            t = self.aligned_text
            p = self.aligned_pattern
            align_len = len(t)
            if len(t) > 30:
                t = t[:30] + "..."
                p = p[:30] + "..."
            c = self.cigarstring
            if len(c) > 30:
                c = c[:30]
            length = "Length: %d" % len(t)
            return "\n".join([p, t, c, score, length])
        return score

    @property
    def aligned_pattern(self):
        """Returns the pattern sequence aligned by the cigar

        Notes
        -----
        This will return `None` if `suppress_sequences` was True when this
        object was created

        :return: Aligned pattern sequence
        :rtype: str
        """
        if self.pattern:
            return self._get_aligned_sequence(self.pattern,
                                              self.cigartuples,
                                              self.pattern_start, self.pattern_end,
                                              "D")

    @property
    def aligned_text(self):
        """Returns the text sequence aligned by the cigar

        Notes
        -----
        This will return `None` if `suppress_sequences` was True when this
        object was created

        :return: Aligned text sequence
        :rtype: str
        """
        if self.text:
            return self._get_aligned_sequence(self.text,
                                              self.cigartuples,
                                              self.text_start,
                                              self.text_end,
                                              "I")

    @property
    def cigarstring(self):
        """Returns the cigar in str format

        :return: cigar sequence
        :rtype: str
        """
        return cigartuples_to_str(self.cigartuples)
    
    @property
    def pretty(self):
        """Returns a str of the alignment result in 'pretty' format

        :return: pretty string
        :rtype: str
        """
        s = f"{self.cigarstring}      ALIGNMENT\n"
        s += f"{cigartuples_to_str([i for i in self.cigartuples if i[0] !=0 and i[0] !=[8]])}      ALIGNMENT.COMPACT\n"
        p = "      PATTERN    "
        g = "                 "
        t = "      TEXT       "
        pat = self.pattern
        pi = 0
        txt = self.text
        ti = 0
        for opp, l in self.cigartuples:
            if opp in (1, 4, 5):
                t += txt[ti: ti+l]
                ti += l
                p += "-"*l
                g += " "*l
            elif opp in (0, 7):
                t += txt[ti: ti+l]
                ti += l
                p += pat[pi: pi+l]
                pi += l
                g += "|"*l
            elif opp == 2:
                t += "-"*l
                p += pat[pi: pi+l]
                pi += l
                g += " "*l
            elif opp == 8:
                t += txt[ti: ti+l]
                ti += l
                p += pat[pi: pi+l]
                pi += l
                g += "*"*l
            else:
                raise ValueError(f"Cigar operation not available for pretty print - {opp}")
                
        s += p + "\n" + g + "\n" + t + "\n"
        return s


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

    :param align_result: AlignmentResult dataclass
    :type align_result: dataclass
    :param min_aligned_bases_left: Minimum allowed length of matched bases at left flank
    :type min_aligned_bases_left: int
    :param min_aligned_bases_right: Minimum allowed length of matched bases at right flank
    :type min_aligned_bases_right: int
    :return: AlignmentResult dataclass with trimmed cigartuples
    :rtype: dataclass
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
        modified.append((4, text_start))
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
    
    :param cigartuples: list of cigartuples
    :type cigartuples: list
    :return: cigartuples
    :rtype: list
    """
    if not cigartuples:
        return []
    modified = []
    cdef int l
    cdef int opp
    cdef int block = 0
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

    :param cigartuples: cigartuples
    :type cigartuples: list 
    :return: cigar in string format
    :rtype: str
    """
    if not cigartuples:
        return ""
    cdef int l
    str_codes = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "B"]
    cigarstring = ""
    for opp, l in cigartuples:
        cigarstring += f"{l}{str_codes[opp]}"
    return cigarstring

ctypedef struct wildcard_fun_args:
    char* pattern
    char* query
    char wildcard

cdef int wildcard_match_fun(int pattern_pos, int query_pos, void* argsptr) noexcept nogil:
    cdef const wildcard_fun_args* args = <wildcard_fun_args*> argsptr
    return args[0].pattern[pattern_pos] == args[0].wildcard or args[0].query[query_pos] == args[0].wildcard or args[0].pattern[pattern_pos] == args[0].query[query_pos]

cdef class WavefrontAligner:
    """Wrapper class for WFA2-lib. If a pattern is supplied, it will be cached for re-use
    """
    def __init__(self,
                 pattern=None,
                 distance="affine",
                 memory_mode="high",
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
                 wildcard=None,
                 int max_steps=0
                 ):
        self.pattern_len = 0
        self.text_len = 0
        if pattern:
            self._pattern = pattern.upper()
            self._bpattern = self._pattern.encode("ascii")
            self.pattern_len = len(self._bpattern)

        # could get a malloc version working
        # self.attributes = <wfa.wavefront_aligner_attr_t* > malloc(sizeof(wfa.wavefront_aligner_attr_default))
        cdef wfa.wavefront_aligner_attr_t attributes = wfa.wavefront_aligner_attr_default
        self.wildcard = wildcard

        if distance == "indel":
            attributes.distance_metric = wfa.indel
        elif distance == "levenshtein":
            attributes.distance_metric = wfa.edit
        elif distance == "linear":
            attributes.distance_metric = wfa.gap_linear
            attributes.linear_penalties.match = match
            attributes.linear_penalties.mismatch = mismatch
            attributes.linear_penalties.indel = gap_extension
        elif distance == "affine":
            attributes.distance_metric = wfa.gap_affine
            attributes.affine_penalties.match = match
            attributes.affine_penalties.mismatch = mismatch
            attributes.affine_penalties.gap_opening = gap_opening
            attributes.affine_penalties.gap_extension = gap_extension
        elif distance == "affine2p":
            attributes.distance_metric = wfa.gap_affine_2p
            attributes.affine2p_penalties.match = match
            attributes.affine2p_penalties.mismatch = mismatch
            attributes.affine2p_penalties.gap_opening1 = gap_opening
            attributes.affine2p_penalties.gap_extension1 = gap_extension
            attributes.affine2p_penalties.gap_opening2 = gap_opening2
            attributes.affine2p_penalties.gap_extension2 = gap_extension2
        else:
            raise NotImplementedError(f'{distance} distance not implemented')
        if scope == "full":
            attributes.alignment_scope = wfa.compute_alignment
        elif scope == "score":
            attributes.alignment_scope = wfa.compute_score
        else:
            raise ValueError(f'{scope} scope not understood')

        if memory_mode == 'high':
            attributes.memory_mode = wfa.wavefront_memory_high
        elif memory_mode == 'medium':
            attributes.memory_mode = wfa.wavefront_memory_med
        elif memory_mode == 'low':
            attributes.memory_mode = wfa.wavefront_memory_low
        elif memory_mode == 'biwfa':
            attributes.memory_mode = wfa.wavefront_memory_ultralow
        else:
            raise ValueError('memory_mode must be one of \'high\', \'medium\', \'low\', \'biwfa\'')

        attributes.alignment_form.pattern_begin_free = pattern_begin_free
        attributes.alignment_form.pattern_end_free = pattern_end_free
        attributes.alignment_form.text_begin_free = text_begin_free
        attributes.alignment_form.text_end_free = text_end_free
        if span == "ends-free":
            attributes.alignment_form.span = wfa.alignment_endsfree
        elif span == "end-to-end":
            attributes.alignment_form.span = wfa.alignment_end2end
        else:
            raise NotImplementedError(f'{span} span not implemented')

        if heuristic is None:
            attributes.heuristic.strategy = wfa.wf_heuristic_none
        elif heuristic == "adaptive":
            attributes.heuristic.strategy = wfa.wf_heuristic_wfadaptive
            attributes.heuristic.min_wavefront_length = min_wavefront_length
            attributes.heuristic.max_distance_threshold = max_distance_threshold
            attributes.heuristic.steps_between_cutoffs = steps_between_cutoffs
        elif heuristic == "X-drop":
            attributes.heuristic.strategy = wfa.wf_heuristic_xdrop
            attributes.heuristic.xdrop = xdrop
            attributes.heuristic.steps_between_cutoffs = steps_between_cutoffs
        else:
            raise NotImplementedError(f'{heuristic} heuristic not implemented')

        if max_steps <= 0:
            max_steps = INT_MAX
        attributes.system.max_alignment_steps = max_steps

        self.wf_aligner = wfa.wavefront_aligner_new(&attributes)

    def wavefront_align(self, text, pattern=None):
        """Perform wavefront alignment.

        :param text: The text sequence to align in uppercase
        :type text: str
        :param pattern: The pattern sequence to align in uppercase
        :type pattern: str
        :return: Alignment score
        :rtype: int
        """
        if pattern is not None:
            self._pattern = pattern.upper()
            self._bpattern = self._pattern.encode("ascii")
            self.pattern_len = len(self._bpattern)
        cdef bytes t = text.upper().encode('ascii')
        self._text = text
        self.text_len = len(t)
        if not self._wildcard:
            wfa.wavefront_align(self.wf_aligner, self._bpattern, <size_t>len(self._bpattern), t, <size_t>len(text))
        else:
            args = wildcard_fun_args(self._bpattern, t, self._bwildcard)
            wfa.wavefront_align_lambda(self.wf_aligner, wildcard_match_fun, &args, <size_t>len(self._bpattern), <size_t>len(text))
        return self.wf_aligner.cigar.score

    def cigar_print_pretty(self, file_name=None):
        cdef bytes t = self._text.encode('ascii')
        cdef bytes fname_bytes
        cdef char* fname
        cdef FILE * outfile
        if file_name:
            fname_bytes = file_name.encode("UTF-8")
            fname = fname_bytes
            outfile = fopen(fname, "w")
        else:
            outfile = stdout
        wfa.cigar_print_pretty(outfile, self.wf_aligner.cigar, self._bpattern, <size_t>len(self._bpattern), t, <size_t>len(t))

        if file_name:
            fclose(outfile)

    @property
    def status(self):
        return self.wf_aligner.align_status.status

    @property
    def score(self):
        return self.wf_aligner.cigar.score

    @property
    def pattern_begin_free(self):
        return self.wf_aligner.alignment_form.pattern_begin_free

    @pattern_begin_free.setter
    def pattern_begin_free(self, int pattern_begin_free):
        self.wf_aligner.alignment_form.pattern_begin_free = pattern_begin_free

    @property
    def pattern_end_free(self):
        return self.wf_aligner.alignment_form.pattern_end_free

    @pattern_end_free.setter
    def pattern_end_free(self, int pattern_end_free):
        self.wf_aligner.alignment_form.pattern_end_free = pattern_end_free

    @property
    def text_begin_free(self):
        return self.wf_aligner.alignment_form.text_begin_free

    @text_begin_free.setter
    def text_begin_free(self, int text_begin_free):
        self.wf_aligner.alignment_form.text_begin_free = text_begin_free

    @property
    def text_end_free(self):
        return self.wf_aligner.alignment_form.text_end_free

    @text_end_free.setter
    def text_end_free(self, int text_end_free):
        self.wf_aligner.alignment_form.text_end_free = text_end_free

    @property
    def scope(self):
        if self.wf_aligner.alignment_scope == wfa.compute_alignment:
            return "full"
        else:
            return "score"

    @scope.setter
    def scope(self, scope):
        if scope == "full":
            self.wf_aligner.alignment_scope = wfa.compute_alignment
        elif scope == "score":
            self.wf_aligner.alignment_scope = wfa.compute_score
        else:
            raise ValueError(f'{scope} scope not understood')

    @property
    def span(self):
        if self.wf_aligner.alignment_form.span == wfa.alignment_endsfree:
            return "ends-free"
        elif self.wf_aligner.alignment_form.span == wfa.alignment_end2end:
            return "end-to-end"

    @span.setter
    def span(self, span):
        if span == "ends-free":
            self.wf_aligner.alignment_form.span = wfa.alignment_endsfree

        elif span == "end-to-end":
            self.wf_aligner.alignment_form.span = wfa.alignment_end2end
        else:
            raise NotImplementedError(f'{span} span not implemented')

    @property
    def memory_mode(self):
        if self.wf_aligner.memory_mode == wfa.wavefront_memory_high:
            return "high"
        elif self.wf_aligner.memory_mode == wfa.wavefront_memory_med:
            return "medium"
        elif self.wf_aligner.memory_mode == wfa.wavefront_memory_low:
            return "low"
        elif self.wf_aligner.memory_mode == wfa.wavefront_memory_ultralow:
            return "biwfa"

    @memory_mode.setter
    def memory_mode(self, memory_mode):
        if memory_mode is "high":
            self.wf_aligner.memory_mode = wfa.wavefront_memory_high
        elif memory_mode == "med":
            self.wf_aligner.memory_mode = wfa.wavefront_memory_med
        elif memory_mode == "low":
            self.wf_aligner.memory_mode = wfa.wavefront_memory_low
        elif memory_mode == "biwfa":
            self.wf_aligner.memory_mode = wfa.wavefront_memory_ultralow
        else:
            raise NotImplementedError(f'{memory_mode} memory_mode not implemented')

    @property
    def heuristic(self):
        if self.wf_aligner.heuristic.strategy == wfa.wf_heuristic_none:
            return None
        elif self.wf_aligner.heuristic.strategy == wfa.wf_heuristic_wfadaptive:
            return "adaptive"
        elif self.wf_aligner.heuristic.strategy == wfa.wf_heuristic_xdrop:
            return "X-drop"

    @heuristic.setter
    def heuristic(self, heuristic):
        if heuristic is None:
            self.wf_aligner.heuristic.strategy = wfa.wf_heuristic_none
        elif heuristic == "adaptive":
            self.wf_aligner.heuristic.strategy = wfa.wf_heuristic_wfadaptive
        elif heuristic == "X-drop":
            self.wf_aligner.heuristic.strategy = wfa.wf_heuristic_xdrop
        else:
            raise NotImplementedError(f'{heuristic} heuristic not implemented')

    @property
    def min_wavefront_length(self):
        return self.wf_aligner.heuristic.min_wavefront_length

    @min_wavefront_length.setter
    def min_wavefront_length(self, int length):
        self.wf_aligner.heuristic.min_wavefront_length = length

    @property
    def max_distance_threshold(self):
        return self.wf_aligner.heuristic.max_distance_threshold

    @max_distance_threshold.setter
    def max_distance_threshold(self, int thresh):
        self.wf_aligner.heuristic.max_distance_threshold = thresh

    @property
    def steps_between_cutoffs(self):
        return self.wf_aligner.heuristic.steps_between_cutoffs

    @steps_between_cutoffs.setter
    def steps_between_cutoffs(self, int steps):
        self.wf_aligner.heuristic.steps_between_cutoffs = steps

    @property
    def xdrop(self):
        return self.wf_aligner.heuristic.xdrop

    @xdrop.setter
    def xdrop(self, int xdrop):
        self.wf_aligner.heuristic.xdrop = xdrop

    def _edit_penalties(self):
        if self.wf_aligner.penalties.distance_metric == wfa.indel:
            wfa.wavefront_penalties_set_indel(&self.wf_aligner.penalties)
        elif self.wf_aligner.penalties.distance_metric == wfa.edit:
            wfa.wavefront_penalties_set_edit(&self.wf_aligner.penalties)
        elif self.wf_aligner.penalties.distance_metric == wfa.gap_linear:
            wfa.wavefront_penalties_set_linear(&self.wf_aligner.penalties, &self.wf_aligner.penalties.linear_penalties)
        elif self.wf_aligner.penalties.distance_metric == wfa.gap_affine:
            wfa.wavefront_penalties_set_affine(&self.wf_aligner.penalties, &self.wf_aligner.penalties.affine_penalties)
        elif self.wf_aligner.penalties.distance_metric == wfa.gap_affine_2p:
            wfa.wavefront_penalties_set_affine2p(&self.wf_aligner.penalties, &self.wf_aligner.penalties.affine2p_penalties)

    @property
    def distance(self):
        if self.wf_aligner.penalties.distance_metric == wfa.indel:
            return "indel"
        elif self.wf_aligner.penalties.distance_metric == wfa.edit:
            return "levenshtein"
        elif self.wf_aligner.penalties.distance_metric == wfa.gap_linear:
            return "linear"
        elif self.wf_aligner.penalties.distance_metric == wfa.gap_affine:
            return "affine"
        elif self.wf_aligner.penalties.distance_metric == wfa.gap_affine_2p:
            return "affine2p"

    @distance.setter
    def distance(self, distance):
        if distance == "indel":
            self.wf_aligner.penalties.distance_metric = wfa.indel
        elif distance == "levenshtein":
            self.wf_aligner.penalties.distance_metric = wfa.edit
        elif distance == "linear":
            self.wf_aligner.penalties.distance_metric = wfa.gap_linear
        elif distance == "affine":
            self.wf_aligner.penalties.distance_metric = wfa.gap_affine
        elif distance == "affine2p":
            self.wf_aligner.penalties.distance_metric = wfa.gap_affine_2p
        else:
            raise NotImplementedError(f'{distance} distance not implemented')
        self._edit_penalties()

    @property
    def match_score(self):
        return self.wf_aligner.penalties.match

    @match_score.setter
    def match_score(self, int match):
        self.wf_aligner.penalties.linear_penalties.match = self.wf_aligner.penalties.affine_penalties.match = self.wf_aligner.penalties.affine2p_penalties.match = match
        self._edit_penalties()

    @property
    def mismatch_penalty(self):
        return self.wf_aligner.penalties.mismatch

    @mismatch_penalty.setter
    def mismatch_penalty(self, int mismatch):
        self.wf_aligner.penalties.linear_penalties.mismatch = self.wf_aligner.penalties.affine_penalties.mismatch = self.wf_aligner.penalties.affine2p_penalties.mismatch = mismatch
        self._edit_penalties()

    @property
    def gap_opening_penalty(self):
        return self.wf_aligner.penalties.gap_opening1

    @gap_opening_penalty.setter
    def gap_opening_penalty(self, int penalty):
        self.wf_aligner.penalties.linear_penalties.indel = self.wf_aligner.penalties.affine_penalties.gap_opening = self.wf_aligner.penalties.affine2p_penalties.gap_opening1 = penalty
        self._edit_penalties()

    @property
    def gap_extension_penalty(self):
        return self.wf_aligner.penalties.gap_extension1

    @gap_extension_penalty.setter
    def gap_extension_penalty(self, int penalty):
        self.wf_aligner.penalties.linear_penalties.indel = self.wf_aligner.penalties.affine_penalties.gap_extension = self.wf_aligner.penalties.affine2p_penalties.gap_extension1 = penalty
        self._edit_penalties()

    @property
    def gap_opening2_penalty(self):
        return self.wf_aligner.penalties.gap_opening2

    @gap_opening2_penalty.setter
    def gap_opening2_penalty(self, int penalty):
        self.wf_aligner.penalties.affine2p_penalties.gap_opening2 = penalty
        self._edit_penalties()

    @property
    def gap_extension2_penalty(self):
        return self.wf_aligner.penalties.gap_extension2

    @gap_extension2_penalty.setter
    def gap_extension2_penalty(self, int penalty):
        self.wf_aligner.penalties.affine2p_penalties.gap_extension2 = penalty
        self._edit_penalties()

    @property
    def wildcard(self):
        return self._wildcard

    @wildcard.setter
    def wildcard(self, wildcard):
        if wildcard is not None:
            if not isinstance(wildcard, str):
                raise TypeError(f"expected wildcard to be a string, but it is {type(wildcard)}")
            if len(wildcard) > 1:
                raise ValueError(f"wildcard must have length 1, but has length {len(wildcard)}")
            self._wildcard = wildcard
            self._bwildcard = wildcard.upper().encode("ascii")[0]
        else:
            self._wildcard = None

    @property
    def max_steps(self):
        return self.wf_aligner.system.max_alignment_steps

    @max_steps.setter
    def max_steps(self, int steps):
        if steps <= 0:
            steps = INT_MAX
        wfa.wavefront_aligner_set_max_alignment_steps(self.wf_aligner, steps)

    @property
    def cigarstring(self):
        cdef wfa.cigar_t* cigar
        cdef char last_op
        cdef int last_op_length, i, length

        cigar = self.wf_aligner.cigar

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

        cigar = self.wf_aligner.cigar

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
        if self.scope == "score":
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
        """
        Align `text` sequence to `pattern` sequence.
        :param text: The text sequence
        :type text: str
        :param pattern: The pattern sequence
        :type pattern: str
        :param clip_cigar: Converts unaligned bases at the flank to soft-clips. Also applies min_aligned_bases threshold
        :type clip_cigar: bool
        :param min_aligned_bases_left: The minimum number of aligned bases at the left flank during clipping
        :type min_aligned_bases_left: int
        :param min_aligned_bases_right: The minimum number of aligned bases at the right flank during clipping
        :type min_aligned_bases_right: int
        :param elide_mismatches: Convert mismatch cigar opperation (X) to matches (M)
        :type elide_mismatches: bool
        :param supress_sequences: Dont output the text and pattern sequences
        :type supress_sequences: bool
        :return: AlignmentResult dataclass
        :rtype: dataclass
        """
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
        if not self.scope == "full":
            if clip_cigar:
                res = clip_cigartuples(res, min_aligned_bases_left, min_aligned_bases_right)
            if elide_mismatches:
                res.cigartuples = elide_mismatches_from_cigar(res.cigartuples)
        return res

    def __dealloc__(self):
        if self.wf_aligner: # if an exception is raised in the constructor, self.wf_aligner does not exist yet
            wfa.wavefront_aligner_delete(self.wf_aligner)
