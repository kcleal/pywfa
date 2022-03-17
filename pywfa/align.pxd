
import pywfa
from pywfa cimport WFA_wrap as wfa


cdef class WavefrontAligner:

    cdef wfa.wavefront_aligner_attr_t* attributes
    cdef wfa.wavefront_aligner_t* wf_aligner
    cdef str _pattern
    cdef bint score_only
    cdef public int match_score, alignment_score, text_len, pattern_len
    # cpdef int wavefront_align(self, text, pattern=*)
