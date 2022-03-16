
from pywfa cimport WFA_wrap as wfa


cdef class WavefrontAligner:

    cdef wfa.wavefront_aligner_attr_t* attributes
    cdef wfa.wavefront_aligner_t* wf_aligner
    cpdef str _pattern
    cdef bint score_only
    cpdef public int match_score, alignment_score

    cpdef int wavefront_align(self, text, pattern=*)
