#cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False

from pywfa cimport WFA_wrap as wfa


cdef class WavefrontAligner:

    cdef wfa.wavefront_aligner_attr_t* attributes
    cdef wfa.wavefront_aligner_t* wf_aligner
    cdef str _pattern, _text
    cdef bytes _bpattern
    cdef char _wildcard
    cdef bint score_only
    cdef public int match_score, alignment_score, text_len, pattern_len
