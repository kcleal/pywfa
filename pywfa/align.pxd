
from pywfa cimport WFA_wrap as wfa


cdef clip_cigartuples(object align_result, int min_aligned_bases_left=*, int min_aligned_bases_right=*)
cdef elide_mismatches_from_cigar(object cigartuples)
cdef cigartuples_to_str(object cigartuples)


cdef class WavefrontAligner:

    cdef wfa.wavefront_aligner_attr_t* attributes
    cdef wfa.wavefront_aligner_t* wf_aligner
    cdef str _pattern
    cdef bint score_only
    cdef public int match_score, alignment_score

    cdef int wavefront_align(self, text, pattern=*)
