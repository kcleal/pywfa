
from libc.stdint cimport uint64_t
from libc.stdio cimport FILE


cdef extern from "WFA2_lib/wavefront/wavefront_penalties.h" nogil:
    ctypedef struct linear_penalties_t:
        int match          # (Penalty representation usually M <= 0)
        int mismatch       # (Penalty representation usually X > 0)
        int indel          # (Penalty representation usually I > 0)


cdef extern from "WFA2_lib/alignment/affine_penalties.h" nogil:
    # Affine penalties
    ctypedef struct affine_penalties_t:
        int match              # (Penalty representation usually M <= 0)
        int mismatch           # (Penalty representation usually X > 0)
        int gap_opening        # (Penalty representation usually O > 0)
        int gap_extension      # (Penalty representation usually E > 0)

    # Affine matrix-type (for backtrace)
    ctypedef enum affine_matrix_type:
        affine_matrix_M
        affine_matrix_I
        affine_matrix_D


cdef extern from "WFA2_lib/alignment/affine2p_penalties.h" nogil:
    ctypedef struct affine2p_penalties_t:
        int match             # (Penalty representation usually M <= 0)
        int mismatch          # (Penalty representation usually X > 0)
        # Usually concave Q1 + E1 < Q2 + E2 and E1 > E2.
        int gap_opening1      # (Penalty representation usually O1 > 0)
        int gap_extension1    # (Penalty representation usually E1 > 0)
        int gap_opening2      # (Penalty representation usually O2 > 0)
        int gap_extension2    # (Penalty representation usually E2 > 0)

    # Affine 2-piece matrix-type (for bcktrace)
    ctypedef enum affine2p_matrix_type:
        affine2p_matrix_M
        affine2p_matrix_I1
        affine2p_matrix_I2
        affine2p_matrix_D1
        affine2p_matrix_D2


cdef extern from "WFA2_lib/wavefront/wavefront_heuristic.h" nogil:
    # Wavefront ahead definition
    # ctypedef struct _wavefront_aligner_t wavefront_aligner_t

    # Wavefront Heuristics
    ctypedef enum wf_heuristic_strategy:
        wf_heuristic_none            = 0x0000000000000000ul
        wf_heuristic_banded_static   = 0x0000000000000001ul
        wf_heuristic_banded_adaptive = 0x0000000000000002ul
        wf_heuristic_wfadaptive      = 0x0000000000000004ul
        wf_heuristic_xdrop           = 0x0000000000000010ul
        wf_heuristic_zdrop           = 0x0000000000000020ul

    ctypedef struct wavefront_heuristic_t:
        # Heuristic
        wf_heuristic_strategy strategy     # Heuristic strategy
        int steps_between_cutoffs          # Score-steps between heuristic cut-offs
        # Banded
        int min_k                          # Banded: Minimum k to consider in band
        int max_k                          # Banded: Maximum k to consider in band
        # Adaptive
        int min_wavefront_length           # Adaptive: Minimum wavefronts length to cut-off
        int max_distance_threshold         # Adaptive: Maximum distance between offsets allowed
        # Drops
        int xdrop                          # X-drop parameter
        int zdrop                          # Z-drop parameter
        # Internals
        int steps_wait                     # Score-steps until next cut-off
        int max_sw_score                   # Maximum score observed (for x/z drops)
        int max_sw_score_offset            # Offset of the maximum score observed
        int max_sw_score_k                 # Diagonal of the maximum score observed


cdef extern from "WFA2_lib/utils/vector.h" nogil:
    # Data Structures
    ctypedef struct vector_t:
        void* memory
        uint64_t used
        uint64_t element_size
        uint64_t elements_allocated


cdef extern from "WFA2_lib/system/mm_allocator.h" nogil:
    # MM-Allocator
    ctypedef struct mm_allocator_t:
        # Metadata
        uint64_t request_ticker        # Request ticker
        # Memory segments
        uint64_t segment_size          # Memory segment size (bytes)
        vector_t* segments             # Memory segments (mm_allocator_segment_t*)
        vector_t* segments_free        # Completely free segments (mm_allocator_segment_t*)
        uint64_t current_segment_idx   # Current segment being used (serving memory)
        # Malloc memory
        vector_t* malloc_requests      # Malloc requests (mm_malloc_request_t)
        uint64_t malloc_requests_freed # Total malloc request freed and still in vector

    void mm_allocator_free(
        mm_allocator_t* const mm_allocator,
        void* const memory)


cdef extern from "WFA2_lib/wavefront/wavefront_attributes.h" nogil:

    # Alignment scope
    ctypedef enum alignment_scope_t:
        compute_score          # Only distance/score
        compute_alignment     # Full alignment CIGAR

    ctypedef enum alignment_span_t:
        alignment_end2end     # End-to-end alignment (aka global)
        alignment_endsfree    # Ends-free alignment  (semiglobal, glocal, etc)

    ctypedef struct alignment_form_t:
        # Mode
        alignment_span_t span   # Alignment form (End-to-end/Ends-free)
        # Ends-free
        int pattern_begin_free  # Allow free-gap at the beginning of the pattern
        int pattern_end_free    # Allow free-gap at the end of the pattern
        int text_begin_free     # Allow free-gap at the beginning of the text
        int text_end_free       # Allow free-gap at the end of the text
        # Limits
        int max_alignment_score # Maximum score allowed before quit

    # Custom extend-match function
    ctypedef int (*alignment_match_funct_t)(int,int,void*)

     # Alignment system configuration

    ctypedef struct alignment_system_t:
        # Debug
        bint check_alignment_correct  # Verify that the alignment CIGAR output is correct
        # Probing intervals
        int probe_interval_global     # Score-ticks interval to check any limits
        int probe_interval_compact    # Score-ticks interval to check BT-buffer compacting
        # Memory
        uint64_t max_partial_compacts # Maximum partial-compacts before attempting full-compact
        uint64_t max_memory_compact   # Maximum BT-buffer memory allowed before trigger compact
        uint64_t max_memory_resident  # Maximum memory allowed to be buffered before reap
        uint64_t max_memory_abort     # Maximum memory allowed to be used before aborting alignment
        # Verbose
        #  0 - Quiet
        #  1 - Report WFA progress and heavy tasks
        #  2 - Report each sequence aligned (brief)
        #  3 - Report each sequence aligned (very verbose)
        int verbose                   # Verbose (regulates messages during alignment)
        # Profile
        # profiler_timer_t timer        # Time alignment


    # Low-memory modes
    ctypedef enum wavefront_memory_t:
        wavefront_memory_high = 0     # High-memore mode (fastest, stores all WFs explicitly)
        wavefront_memory_med = 1      # Succing-memory mode (medium, offloads half-full BT-blocks)
        wavefront_memory_low = 2      # Succing-memory mode (slow, offloads only full BT-blocks)


    #Wavefront Aligner Attributes
    ctypedef struct wavefront_aligner_attr_t:
        # Distance model
        distance_metric_t distance_metric       # Alignment metric/distance used
        alignment_scope_t alignment_scope       # Alignment scope (score only or full-CIGAR)
        alignment_form_t alignment_form         # Alignment mode (end-to-end/ends-free)
        # Penalties
        linear_penalties_t linear_penalties     # Gap-linear penalties (placeholder)
        affine_penalties_t affine_penalties     # Gap-affine penalties (placeholder)
        affine2p_penalties_t affine2p_penalties # Gap-affine-2p penalties (placeholder)
        # Heuristic strategy
        wavefront_heuristic_t heuristic         # Wavefront heuristic
        # Memory model
        wavefront_memory_t memory_mode          # Wavefront memory strategy (modular wavefronts and piggyback)
        # Custom function to compare sequences
        alignment_match_funct_t match_funct     # Custom matching function (match(v,h,args))
        void* match_funct_arguments             # Generic arguments passed to matching function (args)
        # External MM (instead of allocating one inside)
        mm_allocator_t* mm_allocator            # MM-Allocator
        # Display
        # wavefront_plot_params_t plot_params     # Wavefront plot
        # System
        alignment_system_t system               # System related parameters


cdef extern from "WFA2_lib/wavefront/wavefront_penalties.h" nogil:

    # Distance metrics
    ctypedef enum distance_metric_t:
      indel         = 0, # Longest Common Subsequence - LCS
      edit          = 1, # Levenshtein
      gap_linear    = 2, # Needleman-Wunsch
      gap_affine    = 3, # Smith-Waterman-Gotoh
      gap_affine_2p = 4  # Concave 2-pieces

    # Penalty adaptation strategy
    ctypedef enum wf_penalties_strategy_type:
      wavefronts_penalties_force_zero_match
      wavefronts_penalties_shifted_penalties

    # Wavefront Penalties
    ctypedef struct wavefronts_penalties_t:
      distance_metric_t distance_metric  # Alignment metric/distance used
      # int match          # (M = 0)
      int mismatch          # (X > 0)
      int gap_opening1      # (O1 > 0)
      int gap_extension1    # (E1 > 0)
      int gap_opening2      # (O2 > 0)
      int gap_extension2    # (E2 > 0)


cdef extern from "WFA2_lib/alignment/cigar.h" nogil:
    #CIGAR
    ctypedef struct cigar_t:
        # Operations buffer
        char* operations
        int max_operations
        int begin_offset
        int end_offset
        # Score
        int score
        # MM
        mm_allocator_t* mm_allocator

    # Setup
    # void cigar_allocate(
    #     cigar_t* const cigar,
    #     const int max_operations,
    #     mm_allocator_t* const mm_allocator)
    void cigar_clear(
        cigar_t* const cigar)
    void cigar_resize(
        cigar_t* const cigar,
        const int max_operations)
    void cigar_free(
        cigar_t* const cigar)

    # Accessors
    int cigar_get_matches(
        cigar_t* const cigar)
    void cigar_add_mismatches(
        char* const pattern,
        const int pattern_length,
        char* const text,
        const int text_length,
        cigar_t* const cigar)

    # Score
    int cigar_score_edit(
        cigar_t* const cigar)
    int cigar_score_gap_linear(
        cigar_t* const cigar,
        linear_penalties_t* const penalties)
    int cigar_score_gap_affine(
        cigar_t* const cigar,
        affine_penalties_t* const penalties)
    int cigar_score_gap_affine2p(
        cigar_t* const cigar,
        affine2p_penalties_t* const penalties)

    # Utils
    int cigar_cmp(
        cigar_t* const cigar_a,
        cigar_t* const cigar_b)
    void cigar_copy(
        cigar_t* const cigar_dst,
        cigar_t* const cigar_src)
    bint cigar_check_alignment(
        FILE* const stream,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length,
        cigar_t* const cigar,
        const bint verbose)

    # Display
    void cigar_print(
        FILE* const stream,
        cigar_t* const cigar,
        const bint print_matches)
    int cigar_sprint(
        char* buffer,
        cigar_t* const cigar,
        const bint print_matches)
    void cigar_print_pretty(
        FILE* const stream,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length,
        cigar_t* const cigar,
        mm_allocator_t* const mm_allocator)



cdef extern from "WFA2_lib/wavefront/wavefront_aligner.h" nogil:
    # Error codes & messages
    DEF WF_STATUS_SUCCESSFUL = 0
    DEF WF_STATUS_IN_PROGRESS = 1
    DEF WF_STATUS_HEURISTICALY_DROPPED = -1
    DEF WF_STATUS_MAX_SCORE_REACHED = -2
    DEF WF_STATUS_OOM = -3
    extern char* wf_error_msg[5]
    char* wavefront_align_strerror(const int wf_error_code)

    # Alignment status
    ctypedef struct _wavefront_aligner_t
    ctypedef _wavefront_aligner_t wavefront_aligner_t
    ctypedef struct wavefront_align_status_t:
        # Status
        int status                                                     # Status code
        int score                                                      # Current WF-alignment score
        # Wavefront alignment functions
        void (*wf_align_compute)(wavefront_aligner_t* const,const int) # WF Compute function
        bint (*wf_align_extend)(wavefront_aligner_t* const,const int)  # WF Extend function

    # Wavefront Aligner
    ctypedef struct _wavefront_aligner_t:
        # Status
        wavefront_align_status_t align_status   # Current alignment status
        # Sequences
        # strings_padded_t* sequences             # Padded sequences
        char* pattern                           # Pattern sequence (padded)
        int pattern_length                      # Pattern length
        char* text                              # Text sequence (padded)
        int text_length                         # Text length
        # Alignment Attributes
        alignment_scope_t alignment_scope       # Alignment scope (score only or full-CIGAR)
        alignment_form_t alignment_form         # Alignment form (end-to-end/ends-free)
        wavefronts_penalties_t penalties        # Alignment penalties
        wavefront_heuristic_t heuristic         # Heuristic's parameters
        wavefront_memory_t memory_mode          # Wavefront memory strategy (modular wavefronts and piggyback)
        # Custom function to compare sequences
        # alignment_match_funct_t match_funct     # Custom matching function (match(v,h,args))
        # void* match_funct_arguments             # Generic arguments passed to matching function (args)
        # Wavefront components
        # wavefront_components_t wf_components    # Wavefront components
        # CIGAR
        cigar_t cigar                           # Alignment CIGAR
        # MM
        # bint mm_allocator_own                   # Ownership of MM-Allocator
        # mm_allocator_t* mm_allocator            # MM-Allocator
        # wavefront_slab_t* wavefront_slab        # MM-Wavefront-Slab (Allocates/Reuses the individual wavefronts)
        # Display
        # wavefront_plot_params_t plot_params     # Wavefront plot parameters
        # wavefront_plot_t wf_plot                # Wavefront plot
        # System
        # alignment_system_t system               # System related parameters

    # ctypedef _wavefront_aligner_t wavefront_aligner_t

    # Setup
    wavefront_aligner_t* wavefront_aligner_new(
        wavefront_aligner_attr_t* attributes)
    void wavefront_aligner_resize(
        wavefront_aligner_t* const wf_aligner,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length)
    void wavefront_aligner_reap(
        wavefront_aligner_t* const wf_aligner)
    void wavefront_aligner_delete(
        wavefront_aligner_t* const wf_aligner)

    # Span configuration
    void wavefront_aligner_set_alignment_end_to_end(
        wavefront_aligner_t* const wf_aligner)
    void wavefront_aligner_set_alignment_free_ends(
        wavefront_aligner_t* const wf_aligner,
        const int pattern_begin_free,
        const int pattern_end_free,
        const int text_begin_free,
        const int text_end_free)

    # Heuristic configuration
    void wavefront_aligner_set_heuristic_none(
        wavefront_aligner_t* const wf_aligner)
    void wavefront_aligner_set_heuristic_banded_static(
        wavefront_aligner_t* const wf_aligner,
        const int band_min_k,
        const int band_max_k)
    void wavefront_aligner_set_heuristic_banded_adaptive(
        wavefront_aligner_t* const wf_aligner,
        const int band_min_k,
        const int band_max_k,
        const int score_steps)
    void wavefront_aligner_set_heuristic_wfadaptive(
        wavefront_aligner_t* const wf_aligner,
        const int min_wavefront_length,
        const int max_distance_threshold,
        const int score_steps)
    void wavefront_aligner_set_heuristic_xdrop(
        wavefront_aligner_t* const wf_aligner,
        const int xdrop,
        const int score_steps)
    void wavefront_aligner_set_heuristic_zdrop(
        wavefront_aligner_t* const wf_aligner,
        const int ydrop,
        const int score_steps)

    #  Match-funct configuration
    void wavefront_aligner_set_match_funct(
        wavefront_aligner_t* const wf_aligner,
        int (*match_funct)(int,int,void*),
        void* const match_funct_arguments)

    # System configuration
    void wavefront_aligner_set_max_alignment_score(
        wavefront_aligner_t* const wf_aligner,
        const int max_alignment_score)
    void wavefront_aligner_set_max_memory(
        wavefront_aligner_t* const wf_aligner,
        const uint64_t max_memory_compact,
        const uint64_t max_memory_resident,
        const uint64_t max_memory_abort)

    # Utils
    uint64_t wavefront_aligner_get_size(
        wavefront_aligner_t* const wf_aligner)

    # /*
    #  * Display
    #  */
    # void wavefront_aligner_print_status(
    #     FILE* const stream,
    #     wavefront_aligner_t* const wf_aligner,
    #     const int current_score)


cdef extern from "WFA2_lib/wavefront/wavefront_heuristic.h" nogil:
    # Wavefront ahead definition
    # ctypedef struct _wavefront_aligner_t wavefront_aligner_t

    # Wavefront Heuristics
    # ctypedef enum wf_heuristic_strategy:
    #     wf_heuristic_none            = 0x0000000000000000ul
    #     wf_heuristic_banded_static   = 0x0000000000000001ul
    #     wf_heuristic_banded_adaptive = 0x0000000000000002ul
    #     wf_heuristic_wfadaptive      = 0x0000000000000004ul
    #     wf_heuristic_xdrop           = 0x0000000000000010ul
    #     wf_heuristic_zdrop           = 0x0000000000000020ul
    #
    # ctypedef struct wavefront_heuristic_t:
    #     # Heuristic
    #     wf_heuristic_strategy strategy     # Heuristic strategy
    #     int steps_between_cutoffs          # Score-steps between heuristic cut-offs
    #     # Banded
    #     int min_k                          # Banded: Minimum k to consider in band
    #     int max_k                          # Banded: Maximum k to consider in band
    #     # Adaptive
    #     int min_wavefront_length           # Adaptive: Minimum wavefronts length to cut-off
    #     int max_distance_threshold         # Adaptive: Maximum distance between offsets allowed
    #     # Drops
    #     int xdrop                          # X-drop parameter
    #     int zdrop                          # Z-drop parameter
    #     # Internals
    #     int steps_wait                     # Score-steps until next cut-off
    #     int max_sw_score                   # Maximum score observed (for x/z drops)
    #     int max_sw_score_offset            # Offset of the maximum score observed
    #     int max_sw_score_k                 # Diagonal of the maximum score observed

    # Setup
    void wavefront_heuristic_set_none(
        wavefront_heuristic_t* const wf_heuristic)
    void wavefront_heuristic_set_banded_static(
        wavefront_heuristic_t* const wf_heuristic,
        const int band_min_k,
        const int band_max_k)
    void wavefront_heuristic_set_banded_adaptive(
        wavefront_heuristic_t* const wf_heuristic,
        const int band_min_k,
        const int band_max_k,
        const int steps_between_cutoffs)
    void wavefront_heuristic_set_wfadaptive(
        wavefront_heuristic_t* const wf_heuristic,
        const int min_wavefront_length,
        const int max_distance_threshold,
        const int steps_between_cutoffs)
    void wavefront_heuristic_set_xdrop(
        wavefront_heuristic_t* const wf_heuristic,
        const int xdrop,
        const int steps_between_cutoffs)
    void wavefront_heuristic_set_zdrop(
        wavefront_heuristic_t* const wf_heuristic,
        const int ydrop,
        const int steps_between_cutoffs)

    void wavefront_heuristic_clear(
        wavefront_heuristic_t* const wf_heuristic)

    # Wavefront heuristic cut-off
    bint wavefront_heuristic_cufoff(
        wavefront_aligner_t* const wf_aligner,
        const int score)


    # Default parameters
cdef extern from "WFA2_lib/wavefront/wavefront_attributes.c" nogil:
    cdef extern wavefront_aligner_attr_t wavefront_aligner_attr_default



cdef extern from "WFA2_lib/wavefront/wavefront_penalties.h" nogil:

    # Distance metrics
    # ctypedef enum distance_metric_t:
    #   indel         = 0, # Longest Common Subsequence - LCS
    #   edit          = 1, # Levenshtein
    #   gap_linear    = 2, # Needleman-Wunsch
    #   gap_affine    = 3, # Smith-Waterman-Gotoh
    #   gap_affine_2p = 4  # Concave 2-pieces
    #
    # # Penalty adaptation strategy
    # ctypedef enum wf_penalties_strategy_type:
    #   wavefronts_penalties_force_zero_match
    #   wavefronts_penalties_shifted_penalties
    #
    # # Wavefront Penalties
    # ctypedef struct wavefronts_penalties_t:
    #   distance_metric_t distance_metric  # Alignment metric/distance used
    #   # int match          # (M = 0)
    #   int mismatch          # (X > 0)
    #   int gap_opening1      # (O1 > 0)
    #   int gap_extension1    # (E1 > 0)
    #   int gap_opening2      # (O2 > 0)
    #   int gap_extension2    # (E2 > 0)

    # Penalties adjustment
    void wavefronts_penalties_set_indel(
        wavefronts_penalties_t* const wavefronts_penalties)
    void wavefronts_penalties_set_edit(
        wavefronts_penalties_t* const wavefronts_penalties)
    void wavefronts_penalties_set_linear(
        wavefronts_penalties_t* const wavefronts_penalties,
        linear_penalties_t* const linear_penalties,
        const wf_penalties_strategy_type penalties_strategy)
    void wavefronts_penalties_set_affine(
        wavefronts_penalties_t* const wavefronts_penalties,
        affine_penalties_t* const affine_penalties,
        const wf_penalties_strategy_type penalties_strategy)
    void wavefronts_penalties_set_affine2p(
        wavefronts_penalties_t* const wavefronts_penalties,
        affine2p_penalties_t* const affine2p_penalties,
        const wf_penalties_strategy_type penalties_strategy)

    # Display
    # void wavefronts_penalties_print(
    #     FILE* const stream,
    #     wavefronts_penalties_t* const wavefronts_penalties)

    # Display
    # void wavefront_heuristic_print(
    #     FILE* const stream,
    #     wavefront_heuristic_t* const wf_heuristic)


cdef extern from "WFA2_lib/wavefront/wavefront_align.h" nogil:
    # Wavefront Alignment
    int wavefront_align(
        wavefront_aligner_t* const wf_aligner,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length)
    int wavefront_align_resume(
        wavefront_aligner_t* const wf_aligner)


