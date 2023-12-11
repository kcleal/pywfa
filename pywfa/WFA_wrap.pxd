#cython: language_level=3, boundscheck=False, wraparound=False, nonecheck=False

from libc.stdint cimport uint8_t, int32_t, uint32_t, uint64_t
from libc.stdio cimport FILE
from posix.time cimport timespec


cdef extern from "WFA2_lib/utils/vector.h" nogil:
    # Data Structures
    ctypedef struct vector_t:
        void* memory
        uint64_t used
        uint64_t element_size
        uint64_t elements_allocated


cdef extern from "WFA2_lib/alignment/linear_penalties.h" nogil:
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
    ctypedef _wavefront_aligner_t wavefront_aligner_t

    # Wavefront Heuristics
    ctypedef enum wf_heuristic_strategy:
        wf_heuristic_none            = 0x0000000000000000ul
        wf_heuristic_banded_static   = 0x0000000000000001ul
        wf_heuristic_banded_adaptive = 0x0000000000000002ul
        wf_heuristic_wfadaptive      = 0x0000000000000004ul
        wf_heuristic_xdrop           = 0x0000000000000010ul
        wf_heuristic_zdrop           = 0x0000000000000020ul
        wf_heuristic_wfmash          = 0x0000000000000040ul

    ctypedef struct wavefront_heuristic_t:
        # Heuristic
        wf_heuristic_strategy strategy     # Heuristic strategy
        int steps_between_cutoffs          # Score-steps between heuristic cut-offs
        # Static/Adaptive Banded
        int min_k                          # Banded: Minimum k to consider in band
        int max_k                          # Banded: Maximum k to consider in band
        # WFAdaptive
        int min_wavefront_length           # Adaptive: Minimum wavefronts length to cut-off
        int max_distance_threshold         # Adaptive: Maximum distance between offsets allowed
        # Drops
        int xdrop                          # X-drop parameter
        int zdrop                          # Z-drop parameter
        # Internals
        int steps_wait                     # Score-steps until next cut-off
        int max_sw_score                   # Maximum score observed (for x/z drops)
        int max_wf_score                   # Corresponding WF-score (to max_sw_score)
        int max_sw_score_offset            # Offset of the maximum score observed
        int max_sw_score_k                 # Diagonal of the maximum score observed

    # Setup
    void wavefront_heuristic_set_none(
        wavefront_heuristic_t* const wf_heuristic)

    void wavefront_heuristic_set_wfadaptive(
        wavefront_heuristic_t* const wf_heuristic,
        const int min_wavefront_length,
        const int max_distance_threshold,
        const int steps_between_cutoffs)
    void wavefront_heuristic_set_wfmash(
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

    void wavefront_heuristic_set_banded_static(
        wavefront_heuristic_t* const wf_heuristic,
        const int band_min_k,
        const int band_max_k)
    void wavefront_heuristic_set_banded_adaptive(
        wavefront_heuristic_t* const wf_heuristic,
        const int band_min_k,
        const int band_max_k,
        const int steps_between_cutoffs)

    # Wavefront heuristic cut-off
    bint wavefront_heuristic_cufoff(
        wavefront_aligner_t* const wf_aligner,
        const int score,
        const int score_mod)


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

cdef extern from "WFA2_lib/system/profiler_counter.h" nogil:
    # Counters
    ctypedef struct profiler_counter_t:
        uint64_t total
        uint64_t samples
        uint64_t min
        uint64_t max
        double m_oldM
        double m_newM
        double m_oldS
        double m_newS

    void counter_reset(
        profiler_counter_t* const counter)
    void counter_add(
        profiler_counter_t* const counter,
        const uint64_t amount)

    uint64_t counter_get_total(const profiler_counter_t* const counter)
    uint64_t counter_get_num_samples(const profiler_counter_t* const counter)
    uint64_t counter_get_min(const profiler_counter_t* const counter)
    uint64_t counter_get_max(const profiler_counter_t* const counter)
    double counter_get_mean(const profiler_counter_t* const counter)
    double counter_get_variance(const profiler_counter_t* const counter)
    double counter_get_stddev(const profiler_counter_t* const counter)

    void counter_combine_sum(
        profiler_counter_t* const counter_dst,
        profiler_counter_t* const counter_src)

    void counter_print(
        FILE* const stream,
        const profiler_counter_t* const counter,
        const profiler_counter_t* const ref_counter,
        const char* const units,
        const bint full_report)
    void percentage_print(
        FILE* const stream,
        const profiler_counter_t* const counter,
        const char* const units)

    # Reference Counter (Counts wrt a reference counter. Eg ranks)
    ctypedef struct profiler_rcounter_t:
        uint64_t begin_count       # Counter
        profiler_counter_t counter # Total count & samples taken
        uint64_t accumulated       # Total accumulated

    void rcounter_start(
        profiler_rcounter_t* const rcounter,
        const uint64_t reference)
    void rcounter_stop(
        profiler_rcounter_t* const rcounter,
        const uint64_t reference)
    void rcounter_pause(
        profiler_rcounter_t* const rcounter,
        const uint64_t reference)
    void rcounter_restart(
        profiler_rcounter_t* const rcounter,
        const uint64_t reference)
    void rcounter_reset(
        profiler_rcounter_t* const rcounter)

    uint64_t rcounter_get_total(profiler_rcounter_t* const rcounter)
    uint64_t rcounter_get_num_samples(profiler_rcounter_t* const rcounter)
    uint64_t rcounter_get_min(profiler_rcounter_t* const rcounter)
    uint64_t rcounter_get_max(profiler_rcounter_t* const rcounter)
    uint64_t rcounter_get_mean(profiler_rcounter_t* const rcounter)
    uint64_t rcounter_get_variance(profiler_rcounter_t* const rcounter)
    uint64_t rcounter_get_stddev(profiler_rcounter_t* const rcounter)

cdef extern from "WFA2_lib/system/profiler_timer.h" nogil:
    # System time
    void timer_get_system_time(timespec *ts);

    # Timers
    ctypedef struct profiler_timer_t:
        # Timer
        timespec begin_timer;     # Timer begin
        # Total time & samples taken
        profiler_counter_t time_ns;
        uint64_t accumulated;

    void timer_start(profiler_timer_t* const timer);
    void timer_stop(profiler_timer_t* const timer);
    void timer_pause(profiler_timer_t* const timer);
    void timer_continue(profiler_timer_t* const timer);
    void timer_reset(profiler_timer_t* const timer);

    uint64_t timer_get_current_lap_ns(profiler_timer_t* const timer);
    uint64_t timer_get_current_total_ns(profiler_timer_t* const timer);
    uint64_t timer_get_total_ns(const profiler_timer_t* const timer);
    uint64_t timer_get_num_samples(const profiler_timer_t* const timer);
    uint64_t timer_get_min_ns(const profiler_timer_t* const timer);
    uint64_t timer_get_max_ns(const profiler_timer_t* const timer);
    uint64_t timer_get_mean(const profiler_timer_t* const timer);
    uint64_t timer_get_variance(const profiler_timer_t* const timer);
    uint64_t timer_get_stddev(const profiler_timer_t* const timer);

    void timer_print_total(
        FILE* const stream,
        const profiler_timer_t* const timer);

    void timer_print(
        FILE* const stream,
        const profiler_timer_t* const timer,
        const profiler_timer_t* const ref_timer);

cdef extern from "WFA2_lib/utils/heatmap.h" nogil:
    # Heatmap
    ctypedef enum heatmap_type:
        heatmap_min   # Min value stays
        heatmap_max   # Max value stays
        heatmap_value # Last value set stays
    ctypedef struct heatmap_t:
        # Configuration
        heatmap_type type;
        # Dimensions
        int num_rows;
        int num_columns;
        # Range
        int min_v;
        int max_v;
        int min_h;
        int max_h;
        float binning_factor;
        # Data
        int** values;

    # Setup
    heatmap_t* heatmap_new(
        const heatmap_type type,
        const int min_v,
        const int max_v,
        const int min_h,
        const int max_h,
        const int resolution_points)
    void heatmap_clear(
        heatmap_t* const heatmap)
    void heatmap_delete(
        heatmap_t* const heatmap)

    # Accessors
    void heatmap_set(
        heatmap_t* const heatmap,
        const int v,
        const int h,
        const int value)

    # Display
    void heatmap_print(
        FILE* const stream,
        heatmap_t* const heatmap)

cdef extern from "WFA2_lib/wavefront/wavefront_plot.h" nogil:
    # Wavefront ahead definition
    # ctypedef _wavefront_aligner_t wavefront_aligner_t

    # Wavefront Display
    ctypedef struct wavefront_plot_attr_t:
        bint enabled               # Is plotting enabled
        int resolution_points      # Total resolution points
        int align_level            # Level of recursion to plot (-1 == final)
    ctypedef struct wavefront_plot_t:
        # Configuration
        wavefront_plot_attr_t attributes
        distance_metric_t distance_metric
        int min_v
        int max_v
        int min_h
        int max_h
        # Wavefront Heatmaps
        heatmap_t* m_heatmap
        heatmap_t* i1_heatmap
        heatmap_t* d1_heatmap
        heatmap_t* i2_heatmap
        heatmap_t* d2_heatmap
        heatmap_t* behavior_heatmap

    # Setup
    wavefront_plot_t* wavefront_plot_new(
        const distance_metric_t distance_metric,
        const int pattern_length,
        const int text_length,
        wavefront_plot_attr_t* const attributes)
    void wavefront_plot_resize(
        wavefront_plot_t* const wf_plot,
        const int pattern_length,
        const int text_length)
    void wavefront_plot_delete(
        wavefront_plot_t* const wf_plot)

    # Plot record state
    void wavefront_plot(
        wavefront_aligner_t* const wf_aligner,
        const int score,
        const int align_level)

    # Display/Dump
    void wavefront_plot_print(
        FILE* const stream,
        wavefront_aligner_t* const wf_aligner)


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
        # Extension
        bint extension          # Activate extension-like alignment
        # Ends-free
        int pattern_begin_free  # Allow free-gap at the beginning of the pattern
        int pattern_end_free    # Allow free-gap at the end of the pattern
        int text_begin_free     # Allow free-gap at the beginning of the text
        int text_end_free       # Allow free-gap at the end of the text

     # Alignment system configuration
    ctypedef struct alignment_system_t:
        # Limits
        int max_alignment_steps       # Maximum WFA-steps allowed before quit
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
        #  1 - Report each sequence aligned                     (brief)
        #  2 - Report each sequence/subsequence aligned         (brief)
        #  3 - Report WFA progress (heavy tasks)                (verbose)
        #  4 - Full report of each sequence/subsequence aligned (very verbose)
        int verbose                   # Verbose (regulates messages during alignment)
        # Debug
        bint check_alignment_correct  # Verify that the alignment CIGAR output is correct
        # Profile
        profiler_timer_t timer        # Time alignment
        # OS
        int max_num_threads           # Maximum number of threads to use to compute/extend WFs
        int min_offsets_per_thread    # Minumum amount of offsets to spawn a thread


    # Low-memory modes
    ctypedef enum wavefront_memory_t:
        wavefront_memory_high = 0     # High-memore mode (fastest, stores all WFs explicitly)
        wavefront_memory_med = 1      # Succing-memory mode piggyback-based (medium, offloads half-full BT-blocks)
        wavefront_memory_low = 2      # Succing-memory mode piggyback-based (slow, offloads only full BT-blocks)
        wavefront_memory_ultralow = 3 # Bidirectional WFA


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
        # External MM (instead of allocating one inside)
        mm_allocator_t* mm_allocator            # MM-Allocator
        # Display
        wavefront_plot_attr_t plot              # Plot wavefront
        # System
        alignment_system_t system               # System related parameters

    # Default parameters
    cdef extern wavefront_aligner_attr_t wavefront_aligner_attr_default


cdef extern from "WFA2_lib/wavefront/wavefront_penalties.h" nogil:

    # Distance metrics
    ctypedef enum distance_metric_t:
        indel         = 0, # Longest Common Subsequence - LCS
        edit          = 1, # Levenshtein
        gap_linear    = 2, # Needleman-Wunsch
        gap_affine    = 3, # Smith-Waterman-Gotoh
        gap_affine_2p = 4  # Gap-Affine 2-pieces

    # Wavefront Penalties
    ctypedef struct wavefront_penalties_t:
        distance_metric_t distance_metric  # Alignment metric/distance used
        int match             # (M <= 0)
        int mismatch          # (X > 0)
        int gap_opening1      # (O1 >= 0)
        int gap_extension1    # (E1 > 0)
        int gap_opening2      # (O2 >= 0)
        int gap_extension2    # (E2 > 0)
        # Internals
        linear_penalties_t linear_penalties     # Original gap-linear penalties
        affine_penalties_t affine_penalties     # Original gap-affine penalties
        affine2p_penalties_t affine2p_penalties # Original gap-affine2p penalties
        int internal_gap_e                      # Original gap-extension value (used for z-drop)

    # Penalties adjustment
    void wavefront_penalties_set_indel(
        wavefront_penalties_t* const wf_penalties)
    void wavefront_penalties_set_edit(
        wavefront_penalties_t* const wf_penalties)
    void wavefront_penalties_set_linear(
        wavefront_penalties_t* const wf_penalties,
        linear_penalties_t* const linear_penalties)
    void wavefront_penalties_set_affine(
        wavefront_penalties_t* const wf_penalties,
        affine_penalties_t* const affine_penalties)
    void wavefront_penalties_set_affine2p(
        wavefront_penalties_t* const wf_penalties,
        affine2p_penalties_t* const affine2p_penalties)

    # Display
    void wavefront_penalties_print(
        FILE* const stream,
        wavefront_penalties_t* const wf_penalties)

cdef extern from "WFA2_lib/alignment/cigar.h" nogil:
    #CIGAR
    ctypedef struct cigar_t:
        # Alignment operations
        char* operations          # Raw alignment operations
        int max_operations        # Maximum buffer size
        int begin_offset          # Eegin_offset
        int end_offset            # End offset
        # Score
        int score                 # Computed score
        int end_v                 # Alignment-end vertical coordinate (pattern characters aligned)
        int end_h                 # Alignment-end horizontal coordinate (text characters aligned)
        # CIGAR (SAM compliant)
        bint has_misms            # Show 'X' and '=' instead of  just 'M'
        uint32_t* cigar_buffer    # CIGAR-operations (max operations length)
        int cigar_length          # total CIGAR operations

    # Setup
    cigar_t* cigar_new(
        const int max_operations)
    void cigar_clear(
        cigar_t* const cigar)
    void cigar_resize(
        cigar_t* const cigar,
        const int max_operations)
    void cigar_free(
        cigar_t* const cigar)

    # Accessors
    bint cigar_is_null(
        cigar_t* const cigar)

    int cigar_count_matches(
        cigar_t* const cigar)

    void cigar_append_forward(
        cigar_t* const cigar_dst,
        cigar_t* const cigar_src
    )
    void cigar_append_reverse(
        cigar_t* const cigar_dst,
        cigar_t* const cigar_src
    )

    void cigar_append_deletion(
        cigar_t* const cigar,
        const int length
    )
    void cigar_append_insertion(
        cigar_t* const cigar,
        const int length
    )

    # SAM-compliant CIGAR
    void cigar_get_CIGAR(
        cigar_t* const cigar,
        const bint show_mismatches,
        uint32_t** const cigar_buffer,
        int* const cigar_length
    )

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

    void cigar_discover_mismatches(
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length,
        cigar_t* const cigar
    )

    bint cigar_maxtrim_gap_linear(
        cigar_t* const cigar,
        linear_penalties_t* const penalties
    )
    bint cigar_maxtrim_gap_affine(
        cigar_t* const cigar,
        affine_penalties_t* const penalties
    )
    bint cigar_maxtrim_gap_affine2p(
        cigar_t* const cigar,
        affine2p_penalties_t* const penalties
    )

    # Check
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

    void cigar_print_SAM_CIGAR(
        FILE* const stream,
        cigar_t* const cigar,
        const bint show_mismatches)
    void cigar_sprint_SAM_CIGAR(
        char* const buffer,
        cigar_t* const cigar,
        const bint show_mismatches)

    void cigar_print_pretty(
        FILE* const stream,
        cigar_t* const cigar,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length)


cdef extern from "WFA2_lib/wavefront/wavefront_sequences.h" nogil:
    ctypedef int (*alignment_match_funct_t)(int,int,void*)

    # Wavefront sequences
    ctypedef enum wf_sequences_mode_t:
        wf_sequences_ascii       = 0
        wf_sequences_lambda      = 1
        wf_sequences_packed2bits = 2
    ctypedef struct wavefront_sequences_t:
        # Mode
        wf_sequences_mode_t mode              # Sequences mode
        bint reverse                          # Reverse sequences
        # Current sequences & bounds
        char* pattern                         # Pointer to current pattern sequence (padded)
        char* text                            # Pointer to current text sequence (padded)
        int pattern_begin                     # Pattern begin offset
        int pattern_length                    # Pattern length
        int text_begin                        # Text begin offset
        int text_length                       # Text length
        # Lambda Sequence
        alignment_match_funct_t match_funct   # Custom matching function (match(v,h,args))
        void* match_funct_arguments           # Generic arguments passed to matching function (args)
        # Internal buffers (ASCII encoded)
        char* seq_buffer                      # Internal buffer
        int seq_buffer_allocated              # Internal buffer allocated
        char* pattern_buffer                  # Source pattern sequence
        char* text_buffer                     # Source text sequence
        int pattern_buffer_length             # Source pattern length
        int text_buffer_length                # Source text length
        char pattern_eos                      # Source pattern char at EOS
        char text_eos                         # Source pattern char at EOS

    # Setup
    void wavefront_sequences_allocate(
        wavefront_sequences_t* const wf_sequences);
    void wavefront_sequences_free(
        wavefront_sequences_t* const wf_sequences);

    # Init Sequences
    void wavefront_sequences_init_ascii(
        wavefront_sequences_t* const wf_sequences,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length,
        const bint reverse);
    void wavefront_sequences_init_lambda(
        wavefront_sequences_t* const wf_sequences,
        alignment_match_funct_t match_funct,
        void* match_funct_arguments,
        const int pattern_length,
        const int text_length,
        const bint reverse);
    void wavefront_sequences_init_packed2bits(
        wavefront_sequences_t* const wf_sequences,
        const uint8_t* const pattern,
        const int pattern_length,
        const uint8_t* const text,
        const int text_length,
        const bint reverse);

cdef extern from "WFA2_lib/wavefront/wavefront_offset.h" nogil:
    ctypedef int32_t wf_offset_t
    ctypedef uint32_t wf_unsigned_offset_t

cdef extern from "WFA2_lib/wavefront/wavefront_pcigar.h" nogil:
    ctypedef uint32_t pcigar_t;
    # Accessors
    int pcigar_get_length(
        const pcigar_t pcigar);
    int pcigar_unpack(
        pcigar_t pcigar,
        char* cigar_buffer);

    # PCIGAR unpack
    void pcigar_unpack_linear(
        pcigar_t pcigar,
        wavefront_sequences_t* const sequences,
        int* const v_pos,
        int* const h_pos,
        char* cigar_buffer,
        int* const cigar_length);
    void pcigar_unpack_affine(
        pcigar_t pcigar,
        wavefront_sequences_t* const sequences,
        int* const v_pos,
        int* const h_pos,
        char* cigar_buffer,
        int* const cigar_length,
        affine_matrix_type* const current_matrix_type);

    # Display
    void pcigar_print(
        FILE* const stream,
        pcigar_t pcigar);

cdef extern from "WFA2_lib/utils/bitmap.h" nogil:
    # Bitmap
    ctypedef struct bitmap_block_t:
        uint64_t counter
        uint64_t bitmap
    ctypedef struct bitmap_t:
        # Bitmap
        uint64_t num_blocks
        bitmap_block_t* bitmap_blocks
        # MM
        mm_allocator_t* mm_allocator

    # Setup
    bitmap_t* bitmap_new(
        const uint64_t length,
        mm_allocator_t* const mm_allocator)
    void bitmap_delete(
        bitmap_t* const bitmap)

    # Accessors
    void bitmap_set(
        bitmap_t* const bitmap,
        const uint64_t pos)
    bint bitmap_is_set(
        bitmap_t* const bitmap,
        const uint64_t pos)
    bint bitmap_check__set(
        bitmap_t* const bitmap,
        const uint64_t pos)

    # Rank
    void bitmap_update_counters(
        bitmap_t* const bitmap)
    uint64_t bitmap_erank(
        bitmap_t* const bitmap,
        const uint64_t pos)


cdef extern from "WFA2_lib/wavefront/wavefront_backtrace_buffer.h" nogil:
    # Separated Backtrace Block
    ctypedef uint32_t bt_block_idx_t # Up to 2^31 references (~32GB of not-compactable pCIGARs)

    ctypedef packed struct bt_block_t:
        pcigar_t pcigar            # Packed CIGAR
        bt_block_idx_t prev_idx    # Index of the previous BT-block

    # Backtrace initial positions
    ctypedef struct wf_backtrace_init_pos_t:
        int v
        int h

    # Backtrace Buffer
    ctypedef struct wf_backtrace_buffer_t:
        # Locator
        int segment_idx                     # Current segment idx
        int segment_offset                  # Current free position within segment
        bt_block_t* block_next              # Next BT-block free
        # Buffers                           #
        vector_t* segments                  # Memory segments (bt_block_t*)
        vector_t* alignment_init_pos        # Buffer to store alignment's initial coordinates (h,v) (wf_backtrace_init_pos_t)
        bt_block_idx_t num_compacted_blocks # Total compacted blocks in BT-buffer compacted (dense from 0..num_compacted_blocks-1)
        int num_compactions                 # Total compactions performed
        # Internal buffers                  #
        vector_t* alignment_packed          # Temporal buffer to store final alignment (pcigar_t)
        vector_t* prefetch_blocks_idxs      # Temporal buffer to store blocks_idxs (bt_block_idx_t)
        # MM
        mm_allocator_t* mm_allocator

    # Setup
    wf_backtrace_buffer_t* wf_backtrace_buffer_new(
        mm_allocator_t* const mm_allocator)
    void wf_backtrace_buffer_clear(
        wf_backtrace_buffer_t* const bt_buffer)
    void wf_backtrace_buffer_reap(
        wf_backtrace_buffer_t* const bt_buffer)
    void wf_backtrace_buffer_delete(
        wf_backtrace_buffer_t* const bt_buffer)

    # Accessors
    void wf_backtrace_buffer_add_used(
        wf_backtrace_buffer_t* const bt_buffer,
        const int used)
    bt_block_idx_t wf_backtrace_buffer_get_mem(
        wf_backtrace_buffer_t* const bt_buffer,
        bt_block_t** const bt_block_mem,
        int* const bt_blocks_available)

    # Store blocks
    bt_block_idx_t wf_backtrace_buffer_init_block(
        wf_backtrace_buffer_t* const bt_buffer,
        const int v,
        const int h)

    # Unpack CIGAR
    bt_block_t* wf_backtrace_buffer_traceback_pcigar(
        wf_backtrace_buffer_t* const bt_buffer,
        bt_block_t* bt_block)
    void wf_backtrace_buffer_unpack_cigar_linear(
        wf_backtrace_buffer_t* const bt_buffer,
        wavefront_sequences_t* const sequences,
        const int begin_v,
        const int begin_h,
        const int end_v,
        const int end_h,
        cigar_t* const cigar)
    void wf_backtrace_buffer_unpack_cigar_affine(
        wf_backtrace_buffer_t* const bt_buffer,
        wavefront_sequences_t* const sequences,
        const int begin_v,
        const int begin_h,
        const int end_v,
        const int end_h,
        cigar_t* const cigar)

    # Compact
    void wf_backtrace_buffer_mark_backtrace(
        wf_backtrace_buffer_t* const bt_buffer,
        const bt_block_idx_t bt_block_idx,
        bitmap_t* const bitmap)
    void wf_backtrace_buffer_mark_backtrace_batch(
        wf_backtrace_buffer_t* const bt_buffer,
        wf_offset_t* const offsets,
        bt_block_idx_t* const bt_block_idxs,
        const int num_block_idxs,
        bitmap_t* const bitmap)

    bt_block_idx_t wf_backtrace_buffer_compact_marked(
        wf_backtrace_buffer_t* const bt_buffer,
        bitmap_t* const bitmap,
        const int verbose)

    # Utils
    uint64_t wf_backtrace_buffer_get_used(
        wf_backtrace_buffer_t* const bt_buffer)

    bt_block_idx_t wf_backtrace_buffer_get_num_compacted_blocks(
        wf_backtrace_buffer_t* const bt_buffer)
    void wf_backtrace_buffer_set_num_compacted_blocks(
        wf_backtrace_buffer_t* const bt_buffer,
        const bt_block_idx_t num_compacted_blocks)
    void wf_backtrace_buffer_reset_compaction(
        wf_backtrace_buffer_t* const bt_buffer)

    uint64_t wf_backtrace_buffer_get_size_allocated(
        wf_backtrace_buffer_t* const bt_buffer)
    uint64_t wf_backtrace_buffer_get_size_used(
        wf_backtrace_buffer_t* const bt_buffer)


cdef extern from "WFA2_lib/wavefront/wavefront.h" nogil:
    # Alignment position
    ctypedef struct wavefront_pos_t:
        int score          # Score
        int k              # Diagonal
        wf_offset_t offset # Offset

    # Wavefront
    ctypedef enum wavefront_status_type:
        wavefront_status_free
        wavefront_status_busy
        wavefront_status_deallocated
    ctypedef struct wavefront_t:
        # Dimensions
        bint null                           # Is null interval?
        int lo                              # Lowest diagonal (inclusive)
        int hi                              # Highest diagonal (inclusive)
        # Wavefront elements                #
        wf_offset_t* offsets                # Offsets (k-centered)
        wf_offset_t* offsets_mem            # Offsets base memory (Internal)
        # Piggyback backtrace               #
        int bt_occupancy_max                # Maximum number of pcigar-ops stored on the Backtrace-block
        pcigar_t* bt_pcigar                 # Backtrace-block pcigar (k-centered)
        bt_block_idx_t* bt_prev             # Backtrace-block previous-index (k-centered)
        pcigar_t* bt_pcigar_mem             # Backtrace-block (base memory - Internal)
        bt_block_idx_t* bt_prev_mem         # Backtrace-block previous-index (base memory - Internal)
        # Slab internals                    #
        wavefront_status_type status        # Wavefront status (memory state)
        int wf_elements_allocated           # Total wf-elements allocated (max. wf. size)
        int wf_elements_allocated_min       # Minimum diagonal-element wf-element allocated
        int wf_elements_allocated_max       # Maximum diagonal-element wf-element allocated
        int wf_elements_init_min            # Minimum diagonal-element initialized (inclusive)
        int wf_elements_init_max            # Maximum diagonal-element initialized (inclusive)

    # Wavefront Set
    ctypedef struct wavefront_set_t:
        # In Wavefronts
        wavefront_t* in_mwavefront_misms
        wavefront_t* in_mwavefront_open1
        wavefront_t* in_mwavefront_open2
        wavefront_t* in_i1wavefront_ext
        wavefront_t* in_i2wavefront_ext
        wavefront_t* in_d1wavefront_ext
        wavefront_t* in_d2wavefront_ext
        # Out Wavefronts
        wavefront_t* out_mwavefront
        wavefront_t* out_i1wavefront
        wavefront_t* out_i2wavefront
        wavefront_t* out_d1wavefront
        wavefront_t* out_d2wavefront

    # Setup
    void wavefront_allocate(
        wavefront_t* const wavefront,
        const int wf_elements_allocated,
        const bint allocate_backtrace,
        mm_allocator_t* const mm_allocator)
    void wavefront_resize(
        wavefront_t* const wavefront,
        const int wf_elements_allocated,
        mm_allocator_t* const mm_allocator)
    void wavefront_free(
        wavefront_t* const wavefront,
        mm_allocator_t* const mm_allocator)

    # Initialization
    void wavefront_init(
        wavefront_t* const wavefront,
        const int min_lo,
        const int max_hi)
    void wavefront_init_null(
        wavefront_t* const wavefront,
        const int min_lo,
        const int max_hi)
    void wavefront_init_victim(
        wavefront_t* const wavefront,
        const int min_lo,
        const int max_hi)

    # Accessors
    void wavefront_set_limits(
        wavefront_t* const wavefront,
        const int lo,
        const int hi)

    # Utils
    uint64_t wavefront_get_size(
        wavefront_t* const wavefront)


cdef extern from "WFA2_lib/wavefront/wavefront_components.h" nogil:
    # Wavefront Components
    ctypedef struct wavefront_components_t:
        # Configuration
        bint memory_modular                         # Memory strategy (modular wavefronts)
        bint bt_piggyback                           # Backtrace Piggyback
        # Wavefronts dimensions                     #
        int num_wavefronts                          # Total number of allocated wavefronts
        int max_score_scope                         # Maximum score-difference between dependent wavefronts
        int historic_max_hi                         # Maximum WF hi-limit seen during current alignment
        int historic_min_lo                         # Minimum WF lo-limit seen during current alignment
        # Wavefronts                                #
        wavefront_t** mwavefronts                   # M-wavefronts
        wavefront_t** i1wavefronts                  # I1-wavefronts
        wavefront_t** i2wavefronts                  # I2-wavefronts
        wavefront_t** d1wavefronts                  # D1-wavefronts
        wavefront_t** d2wavefronts                  # D2-wavefronts
        wavefront_t* wavefront_null                 # Null wavefront (orthogonal reading)
        wavefront_t* wavefront_victim               # Dummy wavefront (orthogonal writing)
        # BT-Buffer                                 #
        wf_backtrace_buffer_t* bt_buffer            # Backtrace Buffer
        # MM                                        #
        mm_allocator_t* mm_allocator                # MM-Allocator

    # Setup
    void wavefront_components_allocate(
        wavefront_components_t* const wf_components,
        const int max_pattern_length,
        const int max_text_length,
        wavefront_penalties_t* const penalties,
        const bint memory_modular,
        const bint bt_piggyback,
        mm_allocator_t* const mm_allocator)
    void wavefront_components_reap(
        wavefront_components_t* const wf_components)
    void wavefront_components_clear(
        wavefront_components_t* const wf_components)
    void wavefront_components_free(
        wavefront_components_t* const wf_components)

    # Resize
    void wavefront_components_resize(
        wavefront_components_t* const wf_components,
        const int max_pattern_length,
        const int max_text_length,
        wavefront_penalties_t* const penalties)
    void wavefront_components_resize_null__victim(
        wavefront_components_t* const wf_components,
        const int lo,
        const int hi)

    # Compact
    void wavefront_components_compact_bt_buffer(
        wavefront_components_t* const wf_components,
        const int score,
        const int verbose)


cdef extern from "WFA2_lib/wavefront/wavefront_bialigner.h" nogil:
    ctypedef struct wf_bialign_breakpoint_t:
        # Scores
        int score                      # Score total
        int score_forward              # Score (forward)
        int score_reverse              # Score (reverse)
        # Location                     #
        int k_forward                  # Breakpoint diagonal (forward)
        int k_reverse                  # Breakpoint diagonal (reverse)
        wf_offset_t offset_forward     # Offset (forward)
        wf_offset_t offset_reverse     # Offset (reverse)
        affine2p_matrix_type component # Component (M/I/D)

    ctypedef struct wavefront_bialigner_t:
        # Wavefronts
        wavefront_aligner_t* wf_forward  # Breakpoint Forward aligner
        wavefront_aligner_t* wf_reverse  # Breakpoint Reverse aligner
        wavefront_aligner_t* wf_base     # Base/Subsidiary aligner
        # Operators
        void (*wf_align_compute)(wavefront_aligner_t* const,const int)

    # Setup
    wavefront_bialigner_t* wavefront_bialigner_new(
        wavefront_aligner_attr_t* const attributes,
        wavefront_plot_t* const plot)
    void wavefront_bialigner_reap(
        wavefront_bialigner_t* const wf_bialigner)
    void wavefront_bialigner_delete(
        wavefront_bialigner_t* const wf_bialigner)

    # Sequences
    void wavefront_bialigner_set_sequences_ascii(
        wavefront_bialigner_t* const wf_bialigner,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length)
    void wavefront_bialigner_set_sequences_lambda(
        wavefront_bialigner_t* const wf_bialigner,
        alignment_match_funct_t match_funct,
        void* match_funct_arguments,
        const int pattern_length,
        const int text_length)
    void wavefront_bialigner_set_sequences_packed2bits(
        wavefront_bialigner_t* const wf_bialigner,
        const uint8_t* const pattern,
        const int pattern_length,
        const uint8_t* const text,
        const int text_length)
    void wavefront_bialigner_set_sequences_bounds(
        wavefront_bialigner_t* const wf_bialigner,
        const int pattern_begin,
        const int pattern_end,
        const int text_begin,
        const int text_end)

    # Accessors
    uint64_t wavefront_bialigner_get_size(
        wavefront_bialigner_t* const wf_bialigner)
    void wavefront_bialigner_set_heuristic(
        wavefront_bialigner_t* const wf_bialigner,
        wavefront_heuristic_t* const heuristic)
    void wavefront_bialigner_set_max_alignment_steps(
        wavefront_bialigner_t* const wf_bialigner,
        const int max_alignment_steps)
    void wavefront_bialigner_set_max_memory(
        wavefront_bialigner_t* const wf_bialigner,
        const uint64_t max_memory_resident,
        const uint64_t max_memory_abort)
    void wavefront_bialigner_set_max_num_threads(
        wavefront_bialigner_t* const wf_bialigner,
        const int max_num_threads)
    void wavefront_bialigner_set_min_offsets_per_thread(
        wavefront_bialigner_t* const wf_bialigner,
        const int min_offsets_per_thread)


cdef extern from "WFA2_lib/wavefront/wavefront_slab.h" nogil:
    ctypedef enum wf_slab_mode_t:
        wf_slab_reuse = 1 #  Keep all wavefronts (Reap only by demand)
        wf_slab_tight = 2 #  Reap all if wavefronts are resized
    ctypedef struct wavefront_slab_t:
        # Attributes
        bint allocate_backtrace         # WFs require BT-vector
        wf_slab_mode_t slab_mode        # Slab strategy
        # Wavefront Slabs               #
        int init_wf_length              # Initial wf-elements allocated
        int current_wf_length           # Current wf-elements allocated
        vector_t* wavefronts            # All wavefronts (wavefront_t*)
        vector_t* wavefronts_free       # Free wavefronts (wavefront_t*)
        # Stats                         #
        uint64_t memory_used            # Memory used (Bytes)
        # MM                            #
        mm_allocator_t* mm_allocator    # MM-Allocator

    # Setup
    wavefront_slab_t* wavefront_slab_new(
        const int init_wf_length,
        const bint allocate_backtrace,
        const wf_slab_mode_t slab_mode,
        mm_allocator_t* const mm_allocator)
    void wavefront_slab_reap(
        wavefront_slab_t* const wavefront_slab)
    void wavefront_slab_clear(
        wavefront_slab_t* const wavefront_slab)
    void wavefront_slab_delete(
        wavefront_slab_t* const wavefront_slab)

    # Accessors
    void wavefront_slab_set_mode(
        wavefront_slab_t* const wavefront_slab,
        const wf_slab_mode_t slab_mode)

    # Allocator
    wavefront_t* wavefront_slab_allocate(
        wavefront_slab_t* const wavefront_slab,
        const int min_lo,
        const int max_hi)
    void wavefront_slab_free(
        wavefront_slab_t* const wavefront_slab,
        wavefront_t* const wavefront)

    # Utils
    uint64_t wavefront_slab_get_size(
        wavefront_slab_t* const wavefront_slab)


cdef extern from "WFA2_lib/wavefront/wfa.h" nogil:
    # Error codes & messages
    # [OK]
    DEF WF_STATUS_ALG_COMPLETED     = 0    # Success (Complete alignment found)
    DEF WF_STATUS_ALG_PARTIAL       = 1    # Success (Partial alignment found)
    # [FAIL]
    DEF WF_STATUS_MAX_STEPS_REACHED = -100 # Maximum number of WFA-steps reached
    DEF WF_STATUS_OOM               = -200 # Maximum memory limit reached
    DEF WF_STATUS_UNATTAINABLE      = -300 # Alignment unattainable under configured heuristics
    # [INTERNAL]
    DEF WF_STATUS_OK                = -1   # Computing alignment (in progress)
    DEF WF_STATUS_END_REACHED       = -2   # Alignment end reached
    DEF WF_STATUS_END_UNREACHABLE   = -3   # Alignment end unreachable under current configuration (e.g. Z-drop)

    # error messages
    char* wavefront_align_strerror(const int wf_error_code)
    char* wavefront_align_strerror_short (const int error_code)

    # Alignment status
    ctypedef struct _wavefront_aligner_t
    # ctypedef _wavefront_aligner_t wavefront_aligner_t
    ctypedef struct wavefront_align_status_t:
        # Status
        int status                                                     # Status code
        int score                                                      # Current WF-alignment score
        bint dropped                                                   # Heuristically dropped
        int num_null_steps                                             # Total contiguous null-steps performed
        uint64_t memory_used                                           # Total memory used
        # Wavefront alignment functions
        void (*wf_align_compute)(wavefront_aligner_t* const,const int) # WF Compute function
        bint (*wf_align_extend)(wavefront_aligner_t* const,const int)  # WF Extend function

    # Alignment type
    ctypedef enum wavefront_align_mode_t:
        wf_align_regular                  = 0
        wf_align_biwfa                    = 1
        wf_align_biwfa_breakpoint_forward = 2
        wf_align_biwfa_breakpoint_reverse = 3
        wf_align_biwfa_subsidiary         = 4

    # Wavefront Aligner
    ctypedef struct _wavefront_aligner_t:
        # Mode and Status
        wavefront_align_mode_t align_mode       # WFA alignment mode
        char* align_mode_tag                    # WFA mode tag
        wavefront_align_status_t align_status   # Current alignment status
        # Sequences
        wavefront_sequences_t* sequences        # Input sequences
        # Alignment Attributes
        alignment_scope_t alignment_scope       # Alignment scope (score only or full-CIGAR)
        alignment_form_t alignment_form         # Alignment form (end-to-end/ends-free)
        wavefront_penalties_t penalties         # Alignment penalties
        wavefront_heuristic_t heuristic         # Heuristic's parameters
        wavefront_memory_t memory_mode          # Wavefront memory strategy (modular wavefronts and piggyback)
        # Wavefront components
        wavefront_components_t wf_components    # Wavefront components
        affine2p_matrix_type component_begin    # Alignment begin component
        affine2p_matrix_type component_end      # Alignment end component
        wavefront_pos_t alignment_end_pos       # alignment end position
        # Bidirectional alignment
        wavefront_bialigner_t* bialigner        # BiWFA aligner
        # CIGAR
        cigar_t* cigar                          # Alignment CIGAR
        # MM
        bint mm_allocator_own                   # Ownership of MM-Allocator
        mm_allocator_t* mm_allocator            # MM-Allocator
        wavefront_slab_t* wavefront_slab        # MM-Wavefront-Slab (Allocates/Reuses the individual wavefronts)
        # Display
        wavefront_plot_t plot                   # Wavefront plot
        # System
        alignment_system_t system               # System related parameters

    # Setup
    wavefront_aligner_t* wavefront_aligner_new(
        wavefront_aligner_attr_t* attributes)
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
    void wavefront_aligner_set_alignment_extension(
        wavefront_aligner_t* const wf_aligner)

    # Heuristic configuration
    void wavefront_aligner_set_heuristic_none(
        wavefront_aligner_t* const wf_aligner)
    void wavefront_aligner_set_heuristic_wfadaptive(
        wavefront_aligner_t* const wf_aligner,
        const int min_wavefront_length,
        const int max_distance_threshold,
        const int score_steps)
    void wavefront_aligner_set_heuristic_wfmash(
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
    void wavefront_aligner_set_heuristic_banded_static(
        wavefront_aligner_t* const wf_aligner,
        const int band_min_k,
        const int band_max_k)
    void wavefront_aligner_set_heuristic_banded_adaptive(
        wavefront_aligner_t* const wf_aligner,
        const int band_min_k,
        const int band_max_k,
        const int score_steps)

    # System configuration
    void wavefront_aligner_set_max_alignment_steps(
        wavefront_aligner_t* const wf_aligner,
        const int max_alignment_steps)
    void wavefront_aligner_set_max_memory(
        wavefront_aligner_t* const wf_aligner,
        const uint64_t max_memory_resident,
        const uint64_t max_memory_abort)
    void wavefront_aligner_set_max_num_threads(
        wavefront_aligner_t* const wf_aligner,
        const int max_num_threads)
    void wavefront_aligner_set_min_offsets_per_thread(
        wavefront_aligner_t* const wf_aligner,
        const int min_offsets_per_thread)

    # Wavefront Align
    int wavefront_align(
        wavefront_aligner_t* const wf_aligner,
        const char* const pattern,
        const int pattern_length,
        const char* const text,
        const int text_length)
    int wavefront_align_lambda(
        wavefront_aligner_t* const wf_aligner,
        const alignment_match_funct_t match_funct,
        void* match_funct_arguments,
        const int pattern_length,
        const int text_length)
    int wavefront_align_packed2bits(
        wavefront_aligner_t* const wf_aligner,
        const uint8_t* const pattern,
        const int pattern_length,
        const uint8_t* const text,
        const int text_length)

