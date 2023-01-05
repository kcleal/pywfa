=====
pyWFA
=====

A python wrapper for wavefront alignment using `WFA2-lib
<https://github.com/smarco/WFA2-lib/>`_

Installation
------------

To download from pypi::

    pip install pywfa

Build from source::

    git clone https://github.com/kcleal/pywfa
    cd pywfa
    pip install .

Overview
--------

Alignment of pattern and text strings can be performed by accessing WFA2-lib functions directly:

.. code-block:: python

    from pywfa.align import WavefrontAligner

    pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
    text =    "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
    a = WavefrontAligner(pattern)
    score = a.wavefront_align(text)
    assert a.status == 0  # alignment was successful
    assert a.cigarstring == "3M1X4M1D7M1I9M1X6M"
    assert a.score == -24
    a.cigartuples
    >>> [(0, 3), (8, 1), (0, 4), (2, 1), (0, 7), (1, 1), (0, 9), (8, 1), (0, 6)]
    a.cigar_print_pretty()

.. code-block:: text

    >>> 3M1X4M1D7M1I9M1X6M      ALIGNMENT
        1X1D1I1X      ALIGNMENT.COMPACT
        PATTERN    TCTTTACTCGCGCGTT-GGAGAAATACAATAGT
                   ||| |||| ||||||| ||||||||| ||||||
        TEXT       TCTATACT-GCGCGTTTGGAGAAATAAAATAGT

The output of cigar_pretty_print can be directed to a file, rather than stdout using:

.. code-block:: python

    a.cigar_print_pretty("file.txt")

To obtain a python str of this print out, access the results object (see below).

Cigartuples follow the convention:

.. list-table::
   :widths: 15 15
   :header-rows: 1

   * - Operation
     - Code
   * - M
     - 0
   * - I
     - 1
   * - D
     - 2
   * - N
     - 3
   * - S
     - 4
   * - H
     - 5
   * - =
     - 7
   * - X
     - 8
   * - B
     - 9

For convenience, a results object can be obtained by calling the `WavefrontAligner` with a pattern and text:

.. code-block:: python

    pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
    text =    "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
    a = WavefrontAligner(pattern)
    a(text)
    >>> {'pattern_length': 32, 'text_length': 32, 'pattern_start': 0, 'pattern_end': 32, 'text_start': 0, 'text_end': 32, 'cigartuples': [(0, 3), (8, 1), (0, 4), (2, 1), (0, 7), (1, 1), (0, 9), (8, 1), (0, 6)], 'score': -24, 'pattern': 'TCTTTACTCGCGCGTTGGAGAAATACAATAGT', 'text': 'TCTATACTGCGCGTTTGGAGAAATAAAATAGT', 'status': 0}

    # Alignment can also be called with a pattern like this:
    a(text, pattern)

    # obtain a string in the same format as cigar_print_pretty
    a.pretty
    >>> 3M1X4M1D7M1I9M1X6M      ALIGNMENT
        1X1D1I1X      ALIGNMENT.COMPACT
              PATTERN    TCTTTACTCGCGCGTT-GGAGAAATACAATAGT
                         |||*|||| ||||||| |||||||||*||||||
              TEXT       TCTATACT-GCGCGTTTGGAGAAATAAAATAGT


Configure
---------
To configure the `WaveFrontAligner`, options can be provided during initialization:


.. code-block:: python

    from pywfa.align import WavefrontAligner

    a = WavefrontAligner(scope="score",
                         distance="affine2p",
                         span="end-to-end",
                         heuristic="adaptive")

Supported distance metrics are "affine" (default) and "affine2p". Scope can be "full" (default)
or "score". Span can be "ends-free" (default) or "end-to-end". Heuristic can be None (default),
"adaptive" or "X-drop".

When using heuristic functions it is recommended to check the status attribute:


.. code-block:: python

    pattern = "AAAAACCTTTTTAAAAAA"
    text = "GGCCAAAAACCAAAAAA"
    a = WavefrontAligner(heuristic="adaptive")
    a(pattern, text)
    a.status
    >>> 0   # successful alignment, -1 indicates the alignment was stopped due to the heuristic


Default options
---------------

The `WavefrontAligner` will be initialized with the following default options:

.. list-table::
   :widths: 15 10
   :header-rows: 1

   * - Parameter
     - Default value
   * - pattern
     - None
   * - distance
     - "affine"
   * - match
     - 0
   * - gap_opening
     - 6
   * - gep_extension
     - 2
   * - gap_opening2
     - 24
   * - gap_extension2
     - 1
   * - scope
     - "full"
   * - span
     - "ends-free"
   * - pattern_begin_free
     - 0
   * - pattern_end_free
     - 0
   * - text_begin_free
     - 0
   * - text_end_free
     - 0
   * - heuristic
     - None
   * - min_wavefront_length
     - 10
   * - max_distance_threshold
     - 50
   * - steps_between_cutoffs
     - 1
   * - xdrop
     - 20


Modifying the cigar
-------------------

If desired the cigar can be modified so the end operation is either a soft-clip or a match, this makes the
alignment cigar resemble those produced by bwa, for example:

.. code-block:: python

    pattern = "AAAAACCTTTTTAAAAAA"
    text = "GGCCAAAAACCAAAAAA"
    a = WavefrontAligner(pattern)

    res = a(text, clip_cigar=False)
    print(cigartuples_to_str(res.cigartuples))
    >>> 4I7M5D6M

    res(text, clip_cigar=True)
    print(cigartuples_to_str(res.cigartuples))
    >>> 4S7M5D6M


An experimental feature is to trim short matches at the end of alignments. This results in alignments that approximate local alignments:

.. code-block:: python

    pattern = "AAAAAAAAAAAACCTTTTAAAAAAGAAAAAAA"
    text = "ACCCCCCCCCCCAAAAACCAAAAAAAAAAAAA"
    a = WavefrontAligner(pattern)

    # The unmodified cigar may have short matches at the end:
    res = a(text, clip_cigar=False)
    res.cigartuples
    >>> [(0, 1), (1, 5), (8, 6), (0, 7), (2, 5), (0, 5), (8, 1), (0, 7)]
    res.aligned_text
    >>> ACCCCCCCCCCCAAAAACCAAAAAAAAAAAAA
    res.text_start, res.text_end
    >>> 0, 32

    # The minimum allowed block of matches can be set at e.g. 5 bp, which will trim off short matches
    res = a(text, clip_cigar=True, min_aligned_bases_left=5, min_aligned_bases_right=5)
    res.cigartuples
    >>> [(4, 12), (0, 7), (2, 5), (0, 5), (8, 1), (0, 7)]
    res.aligned_text
    >>> AAAAACCAAAAAAAAAAAAA
    res.text_start, res.text_end
    >>> 12, 32

    # Mismatch operations X can also be elided, note this occurs after the clip_cigar stage
    res = a(text, clip_cigar=True, min_aligned_bases_left=5, min_aligned_bases_right=5, elide_mismatches=True)
    res.cigartuples
    >>> [(4, 12), (0, 7), (2, 5), (0, 13)]
    res.aligned_text
    >>> AAAAACCAAAAAAAAAAAAA

Notes: The alignment score is not modified currently by trimming the cigar, however the pattern_start, pattern_end,
test_start and text_end are modified when the cigar is modified.
