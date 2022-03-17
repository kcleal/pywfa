=====
pyWFA
=====

A python wrapper for wavefront alignment using WFA2-lib

Installation
------------

To build from source::

    git clone --recursive https://github.com/kcleal/pywfa
    cd pywfa
    python setup.py install

Overview
--------

Alignment of pattern and text strings can be performed directly by accessing WFA2-lib functions:

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


Cigartuples follow the pysam convention:

.. list-table::
   :widths: 10 10
   :header-rows: 1

   * - Opperation
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

For convenience, a results object can be obtained by calling with a pattern and text:

.. code-block:: python

    pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
    text =    "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
    a = WavefrontAligner(pattern)
    a(text)
    >>> {'pattern_length': 32, 'text_length': 32, 'pattern_start': 0, 'pattern_end': 32, 'text_start': 0, 'text_end': 32, 'cigartuples': [(0, 3), (8, 1), (0, 4), (2, 1), (0, 7), (1, 1), (0, 9), (8, 1), (0, 6)], 'score': -24, 'pattern': 'TCTTTACTCGCGCGTTGGAGAAATACAATAGT', 'text': 'TCTATACTGCGCGTTTGGAGAAATAAAATAGT', 'status': 0}

    # Alignment can also be called with a pattern like this:
    a(text, pattern)


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
    >>> 8S7M5D6M


Additionally, short matches at the end of alignments can be trimmed off:

.. code-block:: python

    pattern = "AAAAAAAAAAAACCTTTTAAAAAAGAAAAAAA"
    text = "ACCCCCCCCCCCAAAAACCAAAAAAAAAAAAA"
    a = WavefrontAligner(pattern)

    res = a(text, clip_cigar=False)
    res.cigartuples
    >>> [(0, 1), (1, 5), (8, 6), (0, 7), (2, 5), (0, 5), (8, 1), (0, 7)]
    res.aligned_text
    >>> ACCCCCCCCCCCAAAAACCAAAAAAAAAAAAA

    # By default the minimum allowed block of matches at each end is 5 bp
    res = a(text, clip_cigar=True)
    res.cigartuples
    >>> [(4, 12), (0, 7), (2, 5), (0, 5), (8, 1), (0, 7)]
    res.aligned_text
    >>> AAAAACCAAAAAAAAAAAAA

    # The minimum length of aligned bases at the end can be controlled using:
    res = a(text, clip_cigar=True, min_aligned_bases_left=1, min_aligned_bases_right=1)
    res.cigartuples
    >>> [(0, 1), (1, 5), (8, 6), (0, 7), (2, 5), (0, 5), (8, 1), (0, 7)]
    res.aligned_text
    >>> ACCCCCCCCCCCAAAAACCAAAAAAAAAAAAA

    # Mismatch operations X can also be elided, note this occurs after the clip_cigar stage
    res = a(text, clip_cigar=True, elide_mismatches=True)
    res.cigartuples
    >>> [(4, 12), (0, 7), (2, 5), (0, 13)]
    res.aligned_text
    >>> AAAAACCAAAAAAAAAAAAA

Notes: The alignment score is not modified currently by trimming the cigar, however the pattern_start, pattern_end,
test_start and text_end are modfied when the cigar is modified.


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

