=====
pyWFA
=====

A python wrapper for wavefront alignment using `WFA2-lib
<https://github.com/smarco/WFA2-lib/>`_

Installation
------------

To download from pypi::

    pip install pywfa

From conda::

    conda install -c bioconda pywfa

Build from source::

    git clone https://github.com/kcleal/pywfa
    cd pywfa
    pip install .

Overview
--------

Alignment of pattern and text strings can be performed by accessing WFA2-lib functions directly:

.. code-block::

    from pywfa import WavefrontAligner

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

.. code-block::

    >>> ALIGNMENT	3M1X4M1D7M1I9M1X6M
      ALIGNMENT.COMPACT	1X1D1I1X
      PATTERN    TCTTTACTCGCGCGTT-GGAGAAATACAATAGT
                 ||| |||| ||||||| ||||||||| ||||||
      TEXT       TCTATACT-GCGCGTTTGGAGAAATAAAATAGT

The output of cigar_pretty_print can be directed to a file, rather than stdout using:

.. code-block::

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

.. code-block::

    pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
    text =    "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
    a = WavefrontAligner(pattern)
    result = a(text)  # alignment result
    result.__dict__
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


.. code-block::

    from pywfa import WavefrontAligner

    a = WavefrontAligner(scope="score",
                         distance="affine2p",
                         span="end-to-end",
                         heuristic="adaptive")

Supported distance metrics are "affine" (default) and "affine2p". Scope can be "full" (default)
or "score". Span can be "ends-free" (default) or "end-to-end". Heuristics can be None (default),
"adaptive", "X-drop", "Z-drop", "banded-static", "banded-adaptive" or "wfmash". The memory
footprint can be selected with memory_mode: "high" (default), "medium", "low" or "biwfa"
(bidirectional WFA, lowest memory). The maximum number of alignment steps and the memory
budget (in bytes) can be capped with max_steps and max_memory, respectively.

When using heuristic functions it is recommended to check the status attribute:


.. code-block::

    pattern = "AAAAACCTTTTTAAAAAA"
    text = "GGCCAAAAACCAAAAAA"
    a = WavefrontAligner(heuristic="adaptive")
    a(pattern, text)
    a.status
    >>> 0   # successful alignment, -1 indicates the alignment was stopped due to the heuristic
    a.status_message
    >>> [WFA] Alignment completed successfully

To inspect which WFA variant performed the alignment (e.g. when using BiWFA), check the
align_mode attribute:

.. code-block::

    a = WavefrontAligner(pattern, memory_mode="biwfa")
    a(text)
    a.align_mode
    >>> 'biwfa'


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
   * - gap_extension
     - 2
   * - gap_opening2
     - 24
   * - gap_extension2
     - 1
   * - scope
     - "full"
   * - span
     - "ends-free"
   * - memory_mode
     - "high"
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
   * - zdrop
     - 20
   * - min_k
     - -10
   * - max_k
     - 10
   * - wildcard
     - None
   * - max_steps
     - 0
   * - max_memory
     - 0


Modifying the cigar
-------------------

If desired the cigar can be modified so the end operation is either a soft-clip or a match, this makes the
alignment cigar resemble those produced by bwa, for example:

.. code-block::

    pattern = "AAAAACCTTTTTAAAAAA"
    text = "GGCCAAAAACCAAAAAA"
    a = WavefrontAligner(pattern)

    res = a(text, clip_cigar=False)
    print(cigartuples_to_str(res.cigartuples))
    >>> 4I7M5D6M

    res = a(text, clip_cigar=True)
    print(cigartuples_to_str(res.cigartuples))
    >>> 4S7M5D6M
