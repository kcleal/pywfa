=====
pyWFA
=====

A python wrapper for wavefront alignment using WFA2-lib

Installation
------------

To build from source::

    git clone https://github.com/kcleal/pywfa
    cd pywfa/WFA2-lib && make && cd ..
    python setup.py install

Overview
--------

Alignment of pattern and text strings can be performed by directly accessing WFA2-lib functions:

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
   * - 3
     - >
   * - S
     - 4
   * - H
     - 5
   * - ?
     - 7
   * - X
     - 8
   * - ?
     - 9

For convenience, a results object can be obtained by calling with a pattern and text:


.. code-block:: python

    from pywfa.align import WavefrontAligner

    pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
    text =    "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
    a = WavefrontAligner(pattern)
    a(text)
    >>> {'pattern_length': 32, 'text_length': 32, 'pattern_start': 0, 'pattern_end': 32, 'text_start': 0, 'text_end': 32, 'cigartuples': [(0, 8), (2, 1), (0, 7), (1, 1), (0, 16)], 'score': -24, 'pattern': 'TCTTTACTCGCGCGTTGGAGAAATACAATAGT', 'text': 'TCTATACTGCGCGTTTGGAGAAATAAAATAGT', 'status': 0}

    # Alignment can also be called with a pattern like this:
    a(text, pattern)


Supplying a pattern during initialization will cache the pattern


To configure the `WaveFrontAligner`, options can be provided during initialization:


.. code-block:: python

    from pywfa.align import WavefrontAligner

    a = WavefrontAligner(scope="score",
                         distance="affine2p",
                         span="end-to-end",
                         heuristic="adaptive")

Supported distance metrics are "affine" (default) and "affine2p". Scope can be "full" (default)
or "score". Span can be "ends-free" (default) or "end-to-end". Heuristic can be None (default),
"adaptive" or "X-drop"
