import pyximport
pyximport.install(reload_support=True)

from pywfa.align import WavefrontAligner, clip_cigartuples, cigartuples_to_str, elide_mismatches_from_cigar

