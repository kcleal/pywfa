
import unittest
import os
import pysam
import time

from pywfa.align import WavefrontAligner, clip_cigartuples, cigartuples_to_str, elide_mismatches_from_cigar


root = os.path.abspath(os.path.dirname(__file__))


class TestConstruct(unittest.TestCase):
    """ Test construction and alignment"""

    def test_affine(self):
        print("Affine")
        pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
        text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
        a = WavefrontAligner(pattern)
        score = a.wavefront_align(text)
        assert a.status == 0
        assert a.cigarstring == "3M1X4M1D7M1I9M1X6M"
        assert a.score == -24
        assert a.score == score
        a.cigar_print_pretty()

        pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
        text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
        a = WavefrontAligner(pattern)
        a(text)

        assert a.status == 0
        assert a.cigarstring == "3M1X4M1D7M1I9M1X6M"
        assert a.score == -24
        assert a.score == score

        print("Affine2")
        pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
        text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
        a = WavefrontAligner()
        a(text, pattern, clip_cigar=False)

        assert a.status == 0
        assert a.cigarstring == "3M1X4M1D7M1I9M1X6M"
        assert a.score == -24
        print("Affine3")
        pattern = "TCTATACTGCGCGTTTGGAGAAATAAAA"
        text = "TCTCCCCATACTGCGCGTTTGGAGAAATAAAA"
        a = WavefrontAligner()
        a(text, pattern, clip_cigar=False)


    def test_scope(self):
        print("Scope")
        pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
        text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
        a = WavefrontAligner(pattern, scope="score")
        a(text)
        assert a.status == 0
        assert a.cigarstring == ""
        assert a.score == -24
        del a

    def test_supress_seqs(self):
        print("Supress seqs")
        pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
        text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
        a = WavefrontAligner(pattern, scope="score")
        res = a(text, supress_sequences=True)
        assert res.aligned_pattern is None and res.aligned_text is None
        assert a.status == 0
        assert a.cigarstring == ""
        assert a.score == -24

        pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
        text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
        a = WavefrontAligner(pattern, scope="full")
        res = a(text, supress_sequences=True)
        assert res.aligned_pattern is None and res.aligned_text is None
        assert a.status == 0
        assert a.cigarstring == "3M1X4M1D7M1I9M1X6M"
        assert a.score == -24

    def test_many(self):
        print("Construct many")
        pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT"
        text = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT"
        for i in range(1000):
            a = WavefrontAligner(pattern)
            a.wavefront_align(text)
            assert a.score == -24

    def test_end_to_end(self):
        print("End to end")
        pattern = "AATTAATTTAAGTCTAGGCTACTTTCGGTACTTTGTTCTT"
        text = "AATTTAAGTCTAGGCTACTTTCGGTACTTTCTT"
        a = WavefrontAligner(pattern, span="end-to-end", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert a.cigarstring == "4M4D26M3D3M"
        assert res.score == -26
        del a

    def test_ends_free(self):
        print("Ends free")
        pattern = "AATTAATTTAAGTCTAGGCTACTTTCGGTACTTTGTTCTT"
        text = "AATTTAAGTCTAGGCTACTTTCGGTACTTTCTT"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text, clip_cigar=True, elide_mismatches=True, min_aligned_bases_left=5, min_aligned_bases_right=5)
        assert res.aligned_pattern == res.aligned_text
        assert a.cigarstring == "4M4D26M3D3M"
        assert res.score == -26
        del a

    def test_ends_free2(self):
        print("Ends free 2")
        pattern = "AAAAACCTTTTTAAAAAA"
        text = "GGCCAAAAACCAAAAAA"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.text_start == 4 and res.text_end == 17
        del a

        pattern = "AAAAACCTTTTTAAAAAA"
        text = "GGCCAAAAACCGGGGGGG"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text
        assert res.text_start == 4 and res.text_end == 11
        del a

        pattern = "AAAAACCGGGG"
        text = "AAAAACC"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text

        pattern = "AAAAACC"
        text = "AAAAACCGGGG"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text

        pattern = "GGGGAAAAACC"
        text = "AAAAACCGGGG"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text

        pattern = "AAAAACCGGGG"
        text = "GGGGAAAAACC"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text

        pattern = "GGGGAAAAACC"
        text = "AAAAACC"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text

        pattern = "GGGGAAAAACC"
        text = "CCCCCAAAAACC"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text

        pattern = "GGGGAAAAACCGGGGG"
        text = "CCCCCAAAAACCTTTTT"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text

        pattern = "AAAAACC"
        text = "CCCCCAAAAACCTTTTT"
        a = WavefrontAligner(pattern, span="ends-free", mismatch=4, gap_opening=6, gap_extension=2)
        res = a(text)
        assert res.aligned_pattern == res.aligned_text

    def test_heuristic(self):
        print("Heuristic")
        pattern = "AAAAACCAAAAAA"
        text = "GGCCAAAAACCAAAAAA"
        a = WavefrontAligner(pattern, distance="affine", mismatch=4, gap_opening=6, gap_extension=2,
                             heuristic='X-drop')
        res = a(text)
        if res.status == 0:
            assert res.aligned_pattern == res.aligned_text

        a = WavefrontAligner(pattern, distance="affine", mismatch=4, gap_opening=6, gap_extension=2,
                             heuristic='adaptive')
        res = a(text)
        if res.status == 0:
            assert res.aligned_pattern == res.aligned_text

    def test_long(self):
        print("Long")
        with pysam.FastxFile(root + '/long.fa') as reads, pysam.FastxFile(root + '/long.reference.fa') as refs:

            for r in reads:
                text = r.sequence.upper()
                pattern = next(refs).sequence.upper()
                l_text = int(len(text) / 2)
                l_pattern = int(len(pattern) / 2)
                a = WavefrontAligner(distance="affine", mismatch=4, gap_opening=6, gap_extension=2,
                                     pattern_begin_free=l_pattern,
                                     pattern_end_free=l_pattern,
                                     text_begin_free=l_text,
                                     text_end_free=l_text
                                     )

                a(text, pattern, clip_cigar=True)

    def test_short(self):
        print('Short')
        with pysam.FastxFile(root + '/short.fa') as reads, pysam.FastxFile(root + '/short.reference.fa') as refs:
            for idx, r in enumerate(reads):
                text = r.sequence.upper()
                pattern = next(refs).sequence.upper()
                a = WavefrontAligner(mismatch=5, gap_opening=6, gap_extension=2)
                res = a(text, pattern)

    def test_short2p(self):
        print('Short2p')
        with pysam.FastxFile(root + '/short.fa') as reads, pysam.FastxFile(root + '/short.reference.fa') as refs:
            for idx, r in enumerate(reads):
                text = r.sequence.upper()
                pattern = next(refs).sequence.upper()
                a = WavefrontAligner(distance="affine2p", mismatch=5, gap_opening=6, gap_extension=2)
                res = a(text, pattern, clip_cigar=True, elide_mismatches=True)
                if r.name == "read6.loci:chr1:13,853,852-13,854,838":
                    assert res.cigartuples[3] == (2, 175)


def main():
    unittest.main()


if __name__ == "__main__":
    unittest.main()
