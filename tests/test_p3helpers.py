# Copyright (C) 2014-2026. Ben Pruitt & Nick Conway; Wyss Institute
# See LICENSE for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
tests.test_p3helpers
~~~~~~~~~~~~~~~~~~~~~~~~~

Unit tests for the primer3-py p3helpers.

'''

import unittest

from primer3 import p3helpers  # type: ignore

from . import _leakcheck


class TestP3Helpers(unittest.TestCase):

    def test_sanitize_sequence(self):
        # 1. Test upper case
        seq = 'ACGTRYMKSWHDBV'
        out_seq = p3helpers.sanitize_sequence(seq)
        self.assertEqual(out_seq, 'ACGTNNNNNNNNNN')

        # 2. Test lower case
        seq = 'acgtrymkswhdbv'
        out_seq = p3helpers.sanitize_sequence(seq)
        self.assertEqual(out_seq, 'acgtnnnnnnnnnn')

        # 3. Test bytes upper case
        seq = b'ACGTRYMKSWHDBV'
        out_seq = p3helpers.sanitize_sequence_b(seq)
        self.assertEqual(out_seq, b'ACGTNNNNNNNNNN')

        # 4. Test bytes lower case
        seq = b'acgtrymkswhdbv'
        out_seq = p3helpers.sanitize_sequence_b(seq)
        self.assertEqual(out_seq, b'acgtnnnnnnnnnn')

        # 5. Test ValueError for bad base
        self.assertRaises(
            ValueError,
            p3helpers.sanitize_sequence,
            'ZACGT',
        )

    def test_reverse_complement(self):
        # 1. Test upper case no sanitize
        seq = 'ACGTRYMKSWHDBV'
        out_seq = p3helpers.reverse_complement(seq, do_sanitize=False)
        self.assertEqual(out_seq, 'BVHDWSMKRYACGT')

        # 2. Test upper case sanitize
        out_seq = p3helpers.reverse_complement(seq, do_sanitize=True)
        self.assertEqual(out_seq, 'NNNNNNNNNNACGT')

        # 3. Test lower case no sanitize
        seq = 'acgtrymkswhdbv'
        out_seq = p3helpers.reverse_complement(seq, do_sanitize=False)
        self.assertEqual(out_seq, 'bvhdwsmkryacgt')

        # 4. Test lower case sanitize
        out_seq = p3helpers.reverse_complement(seq, do_sanitize=True)
        self.assertEqual(out_seq, 'nnnnnnnnnnacgt')

        # 5. Test bytes upper case no sanitize
        seq = b'ACGTRYMKSWHDBV'
        out_seq = p3helpers.reverse_complement_b(seq, do_sanitize=False)
        self.assertEqual(out_seq, b'BVHDWSMKRYACGT')

        # 5. Test bytes upper case sanitize
        out_seq = p3helpers.reverse_complement_b(seq, do_sanitize=True)
        self.assertEqual(out_seq, b'NNNNNNNNNNACGT')

        # 6. Test bytes lower case no sanitize
        seq = b'acgtrymkswhdbv'
        out_seq = p3helpers.reverse_complement_b(seq, do_sanitize=False)
        self.assertEqual(out_seq, b'bvhdwsmkryacgt')

        # 7. Test bytes lower case sanitize
        out_seq = p3helpers.reverse_complement_b(seq, do_sanitize=True)
        self.assertEqual(out_seq, b'nnnnnnnnnnacgt')

        # 8. Test ValueError for bad base
        self.assertRaises(
            ValueError,
            p3helpers.reverse_complement,
            'ZACGT',
        )

    def test_reverse_complement_odd_middle_base(self):
        # The lone middle base of an odd-length sequence must be validated
        # like every other position. Regression: previously the middle base
        # complement was assigned without checking for INVALID_BASE, so an
        # invalid middle base was silently accepted and returned as '?'.
        for bad in ('Z', 'ACZGT', b'Z', b'ACZGT', b'\x80', b'AC\x80GT'):
            with self.assertRaises(ValueError):
                if isinstance(bad, bytes):
                    p3helpers.reverse_complement_b(bad)
                else:
                    p3helpers.reverse_complement(bad)

    def test_high_bytes_rejected_not_oob(self):
        # Bytes >= 0x80 index the lookup tables, which must have 256 entries.
        # A 128-entry table would read out of bounds here; instead these
        # inputs must be rejected cleanly (invalid base), not crash.
        for bad in (b'\x80\x80', b'\xff\xff', b'\xc3\xa9', b'A\xffCG'):
            with self.assertRaises(ValueError):
                p3helpers.reverse_complement_b(bad)
            with self.assertRaises(ValueError):
                p3helpers.sanitize_sequence_b(bad)
            with self.assertRaises(ValueError):
                p3helpers.ensure_acgt_uppercase_b(bad)

    def test_str_non_ascii_rejected(self):
        # str helpers must operate on the UTF-8 byte length, not the code-point
        # count, and must reject non-ASCII characters cleanly (a multi-byte
        # char encodes to bytes >= 0x80). Regression for the code-point/byte
        # length mismatch and in-place bytes mutation.
        for bad in ('ACGTé', 'éACGT', 'ACéGT'):
            with self.assertRaises(ValueError):
                p3helpers.reverse_complement(bad)
            with self.assertRaises(ValueError):
                p3helpers.sanitize_sequence(bad)
            with self.assertRaises(ValueError):
                p3helpers.ensure_acgt_uppercase(bad)

    def test_no_memory_leak(self):
        # Exercise every p3helpers entry point, including the error paths that
        # allocate a scratch buffer before raising, and assert no Python/PyMem
        # allocation growth (deterministic via tracemalloc).
        seq = 'ACGTacgtNRYKMSWB' * 4
        seqb = seq.encode()

        def swallow(fn, bad):
            # error paths allocate a scratch buffer before raising
            try:
                fn(bad)
            except ValueError:
                return
            raise AssertionError(f'{fn.__name__}({bad!r}) did not raise')

        def work():
            p3helpers.reverse_complement(seq)
            p3helpers.reverse_complement(seq, do_sanitize=True)
            p3helpers.reverse_complement_b(seqb)
            p3helpers.sanitize_sequence(seq)
            p3helpers.sanitize_sequence_b(seqb)
            p3helpers.ensure_acgt_uppercase('acgtACGT')
            p3helpers.ensure_acgt_uppercase_b(b'acgtACGT')
            swallow(p3helpers.reverse_complement, 'ZZZZ')
            swallow(p3helpers.sanitize_sequence, 'ZZZZ')
            swallow(p3helpers.ensure_acgt_uppercase, 'ACGTX')
            swallow(p3helpers.ensure_acgt_uppercase_b, b'ACGTX')

        # Non-leaking growth measures ~0.05 KiB; 64 KiB over 20k iters catches
        # even a few-byte-per-call scratch-buffer leak with a wide margin.
        _leakcheck.assert_no_leak(
            self, work, iters=20000, warmup=2000, max_tracemalloc_kib=64,
        )

    def test_ensure_acgt_uppercase(self):
        # 1. Test already uppercase ACGT
        seq = 'ACGT'
        out_seq = p3helpers.ensure_acgt_uppercase(seq)
        self.assertEqual(out_seq, 'ACGT')
        # Verify it's the same object (no new allocation)
        self.assertIs(out_seq, seq)

        # 2. Test lowercase to uppercase conversion
        seq = 'acgt'
        out_seq = p3helpers.ensure_acgt_uppercase(seq)
        self.assertEqual(out_seq, 'ACGT')
        # Verify it's a new object
        self.assertIsNot(out_seq, seq)

        # 3. Test mixed case
        seq = 'AcGt'
        out_seq = p3helpers.ensure_acgt_uppercase(seq)
        self.assertEqual(out_seq, 'ACGT')

        # 4. Test empty sequence
        seq = ''
        out_seq = p3helpers.ensure_acgt_uppercase(seq)
        self.assertEqual(out_seq, '')
        # Verify empty string returns same object
        self.assertIs(out_seq, seq)

        # 5. Test bytes already uppercase ACGT
        seq = b'ACGT'
        out_seq = p3helpers.ensure_acgt_uppercase_b(seq)
        self.assertEqual(out_seq, b'ACGT')
        # Verify it's the same object
        self.assertIs(out_seq, seq)

        # 6. Test bytes lowercase to uppercase
        seq = b'acgt'
        out_seq = p3helpers.ensure_acgt_uppercase_b(seq)
        self.assertEqual(out_seq, b'ACGT')
        # Verify it's a new object
        self.assertIsNot(out_seq, seq)

        # 7. Test bytes mixed case
        seq = b'AcGt'
        out_seq = p3helpers.ensure_acgt_uppercase_b(seq)
        self.assertEqual(out_seq, b'ACGT')

        # 8. Test bytes empty sequence
        seq = b''
        out_seq = p3helpers.ensure_acgt_uppercase_b(seq)
        self.assertEqual(out_seq, b'')
        # Verify empty bytes returns same object
        self.assertIs(out_seq, seq)

        # 9. Test invalid chars at different positions (string)
        invalid_positions = [
            'NACGT',  # Start
            'ACGTN',  # End
            'ACNGT',  # Middle
        ]
        for seq in invalid_positions:
            with self.assertRaises(ValueError) as cm:
                p3helpers.ensure_acgt_uppercase(seq)
            err_msg = str(cm.exception)
            self.assertIn('position', err_msg)
            self.assertIn(seq, err_msg)
            self.assertIn("'N'", err_msg)

        # 10. Test invalid chars at different positions (bytes)
        invalid_positions = [
            b'NACGT',  # Start
            b'ACGTN',  # End
            b'ACNGT',  # Middle
        ]
        for seq in invalid_positions:
            with self.assertRaises(ValueError) as cm:
                p3helpers.ensure_acgt_uppercase_b(seq)
            err_msg = str(cm.exception)
            self.assertIn('position', err_msg)
            self.assertIn(str(seq), err_msg)
            self.assertIn("'N'", err_msg)

        # 12. Test long sequence
        long_seq = 'ACGT' * 1000
        out_seq = p3helpers.ensure_acgt_uppercase(long_seq)
        self.assertIs(out_seq, long_seq)


def suite():
    '''Define the test suite'''
    suite = unittest.TestSuite()
    suite.addTest(TestP3Helpers())
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    test_suite = suite()
    runner.run(test_suite)
