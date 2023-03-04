# Copyright (C) 2014-2020. Ben Pruitt & Nick Conway; Wyss Institute
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


def suite():
    '''Define the test suite'''
    suite = unittest.TestSuite()
    suite.addTest(TestP3Helpers())
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    test_suite = suite()
    runner.run(test_suite)
