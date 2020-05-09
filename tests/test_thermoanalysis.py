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
test_thermoanalysis
~~~~~~~~~~~~~~~~~~~

Unit tests for the primer3-py low level thermodynamic calculation bindings.

'''

from __future__ import print_function
import unittest
from time import sleep
import random
import sys

try:
    import resource
except: # For Windows compatibility
    resource = None

from primer3 import (
    bindings,
    wrappers
)


def _getMemUsage():
    """ Get current process memory usage in bytes """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


@unittest.skipIf(
    sys.platform=='win32',
    "Windows doesn't support resource module and wrappers")
class TestLowLevelBindings(unittest.TestCase):

    def randArgs(self):
        self.seq1 = ''.join([random.choice('ATGC') for _ in
                             range(random.randint(20, 59))])
        self.seq2 = ''.join([random.choice('ATGC') for _ in
                             range(random.randint(20, 59))])
        self.mv_conc = random.uniform(1, 200)
        self.dv_conc = random.uniform(0, 40)
        self.dntp_conc = random.uniform(0, 20)
        self.dna_conc = random.uniform(0, 200)
        self.temp_c = random.randint(10, 70)
        self.max_loop = random.randint(10, 30)

    def test_calcTm(self):
        for x in range(100):
            self.randArgs()
            binding_tm = bindings.calcTm(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc
            )
            wrapper_tm = wrappers.calcTm(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc
            )
            self.assertEqual(int(binding_tm), int(wrapper_tm))

    def test_calcHairpin(self):
        for _ in range(100):
            self.randArgs()
            binding_res = bindings.calcHairpin(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True
            )
            wrapper_res = wrappers.calcHairpin(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop
            )
            self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))
            self.assertEqual(
                binding_res.ascii_structure,
                wrapper_res.ascii_structure
            )

    def test_calcHomodimer(self):
        for _ in range(100):
            self.randArgs()
            binding_res = bindings.calcHomodimer(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True
            )
            wrapper_res = wrappers.calcHomodimer(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop
            )
            self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))
            self.assertEqual(
                binding_res.ascii_structure,
                wrapper_res.ascii_structure
            )

    def test_calcHeterodimer(self):
        for _ in range(100):
            self.randArgs()
            binding_res = bindings.calcHeterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True
            )
            wrapper_res = wrappers.calcHeterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop
            )
            self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))
            self.assertEqual(
                binding_res.ascii_structure,
                wrapper_res.ascii_structure
            )
            # Ensure that order of sequences does not matter
            # Should be fixed as of Primer3 2.3.7 update
            binding_12_res = bindings.calcHeterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True
            )
            binding_21_res = bindings.calcHeterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop
            )
            self.assertEqual(int(binding_12_res.tm), int(binding_21_res.tm))

    def test_calcEndStability(self):
        for _ in range(100):
            self.randArgs()
            binding_res = bindings.calcEndStability(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop
            )
            wrapper_res = wrappers.calcEndStability(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop
            )
            self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))

    def test_correctionMethods(self):
        self.randArgs()
        for sc_method in ['schildkraut', 'santalucia', 'owczarzy']:
            for tm_method in ['breslauer', 'santalucia']:
                binding_tm = bindings.calcTm(
                    seq=self.seq1,
                    mv_conc=self.mv_conc,
                    dv_conc=self.dv_conc,
                    dntp_conc=self.dntp_conc,
                    dna_conc=self.dna_conc,
                    tm_method=tm_method,
                    salt_corrections_method=sc_method
                )
                wrapper_tm = wrappers.calcTm(
                    seq=self.seq1,
                    mv_conc=self.mv_conc,
                    dv_conc=self.dv_conc,
                    dntp_conc=self.dntp_conc,
                    dna_conc=self.dna_conc,
                    tm_method=tm_method,
                    salt_corrections_method=sc_method
                )
                self.assertEqual(int(binding_tm), int(wrapper_tm))
        self.assertRaises(
            ValueError,
            bindings.calcTm,
            seq=self.seq1,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            tm_method='not_a_tm_method'
        )

    def test_memoryLeaks(self):
        sm = _getMemUsage()
        for x in range(100):
            self.randArgs()
            bindings.calcHeterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True
            )
        sleep(0.1)  # Pause for any GC
        em = _getMemUsage()
        print('\n\tMemory usage before 1k runs of calcHeterodimer: ', sm)
        print('\tMemory usage after 1k runs of calcHeterodimer:  ', em)
        print('\t\t\t\t\tDifference: \t', em-sm)
        if em-sm > 500:
            raise AssertionError('Memory usage increase after 1k runs of \n\t'
                                 'calcHeterodimer > 500 bytes -- potential \n\t'
                                 'memory leak (mem increase: {})'.format(em-sm))


def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestLowLevelBindings())
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    test_suite = suite()
    runner.run(test_suite)
