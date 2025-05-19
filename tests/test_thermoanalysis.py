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
tests.test_thermoanalysis
~~~~~~~~~~~~~~~~~~~~~~~~~

Unit tests for the primer3-py low level thermodynamic calculation bindings.

'''

from __future__ import print_function

import random
import sys
import unittest
from time import sleep

try:
    import resource

    def _get_mem_usage():
        ''' Get current process memory usage in bytes '''
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
except (ImportError, ModuleNotFoundError):  # For Windows compatibility
    resource = None  # type: ignore

    def _get_mem_usage():  # type: ignore
        ''' Dummy memory usage function for Windows '''
        return 0

from primer3 import (
    bindings,
    thermoanalysis,
)

from . import wrappers


class TestLowLevelBindings(unittest.TestCase):

    def randArgs(self) -> None:
        self.seq1 = ''.join([
            random.choice('ATGC') for _ in
            range(random.randint(20, 59))
        ])
        self.seq2 = ''.join([
            random.choice('ATGC') for _ in
            range(random.randint(20, 59))
        ])

        self.mv_conc = random.uniform(1, 200)
        self.dv_conc = random.uniform(0, 40)
        self.dntp_conc = random.uniform(0, 20)
        self.dna_conc = random.uniform(0, 200)
        self.temp_c = random.randint(10, 70)
        self.max_loop = random.randint(10, 30)
        self.dmso_conc = 0.0
        self.dmso_fact = random.uniform(0.5, 0.7)
        self.formamide_conc = random.uniform(0.8, 1.0)
        self.annealing_temp_c = -10.0  # see oligotm_main.c

        # NOTE: commented out but useful for figuring out bugs in code
        # if they exist
        # self.seq1 = 'GTGTCCAATCGACATTTAGGAGCTACGCGGCCAGGGGCAAA'
        # self.seq2 = 'TGCATTGAATGCGTGTCACGTTATGCACGC'

        # self.mv_conc = 1.1
        # self.dv_conc = 20.0
        # self.dntp_conc = 10.0
        # self.dna_conc = 20.0
        # self.temp_c = 50.0
        # self.max_loop = 20
        # self.dmso_conc = 0.0
        # self.dmso_fact= 0.6
        # self.formamide_conc = 0.0
        # self.annealing_temp_c = -10.0  # see oligotm_main.c
        # self.max_nn_length = 60

    def test_calc_tm(self) -> None:
        '''Test basic calc_tm input'''
        for _ in range(100):
            self.randArgs()
            binding_tm = bindings.calc_tm(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                dmso_conc=self.dmso_conc,
                dmso_fact=self.dmso_fact,
                formamide_conc=self.formamide_conc,
                annealing_temp_c=self.annealing_temp_c,
            )
            wrapper_tm = wrappers.calc_tm(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                dmso_conc=self.dmso_conc,
                dmso_fact=self.dmso_fact,
                formamide_conc=self.formamide_conc,
                annealing_temp_c=self.annealing_temp_c,
            )
            self.assertAlmostEqual(binding_tm, wrapper_tm, delta=0.5)

    def test_calc_hairpin(self) -> None:
        '''Test basic hairpin input'''
        for _ in range(1):
            self.randArgs()
            binding_res = bindings.calc_hairpin(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True,
            )
            wrapper_res = wrappers.calc_hairpin(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
            )
            self.assertAlmostEqual(binding_res.tm, wrapper_res.tm, 2)
            print(binding_res.ascii_structure)
            print()
            print(wrapper_res.ascii_structure)
            self.assertEqual(
                binding_res.ascii_structure,
                # Replace line endings for windows compat
                wrapper_res.ascii_structure.replace('\r\n', '\n'),
            )

    def test_calc_homodimer(self) -> None:
        '''Test basic homodimer input'''
        for _ in range(100):
            self.randArgs()
            binding_res = bindings.calc_homodimer(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True,
            )
            wrapper_res = wrappers.calc_homodimer(
                seq=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
            )
            self.assertAlmostEqual(binding_res.tm, wrapper_res.tm, 2)

            print(binding_res.ascii_structure)
            print()
            print(wrapper_res.ascii_structure)

            self.assertEqual(
                binding_res.ascii_structure,
                # Replace line endings for windows compat
                wrapper_res.ascii_structure.replace('\r\n', '\n'),
            )

    def test_calc_heterodimer(self) -> None:
        '''Test basic heterodimer input'''
        for _ in range(100):
            self.randArgs()
            binding_res = bindings.calc_heterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True,
            )
            wrapper_res = wrappers.calc_heterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
            )
            self.assertAlmostEqual(binding_res.tm, wrapper_res.tm, 2)

            print(binding_res.ascii_structure)
            print()
            print(wrapper_res.ascii_structure)

            self.assertEqual(
                binding_res.ascii_structure,
                # Replace line endings for windows compat
                wrapper_res.ascii_structure.replace('\r\n', '\n'),
            )
            # Ensure that order of sequences does not matter
            binding_12_res = bindings.calc_heterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True,
            )
            binding_21_res = bindings.calc_heterodimer(
                seq1=self.seq2,
                seq2=self.seq1,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
            )
            self.assertAlmostEqual(binding_12_res.tm, binding_21_res.tm, 2)

    def test_max_length_heterodimer(self) -> None:
        '''Test longest heterodimer input of 10000 mer per `THAL_MAX_SEQ` '''
        self.randArgs()

        seq_small = ''.join([
            random.choice('ATGC') for _ in
            range(59)  # max size for sequence
        ])
        seq_big = ''.join([
            random.choice('ATGC') for _ in
            range(10000)  # max size for sequence 2
        ])
        binding_res = bindings.calc_heterodimer(
            seq1=seq_small,
            seq2=seq_big,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            temp_c=self.temp_c,
            max_loop=self.max_loop,
            output_structure=True,
        )
        wrapper_res = wrappers.calc_heterodimer(
            seq1=seq_small,
            seq2=seq_big,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            temp_c=self.temp_c,
            max_loop=self.max_loop,
        )
        self.assertAlmostEqual(binding_res.tm, wrapper_res.tm, 2)

        print(binding_res.ascii_structure)
        print()
        print(wrapper_res.ascii_structure)

        self.assertEqual(
            binding_res.ascii_structure,
            # Replace line endings for windows compat
            wrapper_res.ascii_structure.replace('\r\n', '\n'),
        )
        # Ensure that order of sequences does not matter
        binding_12_res = bindings.calc_heterodimer(
            seq1=seq_small,
            seq2=seq_big,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            temp_c=self.temp_c,
            max_loop=self.max_loop,
            output_structure=True,
        )
        binding_21_res = bindings.calc_heterodimer(
            seq1=seq_big,
            seq2=seq_small,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            temp_c=self.temp_c,
            max_loop=self.max_loop,
        )
        self.assertAlmostEqual(binding_12_res.tm, binding_21_res.tm, 2)

    def test_calc_end_stability(self) -> None:
        '''Test calc_end_stability'''
        for _ in range(100):
            self.randArgs()
            binding_res = bindings.calc_end_stability(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
            )
            wrapper_res = wrappers.calc_end_stability(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
            )
            self.assertAlmostEqual(binding_res.tm, wrapper_res.tm, 2)

    def test_correction_methods(self) -> None:
        '''Test different correction_methods'''
        self.randArgs()
        for sc_method in ['schildkraut', 'santalucia', 'owczarzy']:
            for tm_method in ['breslauer', 'santalucia']:
                binding_tm = bindings.calc_tm(
                    seq=self.seq1,
                    mv_conc=self.mv_conc,
                    dv_conc=self.dv_conc,
                    dntp_conc=self.dntp_conc,
                    dna_conc=self.dna_conc,
                    tm_method=tm_method,
                    salt_corrections_method=sc_method,
                )
                wrapper_tm = wrappers.calc_tm(
                    seq=self.seq1,
                    mv_conc=self.mv_conc,
                    dv_conc=self.dv_conc,
                    dntp_conc=self.dntp_conc,
                    dna_conc=self.dna_conc,
                    tm_method=tm_method,
                    salt_corrections_method=sc_method,
                )
                self.assertAlmostEqual(binding_tm, wrapper_tm, delta=0.5)

        self.assertRaises(
            ValueError,
            bindings.calc_tm,
            seq=self.seq1,
            mv_conc=self.mv_conc,
            dv_conc=self.dv_conc,
            dntp_conc=self.dntp_conc,
            dna_conc=self.dna_conc,
            tm_method='not_a_tm_method',
        )

    @unittest.skipIf(
        sys.platform == 'win32',
        'Windows does not support resource module',
    )
    def test_memory_leaks(self) -> None:
        '''Test for memory leaks'''
        sm = _get_mem_usage()
        run_count = 100
        for x in range(run_count):
            self.randArgs()
            bindings.calc_heterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop,
                output_structure=True,
            )
        sleep(0.1)  # Pause for any GC
        em = _get_mem_usage()
        print(
            f'\n\tMemory usage before {run_count} runs of '
            f'calc_heterodimer: {sm}',
        )
        print(
            f'\tMemory usage after {run_count} runs of calc_heterodimer: {em}',
        )
        print(f'\t\t\t\t\tDifference: {em - sm}\t')
        delta_bytes_limit = 500
        if em - sm > delta_bytes_limit:
            raise AssertionError(
                f'Memory usage increase after {run_count} runs of \n\t'
                f'calc_heterodimer > {delta_bytes_limit} bytes -- potential '
                f'\n\t memory leak (mem increase: {em - sm})',
            )

    def test_todict(self) -> None:
        '''Unit test coverage for ``Thermoanalysis.todict``'''
        args = {
            'mv_conc': 51.0,
            'dv_conc': 1.7,
            'dntp_conc': 0.5,
            'dna_conc': 52.0,
            'temp_c': 36.0,
            'max_loop': 31,
            'max_nn_length': 61,
            'tm_method': 'santalucia',
            'salt_correction_method': 'santalucia',
        }
        ta = thermoanalysis.ThermoAnalysis(**args)
        dict_out = ta.todict()
        for k, v in args.items():
            assert v == dict_out[k]

    def test_calc_tm_acgt_validation(self) -> None:
        '''Test that calc_tm properly validates and converts ACGT sequences'''
        # Test uppercase ACGT
        seq = 'ACGT' * 10
        tm1 = bindings.calc_tm(seq)
        self.assertTrue(tm1 > 0)  # Should get a valid Tm

        # Test lowercase conversion
        seq_lower = 'acgt' * 10
        tm2 = bindings.calc_tm(seq_lower)
        self.assertEqual(tm1, tm2)  # Should get same Tm after conversion

        # Test mixed case
        seq_mixed = 'AcGt' * 10
        tm3 = bindings.calc_tm(seq_mixed)
        self.assertEqual(tm1, tm3)  # Should get same Tm after conversion

        # Test invalid bases
        with self.assertRaises(ValueError) as cm:
            bindings.calc_tm('ACGTN')
        self.assertIn("'N'", str(cm.exception))

        # Test bytes input
        seq_bytes = b'ACGT' * 10
        tm4 = bindings.calc_tm(seq_bytes)
        self.assertEqual(tm1, tm4)  # Should handle bytes input correctly

        # Test lowercase bytes
        seq_bytes_lower = b'acgt' * 10
        tm5 = bindings.calc_tm(seq_bytes_lower)
        self.assertEqual(tm1, tm5)  # Should handle lowercase bytes correctly

        # Test invalid bytes
        with self.assertRaises(ValueError) as cm:
            bindings.calc_tm(b'ACGTN')
        self.assertIn("'N'", str(cm.exception))


def suite() -> unittest.TestSuite:
    '''Run all tests'''
    s = unittest.TestSuite()
    s.addTests((
        unittest.makeSuite(TestLowLevelBindings),
    ))
    return s


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    test_suite = suite()
    runner.run(test_suite)
