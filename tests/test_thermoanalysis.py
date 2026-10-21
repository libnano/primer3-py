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
tests.test_thermoanalysis
~~~~~~~~~~~~~~~~~~~~~~~~~

Unit tests for the primer3-py low level thermodynamic calculation bindings.

'''

from __future__ import print_function

import math
import random
import unittest
from typing import Any

from primer3 import (
    bindings,
    thermoanalysis,
)

from . import (
    _leakcheck,
    wrappers,
)


def _owczarzy_eq16_tm(
        tm_1m_c: float,
        seq: str,
        mv_mm: float,
        dv_mm: float,
        dntp_mm: float,
) -> float:
    '''Independent implementation of the Owczarzy et al. (2008) divalent
    cation correction (Eqs. 16 and 18-20), applied to the 1 M Na+ Tm.

    Only valid in the regime where Eq. 16 applies (no monovalent cations, or
    sqrt([Mg2+]) / [Mon+] >= 0.22).
    '''
    n = len(seq)
    fgc = sum(c in 'GC' for c in seq.upper()) / n
    mg = (dv_mm - dntp_mm) / 1000.0
    mon = mv_mm / 1000.0
    a, b, c, d = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5
    e, f, g = -4.82e-4, 5.25e-4, 8.31e-5
    if mon > 0 and math.sqrt(mg) / mon < 6.0:
        ln_mon = math.log(mon)
        a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mon) * ln_mon)
        d = 1.42e-5 * (1.279 - 4.03e-3 * ln_mon - 8.03e-3 * ln_mon ** 2)
        g = 8.31e-5 * (0.486 - 0.258 * ln_mon + 5.25e-3 * ln_mon ** 3)
    ln_mg = math.log(mg)
    correction = (
        a + b * ln_mg + fgc * (c + d * ln_mg) +
        (1.0 / (2.0 * (n - 1))) * (e + f * ln_mg + g * ln_mg ** 2)
    )
    return 1.0 / (1.0 / (tm_1m_c + 273.15) + correction) - 273.15


class TestLowLevelBindings(unittest.TestCase):

    def setUp(self) -> None:
        '''Set up test case'''
        random.seed(42)  # Make random operations deterministic

    def rand_args(self) -> None:
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
            self.rand_args()
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

    def test_calc_tm_short_sequences(self) -> None:
        '''Sequences shorter than 2 nt must return the OLIGOTM_ERROR sentinel
        rather than reading out of bounds (empty) or returning a nonsense Tm
        (1 nt). See docs/included_primer3_modifications.md (oligotm.c).
        '''
        OLIGOTM_ERROR = -999999.9999
        for seq in ('', 'A', 'a', b'', b'C'):
            self.assertAlmostEqual(
                bindings.calc_tm(seq),
                OLIGOTM_ERROR,
                places=3,
            )
        # A 2-nt sequence is still computed normally (not the sentinel).
        self.assertGreater(bindings.calc_tm('AT'), OLIGOTM_ERROR)

    def test_calc_tm_owczarzy_length_term(self) -> None:
        '''The length-dependent term of the Owczarzy et al. (2008) Eq. 16
        divalent cation correction must not be dropped by integer division.
        See docs/included_primer3_modifications.md (oligotm.c) and issue #178.
        '''
        kw: dict[str, Any] = {
            'tm_method': 'santalucia',
            'salt_corrections_method': 'owczarzy',
        }
        # (seq, mv_conc, dv_conc, dntp_conc, dna_conc, expected Tm). The
        # expected values are Eq. 16 including the length term (they agree
        # with upstream primer3 >= c5b092b); the first is also within 0.5 C
        # of the measured value in Owczarzy et al. 2008, Table S1 (27.4 C).
        cases = (
            ('ATCAATCATA', 0.0, 3.0, 0.0, 2000.0, 26.954),
            ('TGCATGCATGCA', 0.0, 10.0, 0.0, 2000.0, 58.458),
            ('AGCGTAAGCTGATCGATCCT', 50.0, 3.0, 0.8, 250.0, 62.215),
        )
        for seq, mv, dv, dntp, dna, expected in cases:
            tm = bindings.calc_tm(
                seq,
                mv_conc=mv,
                dv_conc=dv,
                dntp_conc=dntp,
                dna_conc=dna,
                **kw,
            )
            self.assertAlmostEqual(tm, expected, delta=0.01, msg=seq)

    def test_calc_tm_owczarzy_matches_eq16(self) -> None:
        '''calc_tm with the Owczarzy correction must match an independent
        implementation of Eq. 16 across lengths and salt conditions where
        the divalent correction applies.
        '''
        rng = random.Random(178)
        n_checked = 0
        while n_checked < 300:
            seq = ''.join(
                rng.choice('ACGT') for _ in range(rng.randint(8, 45))
            )
            mv = rng.choice((0.0, 0.0, 10.0, 25.0, 50.0, 100.0))
            dv = rng.uniform(0.5, 20.0)
            dntp = rng.uniform(0.0, dv * 0.8)
            dna = rng.uniform(50.0, 2000.0)
            free_mg = (dv - dntp) / 1000.0
            if mv > 0 and math.sqrt(free_mg) / (mv / 1000.0) < 0.22:
                continue  # monovalent-dominated regime (Eq. 22), not Eq. 16
            kw: dict[str, Any] = {
                'dna_conc': dna,
                'tm_method': 'santalucia',
                'salt_corrections_method': 'owczarzy',
            }
            tm = bindings.calc_tm(
                seq, mv_conc=mv, dv_conc=dv, dntp_conc=dntp, **kw,
            )
            tm_1m = bindings.calc_tm(
                seq, mv_conc=1000.0, dv_conc=0.0, dntp_conc=0.0, **kw,
            )
            expected = _owczarzy_eq16_tm(tm_1m, seq, mv, dv, dntp)
            self.assertAlmostEqual(tm, expected, delta=1e-4, msg=seq)
            n_checked += 1

    def test_thermo_mismatch_stack_constants(self) -> None:
        '''Single-mismatch ("1x1 internal loop") stack enthalpies for AC/TC,
        CT/CA (-10 -> 0 cal/mol) and TA/AA (forbidden -> 4700 cal/mol) must
        match upstream primer3's corrected stackmm.dh. See
        docs/included_primer3_modifications.md (thal_parameters.c).

        Expected values were produced by upstream primer3's ntthal built from
        the commit that fixed those parameters (d2844b5); each case below is
        one whose result differs from primer3 2.6.1.
        '''
        kw: dict[str, Any] = {
            'mv_conc': 50.0,
            'dv_conc': 1.5,
            'dntp_conc': 0.6,
            'dna_conc': 50.0,
            'temp_c': 37.0,
        }
        # (callable, args, tm, dg, dh, ds)
        cases = (
            # TA/AA, heterodimer
            (
                bindings.calc_heterodimer, ('GCGCTAGCGC', 'GCGCAAGCGC'),
                30.4142, -9822.87, -64200.0, -175.325,
            ),
            (
                bindings.calc_heterodimer,
                ('GCCCGGAAAGACACTA', 'TTTGACTTACC'),
                -25.0143, -2693.39, -34100.0, -101.263,
            ),
            # TA/AA, hairpin (2.6.1 finds no hairpin at all)
            (
                bindings.calc_hairpin, ('TTCAATGGGCCATAGAA',),
                27.4000, 824.092, -25800.0, -85.8426,
            ),
            # TA/AA, homodimer
            (
                bindings.calc_homodimer, ('GCTTACCTAAGGCAAGC',),
                16.5862, -6452.83, -67600.0, -197.154,
            ),
            # AC/TC and/or CT/CA: -10 -> 0 cal/mol shifts dH by +10 per stack
            (
                bindings.calc_heterodimer,
                ('GGCCCTAACACCG', 'ATCGCTGACCAGC'),
                -44.525, -1338.19, -27700.0, -84.997,
            ),
        )
        for func, args, tm, dg, dh, ds in cases:
            res = func(*args, **kw)
            msg = '{}{}'.format(func.__name__, args)
            self.assertTrue(res.structure_found, msg=msg)
            self.assertAlmostEqual(res.tm, tm, delta=0.01, msg=msg)
            self.assertAlmostEqual(res.dg, dg, delta=0.1, msg=msg)
            self.assertAlmostEqual(res.dh, dh, delta=1.0, msg=msg)
            self.assertAlmostEqual(res.ds, ds, delta=0.01, msg=msg)

    def test_thermo_triloop_tetraloop_tables(self) -> None:
        '''Hairpins closed by special triloop / tetraloop table entries must
        still be parsed and scored correctly now that readTLoop no longer
        writes a NUL terminator one byte past struct triloop/tetraloop::loop.
        See docs/included_primer3_modifications.md (thal.c).
        '''
        kw: dict[str, Any] = {
            'mv_conc': 50.0,
            'dv_conc': 1.5,
            'dntp_conc': 0.6,
            'dna_conc': 50.0,
            'temp_c': 37.0,
        }
        # (seq, tm, dg, dh, ds) from upstream primer3's ntthal (HAIRPIN).
        cases = (
            ('GCGCCGAAGGCGC', 95.0595, -6338.76, -40200.0, -109.177),  # CGAAG
            ('GCGCCAAAAGGCGC', 66.5099, -3231.96, -37200.0, -109.521),  # CAAAAG
            ('GCGCCGATAGGCGC', 93.3972, -7031.44, -45700.0, -124.677),
        )
        for seq, tm, dg, dh, ds in cases:
            res = bindings.calc_hairpin(seq, **kw)
            self.assertAlmostEqual(res.tm, tm, delta=0.01, msg=seq)
            self.assertAlmostEqual(res.dg, dg, delta=0.1, msg=seq)
            self.assertAlmostEqual(res.dh, dh, delta=1.0, msg=seq)
            self.assertAlmostEqual(res.ds, ds, delta=0.01, msg=seq)

    def test_calc_hairpin(self) -> None:
        '''Test basic hairpin input'''
        for _ in range(1):
            self.rand_args()
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
            self.rand_args()
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
            self.rand_args()
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
        self.rand_args()

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
            self.rand_args()
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
        self.rand_args()
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

    def test_memory_leaks(self) -> None:
        '''Leak check for the thermodynamic calc bindings.

        Exercises the ``output_structure`` path (raw ``malloc`` of the ASCII
        structure buffer) and ``calc_tm`` on one shared ``ThermoAnalysis``
        instance. Tracemalloc bounds the Python/PyMem growth deterministically;
        a warmed-up peak-RSS bound guards the raw-malloc structure buffers.
        '''
        self.rand_args()
        ta = thermoanalysis.ThermoAnalysis()
        s1, s2 = self.seq1, self.seq2

        def _work():
            '''Run one Tm plus each structure-output calculation'''
            ta.calc_tm(s1)
            ta.calc_heterodimer(s1, s2, output_structure=True)
            ta.calc_homodimer(s1, output_structure=True)
            ta.calc_hairpin(s1, output_structure=True)

        _leakcheck.assert_no_leak(
            self,
            _work,
            iters=200,
            warmup=100,
            max_tracemalloc_kib=64,
            max_peak_rss_kib=4096,
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
