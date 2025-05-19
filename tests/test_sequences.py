# Copyright (C) 2020-2025. Ben Pruitt & Nick Conway; Wyss Institute
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
tests.test_sequences
~~~~~~~~~~~~~~~~~~~

Unit tests for thermodynamic calculations using standardized test sequences.
Tests verify both regression values and expected thermodynamic relationships
between different types of structures (hairpins, heterodimers, end stability).

'''

from __future__ import print_function

import json
import unittest
from typing import (
    Any,
    Dict,
)

from primer3 import (
    argdefaults,
    thermoanalysis,
)

# Use the default values from argdefaults.py
STANDARD_CONDITIONS: Dict[str, Any] = {
    'mv_conc': argdefaults.Primer3PyArguments.mv_conc,  # 50.0 mM
    'dv_conc': argdefaults.Primer3PyArguments.dv_conc,  # 1.5 mM
    'dntp_conc': argdefaults.Primer3PyArguments.dntp_conc,  # 0.6 mM
    'dna_conc': argdefaults.Primer3PyArguments.dna_conc,  # 50.0 nM
    'temp_c': argdefaults.Primer3PyArguments.temp_c,  # 37.0 °C
    'max_loop': argdefaults.Primer3PyArguments.max_loop,  # 30
    'dmso_conc': argdefaults.Primer3PyArguments.dmso_conc,  # 0.0 %
    'dmso_fact': argdefaults.Primer3PyArguments.dmso_fact,  # 0.6
    'formamide_conc': argdefaults.Primer3PyArguments.formamide_conc,  # 0.0 %
}

# Test sequences organized by structural properties
TEST_SEQUENCES: Dict[str, Dict[str, Any]] = {
    'homopolymers': {
        'polyA': 'AAAAAAAAAAAAAAAAAAAA',
        'polyT': 'TTTTTTTTTTTTTTTTTTTT',
        'polyG': 'GGGGGGGGGGGGGGGGGGGG',
        'polyC': 'CCCCCCCCCCCCCCCCCCCC',
    },
    'palindromes': {
        'gc_rich': 'GCGCGCGCGCGCGCGCGCGC',
        'at_rich': 'ATATATATATATATATATAT',
        'mixed': 'GATCGATCGATCGATCGATC',
    },
    'hairpins': {
        # Perfect hairpin with 8bp stem and 4nt loop
        'perfect': 'GCGCGCGCAAAAGCGCGCGC',
        # Hairpin with a single mismatch in the stem
        'mismatched': 'GCGTGCGCAAAAGCGTGCGC',
        # Hairpin with a 2nt bulge in the stem
        'bulged': 'GCGTTCGCGCGCGCTTGCGC',
    },
    'heterodimers': {
        'complementary': {
            'seq1': 'AGCCCATACCGCTGTTGTTG',
            'seq2': 'CAACAACAGCGGTATGGGCT',
        },
        'mismatched': {
            'seq1': 'ATCGATCGATCGATCGATCG',
            'seq2': 'TTGGTGTCGGCTGGCCTCGG',
        },
        'overlapping': {
            'seq1': 'AGCCCATACCCGATCGATCG',
            'seq2': 'CAACAACAGCATCGATCGAT',
        },
    },
    'end_stability': {
        'perfect': {  # seq1 3' end, 8nt, complementary to middle of seq2
            'seq1': 'ATAAGCACAATTTTAAAGCC',
            'seq2': 'GAGGGAGGCTTTAACTGCTC',
        },
        'partial': {  # same as "perfect" with 1 nt mismatch at position -3
            'seq1': 'ATAAGCACAATTTTAAACCC',
            'seq2': 'GAGGGAGGCTTTAACTGCTC',
        },
        'mismatched': {  # seq1 3' end, 8nt, not complementary to middle of seq2
            'seq1': 'ATAAGCACAATTTTAAAGCC',
            'seq2': 'GAGGGATAGCCCAACTGCTC',
        },
    },
}


def calculate_thermo_values() -> Dict[str, Any]:
    '''Calculate thermodynamic values for all test sequences under standard
    conditions
    '''
    thermo = thermoanalysis.ThermoAnalysis()
    thermo.set_thermo_args(**STANDARD_CONDITIONS)

    results: Dict[str, Dict[str, Any]] = {}

    # Calculate values for single sequences
    for category, sequences in TEST_SEQUENCES.items():
        if category not in ['heterodimers', 'end_stability']:
            results[category] = {}
            for name, seq in sequences.items():
                results[category][name] = {
                    'tm': thermo.calc_tm(seq),
                    'hairpin': thermo.calc_hairpin(seq).todict(),
                    'homodimer': thermo.calc_homodimer(seq).todict(),
                }

    # Calculate values for sequence pairs
    for category in ['heterodimers', 'end_stability']:
        results[category] = {}
        for name, pair in TEST_SEQUENCES[category].items():
            results[category][name] = {
                'heterodimer': thermo.calc_heterodimer(
                    pair['seq1'], pair['seq2'],
                ).todict(),
                'end_stability': thermo.calc_end_stability(
                    pair['seq1'], pair['seq2'],
                ).todict(),
            }

    return results


def save_thermo_values(
        values: Dict[str, Any],
        filename: str = 'tests/thermo_standard_values.json',
) -> None:
    '''Save thermodynamic values to a JSON file'''
    with open(filename, 'w') as f:
        json.dump(values, f, indent=2)


def load_thermo_values(
        filename: str = 'tests/thermo_standard_values.json',
) -> Dict[str, Any]:
    '''Load thermodynamic values from a JSON file'''
    with open(filename, 'r') as f:
        return json.load(f)


class TestThermodynamicRelationships(unittest.TestCase):
    '''Test suite for verifying expected thermodynamic relationships between
    different types of DNA structures and sequences.
    '''

    def setUp(self) -> None:
        self.conditions = STANDARD_CONDITIONS.copy()
        self.thermo = thermoanalysis.ThermoAnalysis()
        self.thermo.set_thermo_args(**self.conditions)

    def test_gc_content_effect(self) -> None:
        '''Test that GC-rich sequences have higher melting temperatures than
        AT-rich ones
        '''
        gc_rich = TEST_SEQUENCES['palindromes']['gc_rich']
        at_rich = TEST_SEQUENCES['palindromes']['at_rich']

        gc_tm = self.thermo.calc_tm(gc_rich)
        at_tm = self.thermo.calc_tm(at_rich)

        self.assertGreater(gc_tm, at_tm)

    def test_hairpin_stability(self) -> None:
        '''Test that perfect hairpins are more stable than mismatched ones'''
        perfect = TEST_SEQUENCES['hairpins']['perfect']
        mismatched = TEST_SEQUENCES['hairpins']['mismatched']

        perfect_dg = self.thermo.calc_hairpin(perfect).dg
        mismatched_dg = self.thermo.calc_hairpin(mismatched).dg

        self.assertLess(perfect_dg, mismatched_dg)

    def test_hairpin_stability_relationships(self) -> None:
        '''Test that hairpin stabilities follow expected order: perfect <
        mismatched < bulged.
        '''
        perfect = TEST_SEQUENCES['hairpins']['perfect']
        mismatched = TEST_SEQUENCES['hairpins']['mismatched']
        bulged = TEST_SEQUENCES['hairpins']['bulged']

        perfect_dg = self.thermo.calc_hairpin(perfect).dg
        mismatched_dg = self.thermo.calc_hairpin(mismatched).dg
        bulged_dg = self.thermo.calc_hairpin(bulged).dg

        self.assertLess(perfect_dg, mismatched_dg)
        self.assertLess(mismatched_dg, bulged_dg)

    def test_heterodimer_stability_relationships(self) -> None:
        '''Test that heterodimer stabilities follow expected order:
        complementary < overlapping < mismatched.
        '''
        complementary = TEST_SEQUENCES['heterodimers']['complementary']
        overlapping = TEST_SEQUENCES['heterodimers']['overlapping']
        mismatched = TEST_SEQUENCES['heterodimers']['mismatched']

        comp_dg = self.thermo.calc_heterodimer(
            complementary['seq1'], complementary['seq2'],
        ).dg
        overlap_dg = self.thermo.calc_heterodimer(
            overlapping['seq1'], overlapping['seq2'],
        ).dg
        mismatch_dg = self.thermo.calc_heterodimer(
            mismatched['seq1'], mismatched['seq2'],
        ).dg

        self.assertLess(comp_dg, overlap_dg)
        self.assertLess(overlap_dg, mismatch_dg)

    def test_end_stability_relationships(self) -> None:
        '''Test that end stability follows expected order: perfect < partial <
        mismatched.
        '''
        perfect = TEST_SEQUENCES['end_stability']['perfect']
        partial = TEST_SEQUENCES['end_stability']['partial']
        mismatched = TEST_SEQUENCES['end_stability']['mismatched']

        perfect_dg = self.thermo.calc_end_stability(
            perfect['seq1'], perfect['seq2'],
        ).dg
        partial_dg = self.thermo.calc_end_stability(
            partial['seq1'], partial['seq2'],
        ).dg
        mismatched_dg = self.thermo.calc_end_stability(
            mismatched['seq1'], mismatched['seq2'],
        ).dg

        self.assertLess(perfect_dg, partial_dg)
        self.assertLess(partial_dg, mismatched_dg)

    def test_parameter_effects(self) -> None:
        '''Test that changing parameters affects all sequences consistently'''
        base_conditions = self.conditions.copy()
        high_salt_conditions = base_conditions.copy()
        high_salt_conditions['mv_conc'] = 200.0

        # Test a few sequences
        for seq in [
            TEST_SEQUENCES['palindromes']['gc_rich'],
            TEST_SEQUENCES['palindromes']['at_rich'],
        ]:
            self.thermo.set_thermo_args(**base_conditions)
            base_tm = self.thermo.calc_tm(seq)
            self.thermo.set_thermo_args(**high_salt_conditions)
            high_salt_tm = self.thermo.calc_tm(seq)

            # Higher salt should increase Tm
            self.assertGreater(high_salt_tm, base_tm)


class TestThermodynamicRegression(unittest.TestCase):
    '''Test suite for verifying that thermodynamic calculations match expected
    previously calculated across different types of DNA structures.
    '''

    def setUp(self) -> None:
        self.conditions = STANDARD_CONDITIONS.copy()
        self.thermo = thermoanalysis.ThermoAnalysis()
        self.thermo.set_thermo_args(**self.conditions)
        self.standard_values = load_thermo_values()['values']

    def assert_thermo_result_equal(self, result1, result2, places=2) -> None:
        '''Assert that two thermodynamic results are approximately equal'''
        for key, value in result1.items():
            if isinstance(value, (int, float)):
                self.assertAlmostEqual(value, result2[key], places=places)
            elif isinstance(value, dict):
                self.assert_thermo_result_equal(value, result2[key], places)

    def test_hairpin_regression(self) -> None:
        '''Test that hairpin calculations match standard values'''
        for name, seq in TEST_SEQUENCES['hairpins'].items():
            result = self.thermo.calc_hairpin(seq).todict()
            standard = self.standard_values['hairpins'][name]['hairpin']
            self.assert_thermo_result_equal(result, standard)

    def test_homodimer_regression(self) -> None:
        '''Test that homodimer calculations match standard values'''
        for category in ['homopolymers', 'palindromes', 'hairpins']:
            for name, seq in TEST_SEQUENCES[category].items():
                result = self.thermo.calc_homodimer(seq).todict()
                standard = self.standard_values[category][name]['homodimer']
                self.assert_thermo_result_equal(result, standard)

    def test_heterodimer_regression(self) -> None:
        '''Test that heterodimer calculations match standard values'''
        for name, pair in TEST_SEQUENCES['heterodimers'].items():
            result = self.thermo.calc_heterodimer(
                pair['seq1'], pair['seq2'],
            ).todict()
            standard = self.standard_values['heterodimers'][name]['heterodimer']
            self.assert_thermo_result_equal(result, standard)

    def test_end_stability_regression(self) -> None:
        '''Test that end stability calculations match standard values'''
        for name, pair in TEST_SEQUENCES['end_stability'].items():
            result = self.thermo.calc_end_stability(
                pair['seq1'], pair['seq2'],
            ).todict()
            standard = (
                self.standard_values['end_stability'][name]['end_stability']
            )
            self.assert_thermo_result_equal(result, standard)


if __name__ == '__main__':
    unittest.main()
