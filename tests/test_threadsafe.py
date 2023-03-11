# Copyright (C) 2023. Ben Pruitt & Nick Conway;
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
tests.test_threadsafe
~~~~~~~~~~~~~~~~~~~~~~~~~

Unit tests for the :class:`thermoanalysis.ThermoAnalysis` calls are thread-safe
using the new `thalflex` code

NOTE: It's challenging to get a race condition to compromise thread safety
as the thal() calls are very fast

'''
import random
import threading
import time
import unittest
from typing import (
    Any,
    Dict,
)

from primer3 import thermoanalysis  # type: ignore


def make_random_params() -> Dict[str, Any]:
    '''
    Returns:
        dictionary of random parameters for a call  to
        :meth`ThermoAnalysis.set_thermo_args`

    '''
    return dict(
        mv_conc=random.uniform(1, 200),
        dv_conc=random.uniform(0, 40),
        dntp_conc=random.uniform(0, 20),
        dna_conc=random.uniform(0, 200),
        temp_c=random.randint(10, 70),
        max_loop=random.randint(10, 30),
    )


def make_rand_seq_pair() -> Dict[str, str]:
    '''
    Returns:
        dictionary of 2 unique sequences keyed by seq1 and seq2

    '''
    return dict(
        seq1=''.join([
            random.choice('ATGC') for _ in
            range(random.randint(20, 59))
        ]),
        seq2=''.join([
            random.choice('ATGC') for _ in
            range(random.randint(20, 59))
        ]),
    )


class TestThermoAnalysisInThread(unittest.TestCase):

    def test_calc_heterodimer_in_thread(self):
        '''Test basic heterodimer input in a thread

        '''
        def het_thread(
                _results_list,
                i,
        ):
            _taf = thermoanalysis.ThermoAnalysis()
            _params = make_random_params()
            _taf.set_thermo_args(
                **_params,
            )
            _seq_pair = make_rand_seq_pair()
            _res = _taf.calc_heterodimer(
                **_seq_pair,
            )
            _results_list[i] = (_res, _params, _seq_pair)

        thread_count = 10
        results_list = [None] * thread_count
        thread_list = []
        t0 = time.time()
        # 1. Run in threads
        for i in range(thread_count):
            t = threading.Thread(
                target=het_thread,
                args=(results_list, i),
            )
            thread_list.append(t)
            t.start()

        for t in thread_list:
            t.join()

        print(f'total: {(time.time()-t0):0.2}')

        taf = thermoanalysis.ThermoAnalysis()
        t0 = time.time()

        # 2. Compare threaded results to sequential results
        for res, params, seq_pair in results_list:
            taf.set_thermo_args(
                **params,
            )
            new_res = taf.calc_heterodimer(
                **seq_pair,
            )
            self.assertEqual(new_res.tm, res.tm)
            self.assertEqual(new_res.dh, res.dh)
        print(f'total: {(time.time()-t0):0.2}')

        def tdesign_run(_global_args, _results_list, i):
            sequence_template = (
                'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACA'
                'GCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGC'
                'TAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTG'
                'TTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATG'
                'GAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGAT'
                'AACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTT'
                'TCCCTCTAGTAG'
            ) + ''.join([
                random.choice('ATGC') for _ in
                range(60)
            ])
            quality_list = [
                random.randint(20, 90)
                for i in range(len(sequence_template))
            ]
            seq_args = {
                'SEQUENCE_ID': 'MH1000',
                'SEQUENCE_TEMPLATE': sequence_template,
                'SEQUENCE_QUALITY': quality_list,
                'SEQUENCE_INCLUDED_REGION': (36, 342),
            }
            _taf = thermoanalysis.ThermoAnalysis()
            result = _taf.run_design(
                global_args=_global_args,
                seq_args=seq_args,
                misprime_lib=None,
                mishyb_lib=None,
            )
            _results_list[i] = (result, seq_args, i)

        thread_count = 10
        results_list = [None] * thread_count
        thread_list = []
        t0 = time.time()

        global_args = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [
                [75, 100],
                [100, 125],
                [125, 150],
                [150, 175],
                [175, 200],
                [200, 225],
            ],
        }

        # 3. Run in threads
        t0 = time.time()
        for i in range(thread_count):
            t = threading.Thread(
                target=tdesign_run,
                args=(global_args, results_list, i),
            )
            thread_list.append(t)
            t.start()
        print(f'total: {(time.time()-t0):0.2}')

        taf = thermoanalysis.ThermoAnalysis()
        t0 = time.time()
        # 4. Compare threaded results to sequential results
        for res, seq_args, i in results_list:
            new_res = taf.run_design(
                global_args=global_args,
                seq_args=seq_args,
                misprime_lib=None,
                mishyb_lib=None,
            )
            self.assertEqual(new_res, res)
        print(f'total: {(time.time()-t0):0.2}')
        # raise ValueError('')
