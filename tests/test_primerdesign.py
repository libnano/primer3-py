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
tests.test_primerdesign
~~~~~~~~~~~~~~~~~~~~~~~

Unit tests for the primer3-py primer design bindings.

'''

from __future__ import print_function

import os
import random
import sys
import unittest
from time import sleep
from typing import (
    Any,
    Dict,
    List,
    Tuple,
)

try:
    import resource
except (ImportError, ModuleNotFoundError):  # For Windows compatibility
    resource = None  # type: ignore

from primer3 import (
    argdefaults,
    bindings,
)

from . import _simulatedbindings as simulatedbindings

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))


def _get_mem_usage():
    ''' Get current process memory usage in bytes '''
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


class TestDesignBindings(unittest.TestCase):

    def _compare_results(
        self,
        binding_res: Dict[str, Any],
        simulated_binding_res: Dict[str, Any],
        verbose: bool = False,
    ) -> str:
        '''Compare results between binding and simulated binding calls to
        ``design_primers``

        Args:
            binding_res: A dictionary of Primer3 results (should be identical
                to the expected BoulderIO output from primer3_main) from the
                bindings call
            simulated_binding_res: A dictionary of Primer3 results (should be
                identical to the expected BoulderIO output from primer3_main)
                from the simulated bindings call
            verbose: if True, print information

        Returns:
            String of disagreements between the binding and simulated results
        '''
        keys_in_sim = set(simulated_binding_res)
        keys_in_binding = set(binding_res)

        if keys_in_sim - keys_in_binding:
            if verbose:
                print(
                    '\n\n\nIn wrapper simulation result but missing'
                    ' from binding:',
                )
                fmt = '{:<30} {:<50}'
                print(fmt.format('Output Key', 'SimBinding Result'))
                print('-' * 80)
                for k in sorted(keys_in_sim - keys_in_binding):
                    print(fmt.format(k, repr(simulated_binding_res[k])))

        if keys_in_binding - keys_in_sim:
            if verbose:
                print(
                    '\n\n\nIn binding result but missing from wrapper '
                    'simulation:',
                )
                fmt = '{:<30} {:<50}'
                print(fmt.format('Output Key', 'Binding Result'))
                print('-' * 80)
                for k in sorted(keys_in_binding - keys_in_sim):
                    print(fmt.format(k, repr(binding_res[k])))

        allowable_relative_difference = 0.05
        discrepencies: List[str] = [
            k for k in keys_in_binding & keys_in_sim
            if simulated_binding_res[k] != binding_res[k]
        ]
        disagreements_list: List[str] = []
        for ds in discrepencies:
            if (
                isinstance(binding_res[ds], (float, int)) and
                binding_res[ds] != 0
            ):
                percent_diff = abs(
                    (binding_res[ds] - simulated_binding_res[ds]) /
                    binding_res[ds],
                )
                if percent_diff > allowable_relative_difference:
                    if simulated_binding_res[ds] == 0.0 and binding_res[ds] < 0:
                        pass
                    else:
                        disagreements_list.append(ds)

        if len(disagreements_list):
            fmt = '{:<30} {:<25} {:<25}'
            disagreements_str = '\n'.join([
                fmt.format(
                    k,
                    repr(simulated_binding_res[k]),
                    repr(binding_res[k]),
                ) for k in sorted(disagreements_list)
            ])
            if verbose:
                print('\n\n\nResults disagree:')
                print(
                    fmt.format(
                        'Output Key',
                        'SimBinding Result',
                        'Binding Result',
                    ),
                )
                print('-' * 80)
            return disagreements_str
        else:
            if verbose:
                print(
                    '\n\n\nAll the results in common '
                    f'({len(keys_in_binding & keys_in_sim)}) agree to within '
                    f'{allowable_relative_difference:.2%}',
                )
            return ''

    def _convert_boulder_input(
            self,
            boulder_str: str,
    ) -> List[Tuple[Dict, Dict, Dict]]:
        ''' Convert a boulder IO-style input dictionary into bindings /
        simulated-bindings-friendly dictionaries.

        Args:
            boulder_str: Boulder formatted string of multiple records

        Returns:
            List of a tuple of dictionaries of the form::

            [
                (global_args, seq_args, p3_args),
                ...
            ]
        '''
        boulder_dicts = argdefaults.parse_multirecord_boulder_io(boulder_str)
        input_dicts_list = []
        for bd in boulder_dicts:
            converted_input = [(k, v) for k, v in bd.items()]
            global_args = dict(
                filter(
                    lambda arg: 'PRIMER_' == arg[0][:7],
                    converted_input,
                ),
            )
            seq_args = dict(
                filter(
                    lambda arg: 'SEQUENCE_' == arg[0][:9],
                    converted_input,
                ),
            )
            p3_args = dict(
                filter(
                    lambda arg: 'P3_' == arg[0][:3],
                    converted_input,
                ),
            )
            input_dicts_list.append((global_args, seq_args, p3_args))
        return input_dicts_list

    def test_compare_sim(self):
        sequence_template = (
            'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCAT'
            'GCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGG'
            'AGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNA'
            'AAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAA'
            'TGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCA'
            'CACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG'
        )
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
        simulated_binding_res = simulatedbindings.design_primers(
            seq_args,
            global_args,
        )
        binding_res = bindings.design_primers(
            seq_args,
            global_args=global_args,
        )
        self._compare_results(
            binding_res,
            simulated_binding_res,
        )

    def test_file_based(self):
        test_file_roots = [
            'dv_conc_vs_dntp_conc',
            'long_seq',
            'p3-tmpl-mispriming',
            'primer_all_settingsfiles',
            'primer_check',
            'primer_end_pathology',
            'primer_first_base_index',
            'primer_gc_end',
            'primer_high_gc_load_set',
            'primer_high_tm_load_set',
            'primer_internal',
            'primer_must_overlap_point',
            'primer_must_use_th',
            'primer_num_best',
            'primer_ok_regions',
            'primer_overlap_junction',
            'primer_renewed_tasks',
            'primer_start_codon',
            'primer_task_th',
            'primer_task',
            'primer_thal_args',
            'primer_thal_max_seq_error',
            'primer_tm_lc_masking',
            'test_compl_error',
            'test_left_to_right_of_right',
        ]
        print()
        failures = []
        for fn_root in test_file_roots:
            base_fp = os.path.join(
                LOCAL_DIR,
                'input_files',
                fn_root,
            )
            input_fp = f'{base_fp}_input'
            log_path = os.path.join(
                LOCAL_DIR,
                'input_files_log',
            )
            os.makedirs(log_path, exist_ok=True)

            with open(input_fp) as input_fd:
                input_raw = input_fd.read()
            input_dicts = self._convert_boulder_input(input_raw)

            sys.stdout.write(f'->Testing file {fn_root:<40}\r')
            sys.stdout.flush()
            current_global_args = {}

            log_filepath_sim = os.path.join(log_path, f'{fn_root}_sim.log')
            if os.path.exists(log_filepath_sim):
                os.remove(log_filepath_sim)
            log_filepath_bind = os.path.join(log_path, f'{fn_root}_bind.log')
            if os.path.exists(log_filepath_bind):
                os.remove(log_filepath_bind)

            # NOTE: commented out for developement intentionally
            # sim_std_out_fp = 'sim_stdout.txt'
            # if os.path.exists(sim_std_out_fp):
            #     os.remove(sim_std_out_fp)
            # fd = open(sim_std_out_fp, 'wb')
            for global_args, seq_args, p3_args in input_dicts:
                current_global_args.clear()
                global_args['dump'] = 1
                test_id = str(seq_args.get('SEQUENCE_ID', ''))
                current_global_args.update(global_args)
                current_global_args['DO_LOG_SETTINGS'] = 1
                current_global_args['LOG_SETTINGS_PATH'] = log_filepath_sim

                simulated_binding_res = simulatedbindings.design_primers(
                    seq_args,
                    current_global_args,
                    # output_log=fd
                )
                current_global_args['LOG_SETTINGS_PATH'] = log_filepath_bind
                wrapper_error = simulated_binding_res.get('PRIMER_ERROR')
                if wrapper_error is not None:
                    print(wrapper_error)
                    with self.assertRaises((OSError, ValueError)):
                        binding_res = bindings.design_primers(
                            seq_args,
                            current_global_args,
                        )
                else:
                    try:
                        binding_res = bindings.design_primers(
                            seq_args,
                            current_global_args,
                        )
                    except (OSError, TypeError, ValueError):
                        if max([
                            x in p3_args.get('P3_COMMENT', '')
                            for x in ('complain', 'fail')
                        ]):
                            pass
                        raise
                    disagreements_str = self._compare_results(
                        binding_res,
                        simulated_binding_res,
                    )
                    if disagreements_str:
                        failures.append(
                            (fn_root, test_id, disagreements_str),
                        )
            # fd.close()
        print(' ' * 60, end='\r')
        if len(failures):
            err_msg = (
                'Failures occured during file testing:\n' +
                '\n'.join([
                    '->{}\t{}\n{}'.format(*f)
                    for f in failures
                ])
            )
            raise RuntimeError(err_msg)

    @unittest.skipIf(
        sys.platform == 'win32',
        'Windows does not support resource module',
    )
    def test_memory_leaks(self):
        sm = _get_mem_usage()
        run_count = 100
        for x in range(run_count):
            bindings.design_primers(
                {
                    'SEQUENCE_ID': 'MH1000',
                    'SEQUENCE_TEMPLATE': (
                        'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAG'
                        'TGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAA'
                        'GTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTC'
                        'TGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGT'
                        'GTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAG'
                        'GCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTA'
                        'TGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCT'
                        'CTAGTAG'
                    ),
                    'SEQUENCE_INCLUDED_REGION': [36, 342],
                },
                {
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
                },
            )
        sleep(0.1)  # Pause for any GC
        em = _get_mem_usage()
        print(
            f'\n\tMemory usage before {run_count} runs of design_primers: {sm}',
        )
        print(
            f'\tMemory usage after {run_count} runs of design_primers: {em}',
        )
        print(f'\t\t\t\t\tDifference: \t {em - sm}')

        # NOTE: MacOS has a different allocation strategy than Linux
        # Periodically revisit this. 2023.01.09
        delta_bytes_limit = 1000 if sys.platform == 'linux' else 10000

        if em - sm > delta_bytes_limit:
            raise AssertionError(
                f'Memory usage increase after {run_count} runs of \n\t'
                f'design_primers > {delta_bytes_limit} bytes -- potential \n\t'
                f'memory leak (mem increase: {em - sm})',
            )

    def test_misprime_lib_mishyb_lib(self):
        '''Test that the misprime library and mishyb library arguments
        successfully run through the bindings.

        '''
        seq_args = {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': (
                'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT'
                'AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA'
                'ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG'
                'CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG'
                'TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA'
                'TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA'
                'ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA'
                'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA'
                'GATGTTTCCCTCTAGTAG'
            ),
            'SEQUENCE_INCLUDED_REGION': [36, 342],
        }
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
                [75, 100], [100, 125], [125, 150],
                [150, 175], [175, 200], [200, 225],
            ],
        }

        result = bindings.design_primers(
            seq_args=seq_args,
            global_args=global_args,
            misprime_lib={
                'SEQ1': 'CACCATGGAGCTCCTGATATTAAAGGCGAATGCCATT',
            },
        )
        self.assertEqual(result['PRIMER_PAIR_NUM_RETURNED'], 5)
        self.assertEqual(result['PRIMER_LEFT_0'], [46, 21])

        result = bindings.design_primers(
            seq_args=seq_args,
            global_args=global_args,
            mishyb_lib={
                'SEQ1': 'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA',
            },
        )
        self.assertEqual(result['PRIMER_PAIR_NUM_RETURNED'], 5)
        self.assertEqual(result['PRIMER_LEFT_0'], [46, 21])
        self.assertEqual(result['PRIMER_LEFT'][0]['COORDS'], [46, 21])

        bindings.design_primers(
            seq_args=seq_args,
            global_args=global_args,
            misprime_lib={
                'SEQ1': 'CACCATGGAGCTCCTGATATTAAAGGCGAATGCCATT',
            },
            mishyb_lib={
                'SEQ1': 'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA',
            },
        )

        self.assertAlmostEqual(
            result['PRIMER_PAIR_0_PENALTY'],
            1.37323,
            places=4,
        )
        self.assertAlmostEqual(
            result['PRIMER_PAIR'][0]['PENALTY'],
            1.37323,
            places=4,
        )
        self.assertEqual(len(result['PRIMER_PAIR']), 5)

        self.assertEqual(result['PRIMER_PAIR_NUM_RETURNED'], 5)
        self.assertEqual(result['PRIMER_LEFT_0'], [46, 21])
        self.assertEqual(result['PRIMER_LEFT'][0]['COORDS'], [46, 21])
        self.assertEqual(len(result['PRIMER_LEFT']), 5)

        self.assertEqual(result['PRIMER_RIGHT_0'], [132, 20])
        self.assertEqual(result['PRIMER_RIGHT'][0]['COORDS'], [132, 20])
        self.assertEqual(len(result['PRIMER_RIGHT']), 5)

        self.assertEqual(result['PRIMER_INTERNAL_0'], [69, 24])
        self.assertEqual(result['PRIMER_INTERNAL'][0]['COORDS'], [69, 24])
        self.assertEqual(len(result['PRIMER_INTERNAL']), 5)

    def test_PRIMER_SECONDARY_STRUCTURE_ALIGNMENT(self):
        '''Ensure all result pointers are initialized to NULL.
        '''
        seq_args = {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': (
                'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT'
                'AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA'
                'ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG'
                'CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG'
                'TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA'
                'TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA'
                'ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA'
                'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA'
                'GATGTTTCCCTCTAGTAG'
            ),
            'SEQUENCE_INCLUDED_REGION': [36, 342],
        }
        global_args = {
            'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT': 1,  # key parameter for test
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
                [75, 100], [100, 125], [125, 150],
                [150, 175], [175, 200], [200, 225],
            ],
        }
        # This should run without a segmentation fault.
        bindings.design_primers(
            seq_args=seq_args,
            global_args=global_args,
        )


def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestDesignBindings())
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    test_suite = suite()
    runner.run(test_suite)
