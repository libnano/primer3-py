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
test_primerdesign
~~~~~~~~~~~~~~~~~

Unit tests for the primer3-py primer design bindings.

'''

from __future__ import print_function

import os
import random
import sys
import unittest
from time import sleep

try:
    import resource
except: # For Windows compatibility
    resource = None

from primer3 import (
    bindings,
    wrappers
)

from . import _simulatedbindings as simulatedbindings

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))


def _getMemUsage():
    """ Get current process memory usage in bytes """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

@unittest.skipIf(
    sys.platform=='win32',
    "Windows doesn't support resource module and wrappers")
class TestDesignBindings(unittest.TestCase):

    def _compareResults(self, binding_res, simulated_binding_res,
                        verbose=False):
        keys_in_sim = set(simulated_binding_res)
        keys_in_binding = set(binding_res)

        if keys_in_sim - keys_in_binding:
            if verbose:
                print('\n\n\nIn wrapper simulation result but missing'
                      ' from binding:')
                fmt = '{:<30} {:<50}'
                print(fmt.format('Output Key', 'SimBinding Result'))
                print('-'*80)
                for k in sorted(keys_in_sim - keys_in_binding):
                    print(fmt.format(k, repr(simulated_binding_res[k])))

        if keys_in_binding - keys_in_sim:
            if verbose:
                print('\n\n\nIn binding result but missing from wrapper '
                      'simulation:')
                fmt = '{:<30} {:<50}'
                print(fmt.format('Output Key', 'Binding Result'))
                print('-'*80)
                for k in sorted(keys_in_binding - keys_in_sim):
                    print(fmt.format(k, repr(binding_res[k])))

        allowable_relative_difference = 0.05
        discrepencies = [k for k in keys_in_binding & keys_in_sim
                         if simulated_binding_res[k] != binding_res[k]]
        disagreements = []
        for ds in discrepencies:
            if (isinstance(binding_res[ds], (float, int)) and
                    binding_res[ds] != 0):
                percent_diff = abs((binding_res[ds] - simulated_binding_res[ds])
                                    / binding_res[ds])
                if percent_diff > allowable_relative_difference:
                    if simulated_binding_res[ds] == 0.0 and binding_res[ds] < 0:
                        pass
                    else:
                        disagreements.append(ds)

        if len(disagreements):
            fmt = '{:<30} {:<25} {:<25}'
            disagreements = '\n'.join([fmt.format(k,
                                        repr(simulated_binding_res[k]),
                                        repr(binding_res[k])) for k in
                                        sorted(disagreements)])
            if verbose:
                print('\n\n\nResults disagree:')
                print(fmt.format('Output Key', 'SimBinding Result',
                                 'Binding Result'))
                print('-'*80)
            return disagreements
        else:
            if verbose:
                print('\n\n\nAll the results in common ({}) agree to within '
                      '{:.2%}'.format(len(keys_in_binding & keys_in_sim),
                                      allowable_relative_difference))

    def _convertBoulderInput(self, boulder_str):
        ''' Convert a boulder IO-style input dictionary into bindings /
        simulated-bindings-friendly dictionaries.
        '''
        boulder_dicts = wrappers._parseMultiRecordBoulderIO(boulder_str)
        input_dicts = []
        for bd in boulder_dicts:
            converted_input = [simulatedbindings.unwrap(arg) for arg in
                               bd.items()]
            global_args = dict(filter(lambda arg: "PRIMER_" == arg[0][:7],
                                      converted_input))
            seq_args = dict(filter(lambda arg: "SEQUENCE_" == arg[0][:9],
                                    converted_input))
            p3_args = dict(filter(lambda arg: "P3_" == arg[0][:3],
                                    converted_input))
            input_dicts.append((global_args, seq_args, p3_args))
        return input_dicts


    def test_compareSim(self):
        sequence_template = 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG'
        quality_list = [random.randint(20,90) for i in range(len(sequence_template))]
        seq_args = {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': sequence_template,
            'SEQUENCE_QUALITY': quality_list,
            'SEQUENCE_INCLUDED_REGION': [36,342]
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
            'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],
        }
        simulated_binding_res = simulatedbindings.designPrimers(seq_args, global_args)
        binding_res = bindings.designPrimers(seq_args, global_args)
        self._compareResults(binding_res, simulated_binding_res)

    def test_fileBased(self):
        test_file_roots = [
            'primer_must_use_th',
            'primer_task_th',
            'primer_thal_args',
            'primer_thal_max_seq_error',
            'primer_first_base_index',
            'test_compl_error',
            'test_left_to_right_of_right',
            'dv_conc_vs_dntp_conc',
            'primer_internal',
            'primer_tm_lc_masking',
            'primer_ok_regions',
            'primer_start_codon',
            'primer_task',
            'primer_renewed_tasks',
            'primer_must_overlap_point',
            'primer_overlap_junction',
            'primer_all_settingsfiles',
            'primer_high_tm_load_set',
            'primer_high_gc_load_set',
            'primer_gc_end',
            'primer_num_best',
            'primer_check',
            'primer_end_pathology',
            'long_seq',
            'p3-tmpl-mispriming'
        ]
        print()
        failures = []
        for fn_root in test_file_roots:
            base_fp = os.path.join(LOCAL_DIR, 'input_files', fn_root)
            input_fp = base_fp + '_input'

            with open(input_fp) as input_fd:
                input_raw = input_fd.read()
            input_dicts = self._convertBoulderInput(input_raw)

            sys.stdout.write('->Testing file {:<40}\r'.format(fn_root))
            sys.stdout.flush()
            current_global_args = {}
            for global_args, seq_args, p3_args in input_dicts:
                test_id = str(seq_args.get('SEQUENCE_ID', ''))
                current_global_args.update(global_args)
                simulated_binding_res = simulatedbindings.designPrimers(
                                            seq_args, current_global_args)
                wrapper_error = simulated_binding_res.get('PRIMER_ERROR')
                if wrapper_error is not None:
                    with self.assertRaises(IOError):
                        binding_res = bindings.designPrimers(seq_args,
                                                            current_global_args)
                else:
                    try:
                        binding_res = bindings.designPrimers(seq_args,
                                                            current_global_args)
                    except IOError:
                        if max([x in p3_args.get('P3_COMMENT', '') for x in
                                ('complain', 'fail')]):
                            pass
                    disagreements = self._compareResults(binding_res,
                                                         simulated_binding_res)
                    if disagreements is not None:
                        failures.append((fn_root, test_id, disagreements))
        print(' '* 60, end='\r')
        if len(failures):
            err_msg = ('Failures occured during file testing:\n' +
                      '\n'.join(['->{}\t{}\n{}'.format(*f) for f in
                                 failures]))
            raise RuntimeError(err_msg)

    def test_memoryLeaks(self):
        sm = _getMemUsage()
        for x in range(100):
            bindings.designPrimers(
                {
                    'SEQUENCE_ID': 'MH1000',
                    'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG',
                    'SEQUENCE_INCLUDED_REGION': [36,342]
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
                    'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],
                })
        sleep(0.1)  # Pause for any GC
        em = _getMemUsage()
        print('\n\tMemory usage before 1k runs of designPrimers: ', sm)
        print('\tMemory usage after 1k runs of designPrimers:  ', em)
        print('\t\t\t\t\tDifference: \t', em-sm)
        if em-sm > 1000:
            raise AssertionError('Memory usage increase after 1k runs of \n\t'
                                 'designPrimers > 1000 bytes -- potential \n\t'
                                 'memory leak (mem increase: {})'.format(em-sm))

def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestDesignBindings())
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    test_suite = suite()
    runner.run(test_suite)
