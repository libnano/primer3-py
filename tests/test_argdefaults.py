# Copyright (C) 2023. Ben Pruitt & Nick Conway
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
tests.test_argdefaults
~~~~~~~~~~~~~~~~~~~~~~~

Unit tests for the primer3-py argdefaults

'''
import difflib
import os.path as op
import unittest

from primer3 import argdefaults

TEST_DIR = op.dirname(op.realpath(__file__))
TEST_INPUT_DIR = op.join(TEST_DIR, 'input_files')
PACKAGE_DIR = op.dirname(TEST_DIR)


class TestDesignBindings(unittest.TestCase):

    def test_wrap_list_of_quads(self):
        '''
        test wrap_list_of_quads
        '''
        result = argdefaults.wrap_list_of_quads([1, 2, 3, 4])
        self.assertEqual(result, '1,2,3,4')

        result = argdefaults.wrap_list_of_quads([[1, 2, 3, 4]])
        self.assertEqual(result, '1,2,3,4')

        result = argdefaults.wrap_list_of_quads([[1, 2, 3, 4], [5, 6, 7, 8]])
        self.assertEqual(result, '1,2,3,4 ; 5,6,7,8')

        result = argdefaults.wrap_list_of_quads(((1, 2, 3, 4), [5, 6, 7, 8]))
        self.assertEqual(result, '1,2,3,4 ; 5,6,7,8')

        result = argdefaults.wrap_list_of_quads((5, 6, -1, -1))
        print('woof', result)
        self.assertEqual(result, '5,6,,')

        result = argdefaults.wrap_list_of_quads(((1, 2, 3, 4), [5, 6, -1, -1]))
        self.assertEqual(result, '1,2,3,4 ; 5,6,,')

    def test_unwrap(self):
        '''
        Make sure unwrap works correctly
        '''
        _, result = argdefaults.unwrap('', '1,2,3,4 ; 5,6,,')
        self.assertEqual(result, ((1, 2, 3, 4), (5, 6, -1, -1)))

    def test_roundtrip(self):
        '''
        Assert roundtrip conversion to dictionary works
        '''
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
        for fn_root in test_file_roots:
            fp = op.join(TEST_INPUT_DIR, f'{fn_root}_input')
            with open(fp, 'r') as fd:
                boulder_str_init = fd.read()

            data_dict_list = argdefaults.parse_multirecord_boulder_io(
                boulder_str_init,
            )
            boulder_str_out = b''
            for data_dict in data_dict_list:
                boulder_str_out += argdefaults.format_boulder_io(
                    data_dict,
                )
            boulder_str_out = boulder_str_out.decode('utf8')
            # Double check that the differences are not just white space
            # and floating point `0` characters
            if boulder_str_init != boulder_str_out:
                for s in difflib.ndiff(boulder_str_init, boulder_str_out):
                    if s[0] == ' ':
                        continue
                    elif s[0] == '-':
                        if s[-1] in (' ', '0'):
                            continue
                    elif s[0] == '+':
                        if s[-1] in (' ', '0'):
                            continue
                    else:
                        # write to file to help with inspection
                        with open('boulder_str_init', 'w') as fd:
                            fd.write(boulder_str_init)
                        with open('boulder_str_out', 'w') as fd:
                            fd.write(boulder_str_out)
                        raise ValueError(f'File difference for {fn_root}')
