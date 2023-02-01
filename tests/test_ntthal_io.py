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
tests.test_ntthal_io
~~~~~~~~~~~~~~~~~~~~~~~

Test to confirm primer3-py updates to the thal.c library result in the same
output given the input, `thal_input` and expected output, `thal_output` file
provided by the primer3 source library
'''
import io
import os
import os.path as op
import subprocess
import sys
import unittest
from typing import List

TEST_DIR = op.dirname(op.realpath(__file__))
PACKAGE_DIR = op.dirname(TEST_DIR)
PRIMER3_DIR = op.join(PACKAGE_DIR, 'primer3')
INPUT_DIR = op.join(TEST_DIR, 'input_files')
LIBPRIMER3_PATH = op.join(PRIMER3_DIR, 'src', 'libprimer3')
NTTHAL_PATH = op.join(LIBPRIMER3_PATH, 'ntthal')


def split_args_line_to_list(line: str) -> List[str]:
    '''
    Args:
        line: a line of arguments to split by space.  Filter out empty strings
            if there are extra spaces between arguments

    Returns:
        List of strings of split arguments
    '''
    args_line_list = line.split(' ')
    # filter out empty strings to prevent parsing errors
    return [x for x in args_line_list if x != '']


def run_ntthal(line: str, stderr_fd: io.BufferedWriter) -> str:
    '''
    Run ntthal given the argument line

    Args:
        line: A line of arguments for the ntthal command
        stderr_fd: Writter to send stderr to such that it does not go to stdout

    Returns:
        On success, the string of lines of the process stdout printing.
        On failure an empty string
    '''
    args_line_list = split_args_line_to_list(line)
    args = [NTTHAL_PATH] + args_line_list
    try:
        out_b = subprocess.check_output(
            args,
            stderr=stderr_fd,
            env=os.environ,
        )
        # Need to trim trailing white space as thal_output is trimmed by
        # pre-commit
        return out_b.decode('utf8')

    except subprocess.CalledProcessError:
        # NOTE: Intentionally commented out for debugging purposes
        # print(line)
        return ''


@unittest.skipIf(
    sys.platform == 'win32',
    'Windows does not support resource module and wrappers',
)
class TestLowLevelBindings(unittest.TestCase):
    def test_confirm_ntthal_binary(self):
        '''Confirm updates to thal.c do not change expected output

        NOTE: thal_input has several argument lines that are INTENTIONAL
        failures. Failure argument lines result in no output in thal_output and
        therefore there are fewer output line sets that input argument lines
        '''
        compute_str_list: List[str] = []
        with open(op.join(INPUT_DIR, 'thal_input'), 'r') as fd:
            input_str = fd.read()
        with open(op.join(INPUT_DIR, 'thal_output'), 'r') as fd:
            output_str = fd.read()
        with open(os.devnull, 'wb') as dev_null:
            for line in input_str.split('\n'):
                compute_str_list.append(run_ntthal(line, dev_null))
        compute_str = ''.join(compute_str_list)
        self.assertTrue(compute_str == output_str)
        # NOTE: Intentionally commented out for debugging purposes
        # if compute_str != output_str:
        #     # print(compute_str)
        #     with open('thal_input_compute_str', 'w') as fd:
        #             fd.write(compute_str)
        #     raise ValueError('')
