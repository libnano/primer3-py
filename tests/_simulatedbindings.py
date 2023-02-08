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
tests._simulatedbindings
~~~~~~~~~~~~~~~~~~~~~~~
'''
import io
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
    Union,
)

from primer3 import argdefaults

from . import wrappers

P3_ARGS: Dict[str, Any] = {}

Sequence_T = Union[List[Any], Tuple[Any, ...]]
Str_Bytes_T = Union[str, bytes]


def convert_result(result):
    unwrap = argdefaults.unwrap
    return dict(unwrap(k, v) for k, v in result.items())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Design bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def design_primers(
        seq_args: Dict[str, Any],
        global_args: Optional[Dict[str, Any]] = None,
        reset_args: bool = True,
        misprime_lib: Optional[Dict[Str_Bytes_T, Str_Bytes_T]] = None,
        mishyb_lib: Optional[Dict[Str_Bytes_T, Str_Bytes_T]] = None,
        input_log: Optional[io.BufferedWriter] = None,
        output_log: Optional[io.BufferedWriter] = None,
        err_log: Optional[io.BufferedWriter] = None,
) -> Dict[str, Any]:
    ''' Run the Primer3 design process, with the same interface as the bindings,
    using the wrapped subprocess of primer3_core to do the work.

    If the global args have been previously set (either by a pervious
    `design_primers` call or by a `set_globals` call), `design_primers` may be
    called with seqArgs alone (as a means of optimization).

    Args:
        seq_args: Primer3 sequence/design args as per Primer3 docs
        global_args: Primer3 global args as per Primer3 docs
        reset_args:
        misprime_lib: `Sequence name: sequence` dictionary for mispriming
            checks.  NOTE: Unused as of now
        mishyb_lib: `Sequence name: sequence` dictionary for mishybridization
            checks.  NOTE: Unused as of now
        input_log: Optional log input file descriptor
        output_log: Optional log output file descriptor
        err_log: Optional log error file descriptor

    Returns:
        A dictionary of Primer3 results (should be identical to the expected
        BoulderIO output from primer3_main)

    '''
    if reset_args:
        P3_ARGS.clear()

    if global_args:
        P3_ARGS.update(global_args)
    P3_ARGS.update(seq_args)

    result = wrappers.design_primers(
        P3_ARGS,
        input_log=input_log,
        output_log=output_log,
        err_log=err_log,
    )
    return convert_result(result)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
