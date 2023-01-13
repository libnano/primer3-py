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

import primer3.wrappers as wrappers

P3_ARGS: Dict[str, Any] = {}

interval_list_tags = set([
    'SEQUENCE_INCLUDED_REGION',
    'SEQUENCE_TARGET',
    'SEQUENCE_EXCLUDED_REGION',
    'SEQUENCE_INTERNAL_EXCLUDED_REGION',
])

size_range_list_tags = set(['PRIMER_PRODUCT_SIZE_RANGE'])

#  (semicolon separated list of integer "quadruples"; default empty)
semi_quad_tags = ['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']

Sequence_T = Union[List[Any], Tuple[Any, ...]]
Str_Bytes_T = Union[str, bytes]


def wrapListOfQuads(v: Sequence_T) -> str:
    ''' Wrap a list of ordered quads, potentially containing None's
    Produces a string 'list of semicolon separated list of integer "quadruples"'
    Used for SEQUENCE_PRIMER_PAIR_OK_REGION_LIST

    >>> wrapListOfQuads([1,2,3,4])
    '1,2,3,4'

    >>> wrapListOfQuads([(1,2,3,4)])
    '1,2,3,4'

    >>> wrapListOfQuads([[1,2,3,4]])
    '1,2,3,4'

    >>> wrapListOfQuads([[1,2,3,4], [5,6,7,8]])
    '1,2,3,4 ; 5,6,7,8'

    >>> wrapListOfQuads(((1,2,3,4), [5,6,7,8]))
    '1,2,3,4 ; 5,6,7,8'

    Args:
        v: Sequence of items to wrap or Sequence of lists/tuples to wrap

    Returns:
        Specially string formatted sequence
    '''
    def int_to_str(i):
        return str(i) if i > -1 else ''

    try:
        rv = ';'.join(
            ','.join(map(int_to_str, quad))
            for quad in v
        )
    except TypeError:
        rv = ','.join(x and str(x) or '' for x in v)
    return rv


def _wrapListWithFormat(v: Sequence_T, fmt: str, sep: str = ' ') -> str:
    '''
    Specially format a list as a string

    Args:
        v: Sequence of items to format
        fmt: format string
        sep: Optional separator

    Returns:
        string formatted list
    '''
    try:
        rv = fmt % tuple(v)
    except TypeError:
        rv = sep.join(fmt % tuple(x) for x in v)
    return rv


def wrap(t: Tuple[str, Any]) -> Tuple[str, str]:
    '''Convert a primer3 input in python-friendly bindings-style form
    to a string form for use by the process wrapper

    >>> wrap(('SEQUENCE_TEMPLATE', 'ATCG'))
    ('SEQUENCE_TEMPLATE', 'ATCG')

    >>> wrap(('SEQUENCE_QUALITY', range(5)))
    ('SEQUENCE_QUALITY', '0 1 2 3 4')

    >>> wrap(('SEQUENCE_EXCLUDED_REGION', (5,7)))
    ('SEQUENCE_EXCLUDED_REGION', '5,7')

    >>> wrap(('SEQUENCE_EXCLUDED_REGION', [(5,7), (11,13)]))
    ('SEQUENCE_EXCLUDED_REGION', '5,7 11,13')

    >>> wrap(('PRIMER_PRODUCT_SIZE_RANGE', (7,11)))
    ('PRIMER_PRODUCT_SIZE_RANGE', '7-11')

    Args:
        t: Key, Value tuple where the value is correctly typed in python

    Returns
        Key, Value tuple where the value is always a string
    '''
    k, v = t

    if isinstance(v, (list, tuple)):
        if len(v) == 0:
            rv = ''
        else:
            if k in semi_quad_tags:
                rv = wrapListOfQuads(v)
            elif k in interval_list_tags:
                rv = _wrapListWithFormat(v, '%d,%d')
            elif k in size_range_list_tags:
                rv = _wrapListWithFormat(v, '%d-%d')
            elif isinstance(v[0], (list, tuple)):
                rv = _wrapListWithFormat(v, '%d-%d')
            elif isinstance(v[0], int):
                rv = ' '.join(map(str, v))
            else:
                rv = v  # type: ignore
    else:
        rv = v
    return k, rv


def unwrap(t: Tuple[str, str]) -> Tuple[str, Any]:
    '''convert a wrapper result into the intended form of the
    bindings result

    >>> unwrap(('A_FLOAT', '1.23'))
    ('A_FLOAT', 1.23)

    >>> unwrap(('AN_INT', '42'))
    ('AN_INT', 42)

    >>> unwrap(('A_SIZE_RANGE', '2020-2520'))
    ('A_SIZE_RANGE', (2020, 2520))

    >>> unwrap(('AN_INTERVAL', '7,11'))
    ('AN_INTERVAL', (7, 11))

    >>> unwrap(('MULTI_SIZE_RANGE', '1-2 3-5 7-11'))
    ('MULTI_SIZE_RANGE', ((1, 2), (3, 5), (7, 11)))

    >>> unwrap(('MULTIPLE_INTS', '2 3 5 7 11 13'))
    ('MULTIPLE_INTS', (2, 3, 5, 7, 11, 13))

    >>> unwrap(('A_SEQUENCE', 'ATCG'))
    ('A_SEQUENCE', 'ATCG')

    Args:
        t: Key, Value tuple where the value is always a string

    Returns
        Key, Value tuple where the value is correctly typed in python
    '''
    k, v = t
    rv = v

    def condInt(s):
        return int(s) if s != '' else -1

    for lam in [
        lambda x: int(x),
        lambda x: float(x),
        lambda x: tuple(int(s) for s in x.split(' ')),
        lambda x: tuple(int(s) for s in x.split('-')),
        lambda x: tuple(int(s) for s in x.split(',')),
        lambda x: tuple(
            tuple(int(ss) for ss in s.split('-'))
            for s in x.split(' ')
        ),
        lambda x: tuple(
            tuple(condInt(ss) for ss in s.split(','))
            for s in x.split(' ')
        ),
        lambda x: tuple(
            tuple(
                condInt(ss) for ss in
                s.strip().split(',')
            ) for s in x.split(';')
        ),
    ]:
        try:
            rv = lam(v)  # type: ignore
        except BaseException:
            pass
        else:
            break
    return k, rv


def setGlobals(
        global_args: Dict[str, Any],
        misprime_lib: Optional[Dict[Str_Bytes_T, Str_Bytes_T]],  # NOTE: unused
        mishyb_lib: Optional[Dict[Str_Bytes_T, Str_Bytes_T]],   # NOTE: unused
) -> None:
    '''Update global Primer3 global, mispriming, and mishybridization args in
    ``P3_ARGS``

    Args:
        seq_args:  dictionary of Primer3 global args
        misprime_lib: mispriming library dictionary
        mishyb_lib: mishybridization library dictionary

    Returns:
        None
    '''
    P3_ARGS.update(dict(wrap(v) for v in global_args.items()))


def setSeqArgs(seq_args: Dict[str, Any]) -> None:
    '''Update global Primer3 sequence args in ``P3_ARGS``

    Args:
        seq_args: dictionary of Primer3 sequence args

    Returns:
        None
    '''
    P3_ARGS.update(dict(wrap(v) for v in seq_args.items()))


def convertResult(result):
    return dict(unwrap(v) for v in result.items())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Design bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def designPrimers(
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
    `designPrimers` call or by a `setGlobals` call), `designPrimers` may be
    called with seqArgs alone (as a means of optimization).

    Args:
        seq_args: Primer3 sequence/design args as per Primer3 docs
        global_args: Primer3 global args as per Primer3 docs
        reset_args:
        misprime_lib: `Sequence name: sequence` dictionary for mispriming
            checks.
        mishyb_lib: `Sequence name: sequence` dictionary for mishybridization
            checks.
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
        setGlobals(global_args, misprime_lib, mishyb_lib)
    setSeqArgs(seq_args)
    result = wrappers.designPrimers(
        P3_ARGS,
        input_log=input_log,
        output_log=output_log,
        err_log=err_log,
    )
    return convertResult(result)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
