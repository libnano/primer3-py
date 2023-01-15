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
primer3.wrappers
~~~~~~~~~~~~~~~~

Simple subprocess wrappers for Primer3 executables. These functions closely
mirror the functions found in bindings.py, but are much slower and should
only be used for testing / comparison purposes.

'''

from __future__ import print_function

import io
import os
import re
import subprocess
from collections import (
    OrderedDict,
    namedtuple,
)
from os.path import join as pjoin
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
    Union,
)

from .argdefaults import Primer3PyArguments

DEFAULT_P3_ARGS = Primer3PyArguments()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PRIMER3 WRAPPERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~ Check to insure that the environment is properly configured ~~~~~~~ #

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

if not os.environ.get('PRIMER3HOME'):
    try:
        os.environ['PRIMER3HOME'] = pjoin(LOCAL_DIR, 'src/libprimer3')
    except BaseException:
        raise OSError('PRIMER3HOME environmental variable is not set.')
LIBPRIMER3_PATH = os.environ['PRIMER3HOME']
if not os.path.isdir(LIBPRIMER3_PATH):
    raise OSError(f'Path {LIBPRIMER3_PATH} does not exist')

THERMO_PATH = pjoin(LIBPRIMER3_PATH, 'primer3_config/')

DEV_NULL = open(os.devnull, 'wb')

_tm_methods = {
    'breslauer': 0,
    'santalucia': 1,
}

_salt_corrections_methods = {
    'schildkraut': 0,
    'santalucia': 1,
    'owczarzy': 2,
}


def calcTm(
        seq: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        max_nn_length: int = DEFAULT_P3_ARGS.max_nn_length,
        tm_method: str = DEFAULT_P3_ARGS.tm_method,
        salt_corrections_method: str = DEFAULT_P3_ARGS.salt_corrections_method,
) -> float:
    ''' Return the tm of `seq` as a float.

    Args:
        seq: DNA sequence
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        max_nn_length: Maximum length for nearest-neighbor calcs
        tm_method: Tm calculation method (breslauer or santalucia)
        salt_corrections_method: Salt correction method (schildkraut, owczarzy,
            santalucia)

    Returns:
        The melting temperature in degrees Celsius (float).
    '''
    tm_meth = _tm_methods.get(tm_method)
    if tm_meth is None:
        raise ValueError(
            '{} is not a valid tm calculation method'.format(
                tm_method,
            ),
        )
    salt_meth = _salt_corrections_methods.get(salt_corrections_method)
    if salt_meth is None:
        raise ValueError(
            '{} is not a valid salt correction method'.format(
                salt_corrections_method,
            ),
        )
    # For whatever reason mv_conc and dna_conc have to be ints
    args = [
        pjoin(LIBPRIMER3_PATH, 'oligotm'),
        '-mv', str(mv_conc),
        '-dv', str(dv_conc),
        '-n', str(dntp_conc),
        '-d', str(dna_conc),
        '-tp', str(tm_meth),
        '-sc', str(salt_meth),
        seq,
    ]
    tm = subprocess.check_output(
        args, stderr=DEV_NULL,
        env=os.environ,
    )
    return float(tm)


_NTTHAL_RE = re.compile(
    r'dS\s+=\s+(\S+)\s+'
    r'dH\s+=\s+(\S+)\s+'
    r'dG\s+=\s+(\S+)\s+'
    r't\s+=\s+(\S+)'.encode('utf8'),
)

THERMORESULT = namedtuple(  # type: ignore
    'thermoresult', [
        'result',           # True if a structure is present
        'ds',               # Entropy (cal/(K*mol))
        'dh',               # Enthalpy (kcal/mol)
        'dg',               # Gibbs free energy
        'tm',               # Melting temperature (deg. Celsius)
        'ascii_structure',   # ASCII representation of structure
    ],
)

NULLTHERMORESULT = THERMORESULT(False, 0, 0, 0, 0, '')


def _parse_ntthal(ntthal_output: bytes) -> THERMORESULT:
    ''' Helper method that uses regex to parse ntthal output. '''
    parsed_vals = re.search(_NTTHAL_RE, ntthal_output)
    if parsed_vals:
        ascii_structure = (
            ntthal_output[ntthal_output.index(b'\n') + 1:].decode('utf8')
        )
        res = THERMORESULT(
            True,                           # Structure found
            float(parsed_vals.group(1)),    # dS
            float(parsed_vals.group(2)),    # dH
            float(parsed_vals.group(3)),    # dG
            float(parsed_vals.group(4)),    # tm
            ascii_structure,
        )
    else:
        res = NULLTHERMORESULT
    return res


def calcThermo(
        seq1: str,
        seq2: str,
        calc_type: str = DEFAULT_P3_ARGS.calc_type_wrapper,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        temp_only: Union[bool, int] = DEFAULT_P3_ARGS.temp_only,
) -> THERMORESULT:
    '''Main subprocess wrapper for calls to the ntthal executable.

    Args:
        seq1: DNA sequence to analyze for 3' end  hybridization against the
            target sequence
        seq2: Target DNA sequence to analyze for seq1 3' end hybridization
        calc_type: alignment type, END1, END2, ANY and HAIRPIN, by default ANY
            (when duplex)
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop: Maximum size of loops in the structure
        temp_only: causes the alignment NOT to be displayed on stderr,
            _only_ Tm is printed

    Returns:
        a named tuple with tm, ds, dh, and dg values or None if no
        structure / complex could be computed.
    '''
    args = [
        pjoin(LIBPRIMER3_PATH, 'ntthal'),
        '-a', str(calc_type),
        '-mv', str(mv_conc),
        '-dv', str(dv_conc),
        '-n', str(dntp_conc),
        '-d', str(dna_conc),
        '-t', str(temp_c),
        '-maxloop', str(max_loop),
        '-path', THERMO_PATH,
        '-s1', seq1,
        '-s2', seq2,
    ]
    if temp_only:
        args += ['-r']
    out = subprocess.check_output(
        args, stderr=DEV_NULL,
        env=os.environ,
    )
    return _parse_ntthal(out)


def calcHairpin(
        seq: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        temp_only: Union[bool, int] = DEFAULT_P3_ARGS.temp_only,
) -> THERMORESULT:
    ''' Return a namedtuple of the dS, dH, dG, and Tm of any hairpin struct
    present.

    Args:
        seq: DNA sequence to analyze for hairpin formation
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop(int, optional): Maximum size of loops in the structure
        temp_only:
        temp_only: print only temp to stderr

    Returns:
       ``THERMORESULT`` tuple
    '''
    return calcThermo(
        seq,
        seq,
        'HAIRPIN',
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        temp_only,
    )


def calcHeterodimer(
        seq1: str,
        seq2: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        temp_only: Union[bool, int] = DEFAULT_P3_ARGS.temp_only,
) -> THERMORESULT:
    '''Returns a tuple of the dS, dH, dG, and Tm of any predicted heterodimer.

    Args:
        seq1: DNA sequence to analyze for heterodimer formation
        seq2: DNA sequence to analyze for heterodimer formation
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop(int, optional): Maximum size of loops in the structure
        temp_only:
        temp_only: print only temp to stderr

    Returns:
       ``THERMORESULT`` tuple
    '''
    return calcThermo(
        seq1,
        seq2,
        'ANY',
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        temp_only,
    )


def calcHomodimer(
        seq: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        temp_only: Union[bool, int] = DEFAULT_P3_ARGS.temp_only,
) -> THERMORESULT:
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted homodimer.

    Args:
        seq: DNA sequence to analyze for homodimer formation
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop(int, optional): Maximum size of loops in the structure
        temp_only:
        temp_only: print only temp to stderr

    Returns:
       ``THERMORESULT`` tuple
    '''
    return calcThermo(
        seq,
        seq,
        'ANY',
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        temp_only,
    )


def calcEndStability(
        seq1: str,
        seq2: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        temp_only: Union[bool, int] = DEFAULT_P3_ARGS.temp_only,
) -> THERMORESULT:
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted heterodimer
    end stability

    Args:
        seq1: DNA sequence to analyze for end stability
        seq2: DNA sequence to analyze for end stability
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop(int, optional): Maximum size of loops in the structure
        temp_only:
        temp_only: print only temp to stderr

    Returns:
       ``THERMORESULT`` tuple
    '''
    return calcThermo(
        seq1,
        seq2,
        'END1',
        mv_conc,
        dv_conc,
        dntp_conc,
        dna_conc,
        temp_c,
        max_loop,
        temp_only,
    )


def assessOligo(seq: str) -> Tuple[THERMORESULT, THERMORESULT]:
    '''
    Return the thermodynamic characteristics of hairpin/homodimer structures.

    Args:
        seq: DNA sequence to analyze

    Returns:
        A tuple of namedtuples (hairpin data, homodimer data) in which each
        individual tuple is structured (dS, dH, dG, Tm).

    '''
    hairpin_out = calcHairpin(seq)
    homodimer_out = calcHomodimer(seq)
    return (hairpin_out, homodimer_out)


# ~~~~~~~ RUDIMENTARY PRIMER3 MAIN WRAPPER (see Primer3 docs for args) ~~~~~~ #
def _formatBoulderIO(p3_args: Dict[str, Any]) -> bytes:
    '''Convert argument dictionary to boulder formatted bytes

    Args:
        p3_args: primer3 arguments to format boulder style

    Returns:
        Boulder formatted byte string
    '''
    boulder_str = ''.join([
        '{}={}\n'.format(k, v) for k, v in
        p3_args.items()
    ])
    boulder_str = f'{boulder_str}=\n'
    return boulder_str.encode('utf8')


def _parseBoulderIO(boulder_bytes: bytes) -> Dict[str, str]:
    '''Convert boulder info to a key/value dictionary
    Args:
        boulder_bytes: Bytes of boulder formatted information to parse

    Returns:
        Dictionary of key/values
    '''
    data_dict = OrderedDict()
    for line in boulder_bytes.decode('utf8').split('\n'):
        try:
            k, v = line.strip().split('=')
            data_dict[k] = v
        except ValueError:
            pass
    return data_dict


def parseMultiRecordBoulderIO(boulder_str: str) -> List[OrderedDict[str, str]]:
    '''
    Args:
        boulder_str: boulder string to parse with multiple records

    Returns:
        List of OrderedDicts per record
    '''
    data_dicts = []
    for record in re.split('=\r?\n', boulder_str):
        if record == '':
            continue
        data_dict = OrderedDict()
        for line in record.split('\n'):
            try:
                k, v = line.strip().split('=')
                data_dict[k] = v
            except ValueError:
                pass
        data_dicts.append(data_dict)
    return data_dicts


def designPrimers(
        p3_args: Dict[str, Any],
        input_log: Optional[io.BufferedWriter] = None,
        output_log: Optional[io.BufferedWriter] = None,
        err_log: Optional[io.BufferedWriter] = None,
) -> Dict[str, Any]:
    ''' Return the raw primer3_core output for the provided primer3 args.

    Args:
        p3_args: Dictionarty of arguments for the primer3
        input_log: Optional log input file descriptor
        output_log: Optional log output file descriptor
        err_log: Optional log error file descriptor

    Returns:
        an ordered dict of the boulderIO-format primer3 output file
    '''
    sp = subprocess.Popen(
        [pjoin(LIBPRIMER3_PATH, 'primer3_core')],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    p3_args.setdefault(
        'PRIMER_THERMODYNAMIC_PARAMETERS_PATH',
        pjoin(LIBPRIMER3_PATH, 'primer3_config/'),
    )
    in_str = _formatBoulderIO(p3_args)
    if input_log:
        input_log.write(in_str)
        input_log.flush()
    out_str, err_str = sp.communicate(input=in_str)
    if output_log:
        output_log.write(out_str)
        output_log.flush()
    if err_log and err_str is not None:
        err_log.write(err_str)
        err_log.flush()
    return _parseBoulderIO(out_str)
