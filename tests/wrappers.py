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
tests.wrappers
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
from collections import namedtuple
from os.path import join as pjoin
from typing import (
    Any,
    Dict,
    Optional,
    Tuple,
    Union,
)

from primer3.argdefaults import (
    Primer3PyArguments,
    format_boulder_io,
    parse_boulder_io,
)

DEFAULT_P3_ARGS = Primer3PyArguments()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PRIMER3 WRAPPERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# TODO: remove after update to primer3 >= 2.5.0
LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
PACKAGE_DIR = os.path.dirname(LOCAL_DIR)
LIBPRIMER3_PATH = pjoin(PACKAGE_DIR, 'primer3', 'src', 'libprimer3')
THERMO_PATH = pjoin(
    LIBPRIMER3_PATH,
    'primer3_config',
    '',  # Add trailing slash (OS-ind) req'd by primer3 lib
)

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


def calc_tm(
        seq: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        dmso_conc: float = DEFAULT_P3_ARGS.dmso_conc,
        dmso_fact: float = DEFAULT_P3_ARGS.dmso_fact,
        formamide_conc: float = DEFAULT_P3_ARGS.formamide_conc,
        annealing_temp_c: float = DEFAULT_P3_ARGS.annealing_temp_c,
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
        dmso_conc: Concentration of DMSO (%)
        dmso_fact: DMSO correction factor, default 0.6
        formamide_conc: Concentration of formamide (mol/l)
        annealing_temp_c: Actual annealing temperature of the PCR reaction
            in (C)
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
            f'{tm_method} is not a valid tm calculation method',
        )
    salt_meth = _salt_corrections_methods.get(salt_corrections_method)
    if salt_meth is None:
        raise ValueError(
            f'{salt_corrections_method} is not a valid salt correction method',
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
        '-dm', str(dmso_conc),
        '-df', str(dmso_fact),
        '-fo', str(formamide_conc),
        seq,
    ]
    tm = subprocess.check_output(
        args, stderr=DEV_NULL,
        env=os.environ,
    )
    return float(tm)


calcTm = calc_tm


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


def calc_thermo(
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


def calc_hairpin(
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
        temp_only: print only temp to stderr

    Returns:
       ``THERMORESULT`` tuple
    '''
    return calc_thermo(
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


calcHairpin = calc_hairpin


def calc_heterodimer(
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
    return calc_thermo(
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


calcHeterodimer = calc_heterodimer


def calc_homodimer(
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
    return calc_thermo(
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


calcHomodimer = calc_homodimer


def calc_end_stability(
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
    return calc_thermo(
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


calcEndStability = calc_end_stability


def assess_oligo(seq: str) -> Tuple[THERMORESULT, THERMORESULT]:
    '''
    Return the thermodynamic characteristics of hairpin/homodimer structures.

    Args:
        seq: DNA sequence to analyze

    Returns:
        A tuple of namedtuples (hairpin data, homodimer data) in which each
        individual tuple is structured (dS, dH, dG, Tm).

    '''
    hairpin_out = calc_hairpin(seq)
    homodimer_out = calc_homodimer(seq)
    return (hairpin_out, homodimer_out)


def design_primers(
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
    p3_args.setdefault('PRIMER_THERMODYNAMIC_PARAMETERS_PATH', THERMO_PATH)
    in_str = format_boulder_io(p3_args)
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
    return parse_boulder_io(out_str)


designPrimers = design_primers
