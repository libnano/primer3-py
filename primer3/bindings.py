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
primer3.bindings
~~~~~~~~~~~~~~~~

This module provides a simple API for the Primer3 primer design / thermodynamic
calculations library.

These are direct bindings to an optimized version of the Primer3 C library,
as opposed to the more commonly used subprocess-based wrappers (we provide a
set of wrappers for comparison / testing purposes as well).

Note that this module effectively abstracts the C API / Cython bindings for
the primer design and thermodynamic analysis functionality of Primer3. This is
done primarly to provide a clean, consistent interface. For applications with
stringent performance requirments, you should consider using the C API
and/or Cython modules directly. See the docs for more details.

'''
from typing import (
    Any,
    Dict,
    Optional,
    Union,
)

from . import (  # type: ignore
    primerdesign,
    thermoanalysis,
)
from .argdefaults import Primer3PyArguments

DEFAULT_P3_ARGS = Primer3PyArguments()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Low level bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


Str_Bytes_T = Union[str, bytes]

_THERMO_ANALYSIS = thermoanalysis.ThermoAnalysis()


def _setThermoArgs(
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        dmso_conc: float = DEFAULT_P3_ARGS.dmso_conc,
        dmso_fact: float = DEFAULT_P3_ARGS.dmso_fact,
        formamide_conc: float = DEFAULT_P3_ARGS.formamide_conc,
        annealing_temp_c: float = DEFAULT_P3_ARGS.annealing_temp_c,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_nn_length: int = DEFAULT_P3_ARGS.max_nn_length,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        tm_method: str = DEFAULT_P3_ARGS.tm_method,
        salt_corrections_method: str = DEFAULT_P3_ARGS.salt_corrections_method,
        **kwargs,
):
    '''
    Set parameters in global ``ThermoAnalysis`` instance
    Args:
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        dmso_conc: Concentration of DMSO (%)
        dmso_fact: DMSO correction factor, default 0.6
        formamide_conc: Concentration of formamide (mol/l)
        annealing_temp_c: Actual annealing temperature of the PCR reaction
            in (C)
        temp_c: Simulation temperature for dG (Celsius)
        max_nn_length: Maximum length for nearest-neighbor calcs
        tm_method: Tm calculation method (breslauer or santalucia)
        salt_corrections_method: Salt correction method (schildkraut, owczarzy,
            santalucia)
    '''
    _THERMO_ANALYSIS.mv_conc = float(mv_conc)
    _THERMO_ANALYSIS.dv_conc = float(dv_conc)
    _THERMO_ANALYSIS.dntp_conc = float(dntp_conc)
    _THERMO_ANALYSIS.dna_conc = float(dna_conc)

    _THERMO_ANALYSIS.dmso_conc = float(dmso_conc)
    _THERMO_ANALYSIS.dmso_fact = float(dmso_fact)
    _THERMO_ANALYSIS.formamide_conc = float(formamide_conc)
    _THERMO_ANALYSIS.annealing_temp_c = float(annealing_temp_c)

    _THERMO_ANALYSIS.temp = float(temp_c)
    _THERMO_ANALYSIS.max_loop = int(max_loop)
    _THERMO_ANALYSIS.tm_method = tm_method
    _THERMO_ANALYSIS.salt_correction_method = salt_corrections_method

    _THERMO_ANALYSIS.max_nn_length = int(max_nn_length)


def calcHairpin(
        seq: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        output_structure: bool = False,
):
    ''' Calculate the hairpin formation thermodynamics of a DNA sequence.

    **Note that the maximum length of `seq` is 60 bp.** This is a cap suggested
    by the Primer3 team as the longest reasonable sequence length for which
    a two-state NN model produces reliable results
    (see primer3/src/libnano/thal.h:50).

    Args:
        seq: DNA sequence to analyze for hairpin formation
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop(int, optional): Maximum size of loops in the structure
        output_structure (bool) : If `True`, the ASCII dimer structure is saved

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        hairpin formation.

    Raises:
        ``RuntimeError``

    '''
    _THERMO_ANALYSIS.set_thermo_args(**locals())
    return _THERMO_ANALYSIS.calcHairpin(seq, output_structure).checkExc()


def calcHomodimer(
        seq: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        output_structure: bool = False,
):
    ''' Calculate the homodimerization thermodynamics of a DNA sequence.

    **Note that the maximum length of ``seq`` is 60 bp.** This is a cap imposed
    by Primer3 as the longest reasonable sequence length for which
    a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).

    Args:
        seq: DNA sequence to analyze for homodimer formation calculations
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop: Maximum size of loops in the structure
        output_structure: If `True`, the ASCII dimer structure
            is saved

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        homodimer interaction.

    Raises:
        ``RuntimeError``

    '''
    _THERMO_ANALYSIS.set_thermo_args(**locals())
    return _THERMO_ANALYSIS.calcHomodimer(seq, output_structure).checkExc()


def calcHeterodimer(
        seq1: str,
        seq2: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
        output_structure: bool = False,
):
    ''' Calculate the heterodimerization thermodynamics of two DNA sequences.

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).

    Args:
        seq1: First DNA sequence to analyze for heterodimer formation
        seq2: Second DNA sequence to analyze for heterodimer formation
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop: Maximum size of loops in the structure
        output_structure: If `True`, the ASCII dimer structure is saved

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        heterodimer interaction.

    Raises:
        ``RuntimeError``

    '''
    _THERMO_ANALYSIS.set_thermo_args(**locals())
    return _THERMO_ANALYSIS.calcHeterodimer(
        seq1,
        seq2,
        output_structure,
    ).checkExc()


def calcEndStability(
        seq1: str,
        seq2: str,
        mv_conc: Union[float, int] = DEFAULT_P3_ARGS.mv_conc,
        dv_conc: Union[float, int] = DEFAULT_P3_ARGS.dv_conc,
        dntp_conc: Union[float, int] = DEFAULT_P3_ARGS.dntp_conc,
        dna_conc: Union[float, int] = DEFAULT_P3_ARGS.dna_conc,
        temp_c: Union[float, int] = DEFAULT_P3_ARGS.temp_c,
        max_loop: int = DEFAULT_P3_ARGS.max_loop,
) -> thermoanalysis.ThermoResult:
    ''' Calculate the 3' end stability of DNA sequence `seq1` against DNA
    sequence `seq2`.

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).

    Args:
        seq1: DNA sequence to analyze for 3' end  hybridization against the
            target sequence
        seq2: Target DNA sequence to analyze for seq1 3' end hybridization
        mv_conc: Monovalent cation conc. (mM)
        dv_conc: Divalent cation conc. (mM)
        dntp_conc: dNTP conc. (mM)
        dna_conc: DNA conc. (nM)
        temp_c: Simulation temperature for dG (Celsius)
        max_loop: Maximum size of loops in the structure

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        3' hybridization interaction.

    Raises:
        ``RuntimeError``

    '''
    _THERMO_ANALYSIS.set_thermo_args(**locals())
    return _THERMO_ANALYSIS.calcEndStability(seq1, seq2).checkExc()


def calcTm(
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
    ''' Calculate the melting temperature (Tm) of a DNA sequence.

    Note that NN thermodynamics will be used to calculate the Tm of sequences
    up to 60 bp in length, after which point the following formula will be
    used::

        Tm = 81.5 + 16.6(log10([mv_conc])) + 0.41(%GC) - 600/length

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
        temp_c: Simulation temperature for dG (Celsius)
        max_nn_length: Maximum length for nearest-neighbor calcs
        tm_method: Tm calculation method (breslauer or santalucia)
        salt_corrections_method: Salt correction method (schildkraut, owczarzy,
            santalucia)

    Returns:
        The melting temperature in degrees Celsius (float).
    '''
    _THERMO_ANALYSIS.set_thermo_args(**locals())
    return _THERMO_ANALYSIS.calcTm(seq)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tm-only aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def calcHairpinTm(*args, **kwargs) -> float:
    return calcHairpin(*args, **kwargs).tm


def calcHomodimerTm(*args, **kwargs) -> float:
    return calcHomodimer(*args, **kwargs).tm


def calcHeterodimerTm(*args, **kwargs) -> float:
    return calcHeterodimer(*args, **kwargs).tm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Design bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def design_load_thermo_params():
    primerdesign.load_thermo_params()


def designPrimers(
        seq_args: Dict[str, Any],
        global_args: Dict[str, Any],
        misprime_lib: Optional[Dict[Str_Bytes_T, Str_Bytes_T]] = None,
        mishyb_lib: Optional[Dict[Str_Bytes_T, Str_Bytes_T]] = None,
) -> Dict[str, Any]:
    '''Run the Primer3 design process.

    If the global args have been previously set (either by a pervious
    `designPrimers` call or by a `setGlobals` call), `designPrimers` may be
    called with seqArgs alone (as a means of optimization).

    Args:
        seq_args: Primer3 sequence/design args as per Primer3 docs
        global_args: Primer3 global args as per Primer3 docs
        misprime_lib: `Sequence name: sequence` dictionary for mispriming
            checks.
        mishyb_lib: `Sequence name: sequence` dictionary for mishybridization
            checks.

    Returns:
        A dictionary of Primer3 results (should be identical to the expected
        BoulderIO output from primer3_main)

    '''
    primerdesign.set_globals_and_seq_args(
        global_args=global_args,
        seq_args=seq_args,
        misprime_lib=misprime_lib,
        mishyb_lib=mishyb_lib,
    )
    return primerdesign.run_design()
