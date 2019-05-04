# Copyright (C) 2014-2018. Ben Pruitt & Nick Conway; Wyss Institute
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

import os

from os.path import join as pjoin

from . import thermoanalysis
from . import primerdesign


# ~~~~~~~ Check to insure that the environment is properly configured ~~~~~~~ #

PRIMER3_HOME = os.environ['PRIMER3HOME']

# ~~~~~~~~~~~~~~~~ Load thermodynamic parameters into memory ~~~~~~~~~~~~~~~~ #

primerdesign.loadThermoParams(pjoin(PRIMER3_HOME, 'primer3_config/'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Low level bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


_THERMO_ANALYSIS = thermoanalysis.ThermoAnalysis()

def _setThermoArgs(mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
                   temp_c=37, max_loop=30, tm_method='santalucia',
                   salt_corrections_method='santalucia', **kwargs):
    _THERMO_ANALYSIS.mv_conc = float(mv_conc)
    _THERMO_ANALYSIS.dv_conc = float(dv_conc)
    _THERMO_ANALYSIS.dntp_conc = float(dntp_conc)
    _THERMO_ANALYSIS.dna_conc = float(dna_conc)
    _THERMO_ANALYSIS.temp = float(temp_c)
    _THERMO_ANALYSIS.max_loop = float(max_loop)
    _THERMO_ANALYSIS.tm_method = tm_method
    _THERMO_ANALYSIS.salt_correction_method = salt_corrections_method


def calcHairpin(seq, mv_conc=50.0, dv_conc=0.0, dntp_conc=0.8, dna_conc=50.0,
                temp_c=37, max_loop=30, output_structure=False):
    ''' Calculate the hairpin formation thermodynamics of a DNA sequence.

    **Note that the maximum length of `seq` is 60 bp.** This is a cap suggested
    by the Primer3 team as the longest reasonable sequence length for which
    a two-state NN model produces reliable results (see primer3/src/libnano/thal.h:50).

    Args:
        seq (str): DNA sequence to analyze for hairpin formation

        mv_conc (float/int, optional): Monovalent cation conc. (mM)
        dv_conc (float/int, optional): Divalent cation conc. (mM)
        dntp_conc (float/int, optional): dNTP conc. (mM)
        dna_conc (float/int, optional): DNA conc. (nM)
        temp_c (int, optional): Simulation temperature for dG (Celsius)
        max_loop(int, optional): Maximum size of loops in the structure
        output_structure (bool) : If `True`, the ASCII dimer structure is saved

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        hairpin formation.

    Raises:
        ``RuntimeError``

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcHairpin(seq, output_structure).checkExc()


def calcHomodimer(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
                  temp_c=37, max_loop=30, output_structure=False):
    ''' Calculate the homodimerization thermodynamics of a DNA sequence.

    **Note that the maximum length of ``seq`` is 60 bp.** This is a cap imposed
    by Primer3 as the longest reasonable sequence length for which
    a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).

    Args:
        seq (str)                       : DNA sequence to analyze for homodimer
                                          formation calculations

        mv_conc (float/int, optional)   : Monovalent cation conc. (mM)
        dv_conc (float/int, optional)   : Divalent cation conc. (mM)
        dntp_conc (float/int, optional) : dNTP conc. (mM)
        dna_conc (float/int, optional)  : DNA conc. (nM)
        temp_c (int, optional)          : Simulation temperature for dG (C)
        max_loop (int, optional)        : Maximum size of loops in the
                                          structure
        output_structure (bool) : If `True`, the ASCII dimer structure is saved

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        homodimer interaction.

    Raises:
        ``RuntimeError``

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcHomodimer(seq, output_structure).checkExc()


def calcHeterodimer(seq1, seq2, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                    dna_conc=50, temp_c=37, max_loop=30,
                    output_structure=False):
    ''' Calculate the heterodimerization thermodynamics of two DNA sequences.

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).

    Args:
        seq1 (str)              : First DNA sequence to analyze for heterodimer
                                  formation
        seq2 (str)              : Second DNA sequence to analyze for
                                  heterodimer formation

        mv_conc (float/int)     : Monovalent cation conc. (mM)
        dv_conc (float/int)     : Divalent cation conc. (mM)
        dntp_conc (float/int)   : dNTP conc. (mM)
        dna_conc (float/int)    : DNA conc. (nM)
        temp_c (int)            : Simulation temperature for dG (Celsius)
        max_loop(int)           : Maximum size of loops in the structure
        output_structure (bool) : If `True`, the ASCII dimer structure is saved

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        heterodimer interaction.

    Raises:
        ``RuntimeError``

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcHeterodimer(seq1, seq2, output_structure).checkExc()


def calcEndStability(seq1, seq2, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                     dna_conc=50, temp_c=37, max_loop=30):
    ''' Calculate the 3' end stability of DNA sequence `seq1` against DNA
    sequence `seq2`.

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).

    Args:
        seq1 (str)                        : DNA sequence to analyze for 3' end
                                            hybridization against the target
                                            sequence
        seq2 (str)                        : Target DNA sequence to analyze for
                                            seq1 3' end hybridization

        mv_conc (float/int, optional)     : Monovalent cation conc. (mM)
        dv_conc (float/int, optional)     : Divalent cation conc. (mM)
        dntp_conc (float/int, optional)   : dNTP conc. (mM)
        dna_conc (float/int, optional)    : DNA conc. (nM)
        temp_c (int, optional)            : Simulation temperature for dG (C)
        max_loop(int, optional)           : Maximum size of loops in the
                                            structure

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        3' hybridization interaction.

    Raises:
        ``RuntimeError``

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcEndStability(seq1, seq2).checkExc()


def calcTm(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
           max_nn_length=60, tm_method='santalucia',
           salt_corrections_method='santalucia'):
    ''' Calculate the melting temperature (Tm) of a DNA sequence.

    Note that NN thermodynamics will be used to calculate the Tm of sequences
    up to 60 bp in length, after which point the following formula will be
    used::

        Tm = 81.5 + 16.6(log10([mv_conc])) + 0.41(%GC) - 600/length

    Args:
        seq (str)                               : DNA sequence
        mv_conc (float/int, optional)           : Monovalent cation conc. (mM)
        dv_conc (float/int, optional)           : Divalent cation conc. (mM)
        dntp_conc (float/int, optional)         : dNTP conc. (mM)
        dna_conc (float/int, optional)          : DNA conc. (nM)
        max_nn_length (int, optional)           : Maximum length for
                                                  nearest-neighbor calcs
        tm_method (str, optional)               : Tm calculation method
                                                  (breslauer or santalucia)
        salt_corrections_method (str, optional) : Salt correction method
                                                  (schildkraut, owczarzy,
                                                  santalucia)

    Returns:
        The melting temperature in degrees Celsius (float).

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcTm(seq)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tm-only aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

calcHairpinTm = lambda *args, **kwargs: calcHairpin(*args, **kwargs).tm
calcHomodimerTm = lambda *args, **kwargs: calcHomodimer(*args, **kwargs).tm
calcHeterodimerTm = lambda *args, **kwargs: calcHeterodimer(*args, **kwargs).tm

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Design bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def designPrimers(seq_args, global_args=None, misprime_lib=None,
                  mishyb_lib=None, debug=False):
    ''' Run the Primer3 design process.

    If the global args have been previously set (either by a pervious
    `designPrimers` call or by a `setGlobals` call), `designPrimers` may be
    called with seqArgs alone (as a means of optimization).

    Args:
        seq_args (dict)               : Primer3 sequence/design args as per
                                        Primer3 docs

        global_args (dict, optional)  : Primer3 global args as per Primer3 docs
        misprime_lib (dict, optional) : `Sequence name: sequence` dictionary
                                        for mispriming checks.
        mishyb_lib (dict, optional)   : `Sequence name: sequence` dictionary
                                        for mishybridization checks.

    Returns:
        A dictionary of Primer3 results (should be identical to the expected
        BoulderIO output from primer3_main)

    '''
    if global_args:
        primerdesign.setGlobals(global_args, misprime_lib, mishyb_lib)
    primerdesign.setSeqArgs(seq_args)
    return primerdesign.runDesign(debug)


'''
The following functions are the modular subunits of `designPrimers` and may
be used in cases where performance or customiziation are of high priority.
'''

def setP3Globals(global_args, misprime_lib=None, mishyb_lib=None):
    ''' Set the Primer3 global args and misprime/mishyb libraries.

    Args:
        global_args (dict)            : Primer3 global parameters as per
                                        Primer3 docs

        misprime_lib (dict, optional) : ``<Sequence name: sequence>`` dict
                                        for mispriming checks.
        mishyb_lib (dict, optional)   : ``<Sequence name: sequence>`` dict
                                        for mishybridization checks.

    Returns:
        ``None``

    '''
    primerdesign.setGlobals(global_args, misprime_lib, mishyb_lib)


def setP3SeqArgs(seq_args):
    ''' Set the Primer3 sequence / design arguments.

    Args:
        seq_args (dict)     : Primer3 seq/design args as per Primer3 docs

    Returns:
        ``None``

    '''
    primerdesign.setSeqArgs(seq_args)


def runP3Design(debug=False):
    ''' Start the Primer3 design process, return a dict of the Primer3 output.

    The global parameters and seq args must have been previously set prior to
    this call (raises IOError).

    Args:
        debug (bool, optional)  : If ``True``, prints the received design
                                  params to stderr for debugging purposes

    Returns:
        A dictionary of Primer3 results (should be identical to the expected
        BoulderIO output from primer3_main)

    '''
    primerdesign.runDesign(debug)
