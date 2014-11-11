'''
primer3.bindings | bindings.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the main API for the Python C API Primer3 calls.

'''

import os

from os.path import join as pjoin

from . import thermoanalysis
from . import primerdesign 


# ~~~~~~~ Check to insure that the environment is properly configured ~~~~~~~ #

PRIMER3_HOME = os.environ.get('PRIMER3HOME')


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
                temp_c=37, max_loop=30):
    ''' Calculate the hairpin formation thermodynamics of a DNA sequence.

    Args:
        seq (str)               : DNA sequence to analyze for hairpin formation

    Kwargs:
        mv_conc (float/int)     : Monovalent cation concentration (mM)
        dv_conc (float/int)     : Divalent cation concentration (mM)
        dntp_conc (float/int)   : dNTP concentration (mM)
        dna_conc (float/int)    : DNA concentration (nM)
        temp_c (int)            : Simulation temperature for dG (Celsius)
        max_loop(int)           : Maximum size of loops in the structure

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        hairpin formation.

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcHairpin(seq)


def calcHomodimer(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
                  temp_c=37, max_loop=30):
    ''' Calculate the homodimerization thermodynamics of a DNA sequence.

    Args:
        seq (str)               : DNA sequence to analyze for homodimer
                                  formation calculations

    Kwargs:
        mv_conc (float/int)     : Monovalent cation concentration (mM)
        dv_conc (float/int)     : Divalent cation concentration (mM)
        dntp_conc (float/int)   : dNTP concentration (mM)
        dna_conc (float/int)    : DNA concentration (nM)
        temp_c (int)            : Simulation temperature for dG (Celsius)
        max_loop (int)          : Maximum size of loops in the structure

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        homodimer interaction. 

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcHomodimer(seq)


def calcHeterodimer(seq1, seq2, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                    dna_conc=50, temp_c=37, max_loop=30):
    ''' Calculate the heterodimerization thermodynamics of two DNA sequences.

    Args:
        seq1 (str)              : First DNA sequence to analyze for heterodimer
                                  formation
        seq2 (str)              : Second DNA sequence to analyze for
                                  heterodimer formation

    Kwargs:
        mv_conc (float/int)     : Monovalent cation concentration (mM)
        dv_conc (float/int)     : Divalent cation concentration (mM)
        dntp_conc (float/int)   : dNTP concentration (mM)
        dna_conc (float/int)    : DNA concentration (nM)
        temp_c (int)            : Simulation temperature for dG (Celsius)
        max_loop(int)           : Maximum size of loops in the structure

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        heterodimer interaction. 

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcHeterodimer(seq1, seq2)


def calcEndStability(seq1, seq2, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                     dna_conc=50, temp_c=37, max_loop=30):
    ''' Calculate the 3' end stability of DNA sequence `seq1` against DNA 
    sequence `seq2`.

    Args:
        seq1 (str)              : DNA sequence to analyze for 3' end
                                  hybridization against the target sequence
        seq2 (str)              : Target DNA sequence to analyze for
                                  seq1 3' end hybridization

    Kwargs:
        mv_conc (float/int)     : Monovalent cation concentration (mM)
        dv_conc (float/int)     : Divalent cation concentration (mM)
        dntp_conc (float/int)   : dNTP concentration (mM)
        dna_conc (float/int)    : DNA concentration (nM)
        temp_c (int)            : Simulation temperature for dG (Celsius)
        max_loop(int)           : Maximum size of loops in the structure

    Returns:
        A `ThermoResult` object with thermodynamic characteristics of the
        3' hybridization interaction. 

    '''
    _setThermoArgs(**locals())
    return _THERMO_ANALYSIS.calcEndStability(seq1, seq2)


def calcTm(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
           max_nn_length=60, tm_method='santalucia',
           salt_corrections_method='santalucia'):
    ''' Calculate the melting temperature of a DNA sequence.

    Args:
        seq (str)               : DNA sequence

    Kwargs:
        mv_conc (float/int)     : Monovalent cation concentration (mM)
        dv_conc (float/int)     : Divalent cation concentration (mM)
        dntp_conc (float/int)   : dNTP concentration (mM)
        dna_conc (float/int)    : DNA concentration (nM)
        max_nn_length (int)     : Maximum length for nearest-neighbor calcs
        tm_method (str)         : Tm calculation method (breslauer or
                                  santalucia)
        salt_corrections_method
                          (str) : Salt correction method (schildkraut,
                                  owczarzy, santalucia)

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
        seq_args (dict)     :   Primer3 sequence/design args as per Primer3 docs

    Kwargs:
        global_args (dict)  :   Primer3 global args as per Primer3 docs
        misprime_lib (dict) :   `Sequence name: sequence` dictionary for
                                mispriming checks.
        mishyb_lib (dict)   :   `Sequence name: sequence` dictionary for
                                mishybridization checks.

    Returns:
        A dictionary of Primer3 results (should be identical to the expected
        BoulderIO output from primer3_main)

    '''
    if global_args:
        primerdesign.setGlobals(global_args, misprime_lib, mishyb_lib)
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
        global_args (dict)  :   Primer3 global parameters as per Primer3 docs

    Kwargs:
        misprime_lib (dict) :   `Sequence name: sequence` dictionary for
                                mispriming checks.
        mishyb_lib (dict)   :   `Sequence name: sequence` dictionary for
                                mishybridization checks.

    '''
    primerdesign.setGlobals(global_args, misprime_lib, mishyb_lib)


def setP3SeqArgs(seq_args):
    ''' Set the Primer3 sequence / design arguments.

    Args:
        seq_args (dict)     :   Primer3 seq/design args as per Primer3 docs

    '''
    primerdesign.setSeqArgs(seq_args)


def runP3Design(debug=False):
    ''' Start the Primer3 design process, return a dict of the Primer3 output.

    The global parameters and seq args must have been previously set prior to
    this call (raises IOError).

    Returns:
        A dictionary of Primer3 results (should be identical to the expected
        BoulderIO output from primer3_main)

    '''
    primerdesign.runDesign(debug)
