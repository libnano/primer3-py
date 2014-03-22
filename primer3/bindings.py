'''
primer3.bindings | bindings.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the main API for the Python C API Primer3 calls.

'''

import os

from collections import namedtuple
from os.path import join as pjoin

from primer3 import _primer3


# ~~~~~~~ Check to insure that the environment is properly configured ~~~~~~~ #

local_dir = os.path.dirname(os.path.realpath(__file__))

if not os.environ.get('PRIMER3HOME'):
    try:
        os.environ['PRIMER3HOME'] = pjoin(local_dir, 'src/primer3-2.3.6')
    except:
        raise ImportError('PRIMER3HOME environmental variable is not set.')
PRIMER3_HOME = os.environ.get('PRIMER3HOME')


# ~~~~~~~~~~~~~~~~ Load thermodynamic parameters into memory ~~~~~~~~~~~~~~~~ #

_primer3.getThermoParams(pjoin(PRIMER3_HOME, 'src', 'primer3_config/'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Low level bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Named tuple returned by all low-level bindings (if a structure is not
# present, the functions simply return `None`)
THERMORESULT = namedtuple('thermoresult', [
    'tm',               # Melting temperature (deg. Celsius)
    'ds',               # Entropy (cal/(K*mol))
    'dh',               # Enthalpy (kcal/mol)
    'dg',               # Gibbs free energy
    'align_end_1',      # Last alignment position in the 1st sequence.
    'align_end_2']      # Last alignment position in the 2nd sequence.
)


def _calcThermo(seq1, seq2, calc_type=0, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                   dna_conc=50, temp_c=37, max_loop=30):
    res = _primer3.calcThermo(seq1, seq2, calc_type, mv_conc, dv_conc,
                                dntp_conc, dna_conc, temp_c + 273.15,
                                max_loop, 0, 0)
    if res[1] == 1: # No structure found
        if res[0] != '':
            raise IOError('Primer3 error: {}'.format(res[0]))
        else:
            return None
    return THERMORESULT(*res[2:])


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
        A `thermoresult` named tuple of thermodynamic characteristics of the
        hairpin formation. If no hairpin is formed, returns `None`.

    '''

    return _calcThermo(seq, seq, 4, mv_conc, dv_conc, dntp_conc, dna_conc,
                       temp_c, max_loop)


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
        A `thermoresult` named tuple of thermodynamic characteristics of the
        homodimer interaction. If no interaction is present, returns `None`.

    '''
    return _calcThermo(seq, seq, 1, mv_conc, dv_conc, dntp_conc, dna_conc,
                       temp_c, max_loop)


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
        A `thermoresult` named tuple of thermodynamic characteristics of the
        heterodimer interaction. If no interaction is present, returns `None`.

    '''
    return _calcThermo(seq1, seq2, 1, mv_conc, dv_conc, dntp_conc, dna_conc,
                       temp_c, max_loop)


_tm_methods = {
    'breslauer': 0,
    'santalucia': 1
}

_salt_correction_methods = {
    'schildkraut': 0,
    'santalucia': 1,
    'owczarzy': 2
}

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
    tm_meth = _tm_methods.get(tm_method)
    if not tm_meth:
        raise ValueError('{} is not a valid tm calculation method'.format(
                         tm_meth))
    salt_meth = _tm_methods.get(tm_method)
    if not salt_meth:
        raise ValueError('{} is not a valid salt correction method'.format(
                         salt_meth))
    return _primer3.calcTm(seq, mv_conc, dv_conc, dntp_conc, dna_conc,
                                  max_nn_length, tm_meth, salt_meth)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Design bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def designPrimers(seq_args, global_args=None, misprime_lib=None,
                  mishyb_lib=None):
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
        _primer3.setGlobals(global_args, misprime_lib, mishyb_lib)
    _primer3.setGlobals(global_args, misprime_lib, mishyb_lib)
    _primer3.setSeqArgs(seq_args)
    return _primer3.runDesign()


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
    _primer3.setGlobals(global_args, misprime_lib, mishyb_lib)


def setP3SeqArgs(seq_args):
    ''' Set the Primer3 sequence / design arguments.

    Args:
        seq_args (dict)     :   Primer3 seq/design args as per Primer3 docs

    '''
    _primer3.setSeqArgs(seq_args)


def runP3Design():
    ''' Start the Primer3 design process, return a dict of the Primer3 output.

    The global parameters and seq args must have been previously set prior to
    this call (raises IOError).

    Returns:
        A dictionary of Primer3 results (should be identical to the expected
        BoulderIO output from primer3_main)

    '''
    _primer3.runDesign()
