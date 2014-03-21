'''
primer3.wrappers | wrappers.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Simple subprocess wrappers for Primer3 executables. These functions closely 
mirror the functions found in bindings.py, but are much slower and should
only be used for testing / comparative purposes.

'''

from __future__ import print_function

import glob
import os
import re
import subprocess
import sys

from collections import OrderedDict, namedtuple

from os.path import join as pjoin


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PRIMER3 WRAPPERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

DEV_NULL = open(os.devnull, 'wb')

local_dir = os.path.dirname(os.path.realpath(__file__))

if not os.environ.get('PRIMER3HOME'):
    try:
        np_dir = glob.glob(os.path.join(local_dir, 'src/primer3-*.*.*'))
        os.environ['PRIMER3HOME'] = os.path.abspath(np_dir[0])
    except:
        raise ImportError('PRIMER3HOME environmental variable is not set.')
PRIMER3_HOME = os.environ.get('PRIMER3HOME')

PRIMER3_SRC = pjoin(PRIMER3_HOME, 'src')
THERMO_PATH = pjoin(PRIMER3_SRC, 'primer3_config/')


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
    ''' Return the tm of `seq` as a float.
    '''
    tm_meth = _tm_methods.get(tm_method)
    if not tm_meth:
        raise ValueError('{} is not a valid tm calculation method'.format(
                         tm_meth))
    salt_meth = _tm_methods.get(tm_method)
    if not salt_meth:
        raise ValueError('{} is not a valid salt correction method'.format(
                         salt_meth))
    # For whatever reason mv_conc and dna_conc have to be ints
    args = [pjoin(PRIMER3_SRC, 'oligotm'),
            '-mv',  str(int(mv_conc)),
            '-dv',  str(dv_conc),
            '-n',   str(dntp_conc),
            '-d',   str(int(dna_conc)),
            '-tp',  str(tm_meth),
            '-sc',  str(salt_meth),
            seq]
    tm = subprocess.check_output(args, stderr=DEV_NULL,
                                 env=os.environ)
    return float(tm)


_ntthal_re = re.compile(b'dS\s+=\s+(\S+)\s+dH\s+=\s+(\S+)\s+' +
                        b'dG\s+=\s+(\S+)\s+t\s+=\s+(\S+)')

NTTHALRESULT = namedtuple('ntthal_result', ['ds', 'dh', 'dg', 'tm'])

def _parse_ntthal(ntthal_output):
    ''' Helper method that uses regex to parse ntthal output. '''
    parsed_vals = re.search(_ntthal_re, ntthal_output)
    return NTTHALRESULT(
        float(parsed_vals.group(1)),    # dS
        float(parsed_vals.group(2)),    # dH
        float(parsed_vals.group(3)),    # dG
        float(parsed_vals.group(4))     # tm
    ) if parsed_vals else None


def calcThermo(seq1, seq2, calc_type='ANY', mv_conc=50, dv_conc=0,
                 dntp_conc=0.8, dna_conc=50, temp_c=37, max_loop=30,
                 temp_only=False):
    """ Main subprocess wrapper for calls to the ntthal executable.

    Returns a named tuple with tm, ds, dh, and dg values or None if no
    structure / complex could be computed.
    """
    args = [pjoin(PRIMER3_SRC, 'ntthal'),
            '-a',       str(calc_type),
            '-mv',      str(mv_conc),
            '-dv',      str(dv_conc),
            '-n',       str(dntp_conc),
            '-d',       str(dna_conc),
            '-t',       str(temp_c),
            '-maxloop', str(max_loop),
            '-path',    THERMO_PATH,
            '-s1',      seq1,
            '-s2',      seq2]
    if temp_only:
        args += ['-r']
    out = subprocess.check_output(args, stderr=DEV_NULL,
                                  env=os.environ)
    return _parse_ntthal(out)


def calcHairpin(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
                 temp_c=37, max_loop=30, temp_only=False):
    ''' Return a namedtuple of the dS, dH, dG, and Tm of any hairpin struct
    present.

    Returns None if the sequence does not form a hairpin.

    '''
    return calcThermo(seq, seq, 'HAIRPIN', mv_conc, dv_conc, dntp_conc,
                        dna_conc, temp_c, max_loop, temp_only)


def calcHeterodimer(seq1, seq2, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                     dna_conc=50, temp_c=37, max_loop=30, temp_only=False):
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted heterodimer.

    Returns None if the sequences do not form a heterodimer.

    '''
    return calcThermo(seq1, seq2, 'ANY', mv_conc, dv_conc, dntp_conc,
                        dna_conc, temp_c, max_loop, temp_only)


def calcHomodimer(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                   dna_conc=50, temp_c=37, max_loop=30, temp_only=False):
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted homodimer.

    Returns None if the sequence does not form a homodimer.

    '''
    return calcThermo(seq, seq, 'ANY', mv_conc, dv_conc, dntp_conc,
                      dna_conc, temp_c, max_loop, temp_only)


def assessOligo(seq):
    '''
    Return the thermodynamic characteristics of hairpin/homodimer structures.

    Returns a tuple of namedtuples (hairpin data, homodimer data) in which each
    individual tuple is structured (dS, dH, dG, Tm).

    '''
    hairpin_out = calcHairpin(seq)
    homodimer_out = calcHomodimer(seq)
    return (hairpin_out, homodimer_out)


# ~~~~~~~ RUDIMENTARY PRIMER3 MAIN WRAPPER (see Primer3 docs for args) ~~~~~~ #

if sys.version_info[0] > 2:

    def _formatBoulderIO(p3_args):
        boulder_str = ''.join(['{}={}\n'.format(k,v) for k,v in 
                              p3_args.items()])
        boulder_str += '=\n'
        return bytes(boulder_str, 'UTF-8')

    def _parseBoulderIO(boulder_str):
        data_dict = OrderedDict()
        for line in boulder_str[0].decode("utf-8").split('\n'):
            try:
                k,v = line.strip().split('=')
                data_dict[k] = v
            except:
                pass
        return data_dict

else:

    def _formatBoulderIO(p3_args):
        boulder_str = ''.join(['{}={}\n'.format(k,v) for k,v in 
                              p3_args.items()])
        boulder_str += '=\n'
        return boulder_str

    def _parseBoulderIO(boulder_str):
        data_dict = OrderedDict()
        for line in boulder_str[0].split('\n'):
            try:
                k,v = line.strip().split('=')
                data_dict[k] = v
            except:
                pass
        return data_dict 


def designPrimers(p3_args):
    ''' Return the raw primer3_core output for the provided primer3 args.

    Returns an ordered dict of the boulderIO-format primer3 output file
    '''
    sp = subprocess.Popen([pjoin(PRIMER3_SRC, 'primer3_core')], 
                          stdout=subprocess.PIPE, stdin=subprocess.PIPE, 
                          stderr=subprocess.STDOUT)
    p3_args.setdefault('PRIMER_THERMODYNAMIC_PARAMETERS_PATH', 
                       pjoin(PRIMER3_SRC, 'primer3_config/'))
    in_str = _formatBoulderIO(p3_args)
    out_str = sp.communicate(input=in_str)
    return _parseBoulderIO(out_str)
