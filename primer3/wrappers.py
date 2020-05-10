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

import glob
import os
import re
import subprocess
import sys

from collections import OrderedDict, namedtuple

from os.path import join as pjoin


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PRIMER3 WRAPPERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~ Check to insure that the environment is properly configured ~~~~~~~ #

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

if not os.environ.get('PRIMER3HOME'):
    try:
        os.environ['PRIMER3HOME'] = pjoin(LOCAL_DIR, 'src/libprimer3')
    except:
        raise ImportError('PRIMER3HOME environmental variable is not set.')
PRIMER3_HOME = os.environ.get('PRIMER3HOME')

THERMO_PATH = pjoin(PRIMER3_HOME, 'primer3_config/')

DEV_NULL = open(os.devnull, 'wb')

_tm_methods = {
    'breslauer': 0,
    'santalucia': 1
}

_salt_corrections_methods = {
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
    if tm_meth is None:
        raise ValueError('{} is not a valid tm calculation method'.format(
                         tm_method))
    salt_meth = _salt_corrections_methods.get(salt_corrections_method)
    if salt_meth is None:
        raise ValueError('{} is not a valid salt correction method'.format(
                         salt_corrections_method))
    # For whatever reason mv_conc and dna_conc have to be ints
    args = [pjoin(PRIMER3_HOME, 'oligotm'),
            '-mv',  str(mv_conc),
            '-dv',  str(dv_conc),
            '-n',   str(dntp_conc),
            '-d',   str(dna_conc),
            '-tp',  str(tm_meth),
            '-sc',  str(salt_meth),
            seq]
    tm = subprocess.check_output(args, stderr=DEV_NULL,
                                 env=os.environ)
    return float(tm)


_ntthal_re = re.compile(b'dS\s+=\s+(\S+)\s+dH\s+=\s+(\S+)\s+' +
                        b'dG\s+=\s+(\S+)\s+t\s+=\s+(\S+)')

THERMORESULT = namedtuple('thermoresult', [
        'result',           # True if a structure is present
        'ds',               # Entropy (cal/(K*mol))
        'dh',               # Enthalpy (kcal/mol)
        'dg',               # Gibbs free energy
        'tm',               # Melting temperature (deg. Celsius)
        'ascii_structure'   # ASCII representation of structure
    ]
)

NULLTHERMORESULT = THERMORESULT(False, 0, 0, 0, 0, '')

def _parse_ntthal(ntthal_output):
    ''' Helper method that uses regex to parse ntthal output. '''
    parsed_vals = re.search(_ntthal_re, ntthal_output)
    if parsed_vals:
        ascii_structure = (
            ntthal_output[ntthal_output.index(b'\n') + 1:].decode('utf-8')
        )
        res = THERMORESULT(
            True,                           # Structure found
            float(parsed_vals.group(1)),    # dS
            float(parsed_vals.group(2)),    # dH
            float(parsed_vals.group(3)),    # dG
            float(parsed_vals.group(4)),    # tm
            ascii_structure
        )
    else:
        res = NULLTHERMORESULT
    return res


def calcThermo(seq1, seq2, calc_type='ANY', mv_conc=50, dv_conc=0,
                 dntp_conc=0.8, dna_conc=50, temp_c=37, max_loop=30,
                 temp_only=False):
    """ Main subprocess wrapper for calls to the ntthal executable.

    Returns a named tuple with tm, ds, dh, and dg values or None if no
    structure / complex could be computed.
    """
    args = [pjoin(PRIMER3_HOME, 'ntthal'),
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
    '''
    return calcThermo(seq, seq, 'HAIRPIN', mv_conc, dv_conc, dntp_conc,
                      dna_conc, temp_c, max_loop, temp_only)


def calcHeterodimer(seq1, seq2, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                     dna_conc=50, temp_c=37, max_loop=30, temp_only=False):
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted heterodimer.
    '''
    return calcThermo(seq1, seq2, 'ANY', mv_conc, dv_conc, dntp_conc,
                      dna_conc, temp_c, max_loop, temp_only)


def calcHomodimer(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                   dna_conc=50, temp_c=37, max_loop=30, temp_only=False):
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted homodimer.
    '''
    return calcThermo(seq, seq, 'ANY', mv_conc, dv_conc, dntp_conc,
                      dna_conc, temp_c, max_loop, temp_only)


def calcEndStability(seq1, seq2, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                     dna_conc=50, temp_c=37, max_loop=30, temp_only=False):
    ''' Return a tuple of the dS, dH, dG, and Tm of any predicted heterodimer.
    '''
    return calcThermo(seq1, seq2, 'END1', mv_conc, dv_conc, dntp_conc,
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
        for line in boulder_str.decode("utf-8").split('\n'):
            try:
                k,v = line.strip().split('=')
                data_dict[k] = v
            except:
                pass
        return data_dict

    def _parseMultiRecordBoulderIO(boulder_str):
        data_dicts = []
        for record in re.split('=\r?\n', boulder_str):
            if record == '':
                continue
            data_dict = OrderedDict()
            for line in record.split('\n'):
                try:
                    k,v = line.strip().split('=')
                    data_dict[k] = v
                except:
                    pass
            data_dicts.append(data_dict)
        return data_dicts

else:

    def _formatBoulderIO(p3_args):
        boulder_str = ''.join(['{}={}\n'.format(k,v) for k,v in
                              p3_args.items()])
        boulder_str += '=\n'
        return boulder_str

    def _parseBoulderIO(boulder_str):
        data_dict = OrderedDict()
        for line in boulder_str.split('\n'):
            try:
                k,v = line.strip().split('=')
                data_dict[k] = v
            except:
                pass
        return data_dict

    def _parseMultiRecordBoulderIO(boulder_str):
        data_dicts = []
        for record in re.split('=\r?\n', boulder_str):
            if record == '':
                continue
            data_dict = OrderedDict()
            for line in record.split('\n'):
                try:
                    k,v = line.strip().split('=')
                    data_dict[k] = v
                except:
                    pass
            data_dicts.append(data_dict)
        return data_dicts


def designPrimers(p3_args, input_log=None, output_log=None, err_log=None):
    ''' Return the raw primer3_core output for the provided primer3 args.

    Returns an ordered dict of the boulderIO-format primer3 output file
    '''
    sp = subprocess.Popen([pjoin(PRIMER3_HOME, 'primer3_core')],
                          stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
    p3_args.setdefault('PRIMER_THERMODYNAMIC_PARAMETERS_PATH',
                       pjoin(PRIMER3_HOME, 'primer3_config/'))
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
