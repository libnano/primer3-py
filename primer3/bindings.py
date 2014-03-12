import glob
import os

from collections import namedtuple
from os.path import join as pjoin

from primer3 import _primer3


# ~~~~~~~ Check to insure that the environment is properly configured ~~~~~~~ #

local_dir = os.path.dirname(os.path.realpath(__file__))

if not os.environ.get('PRIMER3HOME'):
    try:
        np_dir = glob.glob(os.path.join(local_dir, 'src/primer3-*.*.*'))
        os.environ['PRIMER3HOME'] = os.path.abspath(np_dir[0])
    except:
        raise ImportError('PRIMER3HOME environmental variable is not set.')
PRIMER3_HOME = os.environ.get('PRIMER3HOME')


# ~~~~~~~~~~~~~~~~ Load thermodynamic parameters into memory ~~~~~~~~~~~~~~~~ #

_primer3.getThermoParams(pjoin(PRIMER3_HOME, 'src', 'primer3_config/'))



THALRESULT = namedtuple('thal_result', ['msg', 'no_structure', 'temp_c',
                                        'ds', 'dh', 'dg', 'align_end_1',
                                        'align_end_2'])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def thalFunction(seq1, seq2, calc_type=0, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                 dna_conc=50, temp_c=37, max_loop=30, temp_only=False):
    res = _primer3.calcThermo(seq1, seq2, calc_type, mv_conc, dv_conc,
                                dntp_conc, dna_conc, temp_c + 273.15,
                                max_loop, temp_only, 0)
    return THALRESULT(*res)


def calcHairpin(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
                temp_c=37, max_loop=30, temp_only=False):
    return thalFunction(seq, seq, 4, mv_conc, dv_conc, dntp_conc, dna_conc,
                        temp_c, max_loop, temp_only)


def calcDimer(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, dna_conc=50,
                temp_c=37, max_loop=30, temp_only=False):
    return thalFunction(seq, seq, 1, mv_conc, dv_conc, dntp_conc, dna_conc,
                        temp_c, max_loop, temp_only)


def calcHeterodimer(seq1, seq2, mv_conc=50, dv_conc=0, dntp_conc=0.8,
                    dna_conc=50, temp_c=37, max_loop=30, temp_only=False):
    return thalFunction(seq1, seq2, 1, mv_conc, dv_conc, dntp_conc, dna_conc,
                        temp_c, max_loop, temp_only)


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

def designPrimers(p3_args, s_args, misprime_lib=None, mishyb_lib=None):
    return _primer3.designPrimers(p3_args, s_args, misprime_lib, mishyb_lib)

def testDesignPrimers():
    p3_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
    }
    s_args = {
        'PRIMER_PRODUCT_SIZE_RANGE': [75,100,100,125,125,150,150,175,175,200,200,225],
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG',
        'SEQUENCE_INCLUDED_REGION': [36,342]
    }
    print(designPrimers(p3_args, s_args))
