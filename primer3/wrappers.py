'''
Simple subprocess wrappers for the primer3 library. 

'''
from __future__ import print_function

import glob
import os
import re
import subprocess

from collections import OrderedDict, namedtuple

from os.path import join as pjoin


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ UTILITY FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Underscore imports for TermporaryDirectory (see issue #10188)
import warnings as _warnings
import sys as _sys
import os as _os

import atexit

from tempfile import mkdtemp
from uuid import uuid4


# Although it does not have an underscore for historical reasons, this
# variable is an internal implementation detail (see issue 10354).
template = 'tmp'

class TemporaryDirectory(object):
    '''Create and return a temporary directory. This has the same
    behavior as mkdtemp but can be used as a context manager. For
    example:

        with TemporaryDirectory() as tmpdir:
            ...

    Upon exiting the context, the directory and everything contained
    in it are removed. This is an augmented version of the class included
    in the Python 3 tempfile module and it includes methods related to
    module-persistent temperorary directories that are cleaned up at exit.
    '''

    def __init__(self, suffix='', prefix=template, dir=None, persist=False):
        self._closed = False
        self.name = None # Handle mkdtemp raising an exception
        self.name = mkdtemp(suffix, prefix, dir)
        if persist:
            self.persist_until_exit()

    def __repr__(self):
        return '<{} {!r}>'.format(self.__class__.__name__, self.name)

    def __enter__(self):
        return self.name

    def cleanup(self, _warn=False):
        if self.name and not self._closed:
            try:
                self._rmtree(self.name)
            except (TypeError, AttributeError) as ex:
                # Issue #10188: Emit a warning on stderr
                # if the directory could not be cleaned
                # up due to missing globals
                if 'None' not in str(ex):
                    raise
                print('ERROR: {!r} while cleaning up {!r}'.format(ex, self,),
                      file=_sys.stderr)
                return
            self._closed = True
            if _warn:
                self._warn('Implicitly cleaning up {!r}'.format(self),
                           Warning)

    def persist_until_exit(self):
        atexit.register(self.cleanup)

    def __exit__(self, exc, value, tb):
        self.cleanup()

    def __del__(self):
        # Issue a ResourceWarning if implicit cleanup needed
        self.cleanup(_warn=True)

    # XXX (ncoghlan): The following code attempts to make
    # this class tolerant of the module nulling out process
    # that happens during CPython interpreter shutdown
    # Alas, it doesn't actually manage it. See issue #10188
    _listdir = staticmethod(_os.listdir)
    _path_join = staticmethod(_os.path.join)
    _isdir = staticmethod(_os.path.isdir)
    _islink = staticmethod(_os.path.islink)
    _remove = staticmethod(_os.remove)
    _rmdir = staticmethod(_os.rmdir)
    _warn = _warnings.warn

    def _rmtree(self, path):
        # Essentially a stripped down version of shutil.rmtree.  We can't
        # use globals because they may be None'ed out at shutdown.
        for name in self._listdir(path):
            fullname = self._path_join(path, name)
            try:
                isdir = self._isdir(fullname) and not self._islink(fullname)
            except OSError:
                isdir = False
            if isdir:
                self._rmtree(fullname)
            else:
                try:
                    self._remove(fullname)
                except OSError:
                    pass
        try:
            self._rmdir(path)
        except OSError:
            pass


def unique_fn(directory, file_ext=''):
    ''' Generates a unique filename and returns the full file path

    Insures that the file does not already exist but does not return an
    open file (returns the filepath as a string)'''

    while True:
        fp = _os.path.join(directory, str(uuid4()) + file_ext)
        try:
            with open(fp, 'rb') as fd:
                pass
        except:
            return fp



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PRIMER3 WRAPPERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

TMP_DIR = TemporaryDirectory(persist=True).name
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
    tm = subprocess.check_output(args, cwd=TMP_DIR, stderr=DEV_NULL, 
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
    out = subprocess.check_output(args, cwd=TMP_DIR, stderr=DEV_NULL, 
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

# TODO: modify runP3Main to use STDERR and STDIN for IO (see docs)

def _writeBoulderIO(fp, p3_args):
    with open(fp, 'wb') as fd:
        for k, v in p3_args.items():
            fd.write(k + '=' + v + '\n')
        fd.write('=')


def _parseBoulderIO(fp):
    data_dict = OrderedDict()
    with open(fp, 'rb') as fd:
        for line in fd:
            k,v = line.strip().split('=')
            data_dict[k] = v
    return data_dict


def runP3Main(p3_args):
    ''' Return the raw primer3_core output for the provided primer3 args.

    Returns an ordered dict of the boulderIO-format primer3 output file
    '''
    fp_in = unique_fn(TMP_DIR, file_ext='.in')
    _writeBoulderIO(fp_in, p3_args)
    fp_out = fp_in.rstrip('.in') + '.out'
    subprocess.check_call([pjoin(PRIMER3_SRC, 'primer3_core'), '--output', 
                           fp_out, fp_in], stderr=DEV_NULL, env=os.environ)
    return _parseBoulderIO(fp_out)

