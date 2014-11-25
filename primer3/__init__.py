'''
primer3-py
~~~~~~~~~~

Python bindings / abstractions for the Primer3 primer design /
oligonucleotide thermodynamics library. 

'''

import os

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

# ~~~~~~~~~ Get / set the environ. variable for the libprimer3 path ~~~~~~~~~ #

if not os.environ.get('PRIMER3HOME'):
    PRIMER3HOME = os.path.join(LOCAL_DIR, 'src/libprimer3')
    if not os.path.exists(PRIMER3HOME):
        raise OSError('PRIMER3HOME environmental variable is not set '
                      'and a path to libprimer3 could not be found: '
                      '<%s>' % PRIMER3HOME)
    os.environ['PRIMER3HOME'] = PRIMER3HOME

from .bindings import (calcHairpin, calcHomodimer, calcHeterodimer,
                       calcHairpinTm, calcHomodimerTm, calcHeterodimerTm, 
                       calcTm, setP3Globals, setP3SeqArgs, runP3Design, 
                       designPrimers)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

from . import bindings, wrappers, thermoanalysis, primerdesign


def includes():
    return [LOCAL_DIR, os.environ['PRIMER3HOME']]


__all__ = [
    # Low-level Tm-only bindings
    'calcHairpinTm', 'calcHomodimerTm', 'calcHeterodimerTm',
    # Low-level bindings
    'calcHairpin', 'calcHomodimer', 'calcHeterodimer', 'calcTm',
    # Primer3 design bindings
    'setP3Globals', 'setP3SeqArgs', 'runP3Design', 'designPrimers', 
    # Modules (bindings = C API bindings, wrappers = subprocess wrappers)
    'bindings', 'wrappers', 'thermoanalysis', 'primerdesign'
    ]


