'''
primer3-py
~~~~~~~~~~

Python bindings for Primer3.

'''

import os
LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

if not os.environ.get('PRIMER3HOME'):
    try:
        os.environ['PRIMER3HOME'] = os.path.join(LOCAL_DIR, 'src/libprimer3')
    except:
        raise ImportError('PRIMER3HOME environmental variable is not set.')

from primer3.bindings import (calcHairpin, calcHomodimer, calcHeterodimer,
                              calcHairpinTm, calcHomodimerTm,
                              calcHeterodimerTm, calcTm, setP3Globals, 
                              setP3SeqArgs, runP3Design, designPrimers)

from primer3 import bindings, wrappers, simulatedbindings, tests


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
    'bindings', 'wrappers', 'thermoanalyis', 'simulatedbindings',
    # Unit tests
    'tests'
    ]
