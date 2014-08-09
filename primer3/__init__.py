'''
primer3-py
~~~~~~~~~~

Python bindings for Primer3.

Current Primer3 version included in package: 2.3.6
Support for both Python 2.7.x and Python 3.x.x

'''

from primer3.bindings import (calcHairpin, calcHomodimer, calcHeterodimer,
							  calcHairpinTm, calcHomodimerTm,
							  calcHeterodimerTm, calcTm, setP3Globals, 
							  setP3SeqArgs, runP3Design, designPrimers)

import primer3.bindings as bindings
import primer3.wrappers as wrappers

__all__ = [
	# Low-level Tm-only bindings
	'calcHairpinTm', 'calcHomodimerTm', 'calcHeterodimerTm',
	# Low-level bindings
	'calcHairpin', 'calcHomodimer', 'calcHeterodimer', 'calcTm',
	# Primer3 design bindings
    'setP3Globals', 'setP3SeqArgs', 'runP3Design', 'designPrimers', 
    # Modules (bindings = C API bindings, wrappers = subprocess wrappers)
    'bindings', 'wrappers']
