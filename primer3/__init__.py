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
primer3-py
~~~~~~~~~~

Python bindings / abstractions for the Primer3 primer design /
oligonucleotide thermodynamics library.

'''

import os
from typing import List

from . import (  # type: ignore
    bindings,
    primerdesign,
    thermoanalysis,
    wrappers,
)
from .bindings import (
    calcHairpin,
    calcHairpinTm,
    calcHeterodimer,
    calcHeterodimerTm,
    calcHomodimer,
    calcHomodimerTm,
    calcTm,
    designPrimers,
    runP3Design,
    setP3Globals,
    setP3SeqArgs,
)

LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))


def includes() -> List[str]:
    '''
    Returns:
        List of directories to include in egg/wheel during packaging
    '''
    return [LOCAL_DIR, os.path.join(LOCAL_DIR, 'src', 'libprimer3')]


__author__ = 'Ben Pruitt, Nick Conway'
__copyright__ = 'Copyright 2014-2020, Ben Pruitt & Nick Conway; Wyss Institute'
__license__ = 'GPLv2'
__version__ = '1.0.0-alpha.3'

__all__ = [
    # Low-level Tm-only bindings
    'calcHairpinTm', 'calcHomodimerTm', 'calcHeterodimerTm',
    # Low-level bindings
    'calcHairpin', 'calcHomodimer', 'calcHeterodimer', 'calcTm',
    # Primer3 design bindings
    'setP3Globals', 'setP3SeqArgs', 'runP3Design', 'designPrimers',
    # Modules (bindings = C API bindings, wrappers = subprocess wrappers)
    'bindings', 'wrappers', 'thermoanalysis', 'primerdesign',
]
