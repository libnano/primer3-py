# Copyright (C) 2020-2026. Ben Pruitt & Nick Conway; Wyss Institute
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
tests.test_str_rep_ifc_compat
~~~~~~~~~~~~~~~~~~~

Unit tests covering support for objects that implement the __str__ interface.
'''

import primer3


class SeqObj:
    '''Object with str representation.'''

    def __init__(self, seq):
        self._seq = seq

    def __str__(self):
        return self._seq


def test_calc_tm_with_object():
    '''Test Tm calc support for __str__ objs'''
    obj = SeqObj('GTAAAACGACGGCCAGT')
    tm = primer3.calc_tm(obj)

    assert isinstance(tm, float)


def test_calc_end_stability_with_object():
    '''Test thermo analysis function support for __str__ objs'''
    obj = SeqObj('GTAAAACGACGGCCAGT')

    res = primer3.calc_end_stability(obj, obj)

    assert hasattr(res, 'tm')
    assert isinstance(res.tm, float)
