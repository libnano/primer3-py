# Copyright (C) 2014-2023. Ben Pruitt & Nick Conway; Wyss Institute
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
primer3.argdefaults
~~~~~~~~~~~~~~~~~~~
'''
import dataclasses


@dataclasses.dataclass()
class Primer3PyArguments:
    '''Class containing the defaults values for the system

    NOTE: Goal is to match defaults of Primer3web at https://primer3.ut.ee
    '''
    mv_conc: float = 50.0
    dv_conc: float = 1.5
    dntp_conc: float = 0.6
    dna_conc: float = 50.0
    temp_c: float = 37.
    max_loop: int = 30
    output_structure: bool = False
    salt_corrections_method: str = 'santalucia'
    salt_corrections_method_int: int = 1
    max_nn_length: int = 60
    tm_method: str = 'santalucia'
    tm_method_int: int = 1
    temp_only: int = 0
    calc_type_wrapper = 'ANY'

    def todict(self) -> dict:
        return dataclasses.asdict(self)
