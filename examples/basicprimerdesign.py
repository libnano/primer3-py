# Copyright (C) 2023. Ben Pruitt & Nick Conway
# See LICENSE for full GPLv2 license.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
examples.basicprimerdesign

A simple script to demonstrate running primer design
'''
import os.path as op
import random
import sys
from typing import (
    Any,
    Dict,
)

try:
    import primer3
except ModuleNotFoundError:
    PACKAGE_DIR = op.dirname(op.dirname(op.abspath(__file__)))
    sys.path.append(PACKAGE_DIR)
    import primer3


def basic_primer_design() -> Dict[str, Any]:
    '''Peform basic primer design given the parameters defined in the function
    below

    Returns:
        Dictionary of the design results
    '''
    sequence_template = (
        'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCAT'
        'GCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGG'
        'AGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNA'
        'AAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAA'
        'TGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCA'
        'CACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG'
    )
    quality_list = [
        random.randint(20, 90)
        for i in range(len(sequence_template))
    ]
    seq_args = {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': sequence_template,
        'SEQUENCE_QUALITY': quality_list,
        'SEQUENCE_INCLUDED_REGION': (36, 342),
    }
    global_args = {
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
        'PRIMER_PRODUCT_SIZE_RANGE': [
            [75, 100],
            [100, 125],
            [125, 150],
            [150, 175],
            [175, 200],
            [200, 225],
        ],
    }
    design_result_dict = primer3.design_primers(
        seq_args=seq_args,
        global_args=global_args,
    )
    return design_result_dict


if __name__ == '__main__':
    result_dict = basic_primer_design()
    import pprint
    print('Generated the following result')
    pprint.pprint(result_dict)
