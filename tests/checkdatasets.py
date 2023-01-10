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
tests.checkdatasets
~~~~~~~~~~~~~~~~~~~
'''
import glob
import json
import os
import os.path
import subprocess
import sys
from typing import (
    Any,
    Dict,
)

import click

from primer3 import __version__ as p3version
from primer3.thermoanalysis import (
    ThermoAnalysis,
    ThermoResult,
)

dirname = os.path.dirname
extra = [dirname(dirname(os.path.abspath(__file__)))]
sys.path = extra + sys.path


def precision(x: float) -> float:
    '''Round the precision of floats to make comparisons exact

    Args:
        x: Float to round

    Returns:
        float rounded to 2 decimal places
    '''
    return round(x, 2)


def todict_tr(x: ThermoResult) -> Dict[str, Any]:
    '''Convert ``ThermoResult`` instance to a dictionary.
    Round float values for readability

    Args:
        x: ``ThermoResult`` instance

    Returns:
        Dictionary form of a ``ThermoResult``
    '''
    return {
        'structure_found': x.structure_found,
        'tm': precision(x.tm),
        'dg': precision(x.dg),
        'dh': precision(x.dh),
        'ds': precision(x.ds),
    }


def todict_ta(x: ThermoAnalysis) -> Dict[str, Any]:
    '''Convert ``ThermoAnalysis`` instance to a dictionary.

    Args:
        x: ``ThermoAnalysis`` instance

    Returns:
        Dictionary form of a ``ThermoAnalysis``
    '''
    return {
        'mv_conc': x.mv_conc,
        'dv_conc': x.dv_conc,
        'dntp_conc': x.dntp_conc,
        'dna_conc': x.dna_conc,
        'temp_c': x.temp,
        'max_loop': x.max_loop,
        # 'temp_only': x.temp_only,
        # 'debug': x.thalargs.debug,
        'max_nn_length': x.max_nn_length,
        'tm_method': x.tm_method,
        'salt_correction_method': x.salt_correction_method,
    }


def git_hash() -> str:
    '''
    Returns:
        Git hash of repo or NA if not available
    '''
    try:
        with open(os.devnull, 'wb') as devnull:
            hash_str = subprocess.check_output(
                ['git', 'rev-parse', '--short', 'HEAD'],
                stderr=devnull,
            )
        return hash_str.strip().decode('utf8')
    except BaseException:
        return 'NA'
# end def


def create_data_set() -> Dict[str, Any]:
    '''Create data sets for testing

    Returns:
        Dictionary of test dataset dictionaries
    '''
    seq1 = 'GGGGCCCCCCAAATTTTTT'
    seq2 = 'AAAAAATTTGGCCCCAAA'

    ta = ThermoAnalysis()
    tests = [
        ['calcHeterodimer', seq1, seq2],
        ['calcHairpin', seq1],
        ['calcHomodimer', seq1],
        ['calcTm', seq1],
    ]
    cases = []
    for test in tests:
        f = getattr(ta, test[0])
        res = f(*test[1:])
        if isinstance(res, ThermoResult):
            cases.append([test, todict_tr(res)])
        else:
            cases.append([test, res])
    out = {
        'thermoanalysis': todict_ta(ta),
        'cases': cases,
    }
    return out


def dump_data_set(x: Dict[str, Any]) -> str:
    '''
    Dump data set to a JSON file with the data time and git hash

    Args:
        x: Data set dictionary to dump

    Returns:
        Filename formatted in form::

        `primer3_precision_<date>_<githash>_<OS-platform>_<Python-version>.json`

    '''
    import time
    date_str = time.strftime('%Y-%m-%d-%H-%M-%S', time.gmtime())
    git_hash_str = git_hash()
    py_version = '%s.%s.%s' % sys.version_info[0:3]
    platform = sys.platform
    filename = 'primer3_precision_%s_%s_%s_%s.json' % (
        date_str, git_hash_str, platform, py_version,
    )
    with open(filename, 'w') as fd:
        print('dumping')
        json.dump(x, fd)
    return filename


def read_data_set() -> Dict[str, Any]:
    '''
    Read data set JSON file

    Returns:
        Dictionary of a the first JSON file found in the local directory
    '''
    files = glob.glob('*.json')
    if len(files) > 0:
        the_file = files[0]
        file_name = os.path.splitext(the_file)[0]
        git_hash, platform, py_version = file_name.split('_')[3:6]
        print('\t- Checking reference file: %s' % the_file)
        print(
            '\t\tgit hash: %s\n\t\tPlatform: %s\n\t\tPython: %s' %
            (git_hash, platform, py_version),
        )
        with open(the_file, 'r') as fd:
            out = json.load(fd)
        return out
    else:
        raise OSError('No json file found')


def compare_sets(
        reference: Dict[str, Any],
        test_set: Dict[str, Any],
) -> bool:
    '''
    Compare thermoanalysis datasets

    Args:
        reference: Reference thermoanalysis dataset
        test_set: Test thermoanalysis data set

    Returns:
        True if the same, False otherwise
    '''
    for k1, v in reference.items():
        if k1 not in test_set:
            return False
        test_val = test_set[k1]
        if k1 == 'thermoanalysis':
            for k2, item in v.items():
                if item != test_val[k2]:
                    return False
        elif k1 == 'cases':
            for item_ref, item_test in zip(v, test_val):
                ref_args, ref_result = item_ref
                test_args, test_result = item_test
                if ref_args != test_args:
                    print(ref_args, test_args)
                    return False
                if ref_result != test_result:
                    print(ref_result, test_result)
                    return False

    return True


@click.command()
@click.option('--generate', default=False, help='Generate a set')
def runner(generate):
    create_dict = create_data_set()
    if generate:
        filename = dump_data_set(create_dict)
        print('> Dumped data set: %s' % (filename))
    else:
        print('> Comparing reference and current revisions')
        print('\t- Current git hash: %s' % (git_hash()))
        print('\t- Current primer3 version: %s' % (p3version))
        read_dict = read_data_set()
        print('\t- Is the reference set equal to the current revisions set?')
        are_equal = compare_sets(read_dict, create_dict)
        print('\t\t%s' % (str(are_equal)))


if __name__ == '__main__':
    runner()
