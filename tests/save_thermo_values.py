# Copyright (C) 2020-2025. Ben Pruitt & Nick Conway; Wyss Institute
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
tests.save_thermo_values
~~~~~~~~~~~~~~~~~~~~~~~

Script to generate and save gold standard thermodynamic values for
regression testing. Imports calculation and saving functions from
test_sequences.

'''

import json
import platform
import subprocess
import sys
from datetime import datetime
from importlib.metadata import (
    PackageNotFoundError,
    version,
)
from pathlib import Path

import tomli

from primer3 import __version__ as primer3_py_version
from primer3 import thermoanalysis

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from tests.test_sequences import calculate_thermo_values


def get_git_info():
    '''Get git branch and latest tag information.'''
    try:
        branch = subprocess.check_output(
            ['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
            stderr=subprocess.DEVNULL,
        ).decode('utf-8').strip()

        latest_tag = subprocess.check_output(
            ['git', 'describe', '--tags', '--abbrev=0'],
            stderr=subprocess.DEVNULL,
        ).decode('utf-8').strip()

        return {
            'branch': branch,
            'latest_tag': latest_tag,
        }
    except subprocess.CalledProcessError:
        return {
            'branch': 'unknown',
            'latest_tag': 'unknown',
        }


def get_dependency_versions():
    '''Get versions of dependencies from pyproject.toml.'''
    deps = {}

    # Read pyproject.toml
    pyproject_path = Path(__file__).parent.parent / 'pyproject.toml'
    with open(pyproject_path, 'rb') as f:
        pyproject = tomli.load(f)

    # Get build dependencies
    build_deps = pyproject.get('build-system', {}).get('requires', [])
    for dep in build_deps:
        name = dep.split('>=')[0].split('~=')[0].split('==')[0].strip()
        try:
            deps[name] = version(name)
        except PackageNotFoundError:
            deps[name] = 'not installed'

    # Get dev dependencies
    dev_deps = pyproject.get('project', {}).get(
        'optional-dependencies', {},
    ).get('dev', [])
    for dep in dev_deps:
        name = dep.split('>=')[0].split('~=')[0].split('==')[0].strip()
        try:
            deps[name] = version(name)
        except PackageNotFoundError:
            deps[name] = 'not installed'

    return deps


def save_thermo_values_with_metadata(
        values,
        filename='tests/thermo_standard_values.json',
):
    '''Save thermodynamic values with metadata to a JSON file.'''
    metadata = {
        'generation_timestamp': datetime.now().isoformat(),
        'primer3_py_version': primer3_py_version,
        'primer3_lib_version': thermoanalysis.get_libprimer3_version(),
        'python_version': sys.version,
        'platform': {
            'system': platform.system(),
            'release': platform.release(),
            'version': platform.version(),
            'machine': platform.machine(),
        },
        'dependencies': get_dependency_versions(),
        'git_info': get_git_info(),
    }

    output = {
        'metadata': metadata,
        'values': values,
    }

    with open(filename, 'w') as f:
        json.dump(output, f, indent=2)


if __name__ == '__main__':
    # Calculate standard values using the function from test_sequences
    standard_values = calculate_thermo_values()

    # Save standard values to a JSON file with metadata
    save_thermo_values_with_metadata(
        standard_values, 'tests/thermo_standard_values.json',
    )

    print('Gold standard values saved to tests/thermo_standard_values.json')
