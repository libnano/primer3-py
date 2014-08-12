'''
===========================================================
 primer3-py: fast primer design and thermodynamic analysis
===========================================================

Primer3-py is a package of *fast* and easy-to-use Python bindings for the
Primer3 PCR primer design library.

See README.rst for more information and primer3_test.py for usage examples.

'''

import os
import subprocess

from distutils.core import setup, Extension
from os.path import join as pjoin
from os.path import relpath as rpath


with open('README.rst') as fd:
    LONG_DESCRIPTION = fd.read()

PACKAGE_PATH =          os.path.abspath(os.path.dirname(__file__))
MODULE_PATH =           pjoin(PACKAGE_PATH, 'primer3')
SRC_PATH =              pjoin(MODULE_PATH, 'src')
LIBPRIMER3_PATH =       pjoin(SRC_PATH, 'libprimer3')
KLIB_PATH =             pjoin(LIBPRIMER3_PATH, 'klib')

libprimer3_paths = [pjoin(LIBPRIMER3_PATH, 'thal.c'),
                    pjoin(LIBPRIMER3_PATH, 'oligotm.c'),
                    pjoin(LIBPRIMER3_PATH, 'p3_seq_lib.c'),
                    pjoin(LIBPRIMER3_PATH, 'libprimer3.c'),
                    pjoin(LIBPRIMER3_PATH, 'dpal.c'),
                    pjoin(SRC_PATH, 'primer3_py_helpers.c')]


# Build primer3 for subprocess bindings
p3build = subprocess.Popen(['make'], shell=True, cwd=LIBPRIMER3_PATH)
p3build.wait()

# Find all primer3 data files
p3_files = [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
            os.walk(LIBPRIMER3_PATH) for f in files]


primer3_ext = Extension(
    '_primer3',
    sources=['primer3/src/primer3_py.c'] + libprimer3_paths,
    include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
    extra_compile_args=["-Wno-error=declaration-after-statement"]
)

setup (
    name='primer3-py',
    version='0.2.5',
    license='GPLv2',
    author='Ben Pruitt, Nick Conway',
    author_email='bpruittvt@gmail.com',
    url='https://github.com/benpruitt/primer3-py',
    description='Python bindings for Primer3',
    long_description=LONG_DESCRIPTION,
    classifiers=[
        'Programming Language :: C',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
    ],
    packages=['primer3'],
    ext_modules=[primer3_ext],
    package_data={'primer3': p3_files},
)
