#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
=========================================================
Primer3-py: simple oligo analysis and primer design
=========================================================

**Primer3-py** is a collection of Python bindings for a derivative of the
popular Primer3 C library (version 2.3.7).

See README.rst and http://libnano.github.io/primer3-py for more information
and primer3_test.py for usage examples.

Installation
------------

**Primer3-py** has no external library dependencies and should compile on
most linux and OS X systems that are running Python 2.7, 3.3, 3.4, or 3.5.

To build **Primer3-py** within the package directory run::

  $ python setup.py build_ext --inplace

If you would like to install primer3-py in your local Python environment
you may do so using either pip or the setup.py script::

  $ pip install primer3-py
            or
  $ python setup.py install

Create wheels with wheel installed and tar.gz::

  $ python setup.py bdist_wheel
  $ python setup.py sdist --formats=gztar
'''

import os
import shutil
import subprocess
import sys

try:
    from setuptools import setup, Extension
    from setuptools.command import install_lib, sdist, build_clib
except ImportError:
    from distutils.core import setup, Extension
    from distutils.command import install_lib, sdist, build_clib

from distutils import log as setup_log

from os.path import join as pjoin
from os.path import relpath as rpath

with open('README.rst') as fd:
    LONG_DESCRIPTION = fd.read()


PACKAGE_PATH =          os.path.abspath(os.path.dirname(__file__))
MODULE_PATH =           pjoin(PACKAGE_PATH, 'primer3')
SRC_PATH =              pjoin(MODULE_PATH, 'src')
TESTS_PATH =            pjoin(MODULE_PATH, 'tests')
LIBPRIMER3_PATH =       pjoin(SRC_PATH, 'libprimer3')
THERMO_PARAMS_PATH =    pjoin(LIBPRIMER3_PATH, 'primer3_config')
KLIB_PATH =             pjoin(LIBPRIMER3_PATH, 'klib')

libprimer3_paths = [pjoin(LIBPRIMER3_PATH, 'thal.c'),
                    pjoin(LIBPRIMER3_PATH, 'oligotm.c'),
                    pjoin(LIBPRIMER3_PATH, 'p3_seq_lib.c'),
                    pjoin(LIBPRIMER3_PATH, 'libprimer3.c'),
                    pjoin(LIBPRIMER3_PATH, 'dpal.c'),
                    pjoin(SRC_PATH, 'primerdesign_helpers.c')]


maker = 'nmake' if sys.platform == 'win32' else 'make'

def p3Clean():
    # Clean up any previous primer3 builds
    if sys.platform != 'win32':
        p3clean = subprocess.Popen(['%s clean' % (maker)], shell=True,
                                   cwd=LIBPRIMER3_PATH)
        p3clean.wait()


def p3Build():
    # Build primer3
    if sys.platform != 'win32':
        p3build = subprocess.Popen(['%s clean; %s' % (maker, maker)], shell=True,
                                   cwd=LIBPRIMER3_PATH)
        p3build.wait()

P3_BUILT = False


# Find all necessary primer3 binaries / data files to include with the package
p3_binaries = ['oligotm', 'ntthal', 'primer3_core']
if sys.platform == 'win32':
    p3_binaries = []

p3_binary_fps = [pjoin(LIBPRIMER3_PATH, fn) for fn in p3_binaries]

thermo_files = [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
                os.walk(THERMO_PARAMS_PATH) for f in files]

test_files = [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
                os.walk(TESTS_PATH) for f in files]

p3_files = thermo_files + test_files + ['thermoanalysis.pxd', 'thermoanalysis.pyx']


# Insure that the copied binaries are executable
def makeExecutable(fp):
    ''' Adds the executable bit to the file at filepath `fp`
    '''
    mode = ((os.stat(fp).st_mode) | 0o555) & 0o7777
    setup_log.info("Adding executable bit to %s (mode is now %o)", fp, mode)
    os.chmod(fp, mode)


class CustomInstallLib(install_lib.install_lib):

    def run(self):
        global P3_BUILT
        install_lib.install_lib.run(self)
        # Copy binary files over to build directory and make executable
        if not self.dry_run:
            if not P3_BUILT:
                p3Clean()
                p3Build()
                P3_BUILT = True
            new_p3_binary_fps = [pjoin(self.install_dir, 'primer3', 'src',
                                 'libprimer3', fn) for fn in p3_binaries]
            [shutil.copyfile(o, d) for o, d in zip(p3_binary_fps,
                                                   new_p3_binary_fps)]
            list(map(makeExecutable, new_p3_binary_fps))


class CustomSdist(sdist.sdist):

    def run(self):
        global P3_BUILT
        # Clean up the primer3 build prior to sdist command to remove
        # binaries and object/library files
        p3Clean()
        P3_BUILT = False
        sdist.sdist.run(self)


class CustomBuildClib(build_clib.build_clib):

    def run(self):
        global P3_BUILT
        # Build primer3 prior to building the extension, if not already built
        if not self.dry_run and not P3_BUILT:
            p3Clean()
            p3Build()
            P3_BUILT = True
        build_clib.build_clib.run(self)


# Build the C API and Cython extensions
if sys.platform == 'win32':
    extra_compile_args = ['']
else:
    extra_compile_args = ['-Wno-error=declaration-after-statement',
                        '-Wno-unused-function']

primerdesign_ext = Extension(
    'primer3.primerdesign',
    sources=[pjoin('primer3','src','primerdesign_py.c')] + libprimer3_paths,
    include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
    extra_compile_args=extra_compile_args
)


thermoanalysis_ext = Extension(
    'primer3.thermoanalysis',
    sources=[pjoin('primer3','thermoanalysis.pyx')] + libprimer3_paths,
    include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
    extra_compile_args=extra_compile_args
)



if ('build_ext' in sys.argv or 'install' in sys.argv):
    if not P3_BUILT:
        p3Clean()
        p3Build()
        P3_BUILT = True

# Insure that we don't include the built Cython module in the dist
if 'sdist' in sys.argv:
    os.remove(os.path.join(MODULE_PATH, 'thermoanalysis.c'))


setup(
    name='primer3-py',
    version='0.6.0',
    license='GPLv2',
    author='Ben Pruitt, Nick Conway',
    author_email='bpruittvt@gmail.com',
    url='https://github.com/libnano/primer3-py',
    description='Python bindings for Primer3',
    long_description=LONG_DESCRIPTION,
    classifiers=[
        'Programming Language :: C',
        'Programming Language :: Cython',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
    ],
    packages=['primer3'],
    ext_modules=[primerdesign_ext, thermoanalysis_ext],
    package_data={'primer3': p3_files},
    cmdclass={'install_lib': CustomInstallLib, 'sdist': CustomSdist,
              'build_clib': CustomBuildClib},
    test_suite='tests',
    setup_requires=['Cython', 'setuptools>=18.0'],
    zip_safe=False
)
