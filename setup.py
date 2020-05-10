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

# Platform-dependent binary for `make` commands
MAKE_BIN = 'mingw32-make' if sys.platform == 'win32' else 'make'
# Indicates whether the Primer3 library has been built during install
P3_BUILT = False

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Package paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

PACKAGE_PATH = os.path.abspath(os.path.dirname(__file__))
MODULE_PATH = pjoin(PACKAGE_PATH, 'primer3')
SRC_PATH = pjoin(MODULE_PATH, 'src')
TESTS_PATH = pjoin(MODULE_PATH, 'tests')
LIBPRIMER3_PATH = pjoin(SRC_PATH, 'libprimer3')
THERMO_PARAMS_PATH = pjoin(LIBPRIMER3_PATH, 'primer3_config')
KLIB_PATH = pjoin(LIBPRIMER3_PATH, 'klib')

LIBPRIMER3_FPS = [
    rpath(pjoin(LIBPRIMER3_PATH, 'thal.c')),
    rpath(pjoin(LIBPRIMER3_PATH, 'oligotm.c')),
    rpath(pjoin(LIBPRIMER3_PATH, 'p3_seq_lib.c')),
    rpath(pjoin(LIBPRIMER3_PATH, 'libprimer3.c')),
    rpath(pjoin(LIBPRIMER3_PATH, 'dpal.c')),
    rpath(pjoin(SRC_PATH, 'primerdesign_helpers.c')),
]

LIBPRIMER3_BINARIES = ['oligotm', 'ntthal', 'primer3_core']
if sys.platform == 'win32':
    LIBPRIMER3_BINARIES = [bin_fn + '.exe' for bin_fn in LIBPRIMER3_BINARIES]

LIBPRIMER3_BINARY_FPS = [
    pjoin(LIBPRIMER3_PATH, fn) for fn in LIBPRIMER3_BINARIES
]

LIBPRIMER3_THERMO_FPS = [
    rpath(pjoin(root, fp), MODULE_PATH) for root, _, fps in
    os.walk(THERMO_PARAMS_PATH) for fp in fps
]

LIBPRIMER3_TEST_FPS = [
    rpath(pjoin(root, fp), MODULE_PATH) for root, _, fps in
    os.walk(TESTS_PATH) for fp in fps
]

PACKAGE_FPS = (
    LIBPRIMER3_THERMO_FPS +
    LIBPRIMER3_TEST_FPS +
    ['thermoanalysis.pxd', 'thermoanalysis.pyx']
)

# ~~~~~~~~~~~~~~~~~~~~~ Primer3 C library build helpers ~~~~~~~~~~~~~~~~~~~~~ #

def p3Clean():
    '''
    Run `make clean` for libprimer3 in a platform-dependent manner
    '''
    proc = subprocess.Popen(
        '{} clean'.format(MAKE_BIN),
        shell=True,
        cwd=LIBPRIMER3_PATH,
    )
    proc.wait()


def p3Build():
    '''
    Run `make clean && make` for libprimer3 in a platform-dependent manner
    '''
    proc = subprocess.Popen(
        '{} clean; {}'.format(MAKE_BIN, MAKE_BIN),
        shell=True,
        cwd=LIBPRIMER3_PATH,
    )
    proc.wait()


# Insure that the copied binaries are executable
def makeExecutable(fp):
    '''
    Adds the executable bit to the file at filepath `fp`
    '''
    mode = ((os.stat(fp).st_mode) | 0o555) & 0o7777
    setup_log.info("Adding executable bit to %s (mode is now %o)", fp, mode)
    os.chmod(fp, mode)


class CustomInstallLib(install_lib.install_lib):
    '''
    Custom library installer to ensure that libprimer3 binaries are installed
    and made executable.

    Installed and invoked internally by `setuptools`/`distutils`
    '''
    def run(self):
        global P3_BUILT
        super().run()
        if not P3_BUILT:
            p3Clean()
            p3Build()
            P3_BUILT = True
        binary_dest_fps = [
            pjoin(
                self.install_dir,
                'primer3',
                'src',
                'libprimer3',
                fn,
            ) for fn in LIBPRIMER3_BINARIES
        ]
        for src_fp, dest_fp in zip(LIBPRIMER3_BINARY_FPS, binary_dest_fps):
            shutil.copyfile(src_fp, dest_fp)
            makeExecutable(dest_fp)


class CustomSdist(sdist.sdist):
    '''
    Custom sdist packager, ensures that libprimer3 build artifacts are removed
    prior to packaging.

    Installed and invoked internally by `setuptools`/`distutils`
    '''
    def run(self):
        global P3_BUILT
        # Clean up the primer3 build prior to sdist command to remove
        # binaries and object/library files
        p3Clean()
        P3_BUILT = False
        # Remove the build Cython
        os.remove(os.path.join(MODULE_PATH, 'thermoanalysis.c'))
        super().run()


class CustomBuildClib(build_clib.build_clib):
    '''
    Custom C library builder, ensures that libprimer3 is built prior to C
    library builds.

    Installed and invoked internally by `setuptools`/`distutils`
    '''
    def run(self):
        global P3_BUILT
        # Build primer3 prior to building the extension, if not already built
        if not P3_BUILT:
            p3Clean()
            p3Build()
            P3_BUILT = True
        super().run()


# ~~~~~~~~~~~~~~~~~~~~~~ Cython / C API extension setup ~~~~~~~~~~~~~~~~~~~~~ #

if sys.platform == 'win32':
    EXTRA_COMPILE_ARGS = ['']
else:
    EXTRA_COMPILE_ARGS = [
        '-Wno-error=declaration-after-statement',
        '-Wno-unused-function',
    ]

primerdesign_ext = Extension(
    'primer3.primerdesign',
    sources=[pjoin('primer3','src','primerdesign_py.c')] + LIBPRIMER3_FPS,
    include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
    extra_compile_args=EXTRA_COMPILE_ARGS,
)


thermoanalysis_ext = Extension(
    'primer3.thermoanalysis',
    sources=[pjoin('primer3','thermoanalysis.pyx')] + LIBPRIMER3_FPS,
    include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
    extra_compile_args=EXTRA_COMPILE_ARGS,
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Ensure that libprimer3 is built prior to the `build_ext` or `install` cmds
if ('build_ext' in sys.argv or 'install' in sys.argv):
    if not P3_BUILT:
        p3Clean()
        p3Build()
        P3_BUILT = True

setup(
    name='primer3-py',
    version='0.6.1',
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
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
    ],
    packages=['primer3'],
    ext_modules=[primerdesign_ext, thermoanalysis_ext],
    package_data={'primer3': PACKAGE_FPS},
    cmdclass={
        'install_lib': CustomInstallLib,
        'sdist': CustomSdist,
        'build_clib': CustomBuildClib,
    },
    test_suite='tests',
    setup_requires=['Cython', 'setuptools>=18.0'],
    zip_safe=False,
)
