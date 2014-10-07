'''
===========================================================
 primer3-py: fast primer design and thermodynamic analysis
===========================================================

``primer3-py`` is a collection of Python bindings for a derivative of the
popular Primer3 C library (version 2.3.6).

See README.rst for more information and primer3_test.py for usage examples.

Installation
------------

``primer3-py`` has no external library dependencies and should compile on
most linux and OS X systems that are running Python 2.7, 3.3, or 3.4.

To build ``primer3-py`` within the package directory run::

  $ python setup.py build_ext --inplace

If you would like to install primer3-py in your local Python environment
you may do so using either pip or the setup.py script::

  $ pip install primer3-py
            or
  $ python setup.py install

'''

import os
import shutil
import subprocess
import sys

try:
    from setuptools import setup, Extension
    from setuptools.command import install_lib, sdist, build_ext
    from setuptools import log as setup_log 
except ImportError:
    from distutils.core import setup, Extension
    from distutils.command import install_lib, sdist, build_ext
    from distutils import log as setup_log 


from os.path import join as pjoin
from os.path import relpath as rpath

from Cython.Build import cythonize


with open('README.rst') as fd:
    LONG_DESCRIPTION = fd.read()


PACKAGE_PATH =          os.path.abspath(os.path.dirname(__file__))
MODULE_PATH =           pjoin(PACKAGE_PATH, 'primer3')
SRC_PATH =              pjoin(MODULE_PATH, 'src')
LIBPRIMER3_PATH =       pjoin(SRC_PATH, 'libprimer3')
THERMO_PARAMS_PATH =    pjoin(LIBPRIMER3_PATH, 'primer3_config')
KLIB_PATH =             pjoin(LIBPRIMER3_PATH, 'klib')

libprimer3_paths = [pjoin(LIBPRIMER3_PATH, 'thal.c'),
                    pjoin(LIBPRIMER3_PATH, 'oligotm.c'),
                    pjoin(LIBPRIMER3_PATH, 'p3_seq_lib.c'),
                    pjoin(LIBPRIMER3_PATH, 'libprimer3.c'),
                    pjoin(LIBPRIMER3_PATH, 'dpal.c'),
                    pjoin(SRC_PATH, 'primerdesign_helpers.c')]


def p3Clean():
    # Clean up any previous primer3 builds
    p3clean = subprocess.Popen(['make clean'], shell=True, 
                               cwd=LIBPRIMER3_PATH)
    p3clean.wait()


def p3Build():
    # Build primer3
    p3build = subprocess.Popen(['make clean; make'], shell=True, 
                               cwd=LIBPRIMER3_PATH)
    p3build.wait()

P3_BUILT = False

# Find all necessary primer3 binaries / data files to include with the package
p3_binaries = ['oligotm', 'ntthal', 'primer3_core']
p3_binary_fps = [pjoin(LIBPRIMER3_PATH, fn) for fn in p3_binaries]

thermo_files = [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
                os.walk(THERMO_PARAMS_PATH) for f in files]

p3_files = thermo_files

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


class CustomBuildExt(build_ext.build_ext):

    def run(self):
        global P3_BUILT
        # Build primer3 prior to building the extension, if not already built
        if not self.dry_run and not P3_BUILT:
            p3Clean()
            p3Build()
            P3_BUILT = True
        build_ext.build_ext.run(self)


# Build the C API and Cython extensions
primerdesign_ext = Extension(
    'primer3.primerdesign',
    sources=['primer3/src/primerdesign_py.c'] + libprimer3_paths,
    include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
    extra_compile_args=["-Wno-error=declaration-after-statement"]
)


thermoanalysis_ext = Extension(
    'primer3.thermoanalysis',
    sources=['primer3/src/thermoanalysis.pyx'] + libprimer3_paths,
    include_dirs=[LIBPRIMER3_PATH, KLIB_PATH],
    extra_compile_args=['-Wno-error=declaration-after-statement', 
                        '-Wno-unused-function']
    )

if ('build_ext' in sys.argv or 'install' in sys.argv) and not P3_BUILT:
    p3Clean()
    p3Build()
    P3_BUILT = True

thermoanalysis_ext_list = cythonize([thermoanalysis_ext])

# Include the pyd header filein the install
src_pxd = pjoin(MODULE_PATH, 'src', 'thermoanalysis.pxd')
dest_pxd = pjoin(MODULE_PATH, 'thermoanalysis.pxd')
shutil.copyfile(src_pxd, dest_pxd)
p3_files.append('thermoanalysis.pxd')


setup(
    name='primer3-py',
    version='0.3.1',
    license='GPLv2',
    author='Ben Pruitt, Nick Conway',
    author_email='bpruittvt@gmail.com',
    url='https://github.com/benpruitt/primer3-py',
    description='Python bindings for Primer3',
    long_description=LONG_DESCRIPTION,
    classifiers=[
        'Programming Language :: C',
        'Programming Language :: Cython',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
    ],
    packages=['primer3'],
    ext_modules=[primerdesign_ext] + thermoanalysis_ext_list,
    package_data={'primer3': p3_files},
    cmdclass={'install_lib': CustomInstallLib, 'sdist': CustomSdist,
              'build_ext': CustomBuildExt},
    test_suite= "primer3.tests",
    zip_safe=False
)

# Remove the extra pxd file
os.remove(dest_pxd)
