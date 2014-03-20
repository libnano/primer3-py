r"""
primer3-py
~~~~~~~~~~

Python C API bindings for Primer3.

Current Primer3 version included in package: 2.3.6
Support for both Python 2.7.x and Python 3.x.x

"""

import sys
import os
import subprocess
import shutil
import tarfile

from distutils.core import setup, Extension
from os.path import join as pjoin
from os.path import relpath as rpath

from buildutil import patchCfiles

root_path = os.path.abspath(os.path.dirname(__file__))
package_path = os.path.join(root_path, 'primer3')
src_path = pjoin(package_path, 'src')
primer3_path = pjoin(src_path, 'primer3-2.3.6')
primer3_src = pjoin(primer3_path, 'src')

primer3_srcs = [pjoin(primer3_src, 'thal_mod.c'),
                pjoin(primer3_src, 'oligotm.c'),
                pjoin(primer3_src, 'p3_seq_lib_mod.c'),
                pjoin(primer3_src, 'libprimer3_mod.cpp'),
                pjoin(primer3_src, 'dpal.c'),
                pjoin(src_path, 'primer3_py_helpers.c')]

if not os.path.exists(primer3_path):
    p3fd = tarfile.open(os.path.join(src_path, 'primer3-src-2.3.6.tar.gz'))
    p3fd.extractall(path=src_path)
    p3fd.close()

# Copy libprimer3_mod.c to libprimer3_mod.ccp (avoids issues w/ hash_map #include)
shutil.copy(pjoin(primer3_src, 'libprimer3.c'),
            pjoin(primer3_src, 'libprimer3.cpp'))

# Patch primer3 files w/ code for C api bindings
print('Patch Primer3 Files'.center(80, '*'))
patched_files = patchCfiles(pjoin(primer3_path), pjoin(
                            root_path, 'primer3', 'src', 'primer3_patches.c'))
print('END Patch Primer3 Files'.center(80, '*'))
# Build primer3 for subprocess bindings
print('Make Primer3'.center(80, '*'))
p3build = subprocess.Popen(['make'], shell=True, cwd=primer3_src)
p3build.wait()
print('END Make Primer3'.center(80, '*'))

# Find all primer3 data files
p3_files = [rpath(pjoin(root, f), package_path) for root, _, files in
            os.walk(primer3_path) for f in files]

# Temporary fix for OS X / Python compiler flag incompatibilities
if sys.platform == 'darwin':
    os.environ['CFLAGS'] = '-Qunused-arguments'
    os.environ['CCFLAGS'] = '-Qunused-arguments'


primer3_ext = Extension('primer3._primer3',
                        sources=['primer3/src/primer3_py.c'] + primer3_srcs,
                        include_dirs=[primer3_src]
                        )
print('Module Setup'.center(80, '*'))
setup (
    name='primer3-py',
    license='GPLv2',
    author='Ben Pruitt, Nick Conway',
    author_email='bpruittvt@gmail.com',
    description='Python 2/3 bindings for Primer3',
    long_description=__doc__,
    packages=['primer3'],
    ext_modules=[primer3_ext],
    package_data={'primer3': p3_files},
)
print('END Module Setup'.center(80, '*'))
