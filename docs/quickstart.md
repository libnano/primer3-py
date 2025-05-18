# Getting Started

**Primer3-py** is designed to be simple to install and use.


## Requirements

**Primer3-py** is built and tested on MacOS, Linux and Windows 64-bit systems; we do not provide official Windows support. Python versions 3.8 - 3.13 builds are supported.

Wheels are released for CPython versions following the [EOL model](https://devguide.python.org/versions/). Pre-built wheels are available for:
- MacOS (x86-64, arm64)
- Linux (x86-64)
- Windows (x86-64)

For other platforms (e.g., Linux aarch64), the package can be built from source but is not officially supported.


## Installation

### Standard Installation

To install the latest stable version:

```console
$ pip install primer3-py
```

### Development Installation

For development or to work with the latest code:

1. Clone the repository:
   ```console
   $ git clone https://github.com/libnano/primer3-py
   $ cd primer3-py
   ```

2. Install development dependencies:
   ```console
   $ pip install -r dev-requirements.txt
   ```

3. Install in development mode:
   ```console
   $ pip install -e .
   ```

### Windows-Specific Setup

If building on Windows, you'll need the TDM-GCC MinGW Compiler:
1. Download and install from [TDM-GCC MinGW](https://jmeubank.github.io/tdm-gcc/)
2. Ensure the compiler is in your system PATH
3. Use the installation commands above

## Thermodynamic analysis

The thermodynamic {py:mod}`primer3.bindings` include support for **Tm, homodimer, heterodimer,
hairpin,** and **3' end stability calculations**:

All of these low-level thermodynamic functions share a set of keyword arguments
used to define the parameters of the respective calculation:

```{eval-rst}
    `For all low-level calculations`:
        **mv_conc** (float/int)
            Monovalent cation concentration (mM)
        **dv_conc** (float/int)
            Divalent cation concentration (mM)
        **dntp_conc** (float/int)
            dNTP concentration (mM)
        **dna_conc** (float/int)
            DNA concentration (nM)

    `For homodimer/heterodimer/end stabilty calculation`:
        **temp_c** (int)
            Simulation temperature for dG calcs (C)
        **max_loop** (int)
            Maximum size of loops in the structure

    `For Tm calculations`:
        **dmso_conc** (float)
            Concentration of DMSO (%)
        **dmso_fact** (float)
            DMSO correction factor
        **formamide_conc** (float)
            Concentration of formamide (mol/l)
        **annealing_temp_c** (float)
            Actual annealing temperature of the PCR reaction
        **max_nn_length** (int)
            Maximum length for nearest-neighbor calcs
        **tm_method** (str)
            Tm calculation method (breslauer or santalucia)
        **salt_corrections_method**
            Salt correction method (schildkraut, wczarzy, santalucia)

```

For finer grain control of analysis, use {py:mod}`primer3.thermoanalysis`.
NOTE. camelCase methods are deprecated.  Please write all new code with
{py:class}`primer3.thermoanalysis.ThermoAnalysis` snake case methods

## Primer design

**Primer3-py** includes bindings for the Primer3 primer design pipeline. The
parameters for the design process are provided as Python dictionaries that
mirror the BoulderIO input files required by the Primer3 binaries. There
are numerous examples of how to use the pipeline in the `tests/` directory.

For documentation regarding the input and output parameters of the pipeline,
please see the [Primer3 2.6.1 documentation](https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm) (the underlying library for
this package is a derivative of v2.6.1).

It is worth noting that some of the inputs deviate from the string format
described in the Primer3 documentation, with notable exceptions being related
to index lists and ranges (i.e., ranges are typically provided as lists/tuples,
and lists of ranges as lists of lists or tuples of tuples). Here we highlight
the differences between the typical `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST`
input and the Python binding input:

```
Primer3 BoulderIO input:      100,50,300,50 ; 900,60,,
Primer3-py Python input:      [[100,50,300,50], [900,60,-1,-1]]
```

Similarly, `PRIMER_PRODUCT_SIZE_RANGE` is provided in the following forms:

```
Primer3 BoulderIO input:      75-100 100-125 125-150
Primer3-py Python input:      [[75,100],[100,125],[125,150]]
```

### Workflow

The easiest way to run the primer design pipeline is with
{py:func}`primer3.bindings.design_primers`. Notice that Primer3 parameters prefixed with
"SEQUENCE\_" are provided in a separate dictionary from those prefixed with
"PRIMER\_". For more advanced / modular approaches, see the {doc}`api/bindings`
documentation.

Example usage:

```python
bindings.design_primers(
    seq_args={
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': (
            'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT'
            'AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA'
            'ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG'
            'CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG'
            'TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA'
            'TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA'
            'ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA'
            'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA'
            'GATGTTTCCCTCTAGTAG',
        ),
        'SEQUENCE_INCLUDED_REGION': [36,342]
    },
    global_args={
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
            [75,100], [100,125], [125,150],
            [150,175], [175,200], [200,225]
        ],
    })
```

## Advanced Installation

Users interested in contributing to development may want to work with the
latest development build. To get the latest and greatest code, head over to
[our Github repo](https://github.com/libnano/primer3-py) and clone the
repo or download a tarball. Building from source is easy:

```console
# For regular installation from source
$ pip install .

# For development installation (recommended)
$ pip install -e .
```

We recommend using the development installation (`-e` flag) if you are testing or
developing `primer3-py` enhancements. This creates an "editable" installation that
links to your source code, allowing you to modify the code without reinstalling.

NOTE: If you're attempting to build on Windows, please review the `primer3`
documentation regarding environment requirements. You'll need to install
the latest version of the TDM-GCC MinGW Compiler if building in a
`MinGW / Mingw-w64` environment: [TDM-GCC MinGW Compiler](https://jmeubank.github.io/tdm-gcc/)

## Testing

Every commit pushed to [the primer3-py GitHub repo](https://github.com/libnano/primer3-py)
is tested via GitHub Actions to ensure it builds properly and passes our unit tests.

If you'd like to run the tests yourself:

```console
# Clone the repository
$ git clone https://github.com/libnano/primer3-py
$ cd primer3-py

# Create and activate a virtual environment
$ python -m venv p3p-test-env

# On Unix/macOS:
$ source p3p-test-env/bin/activate

# On Windows (run this instead of the source command above):
# p3p-test-env\Scripts\activate

# Install development dependencies
$ pip install -r dev-requirements.txt

# Build the package in-place
# This builds the Cython extensions and primer3 binaries in the package directory
# --no-build-isolation: Use local environment's packages
# --no-deps: Don't install or upgrade any package dependencies
# -e: Install in "editable" mode, creating links to the source code
$ pip install --no-build-isolation --no-deps -e .

# Run tests from the package directory
$ pytest

# When finished, deactivate the virtual environment
$ deactivate
```

The in-place build is necessary because the tests require access to the primer3 binaries
and other test resources that are part of the package directory structure but not installed
to site-packages.

Using a virtual environment is recommended to ensure a clean, isolated development environment
and avoid conflicts with other Python packages in your system.

## Contributing

Contributions are welcomed via pull requests. Contact the `primer3-py` maintainers prior to beginning your work to make sure it makes sense for the project.

By contributing, you also agree to release your code under the GPLv2.

For detailed contribution guidelines, development setup, workflows, and release process, see our [Development Guide](development.md).
