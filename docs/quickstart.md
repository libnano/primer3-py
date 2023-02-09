# Getting Started

**Primer3-py** is designed to be simple to install and use.

```{contents}
```

## Requirements

**Primer3-py** is built and tested on MacOS, Linux and Windows 64-bit systems; we do not provide official Windows support. Python versions 3.8 - 3.11 builds are supported.

## Installation

If you want to install the latest stable build of **Primer3-py**, you can
install it using `pip`:

```bash
$ pip install primer3-py
```

Or via `conda`:

```bash
$ conda install primer3-py
```

**NOTE**: We support wheel builds for PyPi for the 3 most recent CPython versions. Target platforms for wheels are MacOS `x86-64` `arm64`, Linux `x86-64`, and Windows `x86-64`.

If your Python version and platform fall outside this such as Linux `aarch64` it is confirmed `primer3-py` builds on this platform but it is not supported as our build GitHub actions runners do not run these builds expediently.

## Thermodynamic analysis

The thermodynamic {py:mod}`bindings` include support for **Tm, homodimer, heterodimer,
hairpin,** and **3' end stability calculations**:

```{eval-rst}
.. py:function:: calc_tm(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, \
                        dna_conc=50, dmso_conc=0.0, dmso_fact=0.6, \
                        formamide_conc=0.8, annealing_temp_c=-10.0, \
                        max_nn_length=60, tm_method='santalucia', \
                        salt_corrections_method='santalucia')

    Calculates the melting temperature of a DNA sequence, ``seq``. Returns the
    melting temperature (C) as a float::

        >>> primer3.calc_tm('GTAAAACGACGGCCAGT')
        49.16808228911765

    Note that NN thermodynamics will be used to calculate the Tm of sequences
    up to 60 bp in length, after which point the following formula will be
    used::

        Tm = 81.5 + 16.6(log10([mv_conc])) + 0.41(%GC) - 600/length

```

```{eval-rst}
.. py:function:: calc_hairpin(seq, mv_conc=50.0, dv_conc=0.0, dntp_conc=0.8, \
                             dna_conc=50.0, temp_c=37, max_loop=30, \
                             output_structure=False)

    Calculates the hairpin formation thermodynamics of a DNA sequence, ``seq``.
    Returns a :class:`ThermoResult` object that provides access to the
    thermodynamic characteristics of the result::

        >>> res = primer3.calc_hairpin('CCCCCATCCGATCAGGGGG')
        >>> print(res)
        ThermoResult(structure_found=True, tm=34.15, dg=337.09, dh=-36300.00,
                     ds=-118.13, msg=)
        >>> print(res.tm)
        34.14640571476207
        >>> print('%f, %f, %f' % (res.dg, res.dh, res.ds))
        337.086509, -36300.000000, -118.126992

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    ``primer3/src/libnano/thal.h:59``).

```

```{eval-rst}
.. py:function:: calc_homodimer(seq, mv_conc=50.0, dv_conc=0.0, dntp_conc=0.8, \
                               dna_conc=50.0, temp_c=37, max_loop=30, \
                               output_structure=False)

    Calculates the homodimer formation thermodynamics of a DNA sequence,
    ``seq``. Returns a :class:`ThermoResult` object that provides access to the
    thermodynamic characteristics of the result (see :py:func:`calc_hairpin`
    doc for more information).

    **Note that the maximum length of ``seq`` is 60 bp.** This is a cap imposed
    by Primer3 as the longest reasonable sequence length for which
    a two-state NN model produces reliable results (see
    ``primer3/src/libprimer3/thal.h:59``).

```

```{eval-rst}
.. py:function:: calc_heterodimer(seq1, seq2, mv_conc=50.0, dv_conc=0.0, \
                                 dntp_conc=0.8, dna_conc=50.0, temp_c=37, \
                                 max_loop=30, output_structure=False)

    Calculates the heterodimerization thermodynamics of two DNA sequences,
    ``seq1`` and ``seq2``. Returns a :class:`ThermoResult` object that provides
    access to the thermodynamic characteristics of the result
    (see :py:func:`calc_hairpin` doc for more information).

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    ``primer3/src/libprimer3/thal.h:59``).

```

```{eval-rst}
.. py:function:: calc_end_stability(seq1, seq2, mv_conc=50.0, dv_conc=0.0, \
                                  dntp_conc=0.8, dna_conc=50.0, temp_c=37, \
                                  max_loop=30)

    Calculates the 3' end stability of DNA sequence ``seq1`` against DNA
    sequence ``seq2``. Returns a :class:`ThermoResult` object that provides
    access to the thermodynamic characteristics of the result
    (see :py:func:`calc_hairpin` doc for more information).

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    ``primer3/src/libprimer3/thal.h:59``).

```

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
{py:func}`design_primers`. Notice that Primer3 parameters prefixed with
"SEQUENCE\_" are provided in a separate dictionary from those prefixed with
"PRIMER\_". For more advanced / modular approaches, see the {doc}`api/bindings`
documentation.

```{eval-rst}
.. autofunction:: primer3.bindings.design_primers
```

Example usage:

```
bindings.design_primers(
    seq_args={
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT'
                             'AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA'
                             'ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG'
                             'CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG'
                             'TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA'
                             'TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA'
                             'ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA'
                             'TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA'
                             'GATGTTTCCCTCTAGTAG',
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
            [75,100],[100,125],[125,150],
            [150,175],[175,200],[200,225]
        ],
    })
```

## Advanced Installation

Users interested in contributing to development may want to work with the
latest development build. To get the latest and greatest code, head over
[our Github repo](https://github.com/libnano/primer3-py) and clone the
repo or download a tarball. Building from source is easy.

If you don't install the latest build via pip or conda, you might have to install
`Cython`, prior to running the `setup.py` script:

```
$ pip install Cython
```

Or via `conda`:

```
$ conda install Cython
```

Then run:

```
$ python setup.py install
```

or if you are developing `primer3-py` enhancements:

```
$ python setup.py build_ext --inplace
```

We recommend running `setup.py` with either `build_ext --inplace` or
`develop` rather than `install` if you are testing development builds.
`build_ext --inplace` will build the Cython and C API extensions in the
package directory without copying any files to your local environment
site-packages directory (so you can import and run tests from within the
package) and `develop` will build in place and then put symlinks in your
site packages directory (this will allow Primer3-py)

NOTE: If you're attempting to build on Windows, please review the `primer3`
documentation regarding environment requirements. You'll need to install
the latest version of the TDM-GCC MinGW Compiler if building in a
`MinGW / Mingw-w64` environment: [TDM-GCC MinGW Compiler](https://jmeubank.github.io/tdm-gcc/)

## Testing

Every commit pushed to
[the primer3-py GitHub repo](https://github.com/libnano/primer3-py) is tested to
ensure it builds properly and passes our unit testing framework as a GitHub action

If you'd like to run the tests yourself, we suggest the following workflow:

```
$ git clone https://github.com/libnano/primer3-py
$ cd primer3-py
$ python setup.py build_ext --inplace
$ pytest
```

NOTE: `pip` / `conda` install `pytest` if not in your environment


## Contributing

Contributions are welcomed via pull requests.

Contact the `primer3-py` maintainers prior to beginning your work to make sure
it makes sense for the project.

By contributing, you also agree to release your code under the GPLv2

After a successful PR will be listed under the [contributors](https://github.com/libnano/primer3-py/graphs/contributors).


### Forking

A forking workflow is preferred for all pull requests.

### Branch naming

Branch naming is preferred to use the format:

```
<GitHub user-name>-<short keyword description of change>
```

Keep branch names not too long.  A good example would be for the user `grinner`
for a documentation update for the 1.0.0 staging branch:

```bash
$ git checkout -b grinner-docs-update-1.0.0-pass-01
```

With the trailing 01 indicative of it being part of several potential

Another example pass that focuses on code clarity comments would be:

```bash
$ git checkout -b grinner-code-clarity-and-comments
```

### Development

Development requires the use of C Python 3.8+, [pytest](https://docs.pytest.org) and
[pre-commit](https://pre-commit.com) as they are used to build and run primer3-py
code CI in the GitHub Action.

Install these dependencies in your python development environment
(`virtualenv`, `conda`, etc):

```bash
$ pip install cython pre-commit pytest
# or
$ conda install cython pre-commit pytest
```

Install `pre-commit` in repo the with:

```bash
$ pre-commit install
```

To ensure the git hook is excecuted on every commit.

### Pull Requests

Pull Requests should meet the following requirements:

1. Excellent PR description describing all changes made. Please use markdown syntax highlighting to help readability.
2. If change is code related, have test coverage for the changes implemented.
3. Attempt to make the PR 1 commit only. Multiple are OK if it helps illustrate the change better.
4. Commit messages should describe the changes.
5. Provided you contact the maintainers in advance, theu will code review your PR, provide feedback and squash merge your code on approval.

**TIP**: Interactive `rebase` is helpful to fix old commit messages.
For example, run:

```bash
$ git rebase -i HEAD~2
```

To rebase the last 2 commits. Use `s` to mark the most recent commit(s), save, then
modify the collective commit messages to update poor commit messages.
