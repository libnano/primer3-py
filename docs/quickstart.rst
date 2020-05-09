Getting Started
===============


**Primer3-py** is designed to be simple to install and use.

.. contents::


Requirements
------------
**Primer3-py** is built and tested on Mac OS X and Linux systems; we do not
provide official Windows support. Python versions 2.7, 3.5-3.8 are supported.

If you don't install the latest build via pip, you might have to install
``Cython``, prior to running the ``setup.py`` script::

    $ pip install Cython


If you're attempting to build on Windows, please review the ``primer3``
documentation regarding environment requirements. You'll need to install
the latest version of the TDM-GCC MinGW Compiler if building in a
``MinGW / Mingw-w64`` environment:

https://jmeubank.github.io/tdm-gcc/

Installation
------------
If you want to install the latest stable build of **Primer3-py**, you can
install it using ``pip``::

    $ pip install primer3-py


Thermodynamic analysis
----------------------
The thermodynamic bindings include support for **Tm, homodimer, heterodimer,
hairpin,** and **3' end stability calculations**:

.. py:function:: calcTm(seq, mv_conc=50, dv_conc=0, dntp_conc=0.8, \
                        dna_conc=50, max_nn_length=60, tm_method='santalucia', \
                        salt_corrections_method='santalucia')

    Calculates the melting temperature of a DNA sequence, ``seq``. Returns the
    melting temperature (C) as a float::

        >>> primer3.calcTm('GTAAAACGACGGCCAGT')
        49.16808228911765

    Note that NN thermodynamics will be used to calculate the Tm of sequences
    up to 60 bp in length, after which point the following formula will be
    used::

        Tm = 81.5 + 16.6(log10([mv_conc])) + 0.41(%GC) - 600/length


.. py:function:: calcHairpin(seq[, mv_conc=50.0, dv_conc=0.0, dntp_conc=0.8, \
                             dna_conc=50.0, temp_c=37, max_loop=30])

    Calculates the hairpin formation thermodynamics of a DNA sequence, ``seq``.
    Returns a ``ThermoResult`` object that provides access to the thermodynamic
    characteristics of the result::

        >>> res = primer3.calcHairpin('CCCCCATCCGATCAGGGGG')
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
    primer3/src/libnano/thal.h:50).


.. py:function:: calcHomodimer(seq[, mv_conc=50.0, dv_conc=0.0, dntp_conc=0.8, \
                               dna_conc=50.0, temp_c=37, max_loop=30])

    Calculates the homodimer formation thermodynamics of a DNA sequence,
    ``seq``. Returns a ``ThermoResult`` object that provides access to the
    thermodynamic characteristics of the result (see :py:func:`calcHairpin`
    doc for more information).

    **Note that the maximum length of ``seq`` is 60 bp.** This is a cap imposed
    by Primer3 as the longest reasonable sequence length for which
    a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).


.. py:function:: calcHeterodimer(seq1, seq2[, mv_conc=50.0, dv_conc=0.0, \
                                 dntp_conc=0.8, dna_conc=50.0, temp_c=37, \
                                 max_loop=30])

    Calculates the heterodimerization thermodynamics of two DNA sequences,
    ``seq1`` and ``seq2``. Returns a ``ThermoResult`` object that provides
    access to the thermodynamic characteristics of the result
    (see :py:func:`calcHairpin` doc for more information).

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).


.. py:function:: calcEndStability(seq1, seq2[, mv_conc=50.0, dv_conc=0.0, \
                                  dntp_conc=0.8, dna_conc=50.0, temp_c=37, \
                                  max_loop=30])

    Calculates the 3' end stability of DNA sequence `seq1` against DNA
    sequence `seq2`. Returns a ``ThermoResult`` object that provides
    access to the thermodynamic characteristics of the result
    (see :py:func:`calcHairpin` doc for more information).

    **Note that at least one of the two sequences must by <60 bp in length.**
    This is a cap imposed by Primer3 as the longest reasonable sequence length
    for which a two-state NN model produces reliable results (see
    primer3/src/libnano/thal.h:50).


All of these low-level thermodynamic functions share a set of keyword arguments
used to define the parameters of the respective calculation:

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
        **max_nn_length** (int)
            Maximum length for nearest-neighbor calcs
        **tm_method** (str)
            Tm calculation method (breslauer or santalucia)
        **salt_corrections_method**
            Salt correction method (schildkraut, wczarzy, santalucia)


Primer design
-------------
**Primer3-py** includes bindings for the Primer3 primer design pipeline. The
parameters for the design process are provided as Python dictionaries that
mirror the BoulderIO input files required by the Primer3 binaries. There
are numerous examples of how to use the pipeline in the ``tests/`` directory.

For documentation regarding the input and output parameters of the pipeline,
please see the Primer3 2.3.7 documentation (the underlying library for
this package is a derivative of v2.3.7):

http://primer3.sourceforge.net/primer3_manual.htm

It is worth noting that some of the inputs deviate from the string format
described in the Primer3 documentation, with notable exceptions being related
to index lists and ranges (i.e., ranges are typically provided as lists/tuples,
and lists of ranges as lists of lists or tuples of tuples). Here we highlight
the differences between the typical ``SEQUENCE_PRIMER_PAIR_OK_REGION_LIST``
input and the Python binding input::

  Primer3 BoulderIO input:      100,50,300,50 ; 900,60,,
  Primer3-py Python input:      [[100,50,300,50], [900,60,-1,-1]]

Similarly, ``PRIMER_PRODUCT_SIZE_RANGE`` is provided in the following forms::

  Primer3 BoulderIO input:      75-100 100-125 125-150
  Primer3-py Python input:      [[75,100],[100,125],[125,150]]

Workflow
~~~~~~~~
The easiest way to run the primer design pipeline is with
:py:func:`designPrimers`. Notice that Primer3 parameters prefixed with
"SEQUENCE\_" are provided in a separate dictionary from those prefixed with
"PRIMER\_". For more advanced / modular approaches, see the :doc:`api/bindings`
documentation.

.. autofunction:: primer3.bindings.designPrimers

Example usage::

    bindings.designPrimers(
        {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTT
                                  AGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCA
                                  ACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACG
                                  CACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAG
                                  TTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAA
                                  TGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAA
                                  ATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAA
                                  TTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCA
                                  GATGTTTCCCTCTAGTAG',
            'SEQUENCE_INCLUDED_REGION': [36,342]
        },
        {
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
            'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],
                                          [150,175],[175,200],[200,225]],
        })

Advanced Installation
---------------------
Users interested in contributing to development may want to work with the
latest development build. To get the latest and greatest code, head over
`our Github repo <https://github.com/libnano/primer3-py>`_ and clone the
repo or download a tarball. Building from source is easy::

    $ python setup.py install

We recommend running ``setup.py`` with either ``build_ext --inplace`` or
``develop`` rather than ``install`` if you are testing development builds.
``build_ext --inplace`` will build the Cython and C API extensions in the
package directory without copying any files to your local environment
site-packages directory (so you can import and run tests from within the
package) and ``develop`` will build in place and then put symlinks in your
site packages directory (this will allow Primer3-py)


Testing
-------
Every commit pushed to
`our Github repo <https://github.com/libnano/primer3-py>`_ is tested to
insure it builds properly and passes our unit testing framework on
`Travis CI <https://travis-ci.org/libnano/primer3-py>`_.

If you'd like to run the tests yourself, we suggest the following workflow::

    $ git clone https://github.com/libnano/primer3-py
    $ cd primer3-py
    $ python setup.py build_ext --inplace
    $ nosetests
