=====================================================
 primer3-py: simple oligo analysis and primer design
=====================================================

.. image:: https://secure.travis-ci.org/benpruitt/primer3-py.png
  :target: https://travis-ci.org/benpruitt/primer3-py
.. image:: https://img.shields.io/pypi/l/primer3-py.png
  :target: http://www.gnu.org/licenses/gpl-2.0.html
.. image:: https://img.shields.io/pypi/v/primer3-py.png
  :target: https://pypi.python.org/pypi/primer3-py
.. image:: https://img.shields.io/pypi/dm/primer3-py.png
  :target: https://pypi.python.org/pypi/primer3-py


``primer3-py`` is a collection of Python bindings for a derivative of the 
popular Primer3 version 2.3.6 C library. The package provides a simple API 
for low-level thermodynamic calculations pertinent to oligonucleotide design 
(e.g., melting temperature) as well as a simple interface to the Primer3 design 
engine for a more holistic approach to primer design. All of the bindings
are implemented using either the Python C API or Cython, which means that they 
are highly efficient but do require initial compilation (see ``Installation``,
below).

We do not provide any additional abstraction of the Primer3 design engine, 
so we suggest that you refer to the official Primer3 
`documentation <http://primer3.sourceforge.net/>`_ for assistance.


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


Testing
-------

We have included a comprehensive test suite to compare the output of
the Python bindings with the output of the Primer3 binaries. After
building and (optionally) installing ``primer3-py`` you can run the 
tests using the ``primer3/test_primer3.py`` module::

  $ python primer3/test_primer3.py

or for memory checking with valgrind::

  $ valgrind --tool=memcheck --suppressions=valgrind-python.supp 
    --leak-check=full python primer3_test.py


API - low-level thermodynamics
------------------------------

``primer3-py`` includes a unified API for low-level thermodynamic 
calculations that are useful for routine oligonucleotide design. 

The simplest API function is ``calcTm``, which uses nearest-neighbor
thermodynamics to calculate the melting temperature of a provided DNA
sequence::

  >>> import primer3
  >>> primer3.calcTm('ATTTGGGACCAATTTGGACCAGGTT')
  57.02235387653167

Higher-level thermodynamic functions include ``calcHairpin``, 
``calcHomodimer``, ``calcHeterodimer``, and ``calcEndStability``. These 
functions perform a thermodynamic alignment to determine the characteristics 
(dH, dS, dG, Tm) of a secondary / multi-stranded structure. All four functions 
return a ThermoResult object::

  >>> from primer3 import calcHeterodimer
  >>> res = calcHeterodimer('CCGACCCTATGGGACC', 'TTGGTCCCATAAGGGTCGG')
  >>> print(res)

  ThermoResult(
    structure_found=True,
    tm=39.92795428766294, 
    ds=-370.12644214999796, 
    dh=-127200.0, 
    dg=-12405.28396717814, 
  )

  >>> print res.tm
  39.92795428766294

**Performance Note:**
Under the hood, thermodynamic calculations are performed by calls to instance
methods of an instantiated ``ThermoAnalysis`` object from the 
``thermoanalysis`` module. Users interested in performance optimization should 
examine the way in which the ``bindings.py`` module is structured. For the best
performance, you should instantiate a single ThermoAnalysis object, set
the thermodynamic parameters once, and then call the respective instance
methods with your oligo sequences to avoid the overhead of per-call
parameter parsing/setting.

For more detailed documentation and usage examples, see 
``primer3/bindings.py`` and ``primer3_test.py``.


API - primer design
-------------------

``primer3-py`` also includes C API bindings for the Primer3 design engine.
As mentioned above, we do not provide any additional "Pythonic" abstraction
of the original design process (that's up to you!) so the general 
interface is basically Boulder IO input/output in the form of Python
dictionaries. 

There are few deviations from the formats described in the Primer3 
documentation, with notable exceptions being related to index lists and 
ranges (i.e., ranges are typically provided as lists/tuples, and lists
of ranges as lists of lists or tuples of tuples). Here we highlight the
differences between the typical ``SEQUENCE_PRIMER_PAIR_OK_REGION_LIST`` 
input and the Python binding input::

  Primer3 boulder IO input:   100,50,300,50 ; 900,60,,
  Primer3 python input:       [[100,50,300,50], [900,60,-1,-1]]

Similarly, ``PRIMER_PRODUCT_SIZE_RANGE`` is provided in the following forms::

  Primer3 boulder IO input:   75-100 100-125 125-150
  Primer3 python input:       [[75,100],[100,125],[125,150]]

For more detailed documentation and usage examples, see 
``primer3/bindings.py`` and ``primer3_test.py``.


Contact and contributions
-------------------------

We are very grateful for any bug fixes or suggestions that you may have. If
you would like to report an issue or idea, or if you would like to 
contribute to the project, please visit the project's 
`Github page  (http://github.com/benpruitt/primer3-py) 
<http://github.com/benpruitt/primer3-py>`_


Licensing and citations
-----------------------
Citations should reference the `lastest Primer3 paper 
<http://nar.oxfordjournals.org/content/early/2012/06/21/nar.gks596>`_::

  Untergasser, Andreas, et al. "Primer3—new capabilities and interfaces." 
  Nucleic acids research 40.15 (2012): e115-e115.
  doi: 10.1093/nar/gks596

All project code, including the derivative Primer3 library, is licensed 
under GPLv2. The included Python and Python C API bindings are 
Copyright (c) 2014 Ben Pruitt, Nick Conway; Wyss Institute for 
Biologically Inspired Engineering.
