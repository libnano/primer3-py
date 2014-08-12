===========================================================
 primer3-py: fast primer design and thermodynamic analysis
===========================================================

.. image:: https://secure.travis-ci.org/benpruitt/primer3-py.png
        :target: https://travis-ci.org/benpruitt/primer3-py

``primer3-py`` is a collection of Python bindings for a derivative of the 
popular Primer3 C library. The package provides a simple API for low-level
thermodynamic calculations pertinent to oligonucleotide design (e.g., 
melting temperature) as well as a simple interface to the Primer3 design 
engine for a more holistic approach to primer design. 

We do not provide any additional abstraction of the Primer3 design engine, 
so we suggest that you refer to the official Primer3 
`documentation <http://primer3.sourceforge.net/>`_ for assistence.


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
tests using the ``primer3_tests.py`` module::

  $ python primer3_tests.py


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
``calcHomodimer``, and ``calcHeterodimer``. These functions perform a
thermodynamic alignment to determine the characteristics (dH, dS, dG, Tm)
of a secondary / multi-stranded structure. All three functions return
a namedtuple::

  >>> import primer3
  >>> primer3.calcHeterodimer('CCGACCCTATGGGACC', 'TTGGTCCCATAAGGGTCGG')

  thermoresult(
    tm=39.92795428766294, 
    ds=-370.12644214999796, 
    dh=-127200.0, 
    dg=-12405.28396717814, 
    align_end_1=16, 
    align_end_2=17
  )



