==============================================================================
 primer3-py
==============================================================================

.. image:: https://secure.travis-ci.org/benpruitt/primer3-py.png
        :target: https://travis-ci.org/benpruitt/primer3-py

primer3-py is a package of *fast* and easy-to-use Python bindings for the 
Primer3 PCR primer design library. The bindings interface directly with
the Primer3 library via the Python C API. 

**Built and tested for Python 2.7, 3.3, and 3.4**


API Summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[primer3.bindings]
  Module interface to Primer3 C API calls.
  
[primer3.wrappers]
  Module interface to simple subprocess wrappers for the primer3 executables. Provided for comparison purposes.


primer3.bindings.calcTm
  Calculate DNA sequence melting temperature using nearest-neighbor thermodynamics.

primer3.bindings.calcHairpin
  Calculate DNA sequence hairpin formation thermodynamics (dS, dH, dG, Tm).

primer3.bindings.calcHomodimer
  Calculate DNA sequence homodimer formation thermodynamics (dS, dH, dG, Tm).

primer3.bindings.calcHeterodimer
  Calculate heterodimer thermodynamics (dS, dH, dG, Tm) for two DNA sequences.

primer3.bindings.designPrimers
  Run the Primer3 design engine using parameters provided in a Python dictionary. Returns a flat dictionary of the Primer3 output.


**See primer3_test.py for usage examples**
