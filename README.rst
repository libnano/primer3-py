==============================================================================
 Primer3-py
==============================================================================

.. image:: https://secure.travis-ci.org/benpruitt/primer3-py.png
        :target: https://travis-ci.org/benpruitt/primer3-py

Primer3-py is a package of *fast* Python C API bindings for the Primer3
PCR primer design library.

**Built and tested for Python 2.7, 3.2, 3.3, and 3.4**


API Summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

primer3.bindings
 : Module interface to Primer3 C API calls.

primer3.bindings.calcTm
 : Calculate DNA sequence melting temperature using nearest-neighbor
 thermodynamics.

primer3.bindings.calcHairpin
 : Calculate DNA sequence hairpin formation thermodynamics (dS, dH, dG, Tm).

primer3.bindings.calcHomodimer
 : Calculate DNA sequence homodimer formation thermodynamics (dS, dH, dG, Tm).

primer3.bindings.calcHeterodimer
 : Calculate heterdimer thermodynamics (dS, dH, dG, Tm) for two DNA sequences.

primer3.bindings.designPrimers
 : Run the Primer3 design engine using parameters provided in a Python
 dictionary. Returns a flat dictionary of the Primer3 output.

primer3.wrappers
 : Module interface to simple subprocess wrappers for the primer3 executables.
 Provided for comparison purposes.


Development status
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Primer3-py is still undergoing active development. It is highly unlikely that
there will be any major changes to the API; most of our continuing efforts
will be focused on under-the-hood performance enhancements and improvements
to the primer/oligo design I/O.
