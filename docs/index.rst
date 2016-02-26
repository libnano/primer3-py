
**Primer3-py** documentation
============================

**Primer3-py** is a Python-abstracted API for the popular Primer3 library. The 
intention is to provide a simple and reliable interface for automated oligo 
analysis and design.

Routine oligo analysis is extremely simple::

    >>> import primer3 
    >>> primer3.calcTm('GTAAAACGACGGCCAGT')
    49.16808228911765
    >>> primer3.calcHairpin('CCCCCATCCGATCAGGGGG')
    ThermoResult(structure_found=True, tm=34.15, dg=337.09, dh=-36300.00, 
                 ds=-118.13, msg=)

... and `fast` (**1000X** faster than traditional subprocess wrappers)::

    In [1]: import primer3

    In [2]: %timeit primer3.calcTm('GTAAAACGACGGCCAGT')
    100000 loops, best of 3: 4.74 µs per loop

    In [3]: %timeit primer3.wrappers.calcTm('GTAAAACGACGGCCAGT')
    100000 loops, best of 3: 5.78 ms per loop

**Primer3-py** also includes bindings for the Primer3 `primer design engine` 
if you'd prefer to use an established pipeline. The IO parameters mirror those
of the original Primer3, but you don't have to deal with messy and slow file
IO for your automated workflows.

**Please note that while we provide bindings, we do not provide support for 
the Primer3 design engine. Please contact the Primer3 dev team with your
questions: http://primer3.sourceforge.net/.** 

A copy of the Primer3 2.3.7 design parameters manual can be found at:
http://git.io/vnqBx

Contents
--------

.. toctree::
   :maxdepth: 2

   quickstart
   api/index
   miscellany
   changes


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

