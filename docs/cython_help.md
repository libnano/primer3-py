# Building new Cython modules against `thermoanalysis.pxd`

It is of great use to build modules that depend on `primer3-py`.

To do so, first, make sure  `primer3-py` installed in your environment.

Write your code like this below:

```python
# this works with a file named thermoanalysis.pxd in the primer3 repo
cimport primer3.thermoanalysis as thermoanalysis

"""
# but below does not
from primer3 cimport thermoanalysis
as You MUST use absolute package paths in cimports i.e.
cimport A.B and not from A cimport B
"""

from primer3 import thermoanalysis

cdef thermoanalysis.ThermoAnalysis a = thermoanalysis.ThermoAnalysis()
print("MAX NN LENGTH", a.max_nn_length)

cdef char* foo = "ACGTACGT"

print("TM:", a.calc_tm_c(foo))
```


and to build/install using a standard `setup.py` add the lines (fill in the ...s)

```python
    import primer3

    ...
    my_ext = Extension(
        ...
        include_dirs=primer3.includes()
        ...)
```

and it should work.

## Notes
- Many `_ThermoAnalysis` methods (e.g. `calc_heterodimer_c`) have C string argument
`c_ascii_structure` to enable 3rd party use for structures reuse
