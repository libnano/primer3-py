# Buiding against cython pxd headers for primer3-py installed in your environment 

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

cdef char * foo = "ACGTACGT"

print("TM:", a.meltingTemp_c(foo))
```


and to build/install using a standard `setup.py` add the lines (fill in the ...s)

    import primer3

    ...
    my_ext = Extension(
        ...
        include_dirs=primer3.includes()
        ...)

 
 and it should work