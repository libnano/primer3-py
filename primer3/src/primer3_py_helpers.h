/*

primer3_py_helpers.h
~~~~~~~~~~~~~~~~~~~~

This file declares helper functions that facilitate interaction between
Python C API code and primer3 native C code. 

*/

#include    <libprimer3_mod.h>

static p3_global_settings* setGlobalParams(PyObject *self, PyObject *p3_settings_dict, char *err);