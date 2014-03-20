/*

primer3_py_helpers.h
~~~~~~~~~~~~~~~~~~~~

This file declares helper functions that facilitate interaction between
Python C API code and primer3 native C code. 

*/

#include    <libprimer3_mod.h>

/* Set the primer3 global settings according to the key: value pairs in 
 * p3_settings_dict. This almost exactly mirrors the primer3 documentation,
 * with the exception being that any file IO related keys are ignored:
 *      P3_FILE_FLAG
 *      PRIMER_EXPLAIN_FLAG
 *      PRIMER_MISPRIMING_LIBRARY
 *      PRIMER_INTERNAL_MISHYB_LIBRARY
 *      PRIMER_THERMODYNAMIC_PARAMETERS_PATH
 */
p3_global_settings* 
_setGlobals(PyObject *p3_settings_dict);

/* Create a sequence library from a dictionary of key: value pairs in
 * createSeqLib (key = seq name, value = seq)
 */
seq_lib* 
createSeqLib(PyObject *seq_dict);

/* Create sequence args as per primer3 docs / code. Also requires a
 * p3_global_settings pointer to set appropriate flags based on input.
 */
seq_args* 
_setSeqArgs(PyObject *sa_dict, p3_global_settings *pa);

/* Parse the primer3 output to a dictionary. The dictionary will have 
 * a flat structure much like a BoulderIO output file, with the field
 * names as keys and their respective values as values.
 */
PyObject*
p3OutputToDict(const p3_global_settings *pa, const seq_args *sa, 
               const p3retval *retval);