// Copyright (C) 2014. Ben Pruitt & Nick Conway; Wyss Institute
// See LICENSE for full GPLv2 license.
// #
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// #
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// #
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

/*

primer3_py_helpers.h
~~~~~~~~~~~~~~~~~~~~

This file declares helper functions that facilitate interaction between
Python C API code and primer3 native C code.

*/

#include    <libprimer3.h>

/* Set the primer3 global settings according to the key: value pairs in
 * p3_settings_dict. This almost exactly mirrors the primer3 documentation,
 * with the exception being that any file IO related keys are ignored:
 *      P3_FILE_FLAG
 *      PRIMER_EXPLAIN_FLAG
 *      PRIMER_MISPRIMING_LIBRARY
 *      PRIMER_INTERNAL_MISHYB_LIBRARY
 *      PRIMER_THERMODYNAMIC_PARAMETERS_PATH
 */
int
pdh_setGlobals(p3_global_settings *pa, PyObject *p3_settings_dict);

/* Create a sequence library from a dictionary of key: value pairs in
 * createSeqLib (key = seq name, value = seq)
 */
seq_lib*
pdh_createSeqLib(PyObject *seq_dict);

/* Copy sequence args as per primer3 docs / code. Also requires a
 * p3_global_settings pointer to set appropriate flags based on input.
 */
int
pdh_setSeqArgs(PyObject *sa_dict, seq_args *sa);

/* Parse the primer3 output to a dictionary. The dictionary will have
 * a flat structure much like a BoulderIO output file, with the field
 * names as keys and their respective values as values.
 */
PyObject*
pdh_outputToDict(const p3_global_settings *pa, const seq_args *sa,
                 const p3retval *retval);