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
primer3_py_helpers.c
~~~~~~~~~~~~~~~~~~~~

Helper functions to facilitate information passing between the Python C API
and the Primer3 library.

*/

#include    <string.h>
#include    <stdio.h>
#include    <Python.h>
#include    <libprimer3.h>
#include    "p3_seq_lib.h"


#if PY_MAJOR_VERSION < 3
/* see http://python3porting.com/cextensions.html */
    #ifdef PyLong_Check
        #undef PyLong_Check
    #endif
    #ifdef PyLong_AsLong
        #undef PyLong_AsLong
    #endif
    #define PyLong_Check PyInt_Check
    #define PyLong_AsLong PyInt_AsLong
#endif


// Check python dictionary `d` for key `k` and (if it exists) assign the
// value to Py_Object o. Return 1 on success or 0 if the key is not in the dict
#define DICT_GET_OBJ(o, d, k) ((o = PyDict_GetItemString(d, k)) != NULL)

// Wraps DICT_GET_OBJ and takes the dictionary value object, extracts a long,
// casts it to an int and assigns its value to `st`
#define DICT_GET_AND_ASSIGN_INT(o, d, k, st)                                   \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        if (!PyLong_Check(o)) {                                                \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not an integer.", k);              \
            return -1;                                                         \
        }                                                                      \
        st = (int)PyLong_AsLong(o);                                            \
    }

// Wraps DICT_GET_OBJ and takes the dictionary value object, extracts a long,
// casts it to `type` and assigns its value to `st`
#define DICT_GET_AND_ASSIGN_INT_TYPE(o, d, k, st, t)                           \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        if (!PyLong_Check(o)) {                                                \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not of type integer.", k);         \
            return -1;                                                         \
        }                                                                      \
        st = (t)PyLong_AsLong(o);                                              \
    }

// Wraps DICT_GET_OBJ and takes the dictionary value object, extracts a double,
// and assigns its value to `st`
#define DICT_GET_AND_ASSIGN_DOUBLE(o, d, k, st)                                \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        if (!PyFloat_Check(o) && !PyLong_Check(o)) {                           \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not of type float or integer.", k);\
            return -1;                                                         \
        }                                                                      \
        st = PyFloat_AsDouble(o);                                              \
    }

// Wraps DICT_GET_OBJ and takes the dictionary value object, exposes the
// internal string buffer, and copies the string into newly allocated memory
// pointed to by st (note that a pointer to st is passed in, so we need to
// further dereference it to get to the actual string pointer).
// PyString_AsString exposes the internal buffer of the string (null terminated)
// It must not be changed or freed, so we have to malloc new memory for the
// param value.
#if PY_MAJOR_VERSION < 3
    #define DICT_GET_AND_COPY_STR(o, d, k, st, tc, ss)                         \
        if (DICT_GET_OBJ(o, d, k)) {                                           \
            int from_int = 0;                                                  \
            if (!PyString_Check(o)){                                           \
                if (PyLong_Check(o)) {                                         \
                    *st = (char *) malloc(20 * sizeof(char));                  \
                    sprintf(*st, "%d", (int)PyLong_AsLong(o));                 \
                    from_int = 1;                                              \
                } else {                                                       \
                    PyErr_Format(PyExc_TypeError,                              \
                                 "Value of %s is not of type string", k);      \
                return -1;}                                                    \
            }                                                                  \
            if (from_int == 0) {                                               \
                if (PyString_AsStringAndSize(o, &tc, &ss) == -1) {             \
                    return -1;}                                                \
                *st = (char *) malloc((ss + 1) * sizeof(char));                \
                if (*st == NULL) {                                             \
                    PyErr_Format(PyExc_IOError,                                \
                        "Could not allocate memory while copying %s", k);      \
                    return -1;}                                                \
                memcpy(*st, tc, (int)(ss + 1));                                \
            }                                                                  \
        }
#else
    #define DICT_GET_AND_COPY_STR(o, d, k, st, tc, ss)                         \
        if (DICT_GET_OBJ(o, d, k)) {                                           \
            int from_int = 0;                                                  \
            if (PyUnicode_Check(o)) {                                          \
                tc = (char *)PyUnicode_AsUTF8AndSize(o, &ss);                  \
            } else if (PyBytes_Check(o)){                                      \
                if (PyBytes_AsStringAndSize(o, &tc, &ss) == -1) {              \
                    return -1;}                                                \
            } else {                                                           \
                if (PyLong_Check(o)) {                                         \
                    *st = (char *) malloc(20 * sizeof(char));                  \
                    sprintf(*st, "%d", (int)PyLong_AsLong(o));                 \
                    from_int = 1;                                              \
                } else {                                                       \
                    PyErr_Format(PyExc_TypeError,                              \
                        "Value of %s is not of type unicode or bytes", k);     \
                    return -1;}                                                \
            }                                                                  \
            if (tc == NULL){                                                   \
                    PyErr_Format(PyExc_TypeError,                              \
                            "Error processing string in %s", k);               \
                    return -1;                                                 \
                }                                                              \
            if (from_int == 0) {                                               \
                *st = (char *) malloc((ss + 1 ) * sizeof(char));               \
                if (*st == NULL) {                                             \
                    PyErr_Format(PyExc_IOError,                                \
                        "Could not allocate memory while copying %s", k);      \
                    return -1;}                                                \
                memcpy(*st, tc, (int)(ss + 1));                                \
            }                                                                  \
        }
#endif

#define DICT_GET_AND_COPY_ARRAY(o, d, k, st, arr_len)                          \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        int i;                                                                 \
        if (!PySequence_Check(o)){                                             \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not a sequence object", k);        \
            return -1;                                                         \
        } else {                                                               \
            int *arr = NULL;                                                   \
            PyObject *seq_item;                                                \
            *arr_len = (int)PySequence_Size(o);                                \
            arr = (int*) malloc(*arr_len*sizeof(int));                         \
            if (arr == NULL) {                                                 \
                PyErr_Format(PyExc_IOError,                                    \
                            "Could not allocate memory while copying %s", k);  \
                return -1;}                                                    \
            for (i=0; i < *arr_len; i++) {                                     \
                seq_item = PySequence_GetItem(o, i);                           \
                arr[i] = (int)PyLong_AsLong(seq_item);                         \
                Py_DECREF(seq_item);                                           \
            }                                                                  \
            *st = arr;                                                         \
        }                                                                      \
    }                                                                          \

#define DICT_GET_AND_COPY_ARRAY_INTO_ARRAY(o, d, k, st, arr_len)               \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        int i;                                                                 \
        PyObject *seq_item;                                                    \
        if (!PySequence_Check(o)){                                             \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not a sequence object", k);        \
            return -1;                                                         \
        *arr_len = (int)PySequence_Size(o);                                    \
        for (i=0; i < *arr_len; i++) {                                         \
            seq_item = PySequence_GetItem(o, i);                               \
            *st[i] = (int)PyLong_AsLong(seq_item);                             \
            Py_DECREF(seq_item);                                               \
        }                                                                      \
    }                                                                          \

#define DICT_GET_AND_COPY_TO_INTERVAL_ARRAY(o, d, k, st)                       \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        int i, too_long, arr_len, flat_list=0;                                 \
        PyObject *sub_seq, *seq_item1, *seq_item2;                             \
        if (!PySequence_Check(o)){                                             \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not a sequence object", k);        \
            return -1;                                                         \
        }                                                                      \
        st.count = 0;                                                          \
        arr_len = (int)PySequence_Size(o);                                     \
        if (arr_len == 2) {                                                    \
            seq_item1 = PySequence_GetItem(o, 0);                              \
            seq_item2 = PySequence_GetItem(o, 1);                              \
            if (PyLong_Check(seq_item1) && PyLong_Check(seq_item2)) {          \
                p3_add_to_interval_array(&st,                                  \
                     (int)PyLong_AsLong(seq_item1),                            \
                     (int)PyLong_AsLong(seq_item2));                           \
                flat_list = 1;                                                 \
            } else if ((PyLong_Check(seq_item1) && !PyLong_Check(seq_item2))|| \
                       (!PyLong_Check(seq_item1) && PyLong_Check(seq_item2))){ \
                PyErr_Format(PyExc_TypeError,                                  \
                                "Value of %s is a mixture of integers"         \
                                " and other objects", k);                      \
                Py_DECREF(seq_item1);                                          \
                Py_DECREF(seq_item1);                                          \
                return -1;                                                     \
            }                                                                  \
        }                                                                      \
        if (!flat_list) {                                                      \
            for (i = 0; i < arr_len; i++) {                                    \
                sub_seq = PySequence_GetItem(o, i);                            \
                seq_item1 = PySequence_GetItem(sub_seq, 0);                    \
                seq_item2 = PySequence_GetItem(sub_seq, 1);                    \
                if (!PyLong_Check(seq_item1) || !PyLong_Check(seq_item2)) {    \
                    PyErr_Format(PyExc_TypeError,                              \
                                    "Value of %s is not a sequence object"     \
                                    " comprised of two integers", k);          \
                    return -1;                                                 \
                }                                                              \
                too_long = p3_add_to_interval_array(&st,                       \
                     (int)PyLong_AsLong(seq_item1),                            \
                     (int)PyLong_AsLong(seq_item2));                           \
                if (too_long) {                                                \
                    PyErr_Format(PyExc_IOError, "Too many elements for tag "   \
                                                "%s", k);                      \
                    return -1;                                                 \
                }                                                              \
                Py_DECREF(seq_item1);                                          \
                Py_DECREF(seq_item2);                                          \
                Py_DECREF(sub_seq);                                            \
            }                                                                  \
        }                                                                      \
    }                                                                          \


int
pdh_setGlobals(p3_global_settings *pa, PyObject *p3s_dict) {
    /* Creates a new p3_global_settings struct and initializes it with
     * defaults using p3_create_global_settings() from libprimer3.c.
     * Parses the user-provided settings from p3_settings_dict and
     * overwrites the defaults (note that minimal error checking is
     * performed in this function). If there is an error during the process
     * (e.g., a param is not of the correct type), the python error string will
     * be set and the function will return NULL.
     */

    // p3_global_settings      *pa;
    PyObject                *p_obj, *p_obj2, *p_obj3, *p_obj4;
    int                     i;
    Py_ssize_t              str_size;
    char                    *temp_char=NULL, *task_tmp=NULL;


    // if (!(pa = p3_create_global_settings())) {
    //     PyErr_SetString(PyExc_IOError, "Could not allocate memory for p3 globals");
    //     return NULL;
    // }

    /* Note that some of the documented primer3 parameters are ignored in
     * this function. Specifically, any parameters related to file IO (as
     * well as thermodynamic parameter files) are ignored:
     *
     *      P3_FILE_FLAG
     *      PRIMER_EXPLAIN_FLAG
     *      PRIMER_MISPRIMING_LIBRARY
     *      PRIMER_INTERNAL_MISHYB_LIBRARY
     *      PRIMER_THERMODYNAMIC_PARAMETERS_PATH
     *
     * Otherwise, all parameters are generally parsed in the same order as
     * in read_boulder.c in the primer3 source. This code will permit some
     * edge cases that are more directly addressed in read_boulder.c (for
     * example, if you provide both PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE /
     * PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE and
     * PRIMER_MIN_THREE_PRIME_DISTANCE, the latter value will be used for both
     * the left and right three prime distance parameters.
     */
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_OPT_SIZE", pa->p_args.opt_size);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_SIZE", pa->p_args.min_size);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MAX_SIZE", pa->p_args.max_size);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MAX_POLY_X", pa->p_args.max_poly_x);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_OPT_TM", pa->p_args.opt_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_OPT_GC_PERCENT", pa->p_args.opt_gc_content);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MIN_TM", pa->p_args.min_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_TM", pa->p_args.max_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_DIFF_TM", pa->max_diff_tm);
    DICT_GET_AND_ASSIGN_INT_TYPE(p_obj, p3s_dict, "PRIMER_TM_FORMULA", pa->tm_santalucia, tm_method_type);
    DICT_GET_AND_ASSIGN_INT_TYPE(p_obj, p3s_dict, "PRIMER_SALT_CORRECTIONS", pa->salt_corrections, salt_correction_type);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MIN_GC", pa->p_args.min_gc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_GC", pa->p_args.max_gc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_SALT_MONOVALENT", pa->p_args.salt_conc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_SALT_DIVALENT", pa->p_args.divalent_conc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_DNTP_CONC", pa->p_args.dntp_conc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_DNA_CONC", pa->p_args.dna_conc);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MAX_NS_ACCEPTED", pa->p_args.num_ns_accepted);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_PRODUCT_OPT_SIZE", pa->product_opt_size);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_SELF_ANY", pa->p_args.max_self_any);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_SELF_END", pa->p_args.max_self_end);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_SELF_ANY_TH", pa->p_args.max_self_any_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_SELF_END_TH", pa->p_args.max_self_end_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_HAIRPIN_TH", pa->p_args.max_hairpin_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_COMPL_ANY", pa->pair_compl_any);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_COMPL_END", pa->pair_compl_end);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_COMPL_ANY_TH", pa->pair_compl_any_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_COMPL_END_TH", pa->pair_compl_end_th);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_PICK_ANYWAY", pa->pick_anyway);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_GC_CLAMP", pa->gc_clamp);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MAX_END_GC", pa->max_end_gc);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_LIBERAL_BASE", pa->liberal_base);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_FIRST_BASE_INDEX", pa->first_base_index);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_NUM_RETURN", pa->num_return);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_QUALITY", pa->p_args.min_quality);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_END_QUALITY", pa->p_args.min_end_quality);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE", pa->min_left_three_prime_distance);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE", pa->min_right_three_prime_distance);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_THREE_PRIME_DISTANCE", pa->min_left_three_prime_distance);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_THREE_PRIME_DISTANCE", pa->min_right_three_prime_distance);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_QUALITY_RANGE_MIN", pa->quality_range_min);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_QUALITY_RANGE_MAX", pa->quality_range_max);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PRODUCT_MAX_TM", pa->product_max_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PRODUCT_MIN_TM", pa->product_min_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PRODUCT_OPT_TM", pa->product_opt_tm);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_SEQUENCING_LEAD", pa->sequencing.lead);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_SEQUENCING_SPACING", pa->sequencing.spacing);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_SEQUENCING_INTERVAL", pa->sequencing.interval);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_SEQUENCING_ACCURACY", pa->sequencing.accuracy);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION", pa->min_5_prime_overlap_of_junction);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION", pa->min_3_prime_overlap_of_junction);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_PICK_RIGHT_PRIMER", pa->pick_right_primer);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_PICK_INTERNAL_OLIGO", pa->pick_internal_oligo);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_PICK_LEFT_PRIMER", pa->pick_left_primer);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_INTERNAL_OPT_SIZE", pa->o_args.opt_size);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_SIZE", pa->o_args.max_size);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_INTERNAL_MIN_SIZE", pa->o_args.min_size);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_POLY_X", pa->o_args.max_poly_x);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_OPT_TM", pa->o_args.opt_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_OPT_GC_PERCENT", pa->o_args.opt_gc_content);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_TM", pa->o_args.max_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MIN_TM", pa->o_args.min_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MIN_GC", pa->o_args.min_gc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_GC", pa->o_args.max_gc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_SALT_MONOVALENT", pa->o_args.salt_conc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_SALT_DIVALENT", pa->o_args.divalent_conc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_DNTP_CONC", pa->o_args.dntp_conc);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_DNA_CONC", pa->o_args.dna_conc);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_NS_ACCEPTED", pa->o_args.num_ns_accepted);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_INTERNAL_MIN_QUALITY", pa->o_args.min_quality);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_SELF_ANY", pa->o_args.max_self_any);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_SELF_END", pa->p_args.max_self_end);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_SELF_ANY_TH", pa->o_args.max_self_any_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_SELF_END_TH", pa->o_args.max_self_end_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_HAIRPIN_TH", pa->o_args.max_hairpin_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_LIBRARY_MISPRIMING", pa->p_args.max_repeat_compl);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_MAX_LIBRARY_MISHYB", pa->o_args.max_repeat_compl);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING", pa->pair_repeat_compl);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_TEMPLATE_MISPRIMING", pa->p_args.max_template_mispriming);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_TEMPLATE_MISPRIMING_TH", pa->p_args.max_template_mispriming_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING", pa->pair_max_template_mispriming);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH", pa->pair_max_template_mispriming_th);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS", pa->lib_ambiguity_codes_consensus);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INSIDE_PENALTY", pa->inside_penalty);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_OUTSIDE_PENALTY", pa->outside_penalty);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_MAX_END_STABILITY", pa->max_end_stability);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_LOWERCASE_MASKING", pa->lowercase_masking);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT", pa->thermodynamic_oligo_alignment);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT", pa->thermodynamic_template_alignment);
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_MUST_MATCH_FIVE_PRIME", &pa->p_args.must_match_five_prime, temp_char, str_size);
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_MUST_MATCH_THREE_PRIME", &pa->p_args.must_match_three_prime, temp_char, str_size);
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME", &pa->o_args.must_match_five_prime, temp_char, str_size);
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME", &pa->o_args.must_match_three_prime, temp_char, str_size);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_TM_GT", pa->p_args.weights.temp_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_TM_LT", pa->p_args.weights.temp_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_GC_PERCENT_GT", pa->p_args.weights.gc_content_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_GC_PERCENT_LT", pa->p_args.weights.gc_content_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SIZE_LT", pa->p_args.weights.length_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SIZE_GT", pa->p_args.weights.length_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SELF_ANY", pa->p_args.weights.compl_any);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SELF_END", pa->p_args.weights.compl_end);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SELF_ANY_TH", pa->p_args.weights.compl_any_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SELF_END_TH", pa->p_args.weights.compl_end_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_HAIRPIN_TH", pa->p_args.weights.hairpin_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_NUM_NS", pa->p_args.weights.num_ns);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_LIBRARY_MISPRIMING", pa->p_args.weights.repeat_sim);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SEQ_QUAL", pa->p_args.weights.seq_quality);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_END_QUAL", pa->p_args.weights.end_quality);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_POS_PENALTY", pa->p_args.weights.pos_penalty);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_END_STABILITY", pa->p_args.weights.end_stability);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_TEMPLATE_MISPRIMING", pa->p_args.weights.template_mispriming);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_TEMPLATE_MISPRIMING_TH", pa->p_args.weights.template_mispriming_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_TM_GT", pa->o_args.weights.temp_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_TM_LT", pa->o_args.weights.temp_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_GC_PERCENT_GT", pa->o_args.weights.gc_content_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_GC_PERCENT_LT", pa->o_args.weights.gc_content_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_SIZE_LT", pa->o_args.weights.length_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_SIZE_GT", pa->o_args.weights.length_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_SELF_ANY", pa->o_args.weights.compl_any);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_SELF_END", pa->o_args.weights.compl_end);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_SELF_ANY_TH", pa->o_args.weights.compl_any_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_SELF_END_TH", pa->o_args.weights.compl_end_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_HAIRPIN_TH", pa->o_args.weights.hairpin_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_NUM_NS", pa->o_args.weights.num_ns);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_LIBRARY_MISHYB", pa->o_args.weights.repeat_sim);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_SEQ_QUAL", pa->o_args.weights.seq_quality);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_INTERNAL_WT_END_QUAL", pa->o_args.weights.end_quality);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_TEMPLATE_MISPRIMING_TH", pa->o_args.weights.template_mispriming_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_PR_PENALTY", pa->pr_pair_weights.primer_quality);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_IO_PENALTY", pa->pr_pair_weights.io_quality);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_DIFF_TM", pa->pr_pair_weights.diff_tm);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_COMPL_ANY", pa->pr_pair_weights.compl_any);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_COMPL_END", pa->pr_pair_weights.compl_end);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_COMPL_ANY_TH", pa->pr_pair_weights.compl_any_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_COMPL_END_TH", pa->pr_pair_weights.compl_end_th);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_PRODUCT_TM_LT", pa->pr_pair_weights.product_tm_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_PRODUCT_TM_GT", pa->pr_pair_weights.product_tm_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_PRODUCT_SIZE_GT", pa->pr_pair_weights.product_size_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_PRODUCT_SIZE_LT", pa->pr_pair_weights.product_size_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_LIBRARY_MISPRIMING", pa->pr_pair_weights.repeat_sim);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING", pa->pr_pair_weights.template_mispriming);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH", pa->pr_pair_weights.template_mispriming_th);

    if DICT_GET_OBJ(p_obj, p3s_dict, "PRIMER_PRODUCT_SIZE_RANGE") {
        int flat_list = 0;
        int product_size_range_list_length = (int)PySequence_Length(p_obj);
        if (!PySequence_Check(p_obj)) {
            PyErr_SetString(PyExc_TypeError,\
                "Value of \"PRIMER_PRODUCT_SIZE_RANGE\" is not a list or tuple");
            return -1;
        }
        if (product_size_range_list_length == 2) {
            p_obj2 = PySequence_GetItem(p_obj, 0);
            p_obj3 = PySequence_GetItem(p_obj, 1);
            if (PyLong_Check(p_obj2)) {
                pa->pr_min[0] = (int)PyLong_AsLong(p_obj2);
                if (PyLong_Check(p_obj3)) {
                    pa->pr_max[0] = (int)PyLong_AsLong(p_obj3);
                    flat_list = 1;
                } else {
                    PyErr_Format(PyExc_TypeError,\
                        "\"PRIMER_PRODUCT_SIZE_RANGE\" contains mixed sequence objects and integers");
                    Py_DECREF(p_obj2);
                    Py_DECREF(p_obj3);
                    return -1;
                }
            }
            Py_DECREF(p_obj2);
            Py_DECREF(p_obj3);
        }
        if (!flat_list) {
            for (i=0; i < product_size_range_list_length; i++){
                p_obj2 = PySequence_GetItem(p_obj, i);
                if (!PySequence_Check(p_obj2)){
                    PyErr_Format(PyExc_TypeError,\
                        "Object at index %d of \"PRIMER_PRODUCT_SIZE_RANGE\" is not a list or tuple", i);
                    return -1;
                }
                if (PySequence_Length(p_obj2) != 2) {
                    PyErr_Format(PyExc_TypeError,\
                        "Object at index %d of \"PRIMER_PRODUCT_SIZE_RANGE\" is not of length 2", i);
                    return -1;
                }
                p_obj3 = PySequence_GetItem(p_obj2, 0);
                p_obj4 = PySequence_GetItem(p_obj2, 1);
                if ((pa->pr_min[i] = (int)PyLong_AsLong(p_obj3)) == -1) {
                    PyErr_Format(PyExc_TypeError,\
                        "Object 1 at index %d of \"PRIMER_PRODUCT_SIZE_RANGE\" is not an integer", i);
                    Py_DECREF(p_obj2);
                    Py_DECREF(p_obj3);
                    Py_DECREF(p_obj4);
                    return -1;
                }
                if ((pa->pr_max[i] = (int)PyLong_AsLong(p_obj4)) == -1) {
                    PyErr_Format(PyExc_TypeError,\
                        "Object 2 at index %d of \"PRIMER_PRODUCT_SIZE_RANGE\" is not an integer", i);
                    Py_DECREF(p_obj2);
                    Py_DECREF(p_obj3);
                    Py_DECREF(p_obj4);
                    return -1;
                }
                Py_DECREF(p_obj2);
                Py_DECREF(p_obj3);
                Py_DECREF(p_obj4);
            }
            pa->num_intervals = i;
        }
    }

    // Handler primer task
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_TASK", &task_tmp, temp_char, str_size);

    // Directly from read_boulder.c
    if (task_tmp != NULL) {
        if (!strcmp_nocase(task_tmp, "pick_pcr_primers")) {
          pa->primer_task = generic;
          pa->pick_left_primer = 1;
          pa->pick_right_primer = 1;
          pa->pick_internal_oligo = 0;
        } else if (!strcmp_nocase(task_tmp, "pick_pcr_primers_and_hyb_probe")) {
          pa->primer_task = generic;
          pa->pick_left_primer = 1;
          pa->pick_right_primer = 1;
          pa->pick_internal_oligo = 1;
        } else if (!strcmp_nocase(task_tmp, "pick_left_only")) {
          pa->primer_task = generic;
          pa->pick_left_primer = 1;
          pa->pick_right_primer = 0;
          pa->pick_internal_oligo = 0;
        } else if (!strcmp_nocase(task_tmp, "pick_right_only")) {
          pa->primer_task = generic;
          pa->pick_left_primer = 0;
          pa->pick_right_primer = 1;
          pa->pick_internal_oligo = 0;
        } else if (!strcmp_nocase(task_tmp, "pick_hyb_probe_only")) {
          pa->primer_task = generic;
          pa->pick_left_primer = 0;
          pa->pick_right_primer = 0;
          pa->pick_internal_oligo = 1;
        } else if (!strcmp_nocase(task_tmp, "generic")) {
          pa->primer_task = generic;
        } else if (!strcmp_nocase(task_tmp, "pick_detection_primers")) {
          pa->primer_task = generic; /* Deliberate duplication for
                        backward compatibility. */
        } else if (!strcmp_nocase(task_tmp, "pick_cloning_primers")) {
          pa->primer_task = pick_cloning_primers;
        } else if (!strcmp_nocase(task_tmp, "pick_discriminative_primers")) {
          pa->primer_task = pick_discriminative_primers;
        } else if (!strcmp_nocase(task_tmp, "pick_sequencing_primers")) {
          pa->primer_task = pick_sequencing_primers;
        } else if (!strcmp_nocase(task_tmp, "pick_primer_list")) {
          pa->primer_task = pick_primer_list;
        } else if (!strcmp_nocase(task_tmp, "check_primers")) {
          pa->primer_task = check_primers;
          /* check_primers sets the picking flags itself */
        } else {
            PyErr_Format(PyExc_ValueError, "%s is not a valid PRIMER_TASK",\
                         task_tmp);
            free(task_tmp);
            return -1;
        }
        free(task_tmp);
    }

    return 0;
}

seq_lib*
pdh_createSeqLib(PyObject *seq_dict){
    /* Generates a library of sequences for mispriming checks.
     * Input is a Python dictionary with <seq name: sequence> key value
     * pairs. Returns NULL and sets the Python error string on failure.
     */

    seq_lib                 *sl;
    PyObject                *py_seq_name, *py_seq;
    Py_ssize_t              pos;
    char                    *seq_name=NULL, *seq=NULL, *errfrag=NULL;

    if (!(sl = create_empty_seq_lib())) {
        PyErr_SetString(PyExc_IOError, "Could not allocate memory for seq_lib");
        return NULL;
    }

    pos = 0;
    while (PyDict_Next(seq_dict, &pos, &py_seq_name, &py_seq)) {
#if PY_MAJOR_VERSION < 3
            if (PyString_Check(py_seq_name)) {
                seq_name = PyString_AsString(py_seq_name);
            } else {
                PyErr_SetString(PyExc_TypeError,
                    "Cannot add seq name with non-String type to seq_lib");
                goto err_create_seq_lib;
            }
            if (PyString_Check(py_seq)) {
                seq = PyString_AsString(py_seq);
            } else {
                PyErr_SetString(PyExc_TypeError,
                    "Cannot add seq with non-String type to seq_lib");
                goto err_create_seq_lib;
            }
            if (add_seq_to_seq_lib(sl, seq, seq_name, errfrag)) {
                PyErr_SetString(PyExc_IOError, errfrag);
                goto err_create_seq_lib;
            }
#else
            if (PyUnicode_Check(py_seq_name)) {
                seq_name = (char *)PyUnicode_AsUTF8(py_seq_name);
            } else if (PyBytes_Check(py_seq_name)){
                seq_name = PyBytes_AsString(py_seq_name);
            } else {
                PyErr_SetString(PyExc_TypeError,
                    "Cannot add seq name with non-Unicode/Bytes type to seq_lib");
                goto err_create_seq_lib;
            }
            if (PyUnicode_Check(py_seq)) {
                seq = (char *)PyUnicode_AsUTF8(py_seq);
            } else if (PyBytes_Check(py_seq)){
                seq = PyBytes_AsString(py_seq);
            } else {
                PyErr_SetString(PyExc_TypeError,
                    "Cannot add seq with non-Unicode/Bytes type to seq_lib");
                goto err_create_seq_lib;
            }
            if (add_seq_to_seq_lib(sl, seq, seq_name, errfrag)) {
                PyErr_SetString(PyExc_IOError, errfrag);
                goto err_create_seq_lib;
            }
#endif
    }
    reverse_complement_seq_lib(sl);
    return sl;
err_create_seq_lib:
    destroy_seq_lib(sl);
    return NULL;
}


int
pdh_setSeqArgs(PyObject *sa_dict, seq_args *sa) {
    /* Creates a sequence args object that defines a DNA/RNA sequence for
     * which you want to design primers / oligos. Returns NULL and sets the
     * Python error string on failure.
     */

    PyObject                *p_obj, *p_obj2, *p_obj3;
    char                    *temp_char=NULL;
    int                     i, j, len1, len2;
    Py_ssize_t              str_size;
    int overlap_junction_arr_len = 0;

    // only allocate the seq args once since seq_args are a global parameter
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_TEMPLATE", &sa->sequence, temp_char, str_size);
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_ID", &sa->sequence_name, temp_char, str_size);
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_PRIMER", &sa->left_input, temp_char, str_size);
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_PRIMER_REVCOMP", &sa->right_input, temp_char, str_size);
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_INTERNAL_OLIGO", &sa->internal_input, temp_char, str_size);
    DICT_GET_AND_COPY_ARRAY(p_obj, sa_dict, "SEQUENCE_QUALITY", &sa->quality, &sa->n_quality);
    if (DICT_GET_OBJ(p_obj, sa_dict, "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST")){
        int ii[4], flat_list = 0;
        if (!PySequence_Check(p_obj)){
            PyErr_SetString(PyExc_IOError, "Value of 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST' "
                            "must support the seqeunce protocol.");
            return -1;
        }
        sa->ok_regions.count = 0;
        sa->ok_regions.any_pair = 0;
        sa->ok_regions.any_left = sa->ok_regions.any_right = 0;
        len1 = (int)PySequence_Size(p_obj);
        if (len1 == 4){
            p_obj2 = PySequence_GetItem(p_obj, 0);
            if (PyLong_Check(p_obj2)) {
                ii[0] = (int)PyLong_AsLong(p_obj2);
                Py_DECREF(p_obj2);
                for (j = 1; j < 4; j++) {
                    p_obj2 = PySequence_GetItem(p_obj, j);
                    if (PyLong_Check(p_obj2)) {
                        ii[j] = (int)PyLong_AsLong(p_obj2);
                    } else {
                        PyErr_SetString(PyExc_IOError, "'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST' must"
                                                       " be a sequence object of four integers or"
                                                       " must be comprised of sequence objects"
                                                       " comprised of four integers.");
                        Py_DECREF(p_obj2);
                        return -1;
                    }
                }
                p3_add_to_2_interval_array(&sa->ok_regions, ii[0], ii[1], ii[2], ii[3]);
                flat_list = 1;
            }
        }
        if (!flat_list) {
            for (i = 0; i < len1; i++) {
                p_obj2 = PySequence_GetItem(p_obj, i);
                if (!PySequence_Check(p_obj2)) {
                    PyErr_SetString(PyExc_IOError, "'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST' must support "
                                                      "the sequence protocol and be comprised of items "
                                                      "that support the sequence protocol (e.g., it must be "
                                                      "a list of lists, tuple of tuples or some combination "
                                                      "of the two).");
                    return -1;
                }
                len2 = (int)PySequence_Size(p_obj2);
                if (!(len2 == 4)) {
                    PyErr_Format(PyExc_TypeError, "Sub-list/tuple #%d of "
                                 "'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST' must "
                                 "of length 4", i);
                    Py_DECREF(p_obj2);
                    return -1;
                }
                for (j = 0; j < 4; j++) {
                    p_obj3 = PySequence_GetItem(p_obj2, j);
                    if (!PyLong_Check(p_obj3)) {
                        PyErr_Format(PyExc_TypeError, "Object #%d of sub-list/tuple %d"
                             "'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST' must "
                             "be an integer or long", j, i);
                        Py_DECREF(p_obj3);
                        return -1;
                    }
                    ii[j] = (int)PyLong_AsLong(p_obj3);
                    Py_DECREF(p_obj3);
                }
                p3_add_to_2_interval_array(&sa->ok_regions, ii[0], ii[1], ii[2], ii[3]);
            }
        }
    }
    DICT_GET_AND_COPY_TO_INTERVAL_ARRAY(p_obj, sa_dict, "SEQUENCE_TARGET", sa->tar2);
    DICT_GET_AND_COPY_TO_INTERVAL_ARRAY(p_obj, sa_dict, "SEQUENCE_EXCLUDED_REGION", sa->excl2);
    DICT_GET_AND_COPY_TO_INTERVAL_ARRAY(p_obj, sa_dict, "SEQUENCE_INTERNAL_EXCLUDED_REGION", sa->excl_internal2);
    // DICT_GET_AND_COPY_ARRAY_INTO_ARRAY(p_obj, sa_dict,  "SEQUENCE_OVERLAP_JUNCTION_LIST", &sa->primer_overlap_junctions, overlap_junction_len);

    if (DICT_GET_OBJ(p_obj, sa_dict, "SEQUENCE_OVERLAP_JUNCTION_LIST")) {
        int *poj_arr, single_value=0;
        PyObject *arr_item;
        if (!PySequence_Check(p_obj)){
            if (PyLong_Check(p_obj)) {
                sa->primer_overlap_junctions[0] = (int)PyLong_AsLong(p_obj);
                sa->primer_overlap_junctions_count++;
                single_value = 1;
            } else {
            PyErr_Format(PyExc_TypeError,
                            "Value of 'SEQUENCE_OVERLAP_JUNCTION_LIST' is not a sequence object");
            return -1;}
        }
        if (!single_value) {
            overlap_junction_arr_len = (int)PySequence_Size(p_obj);
            if (overlap_junction_arr_len > 200) {
                PyErr_Format(PyExc_TypeError,
                                "'SEQUENCE_OVERLAP_JUNCTION_LIST' cannot have over 200 values");
                return -1;
            }
            sa->primer_overlap_junctions_count = overlap_junction_arr_len;
            poj_arr = &sa->primer_overlap_junctions[0];
            for (i=0; i < overlap_junction_arr_len; i++, poj_arr++) {
                arr_item = PySequence_GetItem(p_obj, i);
                if (!PyLong_Check(arr_item)) {
                    PyErr_Format(PyExc_TypeError,
                                "'SEQUENCE_OVERLAP_JUNCTION_LIST' must contain only integers");
                }
                *poj_arr = (int)PyLong_AsLong(arr_item);
                Py_DECREF(arr_item);
            }
        }
    }
    if DICT_GET_OBJ(p_obj, sa_dict, "SEQUENCE_INCLUDED_REGION") {
        PyObject *seq_item1, *seq_item2;
        if (!PySequence_Check(p_obj)) {
            PyErr_SetString(PyExc_TypeError,\
                "Value of \"SEQUENCE_INCLUDED_REGION\" is not a sequence object");
            return -1;
        } else if (PySequence_Size(p_obj) != 2) {
            PyErr_SetString(PyExc_ValueError,\
                "Length of \"SEQUENCE_INCLUDED_REGION\" is not of length 2");
            return -1;
        } else {
            seq_item1 = PySequence_GetItem(p_obj, 0);
            seq_item2 = PySequence_GetItem(p_obj, 1);
            if (!PyLong_Check(seq_item1) || !PyLong_Check(seq_item2)) {
                PyErr_SetString(PyExc_TypeError,\
                    "\"SEQUENCE_INCLUDED_REGION\" contains non-int value");
                Py_DECREF(seq_item1);
                Py_DECREF(seq_item2);
                return -1;
            }
            sa->incl_s = (int)PyLong_AsLong(seq_item1);
            sa->incl_l = (int)PyLong_AsLong(seq_item2);
            Py_DECREF(seq_item1);
            Py_DECREF(seq_item2);
        }

    }

    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_START_CODON_POSITION", sa->start_codon_pos);
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_FORCE_LEFT_START", sa->force_left_start);
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_FORCE_LEFT_END", sa->force_left_end);
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_FORCE_RIGHT_START", sa->force_right_start);
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_FORCE_RIGHT_END", sa->force_right_end);

    return 0;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* Code for dumping all primer3 output to a python dictionary. This is
 * a modified version of print_boulder and its helper methods from
 * print_boulder.c/h
 */

#define SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)                                \
    if (PyDict_SetItemString(dict, key, obj_ptr) == -1) {                      \
        PyErr_Format(PyExc_IOError, "Could not set dictionary value for %s",   \
                     key);                                                     \
        Py_DECREF(obj_ptr);                                                    \
        return NULL;                                                           \
    }                                                                          \
    Py_DECREF(obj_ptr);                                                        \


#define SET_DICT_KEY_TO_LONG(dict, key, num, obj_ptr)                          \
    if ((obj_ptr = PyLong_FromLong(num)) == NULL){                             \
        PyErr_Format(PyExc_IOError, "Could not convert value for %s to PyLong",\
                     key);                                                     \
        return NULL;                                                           \
    }                                                                          \
    SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)                                    \


#define SET_DICT_KEY_TO_DOUBLE(dict, key, num, obj_ptr)                        \
    if ((obj_ptr = PyFloat_FromDouble(num)) == NULL){                          \
        PyErr_Format(PyExc_IOError, "Could not convert value for %s to PyFloat",\
                     key);                                                     \
        return NULL;                                                           \
    }                                                                          \
    SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)                                    \


#if PY_MAJOR_VERSION < 3
    #define SET_DICT_KEY_TO_STR(dict, key, str, obj_ptr)                       \
        if ((obj_ptr = PyString_FromString(str)) == NULL){                     \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not convert value for %s to PyString", key);   \
            return NULL;                                                       \
        }                                                                      \
        SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)
#else
    #define SET_DICT_KEY_TO_STR(dict, key, str, obj_ptr)                       \
        if ((obj_ptr = PyUnicode_FromString(str)) == NULL){                    \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not convert value for %s to PyUnicode", key);  \
            return NULL;                                                       \
        }                                                                      \
        SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)
#endif

// Note that PyTuple_SetItem steals a ref so the PyLong objs don't need to be
// Py_DECREF'd
#define SET_DICT_KEY_TO_TUPLE_OF_LONGS(dict, key, num1, num2, obj_ptr)         \
    if ((obj_ptr = PyTuple_New(2)) == NULL) {                                  \
        PyErr_Format(PyExc_IOError, "Could not create tuple for %s ", key);    \
        return NULL;                                                           \
    }                                                                          \
    if (PyTuple_SetItem(obj_ptr, 0, PyLong_FromLong(num1))) {                  \
        PyErr_Format(PyExc_IOError, "Could not pack value 1 for %s into tuple",\
                     key);                                                     \
        return NULL;                                                           \
    }                                                                          \
    if (PyTuple_SetItem(obj_ptr, 1, PyLong_FromLong(num2))) {                  \
        PyErr_Format(PyExc_IOError, "Could not pack value 2 for %s into tuple",\
                     key);                                                     \
        return NULL;                                                           \
    }                                                                          \
    SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)

#if PY_MAJOR_VERSION < 3
    #define SET_DICT_KEY_TO_TUPLE_OF_LONG_AND_STR(dict, key, num, str, obj_ptr)\
        if ((obj_ptr = PyTuple_New(2)) == NULL) {                              \
            PyErr_Format(PyExc_IOError, "Could not create tuple for %s ", key);\
            return NULL;                                                       \
        }                                                                      \
        if (PyTuple_SetItem(obj_ptr, 0, PyLong_FromLong(num))) {               \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not pack value 1 for %s into tuple", key);     \
            return NULL;                                                       \
        }                                                                      \
        if (PyTuple_SetItem(obj_ptr, 1, PyString_FromString(str))) {           \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not pack value 2 for %s into tuple",key);      \
            return NULL;                                                       \
        }                                                                      \
        SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)
#else
    #define SET_DICT_KEY_TO_TUPLE_OF_LONG_AND_STR(dict, key, num, str, obj_ptr)\
        if ((obj_ptr = PyTuple_New(2)) == NULL) {                              \
            PyErr_Format(PyExc_IOError, "Could not create tuple for %s ", key);\
            return NULL;                                                       \
        }                                                                      \
        if (PyTuple_SetItem(obj_ptr, 0, PyLong_FromLong(num))) {               \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not pack value 1 for %s into tuple", key);     \
            return NULL;                                                       \
        }                                                                      \
        if (PyTuple_SetItem(obj_ptr, 1, PyUnicode_FromString(str))) {          \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not pack value 2 for %s into tuple",key);      \
            return NULL;                                                       \
        }                                                                      \
        SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)
#endif

#if PY_MAJOR_VERSION < 3
    #define SET_DICT_KEY_TO_TUPLE_OF_FLOAT_AND_STR(dict, key, num, str, obj_ptr)\
        if ((obj_ptr = PyTuple_New(2)) == NULL) {                              \
            PyErr_Format(PyExc_IOError, "Could not create tuple for %s ", key);\
            return NULL;                                                       \
        }                                                                      \
        if (PyTuple_SetItem(obj_ptr, 0, PyFloat_FromDouble(num))) {            \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not pack value 1 for %s into tuple", key);     \
            return NULL;                                                       \
        }                                                                      \
        if (PyTuple_SetItem(obj_ptr, 1, PyString_FromString(str))) {           \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not pack value 2 for %s into tuple",key);      \
            return NULL;                                                       \
        }                                                                      \
        SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)
#else
    #define SET_DICT_KEY_TO_TUPLE_OF_FLOAT_AND_STR(dict, key, num, str, obj_ptr)\
        if ((obj_ptr = PyTuple_New(2)) == NULL) {                              \
            PyErr_Format(PyExc_IOError, "Could not create tuple for %s ", key);\
            return NULL;                                                       \
        }                                                                      \
        if (PyTuple_SetItem(obj_ptr, 0, PyFloat_FromDouble(num))) {            \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not pack value 1 for %s into tuple", key);     \
            return NULL;                                                       \
        }                                                                      \
        if (PyTuple_SetItem(obj_ptr, 1, PyUnicode_FromString(str))) {          \
            PyErr_Format(PyExc_IOError,                                        \
                         "Could not pack value 2 for %s into tuple",key);      \
            return NULL;                                                       \
        }                                                                      \
        SET_DICT_KEY_TO_OBJ(dict, key, obj_ptr)
#endif

PyObject*
pdh_outputToDict(const p3_global_settings *pa, const seq_args *sa,
               const p3retval *retval) {
    PyObject *output_dict;
    PyObject *obj_ptr = NULL;

    /* The pointers to warning tag */
    char *warning;

    char outbuff[100]; // Output buffer for sprintf key string formatting

    /* A place to put a string containing all error messages */
    pr_append_str *combined_retval_err = NULL;

    /* A small spacer; WARNING this is a fixed size
     buffer, but plenty bigger than
     log(2^64, 10), the longest character
     string that is needed for a 64 bit integer. */
    char suffix [100];

    /* Pointers for the primer set just printing */
    primer_rec *fwd, *rev, *intl;

    /* Variables only used for Primer Lists */
    int num_fwd, num_rev, num_int, num_pair, num_print;
    int print_fwd = 0;
    int print_rev = 0;
    int print_int = 0;

    /* Switches for printing this primer */
    int go_fwd = 0;
    int go_rev = 0;
    int go_int = 0;

    /* The number of loop cycles */
    int loop_max;

    /* That links to the included region */
    int i, incl_s = sa->incl_s;

    /* This deals with the renaming of the internal oligo */
    const char *new_oligo_name = "INTERNAL";
    char *int_oligo = (char*) new_oligo_name;

    output_dict = PyDict_New();
    if (output_dict == NULL) {
        PyErr_SetString(PyExc_IOError, "Could not create Primer3 output dict");
        return NULL;
    }

    /* Check if there are warnings and print them */
    if ((warning = p3_get_rv_and_gs_warnings(retval, pa)) != NULL) {
        SET_DICT_KEY_TO_STR(output_dict, "PRIMER_WARNING", warning, obj_ptr);
        // printf("PRIMER_WARNING=%s\n", warning);
        free(warning);
    }

    combined_retval_err = create_pr_append_str();
    if (NULL == combined_retval_err) {
        PyErr_SetString(PyExc_IOError, "Primer3 ran out of memory.");
        return NULL;
    }

    if (pr_append_new_chunk_external(combined_retval_err,
                                   retval->glob_err.data)){
        PyErr_SetString(PyExc_IOError, "Primer3 ran out of memory.");
        return NULL;
    }

    if (pr_append_new_chunk_external(combined_retval_err,
                                   retval->per_sequence_err.data)){
        PyErr_SetString(PyExc_IOError, "Primer3 ran out of memory.");
        return NULL;
    }


    /* Check if there are errors, print and return */
    if (!pr_is_empty(combined_retval_err)) {
        PyErr_SetString(PyExc_IOError, \
            pr_append_str_chars(combined_retval_err));
        destroy_pr_append_str(combined_retval_err);
        return NULL;
    }
    destroy_pr_append_str(combined_retval_err);

    /* Get how many primers are in the array */
    num_fwd = retval->fwd.num_elem;
    num_rev = retval->rev.num_elem;
    num_int = retval->intl.num_elem;
    num_pair = retval->best_pairs.num_pairs;

    /* Prints out selection statistics about the primers */
    if (pa->pick_left_primer == 1
      && !(pa->pick_anyway && sa->left_input)) {
        SET_DICT_KEY_TO_STR(output_dict, "PRIMER_LEFT_EXPLAIN", \
            p3_get_oligo_array_explain_string(p3_get_rv_fwd(retval)), obj_ptr);
    }
    if (pa->pick_right_primer == 1
      && !(pa->pick_anyway && sa->right_input)) {
        SET_DICT_KEY_TO_STR(output_dict, "PRIMER_RIGHT_EXPLAIN", \
            p3_get_oligo_array_explain_string(p3_get_rv_rev(retval)), obj_ptr);
    }

    if ( pa->pick_internal_oligo == 1
      && !(pa->pick_anyway && sa->internal_input)){
        SET_DICT_KEY_TO_STR(output_dict, "PRIMER_INTERNAL_EXPLAIN", \
            p3_get_oligo_array_explain_string(p3_get_rv_intl(retval)), obj_ptr);
    }
    if (pa->pick_right_primer == 1
      && pa->pick_left_primer == 1) {
        SET_DICT_KEY_TO_STR(output_dict, "PRIMER_PAIR_EXPLAIN", \
            p3_get_pair_array_explain_string(p3_get_rv_best_pairs(retval)), \
            obj_ptr);
    }

    /* Print out the stop codon if a reading frame was specified */
    if (!PR_START_CODON_POS_IS_NULL(sa)) {
         sprintf(outbuff, "PRIMER_STOP_CODON_POSITION=%d", \
                 retval->stop_codon_pos);
         SET_DICT_KEY_TO_LONG(output_dict, outbuff, retval->stop_codon_pos, \
                              obj_ptr);
    }

    /* How often has the loop to be done? */
    if (retval->output_type == primer_list) {
        /* For Primer Lists: Figure out how many primers are in
         * the array that can be printed. If more than needed,
         * set it to the number requested. */

        /* Get how may primers should be printed */
        num_print = pa->num_return;
        /* Set how many primers will be printed */
        print_fwd = (num_print < num_fwd) ? num_print : num_fwd;
        print_rev = (num_print < num_rev) ? num_print : num_rev;
        print_int = (num_print < num_int) ? num_print : num_int;
        /* Get which list has to print most primers */
        loop_max = 0;
        if (loop_max < print_fwd) {
            loop_max = print_fwd;
        }
        if (loop_max < print_rev) {
            loop_max = print_rev;
        }
        if (loop_max < print_int) {
            loop_max = print_int;
        }
        /* Now the vars are there how often we have to go
         * through the loop and how many of each primer can
         * be printed. */
        num_pair = 0;
    } else {
        loop_max = num_pair;
        /* Set how many primers will be printed */
        print_fwd = num_pair;
        print_rev = num_pair;
        if (num_int != 0) {
            print_int = num_pair;
        }
    }

    // Save the number of each type of oligo that was found
    SET_DICT_KEY_TO_LONG(output_dict, "PRIMER_LEFT_NUM_RETURNED", print_fwd, obj_ptr);
    SET_DICT_KEY_TO_LONG(output_dict, "PRIMER_RIGHT_NUM_RETURNED", print_rev, obj_ptr);
    sprintf(outbuff, "PRIMER_%s_NUM_RETURNED", int_oligo);
    SET_DICT_KEY_TO_LONG(output_dict, outbuff, print_int, obj_ptr);
    SET_DICT_KEY_TO_LONG(output_dict, "PRIMER_PAIR_NUM_RETURNED", num_pair, obj_ptr);

    /* --------------------------------------- */
    /* Start of the loop printing all pairs or primers or oligos */
    for(i=0; i<loop_max; i++) {
    /* What needs to be printed */
    /* The conditions for primer lists */

    if (retval->output_type == primer_list) {
      /* Attach the selected primers to the pointers */
      fwd = &retval->fwd.oligo[i];
      rev = &retval->rev.oligo[i];
      intl = &retval->intl.oligo[i];
      /* Do fwd oligos have to be printed? */
      if ((pa->pick_left_primer) && (i < print_fwd)) {
        go_fwd = 1;
      } else {
        go_fwd = 0;
      }
      /* Do rev oligos have to be printed? */
      if ((pa->pick_right_primer) && (i < print_rev)) {
        go_rev = 1;
      } else {
        go_rev = 0;
      }
      /* Do int oligos have to be printed? */
      if ((pa->pick_internal_oligo) && (i < print_int)) {
        go_int = 1;
      } else {
        go_int = 0;
      }
    }  else {
      /* We will print primer pairs or pairs plus internal oligos */
      /* Get pointers to the primer_rec's that we will print */
      fwd  = retval->best_pairs.pairs[i].left;
      rev  = retval->best_pairs.pairs[i].right;
      intl = retval->best_pairs.pairs[i].intl;
      /* Pairs must have fwd and rev primers */
      go_fwd = 1;
      go_rev = 1;
      /* Do hyb oligos have to be printed? */
      if (pa->pick_internal_oligo == 1) {
        go_int = 1;
      } else {
        go_int = 0;
      }
    }

    /* Get the number for pimer counting in suffix[0] */
    sprintf(suffix, "_%d", i);

    /* Print out the Pair Penalties */
    if (retval->output_type == primer_pairs) {
      sprintf(outbuff, "PRIMER_PAIR%s_PENALTY", suffix);
      SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
        retval->best_pairs.pairs[i].pair_quality, obj_ptr);
    }

    /* Print single primer penalty */
    if (go_fwd == 1) {
        sprintf(outbuff, "PRIMER_LEFT%s_PENALTY", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->quality, obj_ptr);
    }
    if (go_rev == 1) {
        sprintf(outbuff, "PRIMER_RIGHT%s_PENALTY", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->quality, obj_ptr);
    }
    if (go_int == 1) {
        sprintf(outbuff, "PRIMER_%s%s_PENALTY", int_oligo, suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, intl->quality, obj_ptr);
    }

    /* Print the oligo_problems */

    if (go_fwd == 1 && p3_ol_has_any_problem(fwd)) {
        sprintf(outbuff, "PRIMER_LEFT%s_PROBLEMS", suffix);
        SET_DICT_KEY_TO_STR(output_dict, outbuff, \
            p3_get_ol_problem_string(fwd), obj_ptr);
    }
    if (go_rev == 1 && p3_ol_has_any_problem(rev)) {
        sprintf(outbuff, "PRIMER_RIGHT%s_PROBLEMS", suffix);
        SET_DICT_KEY_TO_STR(output_dict, outbuff, \
            p3_get_ol_problem_string(rev), obj_ptr);
    }
    if (go_int == 1 && p3_ol_has_any_problem(intl)) {
        sprintf(outbuff, "PRIMER_%s%s_PROBLEMS", int_oligo, suffix);
        SET_DICT_KEY_TO_STR(output_dict, outbuff, \
            p3_get_ol_problem_string(intl), obj_ptr);
    }

    /* Print primer sequences. */

    if (go_fwd == 1) {
        sprintf(outbuff, "PRIMER_LEFT%s_SEQUENCE", suffix);
        SET_DICT_KEY_TO_STR(output_dict, outbuff, \
            pr_oligo_sequence(sa, fwd), obj_ptr);
    }
    if (go_rev == 1) {
        sprintf(outbuff, "PRIMER_RIGHT%s_SEQUENCE", suffix);
        SET_DICT_KEY_TO_STR(output_dict, outbuff, \
            pr_oligo_rev_c_sequence(sa, rev), obj_ptr);
    }
    if (go_int == 1) {
        sprintf(outbuff, "PRIMER_%s%s_SEQUENCE", int_oligo, suffix);
        SET_DICT_KEY_TO_STR(output_dict, outbuff, \
            pr_oligo_sequence(sa, intl), obj_ptr);
    }

    /* Print primer start and length */
    if (go_fwd == 1) {
      sprintf(outbuff, "PRIMER_LEFT%s", suffix);
      SET_DICT_KEY_TO_TUPLE_OF_LONGS(output_dict, outbuff, \
        fwd->start + incl_s + pa->first_base_index, fwd->length, obj_ptr);
    }
    if (go_rev == 1) {
      sprintf(outbuff, "PRIMER_RIGHT%s", suffix);
      SET_DICT_KEY_TO_TUPLE_OF_LONGS(output_dict, outbuff, \
        rev->start + incl_s + pa->first_base_index, rev->length, obj_ptr);
    }
    if (go_int == 1) {
      sprintf(outbuff, "PRIMER_%s%s", int_oligo, suffix);
      SET_DICT_KEY_TO_TUPLE_OF_LONGS(output_dict, outbuff, \
        intl->start + incl_s + pa->first_base_index, intl->length, obj_ptr);
    }

    /* Print primer Tm */
    if (go_fwd == 1) {
        sprintf(outbuff, "PRIMER_LEFT%s_TM", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->temp, obj_ptr);
    }
    if (go_rev == 1) {
        sprintf(outbuff, "PRIMER_RIGHT%s_TM", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->temp, obj_ptr);
    }
    if (go_int == 1) {
        sprintf(outbuff, "PRIMER_%s%s_TM", int_oligo, suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, intl->temp, obj_ptr);
    }

    /* Print primer GC content */
    if (go_fwd == 1) {
        sprintf(outbuff, "PRIMER_LEFT%s_GC_PERCENT", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->gc_content, obj_ptr);
    }
    if (go_rev == 1) {
        sprintf(outbuff, "PRIMER_RIGHT%s_GC_PERCENT", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->gc_content, obj_ptr);
    }
    if (go_int == 1) {
        sprintf(outbuff, "PRIMER_%s%s_GC_PERCENT", int_oligo, suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, intl->gc_content, obj_ptr);
    }

    /* Print primer self_any */
    if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==0) {
        sprintf(outbuff, "PRIMER_LEFT%s_SELF_ANY", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->self_any, obj_ptr);
    }
    if (go_rev == 1 && pa->thermodynamic_oligo_alignment==0) {
        sprintf(outbuff, "PRIMER_RIGHT%s_SELF_ANY", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->self_any, obj_ptr);
    }
    if (go_int == 1 && pa->thermodynamic_oligo_alignment==0) {
        sprintf(outbuff, "PRIMER_%s%s_SELF_ANY", int_oligo, suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, intl->self_any, obj_ptr);
    }
    /* Print primer self_any thermodynamical approach */
    if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_LEFT%s_SELF_ANY_TH", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->self_any, obj_ptr);
    }
    if (go_rev == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_RIGHT%s_SELF_ANY_TH", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->self_any, obj_ptr);
    }
    if (go_int == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_%s%s_SELF_ANY_TH", int_oligo, suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, intl->self_any, obj_ptr);
    }
    /* Print primer self_end */
    if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==0) {
        sprintf(outbuff, "PRIMER_LEFT%s_SELF_END", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->self_end, obj_ptr);
    }
    if (go_rev == 1 && pa->thermodynamic_oligo_alignment==0) {
        sprintf(outbuff, "PRIMER_RIGHT%s_SELF_END", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->self_end, obj_ptr);
    }
    if (go_int == 1 && pa->thermodynamic_oligo_alignment==0) {
        sprintf(outbuff, "PRIMER_%s%s_SELF_END", int_oligo, suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, intl->self_end, obj_ptr);
    }
    /* Print primer self_end thermodynamical approach */
    if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_LEFT%s_SELF_END_TH", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->self_end, obj_ptr);
    }
    if (go_rev == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_RIGHT%s_SELF_END_TH", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->self_end, obj_ptr);
    }
    if (go_int == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_%s%s_SELF_END_TH", int_oligo, suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, intl->self_end, obj_ptr);
    }

     /* Print primer hairpin */
    if (go_fwd == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_LEFT%s_HAIRPIN_TH", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->hairpin_th, obj_ptr);
    }
    if (go_rev == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_RIGHT%s_HAIRPIN_TH", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->hairpin_th, obj_ptr);
    }
    if (go_int == 1 && pa->thermodynamic_oligo_alignment==1) {
        sprintf(outbuff, "PRIMER_%s%s_HAIRPIN_TH", int_oligo, suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, intl->hairpin_th, obj_ptr);
    }

     /*Print out primer mispriming scores */
    if (seq_lib_num_seq(pa->p_args.repeat_lib) > 0) {
        if (go_fwd == 1){
            sprintf(outbuff, "PRIMER_LEFT%s_LIBRARY_MISPRIMING", suffix);
            SET_DICT_KEY_TO_TUPLE_OF_FLOAT_AND_STR(output_dict, outbuff, \
                fwd->repeat_sim.score[fwd->repeat_sim.max], \
                fwd->repeat_sim.name, obj_ptr);
        }
        if (go_rev == 1){
            sprintf(outbuff, "PRIMER_RIGHT%s_LIBRARY_MISPRIMING", suffix);
            SET_DICT_KEY_TO_TUPLE_OF_FLOAT_AND_STR(output_dict, outbuff, \
                rev->repeat_sim.score[rev->repeat_sim.max], \
                rev->repeat_sim.name, obj_ptr);
        }
        if (retval->output_type == primer_pairs){
            sprintf(outbuff, "PRIMER_PAIR%s_LIBRARY_MISPRIMING", suffix);
            SET_DICT_KEY_TO_TUPLE_OF_FLOAT_AND_STR(output_dict, outbuff, \
                retval->best_pairs.pairs[i].repeat_sim, \
                retval->best_pairs.pairs[i].rep_name, obj_ptr);
        }
    }
    /* Print out internal oligo mispriming scores */
    if (go_int == 1 && seq_lib_num_seq(pa->o_args.repeat_lib) > 0) {
        sprintf(outbuff, "PRIMER_%s%s_LIBRARY_MISPRIMING", int_oligo, suffix);
        SET_DICT_KEY_TO_TUPLE_OF_FLOAT_AND_STR(output_dict, outbuff, \
            intl->repeat_sim.score[intl->repeat_sim.max], \
            intl->repeat_sim.name, obj_ptr);
    }

    /* If a sequence quality was provided, print it*/
    if (NULL != sa->quality){
        if (go_fwd == 1) {
            sprintf(outbuff, "PRIMER_LEFT%s_MIN_SEQ_QUALITY", suffix);
            SET_DICT_KEY_TO_LONG(output_dict, outbuff, fwd->seq_quality, \
                                 obj_ptr);
        }
        if (go_rev == 1) {
            sprintf(outbuff, "PRIMER_RIGHT%s_MIN_SEQ_QUALITY", suffix);
            SET_DICT_KEY_TO_LONG(output_dict, outbuff, rev->seq_quality, \
                                 obj_ptr);
        }
        if (go_int == 1) {
            sprintf(outbuff, "PRIMER_%s%s_MIN_SEQ_QUALITY", int_oligo, suffix);
            SET_DICT_KEY_TO_LONG(output_dict, outbuff, intl->seq_quality, \
                                 obj_ptr);
        }
      /* Has to be here and in primer pairs for backward compatibility */
    }

    /* Print position penalty, this is for backward compatibility */
    if (!_PR_DEFAULT_POSITION_PENALTIES(pa) || !PR_START_CODON_POS_IS_NULL(sa)){
        sprintf(outbuff, "PRIMER_LEFT%s_POSITION_PENALTY", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->position_penalty, \
                               obj_ptr);
        sprintf(outbuff, "PRIMER_RIGHT%s_POSITION_PENALTY", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->position_penalty, \
                               obj_ptr);
    }

    /* Print primer end stability */
    if (go_fwd == 1) {
        sprintf(outbuff, "PRIMER_LEFT%s_END_STABILITY", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, fwd->end_stability, \
                               obj_ptr);
    }

    if (go_rev == 1){
        sprintf(outbuff, "PRIMER_RIGHT%s_END_STABILITY", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, rev->end_stability, \
                               obj_ptr);
    }

    /* Print primer template mispriming */
    if ((pa->thermodynamic_template_alignment == 0) && (go_fwd == 1) &&
         (oligo_max_template_mispriming(fwd) != ALIGN_SCORE_UNDEF)) {
        sprintf(outbuff, "PRIMER_LEFT%s_TEMPLATE_MISPRIMING", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff,
            oligo_max_template_mispriming(fwd), obj_ptr);
    }
    if ( (pa->thermodynamic_template_alignment == 0) && (go_rev == 1) &&
         (oligo_max_template_mispriming(rev) != ALIGN_SCORE_UNDEF)) {
        sprintf(outbuff, "PRIMER_RIGHT%s_TEMPLATE_MISPRIMING", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
            oligo_max_template_mispriming(rev), obj_ptr);
    }

     /* Print primer template mispriming, thermodynamical approach*/
    if ((pa->thermodynamic_template_alignment == 0) && (go_fwd == 1) &&
         (oligo_max_template_mispriming(fwd) != ALIGN_SCORE_UNDEF)) {
        sprintf(outbuff, "PRIMER_LEFT%s_TEMPLATE_MISPRIMING_TH", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
            oligo_max_template_mispriming_thermod(fwd), obj_ptr);
    }
    if ( (pa->thermodynamic_template_alignment == 0) && (go_rev == 1) &&
         (oligo_max_template_mispriming(rev) != ALIGN_SCORE_UNDEF)) {
        sprintf(outbuff, "PRIMER_RIGHT%s_TEMPLATE_MISPRIMING_TH", suffix);
        SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
            oligo_max_template_mispriming_thermod(rev), obj_ptr);

    }

     /* Print the pair parameters*/
    if (retval->output_type == primer_pairs) {
        if (go_int == 1 && NULL != sa->quality) {
            sprintf(outbuff, "PRIMER_%s%s_MIN_SEQ_QUALITY", int_oligo, suffix);
            SET_DICT_KEY_TO_LONG(output_dict, outbuff, intl->seq_quality, \
                                 obj_ptr);
        }
        /* Print pair comp_any */
        if (pa->thermodynamic_oligo_alignment==0){
            sprintf(outbuff, "PRIMER_PAIR%s_COMPL_ANY", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].compl_any, obj_ptr);
        }
        if (pa->thermodynamic_oligo_alignment==1) {
            sprintf(outbuff, "PRIMER_PAIR%s_COMPL_ANY_TH", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].compl_any, obj_ptr);
        }
        /* Print pair comp_end */
        if (pa->thermodynamic_oligo_alignment==0) {
            sprintf(outbuff, "PRIMER_PAIR%s_COMPL_END", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].compl_end, obj_ptr);
        }
        if (pa->thermodynamic_oligo_alignment==1) {
            sprintf(outbuff, "PRIMER_PAIR%s_COMPL_END_TH", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].compl_end, obj_ptr);
        }
        /* Print product size */
        sprintf(outbuff, "PRIMER_PAIR%s_PRODUCT_SIZE", suffix);
        SET_DICT_KEY_TO_LONG(output_dict, outbuff, \
            retval->best_pairs.pairs[i].product_size, obj_ptr);
        /* Print the product Tm if a Tm range is defined */
        if (pa->product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM ||
            pa->product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM) {
            sprintf(outbuff, "PRIMER_PAIR%s_PRODUCT_TM", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].product_tm, obj_ptr);

            sprintf(outbuff, "PRIMER_PAIR%s_PRODUCT_TM_OLIGO_TM_DIFF", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].product_tm_oligo_tm_diff, obj_ptr);

            sprintf(outbuff, "PRIMER_PAIR%s_T_OPT_A=", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].t_opt_a, obj_ptr);
        }

        /* Print the primer pair template mispriming */
        if ((pa->thermodynamic_template_alignment == 0) &&
            (retval->best_pairs.pairs[i].template_mispriming !=\
             ALIGN_SCORE_UNDEF)) {
            sprintf(outbuff, "PRIMER_PAIR%s_TEMPLATE_MISPRIMING", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].template_mispriming, obj_ptr);
        }
       /* Print the primer pair template mispriming. Thermodynamic approach.  */
       if ((pa->thermodynamic_template_alignment == 1) &&
          (retval->best_pairs.pairs[i].template_mispriming != \
           ALIGN_SCORE_UNDEF)) {
            sprintf(outbuff, "PRIMER_PAIR%s_TEMPLATE_MISPRIMING_TH", suffix);
            SET_DICT_KEY_TO_DOUBLE(output_dict, outbuff, \
                retval->best_pairs.pairs[i].template_mispriming, obj_ptr);
        }
    } /* End of print parameters of primer pairs */

    } /* End of the big loop printing all data */

    return output_dict;
}
