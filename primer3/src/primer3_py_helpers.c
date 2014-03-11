/*

primer3_py_helpers.c
~~~~~~~~~~~~~~~~~~~~

This file defines macros and helper functions that facilitate interaction 
between Python C API code and primer3 native C code. 

*/

#include    <string.h>
#include    <stdio.h>
#include    <Python.h>
#include    <libprimer3_mod.h>


#if PY_MAJOR_VERSION < 3
/* see http://python3porting.com/cextensions.html */
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
            PyErr_Format(PyExc_TypeError,                                   \
                            "Value of %s is not an integer.", k);              \
            return NULL;                                                       \
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
            return NULL;                                                       \
        }                                                                      \
        st = (t)PyLong_AsLong(o);                                              \
    }

// Wraps DICT_GET_OBJ and takes the dictionary value object, extracts a double,
// and assigns its value to `st`
#define DICT_GET_AND_ASSIGN_DOUBLE(o, d, k, st)                                \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        if (!PyFloat_Check(o) || !PyLong_Check(o)) {                           \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not of type float or integer.", k);\
            return NULL;                                                       \
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
    #define DICT_GET_AND_COPY_STR(o, d, k, st)                                 \
        if (DICT_GET_OBJ(o, d, k)) {                                           \
            if (!PyString_Check(o)){                                           \
                PyErr_Format(PyExc_TypeError,                                  \
                            sprintf("Value of %s is not of type string", k));  \
                return NULL;}                                                  \
            *st = (char *) malloc(PyBytes_Size(o));                            \
            if (st == NULL) {                                                  \
                strcpy(err, "Primer3 out of memory");                          \
                return NULL;}                                                  \
            strcpy(*st, PyString_AsString(o));                                 \
        }
#else
    #define DICT_GET_AND_COPY_STR(o, d, k, st)                                 \
        if (DICT_GET_OBJ(o, d, k)) {                                           \
            if (PyUnicode_Check(o)) {                                          \
                o = PyUnicode_AsASCIIString(o);                                \
            } else if (!PyBytes_Check(o)){                                     \
                PyErr_Format(PyExc_TypeError,                                  \
                            "Value of %s is not of type unicode or bytes", k); \
                return NULL;}                                                  \
            *st = (char *) malloc(PyBytes_Size(o));                            \
            if (st == NULL) {                                                  \
                PyErr_Format(PyExc_IOError,                                    \
                            "Could not allocate memory while copying %s", k);  \
                return NULL;}                                                  \
            strcpy(*st, PyBytes_AsString(o));                                  \
        }        
#endif

#define DICT_GET_AND_COPY_ARRAY(o, d, k, st, arr_len)                          \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        if (!PyList_Check(o)){                                                 \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not of type list", k);             \
            return NULL;}                                                      \
        *arr_len = PyList_Size(o);                                             \
        int arr[*arr_len];                                                     \
        int i;                                                                 \
        for (i=0; i < *arr_len; i++) {                                         \
            arr[i] = (int)PyLong_AsLong(PyList_GetItem(o, i));                 \
        }                                                                      \
        *st = arr;                                                             \
    }                                                                          

#define DICT_GET_AND_COPY_ARRAY_INTO_ARRAY(o, d, k, st, arr_len)               \
    if (DICT_GET_OBJ(o, d, k)) {                                               \
        if (!PyList_Check(o)){                                                 \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not of type list", k);             \
            return NULL;}                                                      \
        *arr_len = PyList_Size(o);                                             \
        int i;                                                                 \
        for (i=0; i < *arr_len; i++) {                                         \
            *st[i] = (int)PyLong_AsLong(PyList_GetItem(o, i));                 \
        }                                                                      \
    }    

#define DICT_GET_AND_COPY_TO_2_INTERVAL_ARRAY(o, d, k, st)                     \
    if (DICT_GET_OBJ(o, d, k)){                                                \
        if (!PyList_Check(o)){                                                 \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not of type list", k);             \
            return NULL;}                                                      \
        st.count = 0;                                                          \
        st.any_pair = 0;                                                       \
        st.any_left = st.any_right = 0;                                        \
        *arr_len = (int)PyList_Size(p_obj);                                    \
        if (!arr_len % 4) {                                                    \
            PyErr_Format(PyExc_TypeError,                                      \
                            "%s must be linear multiple of 4 in length", k);   \
            return NULL;                                                       \
        }                                                                      \
        for (i = 0; i < *arr_len / 4; i++) {                                   \
            p3_add_to_2_interval_array(&st,                                    \
                   (int)PyLong_AsLong(PyList_GetItem(p_obj, i)),               \
                   (int)PyLong_AsLong(PyList_GetItem(p_obj, i+1)),             \
                   (int)PyLong_AsLong(PyList_GetItem(p_obj, i+2)),             \
                   (int)PyLong_AsLong(PyList_GetItem(p_obj, i+3)));            \
        }                                                                      \
    }                                                                          \

#define DICT_GET_AND_COPY_TO_INTERVAL_ARRAY(o, d, k, st)                       \
    if (DICT_GET_OBJ(o, d, k)){                                                \
        if (!PyList_Check(o)){                                                 \
            PyErr_Format(PyExc_TypeError,                                      \
                            "Value of %s is not of type list", k);             \
            return NULL;}                                                      \
        st.count = 0;                                                          \
        *arr_len = PyList_Size(p_obj);                                         \
        if (!arr_len % 2) {                                                    \
            PyErr_Format(PyExc_TypeError,                                      \
                            "%s must be linear multiple of 2 in length", k);   \
            return NULL;                                                       \
        }                                                                      \
        for (i = 0; i < *arr_len / 2; i++) {                                   \
            p3_add_to_interval_array(&st,                                      \
                 (int)PyLong_AsLong(PyList_GetItem(p_obj, i)),                 \
                 (int)PyLong_AsLong(PyList_GetItem(p_obj, i+1)));              \
        }                                                                      \
    }                                                                          \


p3_global_settings*
setGlobalParams(PyObject *self, PyObject *p3s_dict) {
    /* Creates a new p3_global_settings struct and initializes it with 
     * defaults using p3_create_global_settings() from libprimer3.c.
     * Parses the user-provided settings from p3_settings_dict and 
     * overwrites the defaults (note that minimal error checking is 
     * performed in this function). If there is an error during the process
     * (e.g., a param is not of the correct time), `err` will be set to the
     * error string and the function will return NULL.
     */

    p3_global_settings      *pa;
    PyObject                *p_obj;

    if (!(pa = p3_create_global_settings())) {
        PyErr_SetString(PyExc_IOError, "Could not allocate memory for p3 globals");
        return NULL;
    }

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
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_OPT_TM", pa->p_args.opt_tm);
    DICT_GET_AND_ASSIGN_INT(p_obj, p3s_dict, "PRIMER_OPT_GC_PERCENT", pa->p_args.opt_gc_content);
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
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_PAIR_MAX_COMPL_END_TH", pa->p_args.divalent_conc);
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
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_SEQUENCING_LEAD", pa->sequencing.lead);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_SEQUENCING_SPACING", pa->sequencing.spacing);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_SEQUENCING_INTERVAL", pa->sequencing.interval);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_SEQUENCING_ACCURACY", pa->sequencing.accuracy);
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
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_MUST_MATCH_FIVE_PRIME", &pa->p_args.must_match_five_prime);
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_MUST_MATCH_THREE_PRIME", &pa->p_args.must_match_three_prime);
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME", &pa->o_args.must_match_five_prime);
    DICT_GET_AND_COPY_STR(p_obj, p3s_dict, "PRIMER_INTERNAL_MUST_MATCH_THREE_PRIME", &pa->o_args.must_match_three_prime);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_TM_GT", pa->p_args.weights.temp_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_TM_LT", pa->p_args.weights.temp_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_GC_PERCENT_GT", pa->p_args.weights.gc_content_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_GC_PERCENT_LT", pa->p_args.weights.gc_content_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SIZE_LT", pa->p_args.weights.gc_content_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SIZE_GT", pa->p_args.weights.length_lt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SELF_ANY", pa->p_args.weights.length_gt);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SELF_END", pa->p_args.weights.compl_any);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SELF_ANY_TH", pa->p_args.weights.compl_end);
    DICT_GET_AND_ASSIGN_DOUBLE(p_obj, p3s_dict, "PRIMER_WT_SELF_END_TH", pa->p_args.weights.compl_any_th);
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

    return pa;
}

seq_lib*
createSeqLib(PyObject *self, PyObject *seq_dict){
    /* Generates a library of sequences for mispriming checks.
     * Input is a Python dictionary with <seq name: sequence> key value
     * pairs. 
     */

    seq_lib                 *sl;
    PyObject                *py_seq_name, *py_seq;
    Py_ssize_t              pos;
    char                    *seq_name, *seq, *errfrag=NULL;

    if (!(sl = create_empty_seq_lib())) {
        PyErr_SetString(PyExc_IOError, "Could not allocate memory for seq_lib");
        return NULL;
    }

    pos = 0;
    while (PyDict_Next(seq_dict, &pos, &py_seq_name, &py_seq)) {
        seq_name = PyBytes_AsString(py_seq_name);
        seq = PyBytes_AsString(py_seq);
        if(add_seq_and_rev_comp_to_seq_lib(sl, seq, seq_name, errfrag)) {
            PyErr_SetString(PyExc_IOError, errfrag);
            return NULL;
        }
    }
    return sl;
}


seq_args*
createSeqArgs(PyObject *self, PyObject *sa_dict){

    seq_args                *sa;
    PyObject                *p_obj;
    int                     i, *arr_len=NULL;   

    if (!(sa = create_seq_arg())) {
        PyErr_SetString(PyExc_IOError, "Could not allocate memory for seq_args");
        return NULL;
    }
    // Sequence strings
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_TEMPLATE", &sa->sequence);
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_ID", &sa->right_input);
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_PRIMER", &sa->left_input);
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_PRIMER_REVCOMP", &sa->right_input);
    DICT_GET_AND_COPY_STR(p_obj, sa_dict, "SEQUENCE_INTERNAL_OLIGO", &sa->internal_input);
    DICT_GET_AND_COPY_ARRAY(p_obj, sa_dict, "SEQUENCE_QUALITY", &sa->quality, arr_len);
    if (sa->quality != NULL) {
        sa->n_quality = *arr_len;
    }
    DICT_GET_AND_COPY_TO_2_INTERVAL_ARRAY(p_obj, sa_dict, "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST", sa->ok_regions);
    DICT_GET_AND_COPY_TO_INTERVAL_ARRAY(p_obj, sa_dict, "SEQUENCE_TARGET", sa->tar2);
    DICT_GET_AND_COPY_TO_INTERVAL_ARRAY(p_obj, sa_dict, "SEQUENCE_EXCLUDED_REGION", sa->excl2);
    DICT_GET_AND_COPY_TO_INTERVAL_ARRAY(p_obj, sa_dict, "SEQUENCE_INTERNAL_EXCLUDED_REGION", sa->excl_internal2);
    DICT_GET_AND_COPY_ARRAY_INTO_ARRAY(p_obj, sa_dict,  "SEQUENCE_OVERLAP_JUNCTION_LIST", &sa->primer_overlap_junctions, arr_len);
    if (sa->primer_overlap_junctions != NULL) {
        sa->primer_overlap_junctions_count = *arr_len;
    }
    if DICT_GET_OBJ(p_obj, sa_dict, "SEQUENCE_INCLUDED_REGION") {
        if (!PyList_Check(p_obj)) {
            PyErr_SetString(PyExc_TypeError,\
                "Value of \"SEQUENCE_INCLUDED_REGION\" is not of type list");
            return NULL;
        } else if (PyList_Size(p_obj) != 2) {
            PyErr_SetString(PyExc_ValueError,\
                "Length of \"SEQUENCE_INCLUDED_REGION\" is not 2");
            return NULL;
        } else if (!PyLong_Check(PyList_GetItem(p_obj, 0)) || \
                   !PyLong_Check(PyList_GetItem(p_obj, 1))) {
            PyErr_SetString(PyExc_TypeError,\
                "\"SEQUENCE_INCLUDED_REGION\" contains non-int value");
            return NULL;                    
        }
        sa->incl_s = PyLong_AsLong(PyList_GetItem(p_obj, 0));
        sa->incl_l = PyLong_AsLong(PyList_GetItem(p_obj, 1));
    }
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_START_CODON_POSITION", sa->start_codon_pos);
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_FORCE_LEFT_START", sa->force_left_start);
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_FORCE_LEFT_END", sa->force_left_end);
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_FORCE_RIGHT_START", sa->force_right_start);
    DICT_GET_AND_ASSIGN_INT(p_obj, sa_dict, "SEQUENCE_FORCE_RIGHT_END", sa->force_right_end);

    return sa;
}
