# cython: language_level=3
# Copyright (C) 2023. Ben Pruitt & Nick Conway;
# See LICENSE for full GPLv2 license.

# Declare external C functions and variables
cdef extern from "p3helpers.h":
    const char COMP_BASE_LUT[128]
    const char SANITIZE_LUT[128]
    const char ACGT_UPPER_LUT[128]
    const char IS_UPPERCASE_ACGT_LUT[128]

# Declare internal C functions
cdef int ss_rev_comp_seq(char* seq, Py_ssize_t length) except -1 nogil
cdef int ss_sanitize_seq(char* seq, Py_ssize_t length) except -1 nogil
cdef int ss_needs_conversion(const char* seq, Py_ssize_t length) except -1 nogil
cdef int ss_ensure_acgt_upper(char* seq, Py_ssize_t length) except -1 nogil

# Declare Python-visible functions - these must match EXACTLY how they're defined in pyx
cpdef str reverse_complement(str seq, bint do_sanitize=*)
cpdef bytes reverse_complement_b(bytes seq, bint do_sanitize=*)
cpdef str sanitize_sequence(str seq)
cpdef bytes sanitize_sequence_b(bytes seq)
cpdef str ensure_acgt_uppercase(str seq)
cpdef bytes ensure_acgt_uppercase_b(bytes seq)
