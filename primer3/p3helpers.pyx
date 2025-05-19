# cython: language_level=3
# Copyright (C) 2023. Ben Pruitt & Nick Conway;
# See LICENSE for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
primer3.p3helpers
~~~~~~~~~~~~~~~~~~~~~~

Contains Cython functions and classes that aid primer design

'''
from cpython.mem cimport (
    PyMem_Free,
    PyMem_Malloc,
)
from cython cimport typeof
from libc.string cimport memcpy


cdef extern from "p3helpers.h":
    const char COMP_BASE_LUT[128]
    const char SANITIZE_LUT[128]
    const char ACGT_UPPER_LUT[128]
    const char IS_UPPERCASE_ACGT_LUT[128]

cdef:
    char INVALID_BASE = 63 # '?' character


cdef inline void copy_c_char_buffer(
    const char* in_buf,
    char* out_buf,
    Py_ssize_t str_length,
) noexcept nogil:
    '''Copy string to buffer without null terminator since we track length.

    Args:
        in_buf: input buffer
        out_buf: output buffer (must be at least str_length bytes)
        str_length: number of bytes to copy
    '''
    memcpy(out_buf, in_buf, str_length)


cdef int ss_rev_comp_seq(
    char* seq,
    Py_ssize_t length,
) except -1 nogil:
    '''Reverse complement sequence C string

    Args:
        seq: sequence to sanitize
        length: length of sequence

    Returns:
        0 on success, 1 on error (Invalid base)
    '''
    cdef:
        char        temp
        char*       end_ptr = NULL

    end_ptr = seq + (length - 1)

    while end_ptr > seq:
        temp = end_ptr[0]
        end_ptr[0] = COMP_BASE_LUT[<unsigned char> seq[0]]
        if end_ptr[0] == INVALID_BASE:
            return 1
        end_ptr -= 1
        seq[0] = COMP_BASE_LUT[<unsigned char> temp]
        if seq[0] == INVALID_BASE:
            return 1
        seq += 1


    if length % 2:
        seq[0] = COMP_BASE_LUT[<unsigned char> seq[0]]

    return 0


cdef int ss_sanitize_seq(
    char* seq,
    Py_ssize_t length,
) except -1 nogil:
    '''Sanitize sequence C string IUPAC non-{A,C,G,T,a,c,g,t} bases with {N, n}

    Args:
        seq: sequence to sanitize
        length: length of sequence

    Returns:
        0 on success, 1 on error (Invalid base)
    '''
    cdef:
        char*       end_ptr = seq + length

    while seq < end_ptr:
        seq[0] = SANITIZE_LUT[<unsigned char> seq[0]]
        if seq[0] == INVALID_BASE:
            return 1
        seq += 1
    return 0


cdef int ss_needs_conversion(
    const char* seq,
    Py_ssize_t length,
) except -1 nogil:
    '''Fast check if sequence needs conversion (contains lowercase or non-ACGT chars)
    using lookup table for O(1) per-character checks.

    Args:
        seq: sequence to check
        length: length of sequence

    Returns:
        1 if sequence needs conversion, 0 if it's already uppercase ACGT
    '''
    cdef:
        const char* end_ptr = seq + length

    while seq < end_ptr:
        if not IS_UPPERCASE_ACGT_LUT[<unsigned char> seq[0]]:
            return 1
        seq += 1
    return 0


cdef int ss_ensure_acgt_upper(
    char* seq,
    Py_ssize_t length,
) except -1 nogil:
    '''Convert sequence to uppercase and validate it contains only A, C, G, T bases

    Args:
        seq: sequence to convert to uppercase and validate
        length: length of sequence

    Returns:
        0 on success, 1 on error (non-ACGT base)
    '''
    cdef:
        char*       end_ptr = seq + length

    while seq < end_ptr:
        seq[0] = ACGT_UPPER_LUT[<unsigned char> seq[0]]
        if seq[0] == INVALID_BASE:
            return 1
        seq += 1
    return 0


cpdef str reverse_complement(
    str seq,
    bint do_sanitize = False,
):
    '''Compute reverse complement of the python string sequence

    Args:
        seq: sequence string
        do_sanitize: If True, convert non-IUPAC characters to N's

    Returns:
        Reverse complement of sequence string

    Raises:
        ValueError: Invalid base in sequence
    '''
    cdef:
        char* seq_c = NULL
        Py_ssize_t seq_len = len(seq)
        int check = 0

    seq_b = seq.encode('utf8')
    # convert to C str
    seq_c = seq_b

    check = ss_rev_comp_seq(seq_c, seq_len)
    if check:
        raise ValueError(f'Invalid base in sequence {seq}')
    if do_sanitize:
        check = ss_sanitize_seq(seq_c, seq_len)
        if check:
            raise ValueError(f'Invalid base in sequence {seq}')
    seq_b = seq_c

    return seq_b.decode('utf8')


cpdef bytes reverse_complement_b(
    bytes seq,
    bint do_sanitize = False,
):
    '''Compute reverse complement of the bytes sequence

    Args:
        seq: sequence in bytes
        do_sanitize: If True, convert non-IUPAC characters to N's

    Returns:
        Reverse complement of sequence in bytes

    Raises:
        ValueError: Invalid base in sequence
        OSError: ``malloc`` failure
    '''
    cdef:
        char* seq_c = NULL
        Py_ssize_t seq_len = len(seq)
        char* seq_operate_c = NULL
        bytes seq_out
        int check = 0

    # Allocate one extra byte for null termination
    seq_operate_c = <char *> PyMem_Malloc((seq_len + 1) * sizeof(char))
    if seq_operate_c == NULL:
        raise OSError('malloc failure')

    # convert to C str
    seq_c = seq

    # MUST COPY as python does not mutate bytes
    copy_c_char_buffer(seq_c, seq_operate_c, seq_len)
    # Ensure null termination
    seq_operate_c[seq_len] = 0

    check = ss_rev_comp_seq(seq_operate_c, seq_len)
    if check:
        PyMem_Free(seq_operate_c)
        seq_operate_c = NULL
        raise ValueError(f'Invalid base in sequence {seq}')
    if do_sanitize:
        check = ss_sanitize_seq(seq_operate_c, seq_len)
        if check:
            PyMem_Free(seq_operate_c)
            seq_operate_c = NULL
            raise ValueError(f'Invalid base in sequence {seq}')
    # Create bytes object of exact length without null terminator
    seq_out = seq_operate_c[:seq_len]
    PyMem_Free(seq_operate_c)
    seq_operate_c = NULL
    return seq_out


cpdef str sanitize_sequence(
    str seq,
):
    '''Sanitize sequence with non-standard bases with `N`s
    IUPAC {R,Y,M,K,S,W,H,D,B,V,r,y,m,k,s,w,h,d,b,v}

    IUPAC non-{A,C,G,T,a,c,g,t} bases become {N, n}

    NOTE: Consider keeping a copy of original ``seq`` argument for record
    keeping

    Args:
        seq: sequence to sanitize nonstandard with `N`s or `n`s

    Returns:
        sanitized version of ``seq``

    Raises:
        ValueError: Invalid base in sequence
    '''
    cdef:
        char* seq_c = NULL
        Py_ssize_t seq_len = len(seq)
        int check = 0

    seq_b = seq.encode('utf8')
    # convert to C str
    seq_c = seq_b

    check = ss_sanitize_seq(seq_c, seq_len)
    if check:
        raise ValueError(f'Invalid base in equence {seq}')
    seq_b = seq_c

    return seq_b.decode('utf8')


cpdef bytes sanitize_sequence_b(
    bytes seq,
):
    '''Sanitize bytes sequence with non-standard bases with `N`s
    IUPAC {R,Y,M,K,S,W,H,D,B,V,r,y,m,k,s,w,h,d,b,v}

    IUPAC non-{A,C,G,T,a,c,g,t} bases become {N, n}

    NOTE: Consider keeping a copy of original ``seq`` argument for record
    keeping

    Args:
        seq: sequence to sanitize nonstandard with `N`s or `n`s

    Returns:
        sanitized version of ``seq``

    Raises:
        ValueError: Invalid base in sequence
        OSError: ``malloc`` failure
    '''
    cdef:
        char* seq_c = NULL
        Py_ssize_t seq_len = len(seq)
        char* seq_operate_c = NULL
        bytes seq_out
        int check = 0

    # Allocate one extra byte for null termination
    seq_operate_c = <char *> PyMem_Malloc((seq_len + 1) * sizeof(char))
    if seq_operate_c == NULL:
        raise OSError('malloc failure')

    # convert to C str
    seq_c = seq

    # MUST COPY as python does not mutate bytes
    copy_c_char_buffer(seq_c, seq_operate_c, seq_len)
    # Ensure null termination
    seq_operate_c[seq_len] = 0

    check = ss_sanitize_seq(seq_operate_c, seq_len)
    if check:
        PyMem_Free(seq_operate_c)
        seq_operate_c = NULL
        raise ValueError(f'Invalid base in sequence {seq}')
    # Create bytes object of exact length without null terminator
    seq_out = seq_operate_c[:seq_len]
    PyMem_Free(seq_operate_c)
    seq_operate_c = NULL
    return seq_out


cpdef str ensure_acgt_uppercase(
    str seq,
):
    '''Convert sequence to uppercase and validate it contains only A, C, G, T bases.
    This is stricter than sanitize_sequence() as it only allows A, C, G, T bases.

    If the sequence is already uppercase ACGT, returns the original string without
    any allocation or copying.

    Args:
        seq: sequence to convert to uppercase and validate

    Returns:
        uppercase version of seq containing only A, C, G, T bases

    Raises:
        ValueError: If sequence contains any characters other than A, C, G, T (case insensitive)
    '''
    cdef:
        char* seq_c = NULL
        Py_ssize_t seq_len = len(seq)
        int needs_conversion = 0

    if seq_len == 0:
        return seq

    seq_b = seq.encode('utf8')
    seq_c = seq_b

    # Fast path: check if sequence needs conversion
    needs_conversion = ss_needs_conversion(seq_c, seq_len)
    if not needs_conversion:
        return seq  # Return original if already uppercase ACGT

    # Slow path: needs conversion
    if ss_ensure_acgt_upper(seq_c, seq_len):
        # Find the problematic character
        for i, c in enumerate(seq):
            if c.upper() not in 'ACGT':
                raise ValueError(
                    f"Sequence contains non-ACGT base '{c}' at position {i}: {seq}"
                )

    return seq_b.decode('utf8')


cpdef bytes ensure_acgt_uppercase_b(
    bytes seq,
):
    '''Convert bytes sequence to uppercase and validate it contains only A, C, G, T bases.
    This is stricter than sanitize_sequence_b() as it only allows A, C, G, T bases.

    If the sequence is already uppercase ACGT, returns the original bytes without
    any allocation or copying.

    Args:
        seq: sequence to convert to uppercase and validate

    Returns:
        uppercase version of seq containing only A, C, G, T bases

    Raises:
        ValueError: If sequence contains any characters other than A, C, G, T (case insensitive)
        OSError: malloc failure
    '''
    cdef:
        char* seq_c = NULL
        Py_ssize_t seq_len = len(seq)
        char* seq_operate_c = NULL
        bytes seq_out
        int needs_conversion = 0

    if seq_len == 0:
        return seq

    # Fast path: check if sequence needs conversion
    seq_c = seq
    needs_conversion = ss_needs_conversion(seq_c, seq_len)
    if not needs_conversion:
        return seq  # Return original if already uppercase ACGT

    # Slow path: needs conversion
    seq_operate_c = <char *> PyMem_Malloc(seq_len * sizeof(char))
    if seq_operate_c == NULL:
        raise OSError('malloc failure')

    try:
        copy_c_char_buffer(seq_c, seq_operate_c, seq_len)

        if ss_ensure_acgt_upper(seq_operate_c, seq_len):
            # Find the problematic character
            for i, c in enumerate(seq.decode('utf8')):
                if c.upper() not in 'ACGT':
                    raise ValueError(
                        f"Sequence contains non-ACGT base '{c}' at position {i}: {seq}"
                    )

        # Create bytes object with the correct length
        seq_out = seq_operate_c[:seq_len]
        return seq_out
    finally:
        PyMem_Free(seq_operate_c)
        seq_operate_c = NULL
