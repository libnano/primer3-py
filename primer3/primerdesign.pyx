# cython: language_level=3
# Copyright (C) 2014-2023. Ben Pruitt & Nick Conway; Wyss Institute
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
primer3.primerdesign
Python C API bindings for the Primer3 design engine
'''
cdef:
    p3_global_settings* global_settings_data = NULL
    seq_args* sequence_args_data = NULL

import atexit
import os.path
import sys
from typing import (
    Any,
    Dict,
    Optional,
    Union,
)

from primer3 import argdefaults

_DID_LOAD_THERM_PARAMS = False
_DEFAULT_WORD_LEN_2 = 16  # see masker.h

def get_dunder_file() -> str:
    return __file__


def load_thermo_params():
    global _DID_LOAD_THERM_PARAMS

    cdef:
        char*           p3_cfg_path_bytes_c
        thal_results    thalres
        thal_parameters thermodynamic_parameters

    if _DID_LOAD_THERM_PARAMS is True:
        return

    p3_cfg_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'src',
        'libprimer3',
        'primer3_config',
        '',  # Add trailing slash (OS-ind) req'd by primer3 lib
    )

    # read default thermodynamic parameters
    p3_cfg_path_bytes = p3_cfg_path.encode('utf-8')
    p3_cfg_path_bytes_c = p3_cfg_path_bytes

    thal_set_null_parameters(&thermodynamic_parameters)
    thal_load_parameters(p3_cfg_path_bytes_c, &thermodynamic_parameters, &thalres)
    # set_default_thal_parameters(&thermodynamic_parameters)
    try:
        if get_thermodynamic_values(&thermodynamic_parameters, &thalres) != 0:
            raise OSError(
                f'Could not load thermodynamic config file {p3_cfg_path}'
            )
    finally:
        thal_free_parameters(&thermodynamic_parameters)
    _DID_LOAD_THERM_PARAMS = True

def _cleanup():
    destroy_thal_structures()


cdef int pr_default_position_penalties(const p3_global_settings* pa):
    if (
        (pa[0].inside_penalty == PR_DEFAULT_INSIDE_PENALTY) and
        (pa[0].outside_penalty == PR_DEFAULT_OUTSIDE_PENALTY)
    ):
        return 1
    return 0


cdef int pdh_wrap_set_seq_args_globals(
        p3_global_settings* global_settings_data,
        seq_args* sequence_args_data,
        object kmer_lists_path,
        char* in_buffer,
) except -1:
    '''
    Creates a new p3_global_settings struct and initializes it with
    defaults using p3_create_global_settings() from libprimer3.c.
    Parses the user-provided settings from p3_settings_dict and
    overwrites the defaults (note that minimal error checking is
    performed in this function). If there is an error during the process
    (e.g., a param is not of the correct type), the python error string will
    be set and the function will return NULL.

    Args:
        global_settings_data: pointer to p3_global_settings data structure
        seq_args: pointer to seq_args data structure
        kmer_lists_path: string path to kmer list directory
        in_buffer: string buffer that is the seq_ar

    Raises:
        ValueError: Error parsing the data
    '''
    cdef:
        # Setup the input data structures handlers
        int strict_tags = 0
        int io_version = 4
        int echo_output = 0
        char*  kmer_lists_path_c = NULL

        read_boulder_record_results read_boulder_record_res
        pr_append_str p3_settings_path
        pr_append_str output_path
        pr_append_str error_path
        pr_append_str fatal_parse_err
        pr_append_str nonfatal_parse_err
        pr_append_str warnings

    read_boulder_record_res.explain_flag = 0
    read_boulder_record_res.file_flag = 0

    init_pr_append_str(&fatal_parse_err)
    init_pr_append_str(&nonfatal_parse_err)
    init_pr_append_str(&warnings)
    init_pr_append_str(&p3_settings_path)
    init_pr_append_str(&output_path)
    init_pr_append_str(&error_path)

    read_boulder_record(
        NULL,
        &strict_tags,
        &io_version,
        echo_output,
        p3_file_type.all_parameters,
        global_settings_data,
        sequence_args_data,
        &fatal_parse_err,
        &nonfatal_parse_err,
        &warnings,
        &read_boulder_record_res,
        in_buffer
    )

    # NOTE: Masking with PRIMER_MASK_KMERLIST_PATH is parsed in a non standard
    # way in primer3_boulder_main.c therefore we need to do validation twice
    # here and in argdefaults.py
    if sys.platform != 'windows':
        if global_settings_data[0].mask_template:
            global_settings_data[0].lowercase_masking = global_settings_data[0].mask_template

        # Check that we found the kmer lists in case masking flag was set to 1.
        if (
            (global_settings_data[0].mask_template == 1) and
            (kmer_lists_path == '')
        ):
            raise ValueError(
                'masking template chosen, but path to '
                'PRIMER_MASK_KMERLIST_PATH not specified'
            )

        # Set up some masking parameters
        if global_settings_data[0].mask_template == 1:
            global_settings_data[0].mp.window_size = _DEFAULT_WORD_LEN_2

            if global_settings_data[0].pick_right_primer == 0:
                global_settings_data[0].mp.mdir = masking_direction.fwd
            elif global_settings_data[0].pick_left_primer == 0:
                global_settings_data[0].mp.mdir = masking_direction.rev
            # Check if masking parameters (k-mer list usage) have changed
            if global_settings_data[0].masking_parameters_changed == 1:
                delete_formula_parameters(
                    global_settings_data[0].mp.fp,
                    global_settings_data[0].mp.nlists,
                )
                if isinstance(kmer_lists_path, str):
                    kmer_lists_path_b = kmer_lists_path.encode('utf8')
                else:
                    kmer_lists_path_b = kmer_lists_path
                kmer_lists_path_c = kmer_lists_path_b
                global_settings_data[0].mp.fp = create_default_formula_parameters(
                    global_settings_data[0].mp.list_prefix,
                    kmer_lists_path_c,
                    &fatal_parse_err,
                )
                global_settings_data[0].masking_parameters_changed = 0

    if (
        (global_settings_data[0].primer_task == task.generic) and
        (global_settings_data[0].pick_internal_oligo == 1)
    ):
        if not global_settings_data[0].pick_internal_oligo:
            raise ValueError(
                'global_settings_data[0].pick_internal_oligo must be set'
            )

    if nonfatal_parse_err.data != NULL:
        err_msg_b = <bytes> nonfatal_parse_err.data
        raise ValueError(err_msg_b.decode('utf8'))
    if fatal_parse_err.data != NULL:
        err_msg_b = <bytes> fatal_parse_err.data
        raise ValueError(err_msg_b.decode('utf8'))
    return 0


cdef seq_lib* pdh_create_seq_lib(object seq_dict) except NULL:
    '''
    Generates a library of sequences for mispriming checks.
    Input is a Python dictionary with <seq name: sequence> key value
    pairs. Returns NULL and sets the Python error string on failure.

    Args:
        seq_dict: Sequence disctionary in the format <seq name: sequence>

    Returns:
        pointer to a generated seq_lib

    Raises:
        OSError: Could not allocate memory for seq_lib
        TypeError: Cannot add seq name with non-Unicode/Bytes type to seq_lib
        OSError: primer3 internal error
    '''

    cdef:
        seq_lib* sl = NULL
        char* seq_name_c = NULL
        char* seq_c = NULL
        char* errfrag = NULL

    if sl == NULL:
        raise OSError('Could not allocate memory for seq_lib')

    for seq_name_str, seq_str in seq_dict.items():
        if isinstance(seq_name_str, str):
            seq_name_b = seq_name_str.encode('utf8')
            seq_name = seq_name_b
        elif isinstance(seq_name_str, bytes):
            seq_name = seq_name
        else:
            destroy_seq_lib(sl)
            raise TypeError(
                'Cannot add seq name with non-Unicode/Bytes type to seq_lib',
            )

        if isinstance(seq_str, str):
            seq_b = seq_str.encode('utf8')
            seq_c = seq_b
        elif isinstance(seq_name_str, bytes):
            seq_c = seq_str
        else:
            destroy_seq_lib(sl)
            raise TypeError(
                'Cannot add seq with non-Unicode/Bytes type to seq_lib',
            )

        if add_seq_to_seq_lib(sl, seq_c, seq_name_c, errfrag) == 1:
            err_msg_b =  <bytes> errfrag
            destroy_seq_lib(sl)
            raise OSError(err_msg_b.decode('utf8'))
    reverse_complement_seq_lib(sl)
    return sl


cdef object pdh_design_output_to_dict(
        const p3_global_settings* global_settings_data,
        const seq_args* sequence_args_data,
        const p3retval *retval,
):
    '''
    Args:
        global_settings_data: primer3 p3_global_settings data pointer
        sequence_args_data: primer3 design seq_args data pointer
        retval: primer3 design return value pointer

    Returns:
        converted Python dictionary of design output created from the ``retval``
        data

    Raises:
        OSError: memory issue
    '''
    cdef:
        # The pointers to warning tag
        char* warning = NULL

        # A place to put a string containing all error messages
        pr_append_str* combined_retval_err = NULL

        # Pointers for the primer set just printing
        primer_rec* fwd = NULL
        primer_rec* rev = NULL
        primer_rec* intl = NULL

        # Variables only used for Primer Lists
        int num_fwd, num_rev, num_int, num_pair
        int num_print = 0
        int print_fwd = 0
        int print_rev = 0
        int print_int = 0

        # Switches for printing this primer
        int go_fwd = 0
        int go_rev = 0
        int go_int = 0

        double temp_double = 0

        # The number of loop cycles
        int loop_max

        # That links to the included region
        int i
        int incl_s = sequence_args_data[0].incl_s

        int product_size = 0

        # This deals with the renaming of the internal oligo
        new_oligo_name = "INTERNAL"
        int_oligo = new_oligo_name

    output_dict: Dict[str, Any] = {}

    # Check if there are warnings and print them
    warning = p3_get_rv_and_gs_warnings(retval, global_settings_data)
    if warning != NULL:
        warning_b = <bytes> warning
        output_dict["PRIMER_WARNING"] = warning_b.decode('utf8')
        free(warning)
        warning = NULL

    combined_retval_err = create_pr_append_str()
    if combined_retval_err == NULL:
        raise OSError("Primer3 ran out of memory.")

    try:
        if pr_append_new_chunk_external(combined_retval_err, retval[0].glob_err.data):
            raise OSError("Primer3 ran out of memory.")

        # NOTE: These are non fatal errors
        if pr_append_new_chunk_external(combined_retval_err, retval[0].per_sequence_err.data):
            raise OSError("Primer3 ran out of memory.")

        # Check if there are errors, print and return
        if not pr_is_empty(combined_retval_err):
            err_msg_b = <bytes> pr_append_str_chars(combined_retval_err)
            raise OSError(err_msg_b.decode('utf8'))

    finally:
        destroy_pr_append_str(combined_retval_err)
        combined_retval_err = NULL

    # Get how many primers are in the array
    num_fwd = retval[0].fwd.num_elem
    num_rev = retval[0].rev.num_elem
    num_int = retval[0].intl.num_elem
    num_pair = retval[0].best_pairs.num_pairs

    # Prints out selection statistics about the primers
    if (
        (global_settings_data[0].pick_left_primer == 1) and
        not (global_settings_data[0].pick_anyway and sequence_args_data[0].left_input)
    ):
        explain_str_b = <bytes> p3_get_oligo_array_explain_string(
            p3_get_rv_fwd(retval),
        )
        output_dict['PRIMER_LEFT_EXPLAIN'] = explain_str_b.decode('utf8')

    if (
        (global_settings_data[0].pick_right_primer == 1) and
        not (global_settings_data[0].pick_anyway and sequence_args_data[0].right_input)
    ):
        explain_str_b = <bytes> p3_get_oligo_array_explain_string(
            p3_get_rv_rev(retval),
        )
        output_dict['PRIMER_RIGHT_EXPLAIN'] = explain_str_b.decode('utf8')

    if (
        (global_settings_data[0].pick_internal_oligo == 1) and
        not (global_settings_data[0].pick_anyway and sequence_args_data[0].internal_input)
    ):
        explain_str_b = <bytes> p3_get_oligo_array_explain_string(
            p3_get_rv_intl(retval),
        )
        output_dict['PRIMER_INTERNAL_EXPLAIN'] = explain_str_b.decode('utf8')

    if (
        (global_settings_data[0].pick_right_primer == 1) and
        (global_settings_data[0].pick_left_primer == 1)
    ):
        explain_str_b = <bytes> p3_get_pair_array_explain_string(
            p3_get_rv_best_pairs(retval),
        )
        output_dict['PRIMER_PAIR_EXPLAIN'] = explain_str_b.decode('utf8')

    # Print out the stop codon if a reading frame was specified
    if not PR_START_CODON_POS_IS_NULL(sequence_args_data):
        stop_codon_pos = retval[0].stop_codon_pos
        output_dict['PRIMER_STOP_CODON_POSITION'] = stop_codon_pos

    # How often has the loop to be done?
    if retval[0].output_type == p3_output_type.primer_list:
        # For Primer Lists: Figure out how many primers are in
        # the array that can be printed. If more than needed,
        #  set it to the number requested.
        #  Get how may primers should be printed
        num_print = global_settings_data[0].num_return
        # Set how many primers will be printed
        print_fwd = num_print if (num_print < num_fwd) else num_fwd
        print_rev = num_print if (num_print < num_rev) else  num_rev
        print_int = num_print if (num_print < num_int) else num_int
        # Get which list has to print most primers
        loop_max = 0
        if loop_max < print_fwd:
            loop_max = print_fwd
        if loop_max < print_rev:
            loop_max = print_rev
        if loop_max < print_int:
            loop_max = print_int

        # Now the vars are there how often we have to go
        # through the loop and how many of each primer can
        #  be printed
        num_pair = 0
    else:
        loop_max = num_pair
        # Set how many primers will be printed
        print_fwd = num_pair
        print_rev = num_pair
        if num_int != 0:
            print_int = num_pair

    # Save the number of each type of oligo that was found
    output_dict['PRIMER_LEFT_NUM_RETURNED'] = print_fwd
    output_dict['PRIMER_RIGHT_NUM_RETURNED'] = print_rev

    output_dict[f'PRIMER_{int_oligo}_NUM_RETURNED'] = print_int
    output_dict['PRIMER_PAIR_NUM_RETURNED'] = num_pair

    # Start of the loop printing all pairs or primers or oligos
    for i in range(loop_max):
        # What needs to be printed the conditions for primer lists
        if retval[0].output_type == p3_output_type.primer_list:
            # Attach the selected primers to the pointers
            fwd = &(retval[0].fwd.oligo[i])
            rev = &(retval[0].rev.oligo[i])
            intl = &(retval[0].intl.oligo[i])

            # Do fwd oligos have to be printed?
            if (global_settings_data[0].pick_left_primer) and (i < print_fwd):
                go_fwd = 1
            else:
                go_fwd = 0

            # Do rev oligos have to be printed?
            if (global_settings_data[0].pick_right_primer) and (i < print_rev):
                go_rev = 1
            else:
                go_rev = 0

            # Do int oligos have to be printed?
            if (global_settings_data[0].pick_internal_oligo) and (i < print_int):
                go_int = 1
            else:
                go_int = 0

        else:
            # We will print primer pairs or pairs plus internal oligos
            #  Get pointers to the primer_rec's that we will print
            # Pairs must have fwd and rev primers
            fwd  = retval[0].best_pairs.pairs[i].left
            rev  = retval[0].best_pairs.pairs[i].right
            intl = retval[0].best_pairs.pairs[i].intl
            go_fwd = 1
            go_rev = 1
            # Do hyb oligos have to be printed?
            if (global_settings_data[0].pick_internal_oligo == 1):
                go_int = 1
            else:
                go_int = 0

        # Print out the Pair Penalties
        if retval[0].output_type == p3_output_type.primer_pairs:
            temp_double = retval[0].best_pairs.pairs[i].pair_quality
            output_dict[f'PRIMER_PAIR_{i}_PENALTY'] = temp_double

        # Print single primer penalty
        if go_fwd == 1:
            temp_double = fwd[0].quality
            output_dict[f'PRIMER_LEFT_{i}_PENALTY'] = temp_double

        if go_rev == 1:
            temp_double = rev[0].quality
            output_dict[f'PRIMER_RIGHT_{i}_PENALTY'] = temp_double

        if go_int == 1:
            temp_double = intl[0].quality
            output_dict[f'PRIMER_RIGHT{int_oligo}_{i}_PENALTY'] = temp_double

        # Print the oligo_problems
        if (go_fwd == 1) and p3_ol_has_any_problem(fwd):
            problem_b = <bytes> p3_get_ol_problem_string(fwd)
            output_dict[f'PRIMER_LEFT_{i}_PROBLEMS'] = problem_b

        if (go_rev == 1) and p3_ol_has_any_problem(rev):
            problem_b = <bytes> p3_get_ol_problem_string(rev)
            output_dict[f'PRIMER_RIGHT_{i}_PROBLEMS'] = problem_b

        if (go_int == 1) and p3_ol_has_any_problem(intl):
            problem_b = <bytes> p3_get_ol_problem_string(intl)
            output_dict[f'PRIMER_RIGHT{int_oligo}_{i}_PROBLEMS'] = problem_b


        # Print primer sequences.
        if go_fwd == 1:
            sqtemp_b = <bytes> pr_oligo_overhang_sequence(
                sequence_args_data,
                fwd,
            )
            output_dict[f'PRIMER_LEFT_{i}_SEQUENCE'] = sqtemp_b.decode('utf8')

        if go_rev == 1:
            sqtemp_b = <bytes> pr_oligo_rev_c_overhang_sequence(
                sequence_args_data,
                rev,
            )
            output_dict[f'PRIMER_RIGHT_{i}_SEQUENCE'] = sqtemp_b.decode('utf8')

        if go_int == 1:
            sqtemp_b = <bytes> pr_oligo_sequence(sequence_args_data, intl)
            output_dict[f'PRIMER_{int_oligo}_{i}_SEQUENCE'] = sqtemp_b.decode(
                'utf8',
            )

        # Print primer start and length
        if go_fwd == 1:
            output_dict[f'PRIMER_LEFT_{i}'] = [
                fwd[0].start + incl_s + global_settings_data[0].first_base_index,
                fwd[0].length,
            ]
        if go_rev == 1:
            output_dict[f'PRIMER_RIGHT_{i}'] = [
                rev[0].start + incl_s + global_settings_data[0].first_base_index,
                rev[0].length,
            ]
        if go_int == 1:
            output_dict[f'PRIMER_{int_oligo}_{i}'] = [
                intl[0].start + incl_s + global_settings_data[0].first_base_index,
                intl[0].length,
            ]

        # Print primer Tm
        if go_fwd == 1:
            output_dict[f'PRIMER_LEFT_{i}_TM'] = fwd[0].temp
        if go_rev == 1:
            output_dict[f'PRIMER_RIGHT_{i}_TM'] = rev[0].temp
        if go_int == 1:
            output_dict[f'PRIMER_{int_oligo}_{i}_TM'] = intl[0].temp

        # Print fraction bound at melting temperature
        if (
            (global_settings_data[0].annealing_temp > 0.0) and
            (global_settings_data[0].salt_corrections != 2)
        ):
            if (go_fwd == 1) and (fwd[0].bound > 0.0):
                output_dict[f'PRIMER_LEFT_{i}_BOUND'] = fwd[0].bound
            if (go_rev == 1) and (rev[0].bound > 0.0):
                output_dict[f'PRIMER_RIGHT_{i}_BOUND'] = rev[0].bound
            if (go_int == 1) and (intl[0].bound > 0.0):
                output_dict[f'PRIMER_{int_oligo}_{i}_BOUND'] = intl[0].bound

        # Print primer GC content
        if go_fwd == 1:
            output_dict[f'PRIMER_LEFT_{i}_GC_PERCENT'] = fwd[0].gc_content
        if go_rev == 1:
            output_dict[f'PRIMER_RIGHT_{i}_GC_PERCENT'] = rev[0].gc_content

        if go_int == 1:
            output_dict[f'PRIMER_{int_oligo}_{i}_GC_PERCENT'] = intl[0].gc_content

        # Print primer self_any
        if (go_fwd == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 0):
            output_dict[f'PRIMER_LEFT_{i}_SELF_ANY'] = fwd[0].self_any
        if (go_rev == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 0):
            output_dict[f'PRIMER_RIGHT_{i}_SELF_ANY'] = rev[0].self_any
        if (go_int == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 0):
            output_dict[f'PRIMER_{int_oligo}_{i}_SELF_ANY'] = intl[0].self_any

        # Print primer self_any thermodynamical approach
        if (go_fwd == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 1):
            output_dict[f'PRIMER_LEFT_{i}_SELF_ANY_TH'] = fwd[0].self_any
        if (go_rev == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 1):
            output_dict[f'PRIMER_RIGHT_{i}_SELF_ANY_TH'] = rev[0].self_any
        if (go_int == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 1):
            output_dict[f'PRIMER_{int_oligo}_{i}_SELF_ANY_TH'] = intl[0].self_any

        # Print primer secondary structures*/
        if (
            (go_fwd == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (fwd[0].self_any_struct != NULL)
        ):
            sqtemp_b = <bytes> fwd[0].self_any_struct
            output_dict[f'PRIMER_LEFT_{i}_SELF_ANY_STUCT'] = sqtemp_b.decode('utf8')
        if (
            (go_rev == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (rev[0].self_any_struct != NULL)
        ):
            sqtemp_b = <bytes> rev[0].self_any_struct
            output_dict[f'PRIMER_RIGHT_{i}_SELF_ANY_STUCT'] = sqtemp_b.decode('utf8')
        if (
            (go_int == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (intl[0].self_any_struct != NULL)
        ):
            sqtemp_b = <bytes> intl[0].self_any_struct
            output_dict[f'PRIMER_{int_oligo}_{i}_SELF_ANY_STUCT'] = sqtemp_b.decode('utf8')


        # Print primer self_end
        if (go_fwd == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 0):
            output_dict[f'PRIMER_LEFT_{i}_SELF_END'] = fwd[0].self_end

        if (go_rev == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 0):
            output_dict[f'PRIMER_RIGHT_{i}_SELF_END'] = rev[0].self_end

        if (go_int == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 0):
            output_dict[f'PRIMER_{int_oligo}_{i}_SELF_END'] = intl[0].self_end

        # Print primer self_end thermodynamical approach
        if (go_fwd == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 1):
            output_dict[f'PRIMER_LEFT_{i}_SELF_END_TH'] = fwd[0].self_end
        if (go_rev == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 1):
            output_dict[f'PRIMER_RIGHT_{i}_SELF_END_TH'] = rev[0].self_end
        if (go_int == 1) and ((global_settings_data[0].thermodynamic_oligo_alignment == 1)):
            output_dict[f'PRIMER_{int_oligo}_{i}_SELF_END_TH'] = intl[0].self_end

        # Print primer secondary structures*/
        if (
            (go_fwd == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (fwd[0].self_end_struct != NULL)
        ):
            sqtemp_b = <bytes> fwd[0].self_end_struct
            output_dict[f'PRIMER_LEFT_{i}_SELF_END_STUCT'] = sqtemp_b.decode('utf8')
        if (
            (go_rev == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (rev[0].self_end_struct != NULL)
        ):
            sqtemp_b = <bytes> rev[0].self_end_struct
            output_dict[f'PRIMER_RIGHT_{i}_SELF_END_STUCT'] = sqtemp_b.decode('utf8')
        if (
            (go_int == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (intl[0].self_end_struct != NULL)
        ):
            sqtemp_b = <bytes> intl[0].self_end_struct
            output_dict[f'PRIMER_{int_oligo}_{i}_SELF_END_STUCT'] = sqtemp_b.decode('utf8')

        # Print primer hairpin
        if (go_fwd == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 1):
            output_dict[f'PRIMER_LEFT_{i}_HAIRPIN_TH'] = fwd[0].hairpin_th
        if (go_rev == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 1):
            output_dict[f'PRIMER_RIGHT_{i}_HAIRPIN_TH'] = rev[0].hairpin_th
        if (go_int == 1) and (global_settings_data[0].thermodynamic_oligo_alignment == 1):
            output_dict[f'PRIMER_{int_oligo}_{i}_HAIRPIN_TH'] = intl[0].hairpin_th

        # Print primer secondary structures*/
        if (
            (go_fwd == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (fwd[0].hairpin_struct != NULL)
        ):
            sqtemp_b = <bytes> fwd[0].hairpin_struct
            output_dict[f'PRIMER_LEFT_{i}_HAIRPIN_STUCT'] = sqtemp_b.decode('utf8')

        if (
            (go_rev == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (rev[0].hairpin_struct != NULL)
        ):
            sqtemp_b = <bytes> rev[0].hairpin_struct
            output_dict[f'PRIMER_RIGHT_{i}_HAIRPIN_STUCT'] = sqtemp_b.decode('utf8')

        if (
            (go_int == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (intl[0].hairpin_struct != NULL)
        ):
            sqtemp_b = <bytes> intl[0].hairpin_struct
            output_dict[f'PRIMER_{int_oligo}_{i}_HAIRPIN_STUCT'] = sqtemp_b.decode('utf8')

        # Print out primer mispriming scores
        if seq_lib_num_seq(global_settings_data[0].p_args.repeat_lib) > 0:
            if go_fwd == 1:
                sqtemp_b = <bytes> fwd[0].repeat_sim.name
                output_dict['PRIMER_LEFT_{i}_LIBRARY_MISPRIMING'] = (
                    fwd[0].repeat_sim.score[fwd[0].repeat_sim.max],
                    sqtemp_b.decode('utf8'),
                )

            if go_rev == 1:
                sqtemp_b = <bytes> rev[0].repeat_sim.name
                output_dict['PRIMER_RIGHT_{i}_LIBRARY_MISPRIMING'] = (
                    rev[0].repeat_sim.score[rev[0].repeat_sim.max],
                    sqtemp_b.decode('utf8'),
                )

            if retval[0].output_type == p3_output_type.primer_pairs:
                sqtemp_b = <bytes> retval[0].best_pairs.pairs[i].rep_name
                output_dict['PRIMER_PAIR_{i}_LIBRARY_MISPRIMING'] = (
                    retval[0].best_pairs.pairs[i].repeat_sim,
                    sqtemp_b.decode('utf8'),
                )

        # Print out internal oligo mispriming scores
        if (
            (go_int == 1) and
            (seq_lib_num_seq(global_settings_data[0].o_args.repeat_lib) > 0)
        ):
            sqtemp_b = <bytes> intl[0].repeat_sim.name
            output_dict[f'PRIMER_{int_oligo}_{i}_LIBRARY_MISPRIMING'] = (
                intl[0].repeat_sim.score[intl[0].repeat_sim.max],
                sqtemp_b.decode('utf8'),
            )


        # If a sequence quality was provided, print it*/
        if sequence_args_data[0].quality != NULL:
            if go_fwd == 1:
                output_dict[f'PRIMER_LEFT_{i}_MIN_SEQ_QUALITY'] = fwd[0].seq_quality
            if go_rev == 1:
                output_dict[f'PRIMER_RIGHT_{i}_MIN_SEQ_QUALITY'] = rev[0].seq_quality
            if go_int == 1:
                output_dict[f'PRIMER_{int_oligo}_{i}_MIN_SEQ_QUALITY'] = intl[0].seq_quality

        # Print position penalty, this is for backward compatibility
        if (
            not pr_default_position_penalties(global_settings_data) or
            not PR_START_CODON_POS_IS_NULL(sequence_args_data)
        ):
            output_dict[f'PRIMER_LEFT_{i}_POSITION_PENALTY'] = fwd[0].position_penalty
            output_dict[f'PRIMER_RIGHT_{i}_POSITION_PENALTY'] = rev[0].position_penalty

        # Print primer end stability
        if go_fwd == 1:
            output_dict[f'PRIMER_LEFT_{i}_END_STABILITY'] = fwd[0].end_stability

        if go_rev == 1:
            output_dict[f'PRIMER_RIGHT_{i}_END_STABILITY'] = rev[0].end_stability


        # Print primer template mispriming
        if (
            (global_settings_data[0].thermodynamic_template_alignment == 0) and
            (go_fwd == 1) and
            (oligo_max_template_mispriming(fwd) != ALIGN_SCORE_UNDEF)
        ):
            output_dict[f'PRIMER_LEFT_{i}_TEMPLATE_MISPRIMING'] = \
                oligo_max_template_mispriming(fwd)
        if (
            (global_settings_data[0].thermodynamic_template_alignment == 0) and
            (go_rev == 1) and
            (oligo_max_template_mispriming(rev) != ALIGN_SCORE_UNDEF)
        ):
            output_dict[f'PRIMER_RIGHT_{i}_TEMPLATE_MISPRIMING'] = \
                oligo_max_template_mispriming(rev)


        # Print primer template mispriming, thermodynamical approach*/
        if (
            (global_settings_data[0].thermodynamic_template_alignment == 0) and
            (go_fwd == 1) and
            (oligo_max_template_mispriming(fwd) != ALIGN_SCORE_UNDEF)
        ):
            output_dict[f'PRIMER_LEFT_{i}_TEMPLATE_MISPRIMING_TH'] = \
                oligo_max_template_mispriming_thermod(fwd)

        if (
            (global_settings_data[0].thermodynamic_template_alignment == 0) and
            (go_rev == 1) and
            (oligo_max_template_mispriming(rev) != ALIGN_SCORE_UNDEF)
        ):
            output_dict[f'PRIMER_RIGHT_{i}_TEMPLATE_MISPRIMING_TH'] = \
                oligo_max_template_mispriming_thermod(rev)

        # Print primer secondary structures*/
        if (
            (go_fwd == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (oligo_max_template_mispriming_struct(fwd) != NULL)
        ):
            sqtemp_b = <bytes> oligo_max_template_mispriming_struct(fwd)
            output_dict[f'PRIMER_LEFT_{i}_TEMPLATE_MISPRIMING_STUCT'] = sqtemp_b.decode('utf8')
        if (
            (go_rev == 1) and
            (global_settings_data[0].show_secondary_structure_alignment == 1) and
            (oligo_max_template_mispriming_struct(rev) != NULL)
        ):
            sqtemp_b = <bytes> oligo_max_template_mispriming_struct(rev)
            output_dict[f'PRIMER_RIGHT_{i}_TEMPLATE_MISPRIMING_STUCT'] = sqtemp_b.decode('utf8')

        # Print the pair parameters*/
        if retval[0].output_type == p3_output_type.primer_pairs:
            if (go_int == 1) and (sequence_args_data[0].quality != NULL):
                output_dict[f'PRIMER_{int_oligo}_{i}_MIN_SEQ_QUALITY'] = intl[0].seq_quality
            # Print pair comp_any
            if global_settings_data[0].thermodynamic_oligo_alignment == 0:
                output_dict[f'PRIMER_PAIR_{i}_COMPL_ANY'] = \
                    retval[0].best_pairs.pairs[i].compl_any

            if global_settings_data[0].thermodynamic_oligo_alignment == 1:
                output_dict[f'PRIMER_PAIR_{i}_COMPL_ANY_TH'] = \
                    retval[0].best_pairs.pairs[i].compl_any

            # Print primer secondary structures */
            if (
                (global_settings_data[0].show_secondary_structure_alignment == 1) and
                (retval[0].best_pairs.pairs[i].compl_any_struct != NULL)
            ):
                sqtemp_b = <bytes> retval[0].best_pairs.pairs[i].compl_any_struct
                output_dict[f'PRIMER_PAIR_{i}_COMPL_ANY_STUCT'] = sqtemp_b.decode('utf8')

            # Print pair comp_end
            if global_settings_data[0].thermodynamic_oligo_alignment == 0:
                output_dict[f'PRIMER_PAIR_{i}_COMPL_END'] = \
                    retval[0].best_pairs.pairs[i].compl_end

            if global_settings_data[0].thermodynamic_oligo_alignment == 1:
                output_dict[f'PRIMER_PAIR_{i}_COMPL_END_TH'] = \
                    retval[0].best_pairs.pairs[i].compl_end

            # Print primer secondary structures*/
            if (
                (global_settings_data[0].show_secondary_structure_alignment == 1) and
                (retval[0].best_pairs.pairs[i].compl_end_struct != NULL)
            ):
                sqtemp_b = <bytes> retval[0].best_pairs.pairs[i].compl_end_struct
                output_dict[f'PRIMER_PAIR_{i}_COMPL_END_STUCT'] = sqtemp_b.decode('utf8')

            # Print product size
            product_size = retval[0].best_pairs.pairs[i].product_size
            if sequence_args_data[0].overhang_left != NULL:
                product_size += strlen(sequence_args_data[0].overhang_left)
            if sequence_args_data[0].overhang_right != NULL:
                product_size += strlen(sequence_args_data[0].overhang_right)
            output_dict[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'] = product_size


            # Print the product Tm if a Tm range is defined
            if (
                (global_settings_data[0].product_max_tm != PR_DEFAULT_PRODUCT_MAX_TM) or
                (global_settings_data[0].product_min_tm != PR_DEFAULT_PRODUCT_MIN_TM)
            ):
                output_dict[f'PRIMER_PAIR_{i}_PRODUCT_TM'] = \
                    retval[0].best_pairs.pairs[i].product_tm
                output_dict[f'PRIMER_PAIR_{i}_PRODUCT_TM_OLIGO_TM_DIFF'] = \
                    retval[0].best_pairs.pairs[i].product_tm_oligo_tm_diff
                output_dict[f'PRIMER_PAIR_{i}_T_OPT_A'] = retval[0].best_pairs.pairs[i].t_opt_a
            else:
                output_dict[f'PRIMER_PAIR_{i}_PRODUCT_TM'] = \
                retval[0].best_pairs.pairs[i].product_tm

            # Print the primer pair template mispriming
            if (
                (global_settings_data[0].thermodynamic_template_alignment == 0) and
                (retval[0].best_pairs.pairs[i].template_mispriming != ALIGN_SCORE_UNDEF)
            ):
                output_dict[f'PRIMER_PAIR_{i}_TEMPLATE_MISPRIMING'] = \
                    retval[0].best_pairs.pairs[i].template_mispriming
            # Print the primer pair template mispriming. Thermodynamic approach.
            if (
                    (global_settings_data[0].thermodynamic_template_alignment == 1) and
                    (retval[0].best_pairs.pairs[i].template_mispriming != ALIGN_SCORE_UNDEF)
            ):
                output_dict[f'PRIMER_PAIR_{i}_TEMPLATE_MISPRIMING_TH'] = \
                    retval[0].best_pairs.pairs[i].template_mispriming

            # Print primer secondary structures*/
            if (
                 (global_settings_data[0].show_secondary_structure_alignment == 1) and
                (retval[0].best_pairs.pairs[i].template_mispriming_struct != NULL)
            ):
                sqtemp_b = <bytes> retval[0].best_pairs.pairs[i].template_mispriming_struct
                output_dict[f'PRIMER_PAIR_{i}_TEMPLATE_MISPRIMING_STUCT'] = \
                    sqtemp_b.decode('utf8')
        # End of print parameters of primer pairs
    # End of the for big loop printing all data
    return output_dict


def set_globals_and_seq_args(
        global_args: Dict[str, Any],
        seq_args: Optional[Dict[str, Any]],
        misprime_lib: Optional[Dict[str, Any]] = None,
        mishyb_lib: Optional[Dict[str, Any]] = None,
) -> None:
    '''
    Sets the Primer3 global settings and sequence settings from a Python
    dictionaries containing `key: value` pairs that correspond to the
    documented Primer3 global and sequence argument parameters.
    Also accepts a mispriming or mishybridization library organized as
    `seq_name`:`seq_value` key:value pairs.

    Args:
        seq_args: Primer3 sequence/design args as per Primer3 docs
        global_args: Primer3 global args as per Primer3 docs
        misprime_lib: `Sequence name: sequence` dictionary for mispriming
            checks.
        mishyb_lib: `Sequence name: sequence` dictionary for mishybridization
            checks.

    Raises:
        OSError: Could not allocate memory
    '''
    global global_settings_data
    global sequence_args_data

    cdef:
        seq_lib* mp_lib = NULL
        seq_lib* mh_lib = NULL
        char* arg_input_buffer = NULL


    err_msg = ''

    if sequence_args_data != NULL:
        # Free memory for previous seq args
        destroy_seq_args(sequence_args_data)
        sequence_args_data = NULL

    if seq_args:
        sequence_args_data = create_seq_arg()

        if sequence_args_data == NULL:
            raise OSError('Could not allocate memory for seq_arg')

        global_args.update(seq_args)

    global_arg_bytes = argdefaults.format_boulder_io(global_args)
    arg_input_buffer = global_arg_bytes
    if arg_input_buffer == NULL:
        raise ValueError(global_arg_bytes)

    if global_settings_data != NULL:
        # Free memory for previous global settings
        p3_destroy_global_settings(global_settings_data)
        global_settings_data = NULL

    # Allocate memory for global settings
    global_settings_data = p3_create_global_settings()
    if global_settings_data == NULL:
        raise OSError('Could not allocate memory for p3 globals')

    kmer_lists_path = global_args.get('PRIMER_MASK_KMERLIST_PATH', '')
    if kmer_lists_path:
        local_dir = os.path.dirname(os.path.realpath(get_dunder_file()))
        libprimer3_dir = os.path.join(local_dir, 'src', 'libprimer3')
        if not os.path.isdir(kmer_lists_path):
            if kmer_lists_path[0:2] == '../':
                kmer_lists_path = os.path.join(
                    libprimer3_dir,
                    kmer_lists_path[3:-1],
                )
            else:
                kmer_lists_path = os.path.join(
                    libprimer3_dir,
                    kmer_lists_path,
                )
        if not os.path.isdir(kmer_lists_path):
            raise ValueError(
                f'PRIMER_MASK_KMERLIST_PATH: path {kmer_lists_path} not found'
            )

    try:
        pdh_wrap_set_seq_args_globals(
            global_settings_data,
            sequence_args_data,
            kmer_lists_path,
            arg_input_buffer,
        )
    except BaseException:
        print(f'Issue setting globals. bytes provided: \n\t{global_arg_bytes}')
        p3_destroy_global_settings(global_settings_data)
        global_settings_data = NULL
        if seq_args:
            destroy_seq_args(sequence_args_data)
            sequence_args_data = NULL
        raise

    err_msg = ''
    try:
        if misprime_lib != None:
            mp_lib = pdh_create_seq_lib(misprime_lib)
            if mp_lib == NULL:
                err_msg = f'Issue creating misprime_lib {misprime_lib}'
                raise ValueError(f'Issue creating misprime_lib {misprime_lib}')
            global_settings_data[0].p_args.repeat_lib = mp_lib

        if mishyb_lib != None:
            mh_lib = pdh_create_seq_lib(mishyb_lib)
            if mh_lib == NULL:
                err_msg = f'Issue creating mishyb_lib: {mishyb_lib}'
                raise ValueError(err_msg)
            global_settings_data[0].o_args.repeat_lib = mh_lib
    except (OSError, TypeError) as exc:
        p3_destroy_global_settings(global_settings_data)
        global_settings_data = NULL
        destroy_seq_args(sequence_args_data)
        sequence_args_data = NULL
        raise OSError(err_msg) from exc


def run_design() -> None:
    '''
    Wraps the primer design functionality of Primer3. Should be called
    after setting the global and sequence-specific Primer3 parameters
    (see setGlobals and setSeqArgs, above)
    '''
    global global_settings_data
    global sequence_args_data

    cdef:
        p3retval* retval = NULL

    results_dict: dict = {}

    if global_settings_data == NULL or sequence_args_data == NULL:
        raise ValueError(
            'Primer3 global args and sequence args must be set prior '
            'to calling runDesign'
        )

    load_thermo_params()

    retval = choose_primers(
        global_settings_data,
        sequence_args_data,
    )
    if retval == NULL:
        raise ValueError('Issue choosing primers')
    try:
        results_dict = pdh_design_output_to_dict(
            global_settings_data,
            sequence_args_data,
            retval,
        )
    finally:
        destroy_secundary_structures(
            global_settings_data,
            retval,
        )
        destroy_p3retval(retval)
        retval = NULL
        destroy_dpal_thal_arg_holder()
    return results_dict

def design_clean_up():
    # Free any remaining global Primer3 objects
    global global_settings_data
    global sequence_args_data

    destroy_thal_structures()
    if global_settings_data != NULL:
        # Free memory for previous global settings
        p3_destroy_global_settings(global_settings_data)
        global_settings_data = NULL

    if sequence_args_data != NULL:
        # Free memory for previous seq args
        destroy_seq_args(sequence_args_data)
        sequence_args_data = NULL

atexit.register(design_clean_up)
