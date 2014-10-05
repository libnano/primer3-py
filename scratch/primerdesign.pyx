
cdef extern from "libprimer3.h":
    ctypedef struct p3_global_settings:
        task   primer_task
        int    pick_left_primer
        int    pick_right_primer
        int    pick_internal_oligo
        int    file_flag
        int    first_base_index  
        int    liberal_base  
        int    num_return 
        int    pick_anyway  
        int    lib_ambiguity_codes_consensus
        int    quality_range_min
        int    quality_range_max
        args_for_one_oligo_or_primer p_args
        args_for_one_oligo_or_primer o_args
        tm_method_type tm_santalucia
        salt_correction_type salt_corrections
        double max_end_stability
        int    max_end_gc
        int    gc_clamp 
        int lowercase_masking
        sequencing_parameters sequencing 
        double outside_penalty 
        double inside_penalty
        int    pr_min[PR_MAX_INTERVAL_ARRAY]
        int    pr_max[PR_MAX_INTERVAL_ARRAY]
        int    num_intervals      
        int    product_opt_size
        double product_max_tm
        double product_min_tm
        double product_opt_tm
        double pair_max_template_mispriming
        double pair_max_template_mispriming_th
        double pair_repeat_compl
        double pair_compl_any
        double pair_compl_any_th
        double pair_compl_end
        double pair_compl_end_th
        int thermodynamic_oligo_alignment
        int thermodynamic_template_alignment
        double max_diff_tm
        pair_weights  pr_pair_weights
        int    min_left_three_prime_distance
        int    min_right_three_prime_distance
        int    min_5_prime_overlap_of_junction   
        int    min_3_prime_overlap_of_junction
        int dump   

cdef extern from "primerdesign_helpers.h":

    p3_global_settings* \
    pdh_setGlobals(p3_global_settings *pa, PyObject *p3_settings_dict)

    seq_lib* \
    pdh_createSeqLib(PyObject *seq_dict)

    seq_args* \
    pdh_setSeqArgs(PyObject *sa_dict, p3_global_settings *pa)

    PyObject* \
    pdh_outputToDict(const p3_global_settings *pa, const seq_args *sa,
                     const p3retval *retval)


cdef class PrimerDesign:
    ''' Provides an entry point for the Primer3 primer design process.
    '''

    cdef p3_global_settings     *p3_global_args
    cdef seq_args               *p3_seq_args

    def __cinit__(self, global_args, seq_args=None, misprime_lib=None,
                  mishyb_lib=None):
        pass

    def setGlobals(self, global_args, misprime_lib=None, mishyb_lib=None):
        if self.p3_global_args != NULL:
            p3_destroy_global_settings(self.p3_global_args)
            self.p3_global_args = p3_create_global_settings()
            if self.p3_global_args == NULL:
                raise IOError('Could not allocate memory for '
                              'p3_global_settings')

        if (pdh_setGlobals(self.p3_global_args, global_args)) == NULL:
            raise

        if misprime_lib:
            self.p3_global_args.p_args.repeat_lib = pdh_createSeqLib(misprime_lib)
            if self.p3_global_args.p_args.repeat_lib == NULL:
                raise

        if mishyb_lib:
            self.p3_global_args.p_args.repeat_lib = pdh_createSeqLib(misprime_lib)
            if self.p3_global_args.p_args.repeat_lib == NULL:
                raise            



