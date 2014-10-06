'''
thermoanalysis.pxd
~~~~~~~~~~~~~~~~~~

Cython header file for thermoanalysis.pyx -- allows for cross-project Cython / 
C integration of the low-level thermodynamic analysis bindings.

'''

cdef extern from "thal.h":
    ctypedef enum thal_alignment_type:
        thal_any = 1
        thal_end1 = 2
        thal_end2 = 3
        thal_hairpin = 4

    ctypedef struct thal_args:
        int debug                 # if non zero, print debugging info to stderr
        thal_alignment_type type  # type of thermodynamic alignment
        int maxLoop               # maximum size of loop to consider in calcs
        double mv                 # [ ] of monovalent cations (mM)
        double dv                 # [ ] of divalent cations (mM)
        double dntp               # [ ] of dNTPs (mM)
        double dna_conc           # [ ] of oligos (nM)
        double temp               # temp at which hairpins will be calculated
        int temponly              # print only temp to stderr
        int dimer                 # if non-zero dimer structure is calculated
    
    ctypedef struct thal_results:
        char msg[255]
        int no_structure # Added no structure (1 if no structure found)
        double temp
        double ds # Added entropy value
        double dh # Added enthalpy value
        double dg # Added gibbs free energy value
        int align_end_1
        int align_end_2



    void thal(  const unsigned char*,
                const unsigned char*,
                const thal_args*,
                thal_results*,
                const int)

    int get_thermodynamic_values(const char*, thal_results *)


cdef class ThermoResult:
    cdef thal_results thalres

cdef class ThermoAnalysis:
    cdef thal_args thalargs
    cdef public int max_nn_length
    cdef public int _tm_method
    cdef public int _salt_correction_method

    cdef inline ThermoResult calcHeterodimer_c(ThermoAnalysis self,
                                               unsigned char*s1,
                                               unsigned char* s2)

    cdef inline ThermoResult calcHomodimer_c(ThermoAnalysis self, unsigned char*s1)

    cdef inline ThermoResult calcHairpin_c(ThermoAnalysis self, unsigned char*s1)

    cdef inline ThermoResult calcEndStability_c(ThermoAnalysis self,
                                                unsigned char*s1,
                                                unsigned char* s2)

    cdef inline double calcTm_c(ThermoAnalysis self, char* s1)
