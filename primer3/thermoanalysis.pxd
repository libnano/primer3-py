# Copyright (C) 2014-2020. Ben Pruitt & Nick Conway; Wyss Institute
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
        double dmso_conc          # [ ] of DMSO (%)
        double dmso_fact          # [ ] of DMSO (correction factor, default 0.6)
        double formamide_conc     # [ ] of formamide (mol/l) */
        double dna_conc           # [ ] of oligos (nM)
        double annealing_temp     # [ ] Actual annealing temperature of the PCR reaction
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
        char *sec_struct
    
    ctypedef struct thal_parameters:
        char *dangle_dh
        char *dangle_ds
        char *loops_dh
        char *loops_ds
        char *stack_dh
        char *stack_ds
        char *stackmm_dh
        char *stackmm_ds
        char *tetraloop_dh
        char *tetraloop_ds
        char *triloop_dh
        char *triloop_ds
        char *tstack_tm_inf_ds
        char *tstack_dh
        char *tstack2_dh
        char *tstack2_ds

    ctypedef enum thal_mode:
      THL_FAST    = 0
      THL_GENERAL = 1
      THL_DEBUG_F = 2
      THL_DEBUG   = 3
      THL_STRUCT  = 4



    void thal(  const unsigned char*,
                const unsigned char*,
                const thal_args*,
                const thal_mode mode, 
                thal_results*)

    
    int thal_load_parameters(const char*, thal_parameters *, thal_results *);
    
    int thal_set_null_parameters(thal_parameters *);
    
    int get_thermodynamic_values(const thal_parameters *, thal_results *)

    void destroy_thal_structures()


cdef class ThermoResult:
    cdef thal_results thalres
    cdef public object ascii_structure

cdef class ThermoAnalysis:
    cdef thal_args thalargs
    cdef public int max_nn_length
    cdef public int _tm_method
    cdef public int _salt_correction_method

    cdef inline ThermoResult calcHeterodimer_c(
        ThermoAnalysis self,
        unsigned char*s1,
        unsigned char* s2,
        bint output_structure
    )

    cdef inline ThermoResult calcHomodimer_c(
        ThermoAnalysis self,
        unsigned char*s1,
        bint output_structure
    )

    cdef inline ThermoResult calcHairpin_c(
        ThermoAnalysis self,
        unsigned char*s1,
        bint output_structure
    )

    cdef inline ThermoResult calcEndStability_c(ThermoAnalysis self,
                                                unsigned char*s1,
                                                unsigned char* s2)

    cdef inline double calcTm_c(ThermoAnalysis self, char* s1)

    cpdef calcHeterodimer(ThermoAnalysis self, seq1, seq2, output_structure=*)

    cpdef calcHomodimer(ThermoAnalysis self, seq1, output_structure=*)

    cpdef calcHairpin(ThermoAnalysis self, seq1, output_structure=*)

    cpdef misprimingCheck(ThermoAnalysis self, putative_seq, sequences,  double tm_threshold)
