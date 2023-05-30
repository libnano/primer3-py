/* libprimer3/thalflex.c

Source original thal.c

 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2009,2010,
               2011,2012,2016
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
 All rights reserved.

       This file is part of primer3 software suite.

       This software suite is is free software;
       you can redistribute it and/or modify it under the terms
       of the GNU General Public License as published by the Free
       Software Foundation; either version 2 of the License, or (at
       your option) any later version.

       This software is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this software (file gpl-2.0.txt in the source
       distribution); if not, write to the Free Software
       Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Modifications to thal.c to create thalflex.c
Copyright (c) 2023. Ben Pruitt & Nick Conway;

To enable non-global mutable variable data, thal.c has been modified and split
into 3 files, thalflex.c, thalflex.h, and thalflexsignatures.h to enhance code
readability

> thalflex.c
1. Top level macros
2. Core routines for thal() execution

> thalflex.h
1. C compiler support macros for thal() execution
2. Data structures
3. Support functions

> thalflexsignatures.h
1. Data structures
2. Function signatures to thalflex.c

*/
#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <ctype.h>
#include <math.h>

/* MacOS specific for memory investigations commented out for now*/
// #include <malloc/malloc.h>

/* #include <unistd.h> primer3-py COMMENTED OUT as is unused and breaks Windows build */

#if defined(__sun)
#include <ieeefp.h>
#endif

/* Check on which OS we compile */
#if defined(_WIN32) || defined(WIN32) || defined (__WIN32__) || defined(__CYGWIN__) || defined(__MINGW32__)
#define OS_WIN
#endif

/*#define DEBUG*/
#define THAL_EXIT_ON_ERROR 0
#define CHECK_ERROR(COND,MSG) if (COND) { strcpy(o->msg, MSG); errno = 0; longjmp(_jmp_buf, 1); }
#define THAL_OOM_ERROR { strcpy(o->msg, "Out of memory"); errno = ENOMEM; longjmp(_jmp_buf, 1); }
#define THAL_IO_ERROR(f) { snprintf(o->msg, 255, "Unable to open file %s", f); longjmp(_jmp_buf, 1); }

#define STR(X) #X
#define LONG_SEQ_ERR_STR(MAX_LEN) "Target sequence length > maximum allowed (" STR(MAX_LEN) ") in thermodynamic alignment"
#define XSTR(X) STR(X)

#define SMALL_NON_ZERO 0.000001
#define DBL_EQ(X,Y) (((X) - (Y)) < (SMALL_NON_ZERO) ? (1) : (2)) /* 1 when numbers are equal */

#ifdef INTEGER
#define isFinite(x) (x < _INFINITY / 2)
#else
#define isFinite(x) isfinite(x)
#endif

#define isPositive(x) ((x) > 0 ? (1) : (0))

/*** BEGIN CONSTANTS ***/
#define MIN_HRPN_LOOP 3 /*  minimum size of hairpin loop */


#ifdef INTEGER
const double _INFINITY = 999999.0;
#else
#ifdef INFINITY
const double _INFINITY = INFINITY;
#else
const double _INFINITY = 1.0 / 0.0;
#endif
#endif

static const double R = 1.9872; /* cal/Kmol */
static const double ILAS = (-300 / 310.15); /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
static const double ILAH = 0.0; /* Internal Loop EntHalpy Asymmetry correction */
static const double AT_H = 2200.0; /* AT penalty */
static const double AT_S = 6.9; /* AT penalty */
static const double MinEntropyCutoff = -2500.0; /* to filter out non-existing entropies */
static const double MinEntropy = -3224.0; /* initiation */
static const double G2 = 0.0; /* structures w higher G are considered to be unstabile */
static const double ABSOLUTE_ZERO = 273.15;
static const double TEMP_KELVIN = 310.15;
static const int MAX_LOOP = 30; /* the maximum size of loop that can be calculated; for larger loops formula must be implemented */
// static const int MIN_LOOP = 0;
// static const char BASES[5] = {'A', 'C', 'G', 'T', 'N'}; /* bases to be considered - N is every symbol that is not A, G, C,$
// static const char BASE_PAIRS[4][4] = {"A-T", "C-G", "G-C", "T-A" }; /* allowed basepairs */
/* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
static const int BPI[5][5] =  {
    {0, 0, 0, 1, 0}, /* A, C, G, T, N; */
    {0, 0, 1, 0, 0},
    {0, 1, 0, 0, 0},
    {1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0}
};

/*** END OF CONSTANTS ***/

/*** INCLUDE FUNCTION SIGNATURES ***/
#include "thalflexsignatures.h"

/**********************************************************************/
/*             Begin definitions for thermodynamic data               */
/* This data remains constant once configuration files data is loaded */
/**********************************************************************/
/* hairpin triloop penalties */
static int numTriloops;

/* hairpin tetraloop penalties */
static int numTetraloops;

/* AT penalty */
static double atpS[5][5];
static double atpH[5][5];

/* thermodynamic parameters for 3' dangling ends */
static double dangleEntropies3[5][5][5];
static double dangleEnthalpies3[5][5][5]; /* ther params for 3' dangling ends */

/* thermodynamic parameters for 5' dangling ends */
static double dangleEntropies5[5][5][5];
static double dangleEnthalpies5[5][5][5];

/* thermodynamic parameters for perfect match pairs */
static double stackEntropies[5][5][5][5];
static double stackEnthalpies[5][5][5][5]; /* ther params for perfect match pairs */

/* thermodynamic parameters for perfect match and internal mm */
static double stackint2Entropies[5][5][5][5];
static double stackint2Enthalpies[5][5][5][5];

/* interior loop params according to length of the loop */
static double interiorLoopEntropies[30];
static double interiorLoopEnthalpies[30];

/* bulge loop params according to length of the loop */
static double bulgeLoopEntropies[30];
static double bulgeLoopEnthalpies[30];

/* hairpin loop params accordint to length of the loop */
static double hairpinLoopEntropies[30];
static double hairpinLoopEnthalpies[30];

/* ther params for terminal mismatches */
static double tstackEntropies[5][5][5][5];
static double tstackEnthalpies[5][5][5][5];

/* ther params for internal terminal mismatches */
static double tstack2Entropies[5][5][5][5];
static double tstack2Enthalpies[5][5][5][5];

/* ther penalties for given triloop seq-s */
// Inspection show 256 bytes allocated, so that corresponds to 19 elements of
// 5 chars + 1 double (13 bytes) = 247
static triloop_t triloopEntropies[20];
static triloop_t triloopEnthalpies[20];

/* ther penalties for given tetraloop seq-s */
// Inspection show 2048 bytes allocated, so that corresponds to 146 elements of
// 6 chars + 1 double (14 bytes) = 2044
static tetraloop_t tetraloopEntropies[146];
static tetraloop_t tetraloopEnthalpies[146];

/**********************************************************************/
/*             End definitions for thermodynamic data                 */
/**********************************************************************/

static jmp_buf _jmp_buf;


/*** INCLUDE THAL FLEX HEADER CODE ***/
#include "thalflex.h"

/* Initialize the thermodynamic values (parameters) */
int thal_set_null_parameters(thal_parameters *a) {
  a->dangle_dh = NULL;
  a->dangle_ds = NULL;
  a->loops_dh = NULL;
  a->loops_ds = NULL;
  a->stack_dh = NULL;
  a->stack_ds = NULL;
  a->stackmm_dh = NULL;
  a->stackmm_ds = NULL;
  a->tetraloop_dh = NULL;
  a->tetraloop_ds = NULL;
  a->triloop_dh = NULL;
  a->triloop_ds = NULL;
  a->tstack_tm_inf_ds = NULL;
  a->tstack_dh = NULL;
  a->tstack2_dh = NULL;
  a->tstack2_ds = NULL;
  return 0;
}

/* Free the thermodynamic values (parameters) */
int thal_free_parameters(thal_parameters *a) {
  if (NULL != a->dangle_dh) {
    free(a->dangle_dh);
    a->dangle_dh = NULL;
  }
  if (NULL != a->dangle_ds) {
    free(a->dangle_ds);
    a->dangle_ds = NULL;
  }
  if (NULL != a->loops_dh) {
    free(a->loops_dh);
    a->loops_dh = NULL;
  }
  if (NULL != a->loops_ds) {
    free(a->loops_ds);
    a->loops_ds = NULL;
  }
  if (NULL != a->stack_dh) {
    free(a->stack_dh);
    a->stack_dh = NULL;
  }
  if (NULL != a->stack_ds) {
    free(a->stack_ds);
    a->stack_ds = NULL;
  }
  if (NULL != a->stackmm_dh) {
    free(a->stackmm_dh);
    a->stackmm_dh = NULL;
  }
  if (NULL != a->stackmm_ds) {
    free(a->stackmm_ds);
    a->stackmm_ds = NULL;
  }
  if (NULL != a->tetraloop_dh) {
    free(a->tetraloop_dh);
    a->tetraloop_dh = NULL;
  }
  if (NULL != a->tetraloop_ds) {
    free(a->tetraloop_ds);
    a->tetraloop_ds = NULL;
  }
  if (NULL != a->triloop_dh) {
    free(a->triloop_dh);
    a->triloop_dh = NULL;
  }
  if (NULL != a->triloop_ds) {
    free(a->triloop_ds);
    a->triloop_ds = NULL;
  }
  if (NULL != a->tstack_tm_inf_ds) {
    free(a->tstack_tm_inf_ds);
    a->tstack_tm_inf_ds = NULL;
  }
  if (NULL != a->tstack_dh) {
    free(a->tstack_dh);
    a->tstack_dh = NULL;
  }
  if (NULL != a->tstack2_dh) {
    free(a->tstack2_dh);
    a->tstack2_dh = NULL;
  }
  if (NULL != a->tstack2_ds) {
    free(a->tstack2_ds);
    a->tstack2_ds = NULL;
  }
  return 0;
}

/* Read the thermodynamic values (parameters) from the parameter files
   in the directory specified by 'path'.  Return 0 on success and -1
   on error. The thermodynamic values are stored in multiple static
   variables. */
int
get_thermodynamic_values(const thal_parameters *tp, thal_results *o)
{

  if (setjmp(_jmp_buf) != 0) {
     return -1;
  }

  getStack(stackEntropies, stackEnthalpies, tp, o);
  getStackint2(stackint2Entropies, stackint2Enthalpies, tp, o);
  getDangle(
    dangleEntropies3,
    dangleEnthalpies3,
    dangleEntropies5,
    dangleEnthalpies5,
    tp,
    o
  );
  getLoop(
    hairpinLoopEntropies,
    interiorLoopEntropies,
    bulgeLoopEntropies,
    hairpinLoopEnthalpies,
    interiorLoopEnthalpies,
    bulgeLoopEnthalpies,
    tp,
    o
  );
  getTstack(tstackEntropies, tstackEnthalpies, tp, o);
  getTstack2(tstack2Entropies, tstack2Enthalpies, tp, o);
  getTriloop(
    triloopEntropies,
    triloopEnthalpies,
    &numTriloops,
    tp,
    o
  );
  getTetraloop(
    tetraloopEntropies,
    tetraloopEnthalpies,
    &numTetraloops,
    tp,
    o
  );
  /* getting the AT-penalties */
  tableStartATS(AT_S, atpS);
  tableStartATH(AT_H, atpH);

  // p3_read_line((FILE*) NULL, (thal_results*) NULL); // free static memory in function
  return 0;
}

void
destroy_thal_structures(void) {
  /* Empty fucntion since these data stuctures are no longer malloc'd */
}

/* Central method: execute all sub-methods for calculating secondary
   structure for dimer or for monomer */
void
thal(
    const unsigned char *oligo_f,
    const unsigned char *oligo_r,
    const thal_args *a,
    const thal_mode mode,
    thal_results *o,
    const int print_output  /* primer3-py modification argument */
) {
  double* SH = NULL;
  int i, j;
  int len_f, len_r;
  // double T1; primer3-py commented out variable
  int k;
  int *bp = NULL;
  unsigned char *oligo2_rev = NULL;
  double mh, ms;
  double G1, bestG;

  /***********************************************/
  /* Begin movement of globals to function local */
  /***********************************************/
  double* send5 = NULL;   /* calc 5' */
  double* hend5 = NULL;   /* calc 5' */
  double dplx_init_H;     /* initiation enthalpy; for duplex 200, for unimolecular structure 0 */
  double dplx_init_S;     /* initiation entropy; for duplex -5.7, for unimoleculat structure 0 */

  /* value calculated by saltCorrectS, includes correction for monovalent and divalent cations */
  double saltCorrection;

  double RC;              /* universal gas constant multiplied w DNA conc - for melting temperature */
  double SHleft = 0.0;          /* var that helps to find str w highest melting temperature */

  int bestI, bestJ;       /* starting position of most stable str */

  double* enthalpyDPT = NULL;   /* matrix for values of enthalpy */
  double* entropyDPT = NULL;    /* matrix for values of entropy */

  unsigned char* oligo1 = NULL;   /* inserted oligo sequenced */
  unsigned char* oligo2 = NULL;   /* inserted oligo sequenced */
  unsigned char* numSeq1 = NULL;
  unsigned char* numSeq2 = NULL;  /* same as oligo1 and oligo2 but converted to numbers */
  int len1, len2, len3;           /* length of sequense 1 and 2 */
  /***********************************************/
  /* End movement of globals to function local */
  /***********************************************/

  strcpy(o->msg, "");
  o->temp = THAL_ERROR_SCORE;
  errno = 0;

  if (setjmp(_jmp_buf) != 0) {
    o->temp = THAL_ERROR_SCORE;
    return;  /* If we get here, that means we returned via a
                longjmp.  In this case errno might be ENOMEM,
                but not necessarily. */
  }

  CHECK_ERROR(NULL == oligo_f, "NULL first sequence");
  CHECK_ERROR(NULL == oligo_r, "NULL second sequence");
  len_f = length_unsig_char(oligo_f);
  len_r = length_unsig_char(oligo_r);

  /* The following error messages will be seen by end users and will
     not be easy to understand. */
  CHECK_ERROR((len_f > THAL_MAX_ALIGN) && (len_r > THAL_MAX_ALIGN),
               "At least one sequence must be equal to or shorter than "
               XSTR(THAL_MAX_ALIGN) "bp for thermodynamic calculations");
  CHECK_ERROR((len_f > THAL_MAX_SEQ),
              LONG_SEQ_ERR_STR(THAL_MAX_SEQ) " (1)");
  CHECK_ERROR((len_r > THAL_MAX_SEQ),
              LONG_SEQ_ERR_STR(THAL_MAX_SEQ) " (2)");

  CHECK_ERROR(NULL == a,  "NULL 'in' pointer");
  if (NULL == o) { return; } /* Leave it to the caller to crash */
  CHECK_ERROR(
    (a->type != thal_any) &&
    (a->type != thal_end1) &&
    (a->type != thal_end2) &&
    (a->type != thal_hairpin),
    "Illegal type"
  );
  o->align_end_1 = -1;
  o->align_end_2 = -1;
  if (oligo_f != NULL && '\0' == *oligo_f) { /* primer3-py bug fix for type */
    strcpy(o->msg, "Empty first sequence");
    o->temp = 0.0;
    return;
  }
  if (oligo_r != NULL && '\0' == *oligo_r) { /* primer3-py bug fix for type */
    strcpy(o->msg, "Empty second sequence");
    o->temp = 0.0;
    return;
  }
  if (0 == len_f) {
    o->temp = 0.0;
    return;
  }
  if (0 == len_r) {
    o->temp = 0.0;
    return;
  }
  if (a->type != 3) {
    oligo1 = (unsigned char*) safe_malloc(
      (len_f + 1) * sizeof(unsigned char),
      o
    );
    oligo2 = (unsigned char*) safe_malloc(
      (len_r + 1) * sizeof(unsigned char),
      o
    );
    strcpy((char*) oligo1, (const char*) oligo_f);
    strcpy((char*) oligo2, (const char*) oligo_r);
  } else {
    oligo1 = (unsigned char*) safe_malloc(
      (len_r + 1) * sizeof(unsigned char),
      o
    );
    oligo2 = (unsigned char*) safe_malloc(
      (len_f + 1) * sizeof(unsigned char),
      o
    );
    strcpy((char*) oligo1, (const char*) oligo_r);
    strcpy((char*) oligo2, (const char*) oligo_f);
  }
  /*** INIT values for unimolecular and bimolecular structures ***/
  if (a->type == 4) { /* unimolecular folding */
    len2 = length_unsig_char(oligo2);
    len3 = len2 -1;
    dplx_init_H = 0.0;
    dplx_init_S = -0.00000000001;
    RC = 0;
  } else if (a->type != 4) {
    /* hybridization of two oligos */
    dplx_init_H = 200;
    dplx_init_S = -5.7;
    if(symmetry_thermo(oligo1) && symmetry_thermo(oligo2)) {
      RC = R  * log(a->dna_conc / 1000000000.0);
    } else {
      RC = R  * log(a->dna_conc / 4000000000.0);
    }
    if(a->type != 3) {
      oligo2_rev = (unsigned char*) safe_malloc(
        (length_unsig_char(oligo_r) + 1) * sizeof(unsigned char),
        o
      );
      strcpy((char*) oligo2_rev, (const char*) oligo_r);
    } else {
      oligo2_rev = (unsigned char*) safe_malloc(
        (length_unsig_char(oligo_f) + 1) * sizeof(unsigned char),
        o
      );
      strcpy((char*) oligo2_rev, (const char*) oligo_f);
    }
    /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
    reverse(oligo2_rev);
    free(oligo2);
    oligo2 = NULL;
    oligo2 = &oligo2_rev[0];
  } else {
    strcpy(o->msg, "Wrong alignment type!");
    o->temp = THAL_ERROR_SCORE;
    errno = 0;
#ifdef DEBUG
    fprintf(stderr, o->msg);
#endif
    return;
  }
  len1 = length_unsig_char(oligo1);
  len2 = length_unsig_char(oligo2);

  /* Convert nucleotides to numbers */
  numSeq1 = (unsigned char*) safe_realloc(numSeq1, len1 + 2, o);
  numSeq2 = (unsigned char*) safe_realloc(numSeq2, len2 + 2, o);

  /*** Calc part of the salt correction ***/
  /* salt correction for entropy, must be multiplied with N, which is
  * the total number of phosphates in the duplex divided by 2; 8bp dplx N=7
  */
  saltCorrection = saltCorrectS(a->mv, a->dv, a->dntp);

  if(a->type == 4) { /* monomer */
    /* terminal basepairs */
    send5 = (double*) safe_realloc(send5, (len1 + 1) * sizeof(double), o);
    hend5 = (double*) safe_realloc(hend5, (len1 + 1) * sizeof(double), o);
  }
  for(i = 0; i < len1; i++) { oligo1[i] = toupper(oligo1[i]); }
  for(i = 0; i < len2; i++) { oligo2[i] = toupper(oligo2[i]); }
  for(i = 1; i <= len1; ++i) { numSeq1[i] = str2int(oligo1[i - 1]); }
  for(i = 1; i <= len2; ++i) { numSeq2[i] = str2int(oligo2[i - 1]); }

  /* mark as N-s */
  numSeq1[0] = numSeq1[len1 + 1] = numSeq2[0] = numSeq2[len2 + 1] = 4;

  if (a->type == 4) { /* calculate structure of monomer */
    enthalpyDPT = safe_recalloc(enthalpyDPT, len1, len2, o);
    entropyDPT = safe_recalloc(entropyDPT, len1, len2, o);
    initMatrix2(
      numSeq1, len1,
      numSeq2, len2,
      len3,
      entropyDPT,
      enthalpyDPT
    );
    fillMatrix2(
      a->maxLoop, o,
      numSeq1, len1,
      numSeq2, len2,
      len3,
      entropyDPT,
      enthalpyDPT,
      dplx_init_S,
      dplx_init_H,
      RC
    );
    calc_terminal_bp(
      a->temp,
      send5, hend5,
      entropyDPT, enthalpyDPT,
      numSeq1, len1,
      numSeq2, len2,
      len3,
      dplx_init_S, dplx_init_H, RC
    );
    mh = hend5[len1];
    ms = send5[len1];
    o->align_end_1 = (int) mh;
    o->align_end_2 = (int) ms;
    bp = (int*) safe_calloc(len1, sizeof(int), o);
    for (k = 0; k < len1; ++k) { bp[k] = 0; }
    if(isFinite(mh)) {
      tracebacku(
        bp,
        a->maxLoop,
        o,
        send5, hend5,
        entropyDPT, enthalpyDPT,
        numSeq1, len1,
        numSeq2, len2,
        len3,
        dplx_init_S, dplx_init_H, RC
      );
      /* traceback for unimolecular structure */
      int do_temponly = 0;
      if (mode == THL_FAST) {
        do_temponly = 1;
      }
      if (print_output == 0) {
        calcHairpin(
          bp, mh, ms, do_temponly, a->temp, o,
          len1,
          saltCorrection
        );
      } else {
        /* NOTE: primer3-py This block is enable code running in
         * `libprimer3.choose_primers` to allocate memory for running.  It
         * is that codes resposibility to free the this memory
         * `ThermoAnalysis` manages its own allocations
         */
        if ((o->sec_struct == NULL) && (mode == THL_STRUCT))  {
          int seq_struct_str_len = (len1 + len1 + 1) * 4 + 64;
          char* seq_struct_ptr = (char*) malloc(seq_struct_str_len * sizeof(char));
          if (seq_struct_ptr == NULL) {
            #ifdef DEBUG
              fputs("Error in malloc()\n", stderr);
            #endif
              THAL_OOM_ERROR;
          }
          o->sec_struct = seq_struct_ptr;
          seq_struct_ptr[0] = '\0';
        }
        drawHairpin(
          bp, mh, ms,
          do_temponly, a->temp, o,
          oligo1, len1, saltCorrection
        );
      }

    } else if((mode != THL_FAST) && (mode != THL_DEBUG_F) && (mode != THL_STRUCT)) {
        if (print_output == 1) { /* primer3-py update to supress undesired printing */
          fputs("No secondary structure could be calculated\n", stderr);
        }
        o->no_structure = 1;
    }

    if(o->temp == -_INFINITY && (!strcmp(o->msg, ""))) { o->temp=0.0; }
    free(bp);
    free(enthalpyDPT);
    free(entropyDPT);
    free(numSeq1);
    free(numSeq2);
    free(send5);
    free(hend5);
    free(oligo1);
    free(oligo2);
    return;
  } else if (a->type != 4) { /* Hybridization of two molecules LN 545 */
    int *ps1, *ps2;
    double dH, dS;
    len3 = len2;

    /* dyn. programming table for dS and dH */
    enthalpyDPT = safe_recalloc(enthalpyDPT, len1, len2, o);

    /* enthalpyDPT is 3D array represented as 1D array */
    entropyDPT = safe_recalloc(entropyDPT, len1, len2, o);
    initMatrix(
        numSeq1, len1,
        numSeq2, len2,
        len3,
        entropyDPT,
        enthalpyDPT
    );
    fillMatrix(
      a->maxLoop, o,
      numSeq1, len1,
      numSeq2, len2,
      len3,
      entropyDPT,
      enthalpyDPT,
      dplx_init_S,
      dplx_init_H,
      RC
    );
    /* SHleft = -_INFINITY; primer3-py commented out due to unused */
    SH = (double*) safe_malloc(2 * sizeof(double), o);
    /* calculate terminal basepairs */
    bestI = bestJ = 0;
    G1 = bestG = _INFINITY;
    if (a->type == 1) {
      for (i = 1; i <= len1; i++) {
        for (j = 1; j <= len2; j++) {
          RSH(i, j, SH);

          /* this adding is done for compiler, optimization -O2 vs -O0 */
          SH[0] = SH[0] + SMALL_NON_ZERO;
          SH[1] = SH[1] + SMALL_NON_ZERO;

          G1 = (
            (EnthalpyDPT(i, j) + SH[1] + dplx_init_H) -
            TEMP_KELVIN * (EntropyDPT(i, j) + SH[0] + dplx_init_S)
          );
          if( G1 < bestG) {
            bestG = G1;
            bestI = i;
            bestJ = j;
          }
        }
      }
    }
    ps1 = (int*) safe_calloc(len1, sizeof(int), o);
    ps2 = (int*) safe_calloc(len2, sizeof(int), o);
    for (i = 0; i < len1; ++i) {
      ps1[i] = 0;
    }
    for (j = 0; j < len2; ++j) {
      ps2[j] = 0;
    }
    if (a->type == 2 || a->type == 3) {
      /* THAL_END1 */
      bestI = bestJ = 0;
      bestI = len1;
      i = len1;
      /* SHleft = 0.0 = -_INFINITY; primer3-py commented out due to unused */
      G1 = bestG = _INFINITY;
      for (j = 1; j <= len2; ++j) {
        RSH(i, j, SH);

        /* this adding is done for compiler, optimization -O2 vs -O0
        * that compiler could understand that SH is changed in this cycle
        */
        SH[0] = SH[0] + SMALL_NON_ZERO;

        SH[1] = SH[1] + SMALL_NON_ZERO;
        G1 = (
          (EnthalpyDPT(i, j)+ SH[1] + dplx_init_H) -
          TEMP_KELVIN*(EntropyDPT(i, j) + SH[0] + dplx_init_S)
        );
        if (G1 < bestG) {
          bestG = G1;
          bestJ = j;
        }
      }
    }
    /* if (!isFinite(SHleft = 0.0)) { bestI = bestJ = 1; } primer3-py commented out due to unused DEL */
    if (!isFinite(bestG)) { bestI = bestJ = 1; } // LN 594
    RSH(bestI, bestJ, SH);
    dH = EnthalpyDPT(bestI, bestJ)+ SH[1] + dplx_init_H;
    dS = (EntropyDPT(bestI, bestJ) + SH[0] + dplx_init_S);
    /* tracebacking */
    for (i = 0; i < len1; ++i) {
      ps1[i] = 0;
    }
    for (j = 0; j < len2; ++j) {
      ps2[j] = 0;
    }

    if (isFinite(EnthalpyDPT(bestI, bestJ))) {
      traceback(
        bestI, bestJ,
        ps1, ps2, a->maxLoop, o,
        entropyDPT, enthalpyDPT,
        numSeq1, len1, numSeq2, len2, len3,
        dplx_init_S, dplx_init_H, RC
      );
      // new
      int do_temponly = 0;
      if (mode == THL_FAST) {
        do_temponly = 1;
      }
      if (print_output == 0) {
        calcDimer(
          ps1, ps2,
          SHleft = 0.0, dH, dS,
          do_temponly, a->temp, o,
          len1, len2, saltCorrection, RC
        );
      } else {
        /* NOTE: primer3-py This block is enable code running in
         * `libprimer3.choose_primers` to allocate memory for running.  It
         * is that codes resposibility to free the this memory
         * Thermoanalysis manages its own allocations
         */
        if ((o->sec_struct == NULL) && (mode == THL_STRUCT)) {
          int seq_struct_str_len = (len1 + len2 + 1) * 4 + 64;
          char* seq_struct_ptr = (char*) malloc(seq_struct_str_len * sizeof(char));
          if (seq_struct_ptr == NULL) {
            #ifdef DEBUG
              fputs("Error in malloc()\n", stderr);
            #endif
              THAL_OOM_ERROR;
          }
          o->sec_struct = seq_struct_ptr;
          seq_struct_ptr[0] = '\0';
        }
        drawDimer(
          ps1, ps2,
          SHleft = 0.0, dH, dS,
          do_temponly, a->temp, o,
          oligo1, len1,
          oligo2, len2,
          saltCorrection, RC
        );
      }
      o->align_end_1 = bestI;
      o->align_end_2 = bestJ;
    } else {
      o->no_structure = 1;
      o->temp = 0.0;
    }
    free(ps1);
    free(ps2);
    free(SH);
    free(oligo2_rev);
    free(enthalpyDPT);
    free(entropyDPT);
    free(numSeq1);
    free(numSeq2);
    free(oligo1);
    return;
  }
  return;
}
/*** END thal() ***/

/* Set default args */
void
set_thal_default_args(thal_args *a)
{
  memset(a, 0, sizeof(thal_args));
  a->type = thal_any; /* thal_alignment_type THAL_ANY */
  a->maxLoop = MAX_LOOP;
  a->mv = 50; /* mM */
  a->dv = 0.0; /* mM */
  a->dntp = 0.8; /* mM */
  a->dna_conc = 50; /* nM */
  a->temp = TEMP_KELVIN; /* Kelvin */
  a->dimer = 1; /* by default dimer structure is calculated */
}

/* Set default args for oligo */
void
set_thal_oligo_default_args(thal_args *a)
{
  memset(a, 0, sizeof(thal_args));
  a->type = thal_any; /* thal_alignment_type THAL_ANY */
  a->maxLoop = MAX_LOOP;
  a->mv = 50; /* mM */
  a->dv = 0.0; /* mM */
  a->dntp = 0.0; /* mM */
  a->dna_conc = 50; /* nM */
  a->temp = TEMP_KELVIN; /* Kelvin */
  a->dimer = 1; /* by default dimer structure is calculated */
}


static unsigned char
str2int(char c)
{
  switch (c) {
    case 'A': case '0':
      return 0;
    case 'C': case '1':
      return 1;
    case 'G': case '2':
      return 2;
    case 'T': case '3':
      return 3;
  }
  return 4;
}

/* memory stuff */

static double*
safe_recalloc(double* ptr, int m, int n, thal_results* o)
{
  return (double*) safe_realloc(ptr, m * n * sizeof(double), o);
}

static void*
safe_calloc(size_t m, size_t n, thal_results *o)
{
  void* ptr;
  if ((ptr = calloc(m, n)) == NULL) {
#ifdef DEBUG
    fputs("Error in calloc()\n", stderr);
#endif
    THAL_OOM_ERROR;
  }
  return ptr;
}

static void*
safe_malloc(size_t n, thal_results *o)
{
  void* ptr;
  if ((ptr = malloc(n)) == NULL) {
#ifdef DEBUG
    fputs("Error in malloc()\n", stderr);
#endif
    THAL_OOM_ERROR;
  }
  return ptr;
}

static void*
safe_realloc(void* ptr, size_t n, thal_results *o)
{
  ptr = realloc(ptr, n);
  if (ptr == NULL) {
#ifdef DEBUG
    fputs("Error in realloc()\n", stderr);
#endif
    THAL_OOM_ERROR;
  }
  return ptr;
}

static int
max5(double a, double b, double c, double d, double e)
{
  if(a > b && a > c && a > d && a > e) return 1;
  else if(b > c && b > d && b > e) return 2;
  else if(c > d && c > e) return 3;
  else if(d > e) return 4;
  else return 5;
}

static void
push(struct tracer** stack, int i, int j, int mtrx, thal_results* o)
{
  struct tracer* new_top;
  new_top = (struct tracer*) safe_malloc(sizeof(struct tracer), o);
  new_top->i = i;
  new_top->j = j;
  new_top->mtrx = mtrx;
  new_top->next = *stack;
  *stack = new_top;
}

static void
reverse(unsigned char *s)
{
   int i,j;
   char c;
   for (i = 0, j = length_unsig_char(s)-1; i < j; i++, j--) {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}

#define INIT_BUF_SIZE 1024

static char*
readParamFile(const char* dirname, const char* fname, thal_results* o)
{
  FILE* file;
  char* ret = NULL;
  char* paramdir = NULL;
  paramdir = (char*) safe_malloc(strlen(dirname) + strlen(fname) + 2, o);
  strcpy(paramdir, dirname);
#ifdef OS_WIN
  if (paramdir[strlen(paramdir) - 1] != '\\') {
    strcat(paramdir, "\\\0");
  }
#else
  if (paramdir[strlen(paramdir) - 1] != '/') {
    strcat(paramdir, "/\0");
  }
#endif
  strcat(paramdir, fname);
  if (!(file = fopen(paramdir, "r"))) {
    snprintf(o->msg, 255, "Unable to open file %s", paramdir);
    if (paramdir != NULL) {
      free(paramdir);
      paramdir = NULL;
    }
    longjmp(_jmp_buf, 1);
    return NULL;
  }
  if (paramdir != NULL) {
    free(paramdir);
    paramdir = NULL;
  }
  char c;
  int i = 0;
  size_t ssz = INIT_BUF_SIZE;
  size_t remaining_size;
  remaining_size = ssz;
  ret = (char*) safe_malloc(ssz, o);
  while (1) {
    if (feof(file)) {
      ret[i] = '\0';
      fclose(file);
      return ret;
    }
    c = fgetc(file);
    remaining_size -= sizeof(char);
    if (remaining_size <= 0) {
      if (ssz >= INT_MAX / 2) {
        strcpy(o->msg, "Out of memory");
        free(ret);
        longjmp(_jmp_buf, 1);
        return NULL;
      } else {
        ssz += INIT_BUF_SIZE;
        remaining_size += INIT_BUF_SIZE;
      }
      ret = (char *) safe_realloc(ret, ssz, o);
    }
    ret[i] = c;
    i++;
  }
}

int
thal_load_parameters(const char *path, thal_parameters *a, thal_results* o)
{
  thal_free_parameters(a);
  if (setjmp(_jmp_buf) != 0) {
    printf("longjump\n");
    return -1;
  }
  a->dangle_dh = readParamFile(path, "dangle.dh", o);
  a->dangle_ds = readParamFile(path, "dangle.ds", o);
  a->loops_dh = readParamFile(path, "loops.dh", o);
  a->loops_ds = readParamFile(path, "loops.ds", o);
  a->stack_dh = readParamFile(path, "stack.dh", o);
  a->stack_ds = readParamFile(path, "stack.ds", o);
  a->stackmm_dh = readParamFile(path, "stackmm.dh", o);
  a->stackmm_ds = readParamFile(path, "stackmm.ds", o);
  a->tetraloop_dh = readParamFile(path, "tetraloop.dh", o);
  a->tetraloop_ds = readParamFile(path, "tetraloop.ds", o);
  a->triloop_dh = readParamFile(path, "triloop.dh", o);
  a->triloop_ds = readParamFile(path, "triloop.ds", o);
  a->tstack_tm_inf_ds = readParamFile(path, "tstack_tm_inf.ds", o);
  a->tstack_dh = readParamFile(path, "tstack.dh", o);
  a->tstack2_dh = readParamFile(path, "tstack2.dh", o);
  a->tstack2_ds = readParamFile(path, "tstack2.ds", o);
  return 0;
}

static double
saltCorrectS (double mv, double dv, double dntp)
{
  if(dv <= 0) dntp = dv;
  return 0.368*((log((mv + 120*(sqrt(fmax(0.0, dv - dntp)))) / 1000)));
}

static char*
th_read_str_line(char **str, thal_results* o)
{
  if (*str == NULL) {
    return NULL;
  }
  char *ptr = *str;
  char *ini = *str;
  while(1) {
    if ((*ptr == '\n') || (*ptr == '\0')) {
      char *ret = NULL;
      if (!(ret = (char *) malloc(sizeof(char) * (ptr - ini + 1)))) {
#ifdef DEBUG
        fputs("Error in malloc()\n", stderr);
#endif
        THAL_OOM_ERROR;
      }
      /* copy line */
      strncpy(ret, ini, (ptr - ini + 1));
      ret[ptr - ini] = '\0';

      if (*ptr == '\0') { /* End of String */
        *str = NULL;
      } else {
        ptr++;
        if (*ptr == '\0') { /* End of String */
          *str = NULL;
        } else {
          *str = ptr;
        }
      }
      if (ptr == ini) {
        if (ret != NULL) {
          free(ret);
        }
        return NULL;
      } else {
        return ret;
      }
    }
    ptr++;
  }
}


/* These functions are needed as "inf" cannot be read on Windows directly */
static double
readDouble(char **str, thal_results* o)
{
  double result;
  char *line = th_read_str_line(str, o);
  /* skip any spaces at beginning of the line */
  while (isspace(*line)) line++;
  if (!strncmp(line, "inf", 3)) {
    free(line);
    return _INFINITY;
  }
  sscanf(line, "%lf", &result);
  if (line != NULL) {
    free(line);
  }
  return result;
}

/* Reads a line containing 4 doubles, which can be specified as "inf". */
static void
readLoop(char **str, double *v1, double *v2, double *v3, thal_results *o)
{
  char *line = th_read_str_line(str, o);
  char *p = line, *q;
  /* skip first number on the line */
  while (isspace(*p)) p++;
  while (isdigit(*p)) p++;
  while (isspace(*p)) p++;
  /* read second number */
  q = p;
  while (!isspace(*q)) q++;
  *q = '\0'; q++;
  if (!strcmp(p, "inf")) *v1 = _INFINITY;
  else sscanf(p, "%lf", v1);
  while (isspace(*q)) q++;
  /* read third number */
  p = q;
  while (!isspace(*p)) p++;
  *p = '\0'; p++;
  if (!strcmp(q, "inf")) *v2 = _INFINITY;
  else sscanf(q, "%lf", v2);
  while (isspace(*p)) p++;
  /* read last number */
  q = p;
  while (!isspace(*q) && (*q != '\0')) q++;
  *q = '\0';
  if (!strcmp(p, "inf")) *v3 = _INFINITY;
  else sscanf(p, "%lf", v3);
  if (line != NULL) {
    free(line);
  }
}

/* Reads a line containing a short string and a double, used for reading a triloop or tetraloop. */
static int
readTLoop(char **str, char *s, double *v, int triloop, thal_results *o)
{
  char *line = NULL;
  char *p = NULL;
  char *q = NULL;
  line = th_read_str_line(str, o);
  if (!line) return -1;
  p = line;
  /* skip first spaces */
  while (isspace(*p)) p++;
  /* read the string */
  q = p;
  while (isalpha(*q)) q++;
  *q = '\0'; q++;
  if (triloop) {
    strncpy(s, p, 5);   /*triloop string has 5 characters*/
    s[5] = '\0';
  } else {
    strncpy(s, p, 6);   /*tetraloop string has 6 characters*/
    s[6] = '\0';
  }
  /* skip all spaces */
  while (isspace(*q)) q++;
  p = q;
  while (!isspace(*p) && (*p != '\0')) p++;
  *p = '\0';
  if (!strcmp(q, "inf")) *v = _INFINITY;
  else sscanf(q, "%lg", v);
  if (line != NULL) {
    free(line);
  }
  return 0;
}

static void
getStack(
  double stackEntropies[5][5][5][5],
  double stackEnthalpies[5][5][5][5],
  const thal_parameters *tp, thal_results* o)
{
  int i, j, ii, jj;
  char *pt_ds = tp->stack_ds;
  char *pt_dh = tp->stack_dh;

  for (i = 0; i < 5; ++i) {
    for (ii = 0; ii < 5; ++ii) {
      for (j = 0; j < 5; ++j) {
        for (jj = 0; jj < 5; ++jj) {
          if (i == 4 || j == 4 || ii == 4 || jj == 4) {
            stackEntropies[i][ii][j][jj] = -1.0;
            stackEnthalpies[i][ii][j][jj] = _INFINITY;
          } else {
            stackEntropies[i][ii][j][jj] = readDouble(&pt_ds, o);
            stackEnthalpies[i][ii][j][jj] = readDouble(&pt_dh, o);
            if (!isFinite(stackEntropies[i][ii][j][jj]) || !isFinite(stackEnthalpies[i][ii][j][jj])) {
              stackEntropies[i][ii][j][jj] = -1.0;
              stackEnthalpies[i][ii][j][jj] = _INFINITY;
            }
          }
        } // end for jj
      } // end for j
    } // end for ii
  } // end for i
}

static void
getStackint2(
  double stackint2Entropies[5][5][5][5],
  double stackint2Enthalpies[5][5][5][5],
  const thal_parameters *tp, thal_results* o)
{
  int i, j, ii, jj;
  char *pt_ds = tp->stackmm_ds;
  char *pt_dh = tp->stackmm_dh;

  for (i = 0; i < 5; ++i) {
    for (ii = 0; ii < 5; ++ii) {
      for (j = 0; j < 5; ++j) {
        for (jj = 0; jj < 5; ++jj) {
          if (i == 4 || j == 4 || ii == 4 || jj == 4) {
            stackint2Entropies[i][ii][j][jj] = -1.0;
            stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
          } else {
              stackint2Entropies[i][ii][j][jj] = readDouble(&pt_ds, o);
              stackint2Enthalpies[i][ii][j][jj] = readDouble(&pt_dh, o);
              if (!isFinite(stackint2Entropies[i][ii][j][jj]) || !isFinite(stackint2Enthalpies[i][ii][j][jj])) {
                stackint2Entropies[i][ii][j][jj] = -1.0;
                stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
              }
          }
        } // end for jj
      } // end for j
    } // end for ii
  } // end for i
}

/*
static void
verifyStackTable(double stack[5][5][5][5], char* type)
{

   int i, j, ii, jj;
   for (i = 0; i < 4; ++i)
     for (j = 0; j < 4; ++j)
       for (ii = 0; ii < 4; ++ii)
     for (jj = 0; jj < 4; ++jj)
       if (stack[i][j][ii][jj] != stack[jj][ii][j][i])
#ifdef DEBUG
         fprintf(stderr, "Warning: symmetrical stacks _are_ _not_ equal: %c-%c/%c-%c stack %s is %g; %c-%c/%c-%c stack %s is %g\n",
#endif
             BASES[i], BASES[j], BASES[ii], BASES[jj], type, stack[i][j][ii][jj], BASES[jj],
             BASES[ii], BASES[j], BASES[i], type, stack[jj][ii][j][i]);
}
*/

static void
getDangle(
  double dangleEntropies3[5][5][5],
  double dangleEnthalpies3[5][5][5],
  double dangleEntropies5[5][5][5],
  double dangleEnthalpies5[5][5][5],
  const thal_parameters *tp, thal_results* o)
{
  int i, j, k;
  char *pt_ds = tp->dangle_ds;
  char *pt_dh = tp->dangle_dh;

  for (i = 0; i < 5; ++i) {
    for (j = 0; j < 5; ++j) {
      for (k = 0; k < 5; ++k) {
        if (i == 4 || j == 4) {
          dangleEntropies3[i][k][j] = -1.0;
          dangleEnthalpies3[i][k][j] = _INFINITY;
        } else if (k == 4) {
          dangleEntropies3[i][k][j] = -1.0;
          dangleEnthalpies3[i][k][j] = _INFINITY;
        } else {
          dangleEntropies3[i][k][j] = readDouble(&pt_ds, o);
          dangleEnthalpies3[i][k][j] = readDouble(&pt_dh, o);
          if(!isFinite(dangleEntropies3[i][k][j]) || !isFinite(dangleEnthalpies3[i][k][j])) {
            dangleEntropies3[i][k][j] = -1.0;
            dangleEnthalpies3[i][k][j] = _INFINITY;
          }
        }
      } // end for k
    } // end for j
  } // end for i
  for (i = 0; i < 5; ++i) {
    for (j = 0; j < 5; ++j) {
      for (k = 0; k < 5; ++k) {
        if (i == 4 || j == 4) {
          dangleEntropies5[i][j][k] = -1.0;
          dangleEnthalpies5[i][j][k] = _INFINITY;
        } else if (k == 4) {
          dangleEntropies5[i][j][k] = -1.0;
          dangleEnthalpies5[i][j][k] = _INFINITY;
        } else {
          dangleEntropies5[i][j][k] = readDouble(&pt_ds, o);
          dangleEnthalpies5[i][j][k] = readDouble(&pt_dh, o);
          if(!isFinite(dangleEntropies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k])) {
            dangleEntropies5[i][j][k] = -1.0;
            dangleEnthalpies5[i][j][k] = _INFINITY;
          }
        }
      } // end for k
    } // end for j
  } // end for i
}

static void
getLoop(
  double hairpinLoopEntropies[30],
  double interiorLoopEntropies[30],
  double bulgeLoopEntropies[30],
  double hairpinLoopEnthalpies[30],
  double interiorLoopEnthalpies[30],
  double bulgeLoopEnthalpies[30],
  const thal_parameters *tp, thal_results* o)
{
  int k;
  char *pt_ds = tp->loops_ds;
  char *pt_dh = tp->loops_dh;

  for (k = 0; k < 30; ++k) {
      readLoop(&pt_ds, &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k], o);
      readLoop(&pt_dh, &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k], o);
  }
}

static void
getTstack(
  double tstackEntropies[5][5][5][5],
  double tstackEnthalpies[5][5][5][5],
  const thal_parameters *tp, thal_results* o)
{
  int i1, j1, i2, j2;
  char *pt_ds = tp->tstack_tm_inf_ds;
  char *pt_dh = tp->tstack_dh;

  for (i1 = 0; i1 < 5; ++i1) {
    for (i2 = 0; i2 < 5; ++i2) {
      for (j1 = 0; j1 < 5; ++j1) {
        for (j2 = 0; j2 < 5; ++j2) {
          if (i1 == 4 || j1 == 4) {
            tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
            tstackEntropies[i1][i2][j1][j2] = -1.0;
          } else if (i2 == 4 || j2 == 4) {
            tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
            tstackEnthalpies[i1][i2][j1][j2] = 0.0;
          } else {
            tstackEntropies[i1][i2][j1][j2] = readDouble(&pt_ds, o);
            tstackEnthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, o);
            if (!isFinite(tstackEntropies[i1][i2][j1][j2]) || !isFinite(tstackEnthalpies[i1][i2][j1][j2])) {
              tstackEntropies[i1][i2][j1][j2] = -1.0;
              tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
            }
          }
        } // end for j2
      } // end for j1
    } // end for i2
  } // end for i1
}

static void
getTstack2(
  double tstack2Entropies[5][5][5][5],
  double tstack2Enthalpies[5][5][5][5],
  const thal_parameters *tp, thal_results* o)
{

  int i1, j1, i2, j2;
  char *pt_ds = tp->tstack2_ds;
  char *pt_dh = tp->tstack2_dh;

  for (i1 = 0; i1 < 5; ++i1) {
    for (i2 = 0; i2 < 5; ++i2) {
      for (j1 = 0; j1 < 5; ++j1) {
        for (j2 = 0; j2 < 5; ++j2) {
          if (i1 == 4 || j1 == 4)  {
            tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
            tstack2Entropies[i1][i2][j1][j2] = -1.0;
          } else if (i2 == 4 || j2 == 4) {
            tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
            tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
          } else {
            tstack2Entropies[i1][i2][j1][j2] = readDouble(&pt_ds, o);
            tstack2Enthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, o);
            if (!isFinite(tstack2Entropies[i1][i2][j1][j2]) || !isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
              tstack2Entropies[i1][i2][j1][j2] = -1.0;
              tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
            }
          }
        } // end for j2
      } // end for j1
    } // end for i2
  } // end for i1
}

static void
getTriloop(
  triloop_t* _triloopEntropies,
  triloop_t* _triloopEnthalpies,
  int* num,
  const thal_parameters *tp, thal_results* o)
{
  /* NOTE: old malloc code for data remains as a reference */
  int i, size;
  double value;
  char *pt_ds = tp->triloop_ds;
  *num = 0;
  size = 16;
  // if (*triloopEntropies != NULL) {
  //   free(*triloopEntropies);
  //   *triloopEntropies = NULL;
  // }
  // *triloopEntropies = (triloop_t*) safe_calloc(16, sizeof(triloop_t), o);
  while (readTLoop(&pt_ds, _triloopEntropies[*num].loop, &value, 1, o) != -1) {
    for (i = 0; i < 5; ++i) {
      _triloopEntropies[*num].loop[i] = str2int(_triloopEntropies[*num].loop[i]);
    }
    _triloopEntropies[*num].value = value;
    ++*num;
    if (*num == size) {
      size *= 2;
      // *triloopEntropies = (triloop_t*) safe_realloc(*triloopEntropies, size * sizeof(triloop_t), o);
    }
  }
  // *triloopEntropies = (triloop_t*) safe_realloc(*triloopEntropies, *num * sizeof(triloop_t), o);

  char *pt_dh = tp->triloop_dh;
  *num = 0;
  size = 16;
  // if (*triloopEnthalpies != NULL) {
  //   free(*triloopEnthalpies);
  //   *triloopEnthalpies = NULL;
  // }
  // *triloopEnthalpies = (triloop_t*) safe_calloc(16, sizeof(triloop_t), o);
  while (readTLoop(&pt_dh, _triloopEnthalpies[*num].loop, &value, 1, o) != -1) {
    for (i = 0; i < 5; ++i) {
      triloopEnthalpies[*num].loop[i] = str2int(_triloopEnthalpies[*num].loop[i]);
    }
    _triloopEnthalpies[*num].value = value;
    ++*num;
    if (*num == size) {
      size *= 2;
      // *triloopEnthalpies = (triloop_t*) safe_realloc(*triloopEnthalpies, size * sizeof(triloop_t), o);
    }
  }
  // *triloopEnthalpies = (triloop_t*) safe_realloc(*triloopEnthalpies, *num * sizeof(triloop_t), o);

  /* NOTE: Commented out code to investigate allocation sizes */
  // size_t triloop_entropies_size = malloc_size(*triloopEnthalpies);
  // size_t triloop_enthalpies_size = malloc_size(*triloopEnthalpies);
  // printf("triloop_entropies_size = %zu\n", triloop_entropies_size);
  // printf("triloop_enthalpies_size = %zu\n", triloop_enthalpies_size);
}

static void
getTetraloop(
  tetraloop_t* _tetraloopEntropies,
  tetraloop_t* _tetraloopEnthalpies,
  int* num,
  const thal_parameters *tp, thal_results* o)
{
  /* NOTE: old malloc code for data remains as a reference */
  int i, size;
  double value;
  char *pt_ds = tp->tetraloop_ds;
  *num = 0;
  size = 16;

  // if (*tetraloopEntropies != NULL) {
  //    free(*tetraloopEntropies);
  //    *tetraloopEntropies = NULL;
  // }
  // *tetraloopEntropies = (tetraloop_t*) safe_calloc(16, sizeof(tetraloop_t), o);
  // *tetraloopEntropies = (tetraloop_t*) safe_calloc(1, 2048, o);

  while (readTLoop(&pt_ds, _tetraloopEntropies[*num].loop, &value, 0, o) != -1) {
    for (i = 0; i < 6; ++i) {
      _tetraloopEntropies[*num].loop[i] = str2int(_tetraloopEntropies[*num].loop[i]);
    }
    _tetraloopEntropies[*num].value = value;
    ++*num;
    if (*num == size) {
      size *= 2;
      // *tetraloopEntropies = (tetraloop_t*) safe_realloc(
      //   *tetraloopEntropies,
      //   size * sizeof(tetraloop_t),
      //   o
      // );
    }
  }
  // *tetraloopEntropies = (tetraloop_t*) safe_realloc(
  //    *tetraloopEntropies,
  //    *num * sizeof(tetraloop_t),
  //    o
  // );

  char *pt_dh = tp->tetraloop_dh;
  *num = 0;
  size = 16;

  // if (*tetraloopEnthalpies != NULL) {
  //   free(*tetraloopEnthalpies);
  //   *tetraloopEnthalpies = NULL;
  // }
  // *tetraloopEnthalpies = (tetraloop_t*) safe_calloc(16, sizeof(tetraloop_t), o);
  // *tetraloopEnthalpies = (tetraloop_t*) safe_calloc(1, 2048, o);


  while (readTLoop(&pt_dh, _tetraloopEnthalpies[*num].loop, &value, 0, o) != -1) {
    for (i = 0; i < 6; ++i) {
      _tetraloopEnthalpies[*num].loop[i] = str2int(_tetraloopEnthalpies[*num].loop[i]);
    }
    _tetraloopEnthalpies[*num].value = value;
    ++*num;
    if (*num == size) {
      size *= 2;
      // *tetraloopEnthalpies = (tetraloop_t*) safe_realloc(
      //   *tetraloopEnthalpies,
      //   size * sizeof(tetraloop_t),
      //   o
      // );
    }
  }
  // *tetraloopEnthalpies = (tetraloop_t*) safe_realloc(
  //   *tetraloopEnthalpies,
  //   *num * sizeof(tetraloop_t),
  //   o
  // );

  /* NOTE: Commented out code to investigate allocation sizes */
  // size_t tetraloop_entropies_size = malloc_size(*tetraloopEnthalpies);
  // size_t tetraloop_enthalpies_size = malloc_size(*tetraloopEnthalpies);
  // printf("tetraloop_entropies_size = %zu\n", tetraloop_entropies_size);
  // printf("tetraloop_enthalpies_size = %zu\n", tetraloop_enthalpies_size);
}

static void
tableStartATS(double atp_value, double atpS[5][5])
{

  int i, j;
  for (i = 0; i < 5; ++i) {
    for (j = 0; j < 5; ++j) {
      atpS[i][j] = 0.00000000001;
    }
  }
  atpS[0][3] = atpS[3][0] = atp_value;
}


static void
tableStartATH(double atp_value, double atpH[5][5])
{

  int i, j;
  for (i = 0; i < 5; ++i) {
    for (j = 0; j < 5; ++j) {
      atpH[i][j] = 0.0;
    }
  }
  atpH[0][3] = atpH[3][0] = atp_value;
}

static int
comp3loop(const void* loop1, const void* loop2)
{

  int i;
  const unsigned char* h1 = (const unsigned char*) loop1;
  const triloop_t *h2 = (const triloop_t*) loop2;

  for (i = 0; i < 5; ++i) {
    if (h1[i] < h2->loop[i]) {
      return -1;
    } else if (h1[i] > h2->loop[i]) {
      return 1;
    }
  }
  return 0;
}

static int
comp4loop(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = (const unsigned char*) loop1;
  const tetraloop_t *h2 = (const tetraloop_t*) loop2;

  for (i = 0; i < 6; ++i) {
    if (h1[i] < h2->loop[i]) {
      return -1;
    } else if (h1[i] > h2->loop[i]) {
      return 1;
    }
  }
  return 0;
}


static void
initMatrix(
    unsigned char* numSeq1, int len1,
    unsigned char* numSeq2, int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT
) {
  int i, j;
  for (i = 1; i <= len1; ++i) {
    for (j = 1; j <= len2; ++j) {
      if (bpIndx(numSeq1[i], numSeq2[j]) == 0)  {
        EnthalpyDPT(i, j) = _INFINITY;
        EntropyDPT(i, j) = -1.0;
      } else {
        EnthalpyDPT(i, j) = 0.0;
        EntropyDPT(i, j) = MinEntropy;
      }
    }
  }
}

static void
initMatrix2(
    unsigned char* numSeq1, int len1,
    unsigned char* numSeq2, int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT
) {
  int i, j;
  for (i = 1; i <= len1; ++i) {
    for (j = i; j <= len2; ++j) {
      if (j - i < MIN_HRPN_LOOP + 1 || (bpIndx(numSeq1[i], numSeq1[j]) == 0)) {
        EnthalpyDPT(i, j) = _INFINITY;
        EntropyDPT(i, j) = -1.0;
      } else {
        EnthalpyDPT(i, j) = 0.0;
        EntropyDPT(i, j) = MinEntropy;
      }
    }
  }
}

static void
fillMatrix(
    int maxLoop,
    thal_results *o,
    unsigned char* numSeq1, int len1,
    unsigned char* numSeq2, int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT,
    double dplx_init_S,
    double dplx_init_H,
    double RC
) {
  int d, i, j, ii, jj;
  double* SH;

  SH = (double*) safe_malloc(2 * sizeof(double), o);
  for (i = 1; i <= len1; ++i) {
    for (j = 1; j <= len2; ++j) {
      if (isFinite(EnthalpyDPT(i, j))) { /* if finite */
        SH[0] = -1.0;
        SH[1] = _INFINITY;
        LSH(i, j, SH);
        if(isFinite(SH[1])) {
           EntropyDPT(i,j) = SH[0];
           EnthalpyDPT(i,j) = SH[1];
        }
        if (i > 1 && j > 1) {
          maxTM(
            i, j,
            entropyDPT,
            enthalpyDPT,
            numSeq1, len1,
            numSeq2, len2,
            len3,
            dplx_init_S, dplx_init_H, RC
          ); /* stack: sets EntropyDPT(i, j) and EnthalpyDPT(i, j) */
          for(d = 3; d <= maxLoop + 2; d++) { /* max=30, length over 30 is not allowed */
            ii = i - 1;
            jj = - ii - d + (j + i);
            if (jj < 1) {
              ii -= abs(jj-1);
              jj = 1;
            }
            for (; ii > 0 && jj < j; --ii, ++jj) {
              if (isFinite(EnthalpyDPT(ii, jj))) {
                SH[0] = -1.0;
                SH[1] = _INFINITY;
                calc_bulge_internal(
                    ii, jj,
                    i, j, SH,
                    0, maxLoop,
                    numSeq1, len1,
                    numSeq2,  len2,
                    len3,
                    entropyDPT, enthalpyDPT,
                    dplx_init_S, dplx_init_H, RC
                );
                if(SH[0] < MinEntropyCutoff) {
                  /* to not give dH any value if dS is unreasonable */
                  SH[0] = MinEntropy;
                  SH[1] = 0.0;
                }
                if(isFinite(SH[1])) {
                  EnthalpyDPT(i, j) = SH[1];
                  EntropyDPT(i, j) = SH[0];
                }
              }
            } // end for ii jj
          } // end for d
        } /* if */
      }
    } /* end for j */
  } /*  end for i */
  free(SH);
}

static void
fillMatrix2(
    int maxLoop,
    thal_results* o,
    unsigned char* numSeq1, int len1,
    unsigned char* numSeq2, int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT,
    double dplx_init_S,
    double dplx_init_H,
    double RC
) {
  int i, j;

  thalflex_packet_t tf_packet = make_thalflex_packet(
    NULL, NULL,
    entropyDPT, enthalpyDPT,
    numSeq1, len1,
    numSeq2, len2,
    len3,
    dplx_init_S, dplx_init_H, RC
  );

  double* SH = NULL;
  SH = (double*) safe_malloc(2 * sizeof(double), o);
  for (j = 2; j <= len2; ++j) {
    for (i = j - MIN_HRPN_LOOP - 1; i >= 1; --i) {
      if (isFinite(EnthalpyDPT(i, j))) {
        SH[0] = -1.0;
        SH[1] = _INFINITY;
        maxTM2(
          i, j,
          entropyDPT, enthalpyDPT,
          numSeq1, len1,
          numSeq2, len2,
          len3,
          dplx_init_S,
          dplx_init_H,
          RC
        ); /* calculate stack */

        /* calculate Bulge and Internal loop and stack */
        CBI(i, j, SH, 0, maxLoop, &tf_packet);

        SH[0] = -1.0;
        SH[1] = _INFINITY;
        calc_hairpin(
          i, j, SH, 0,
          numSeq1, len1,
          numSeq2, len2,
          len3,
          entropyDPT, enthalpyDPT,
          dplx_init_S, dplx_init_H, RC
        );
        if(isFinite(SH[1])) {
          /* To not give dH any value if dS is unreasonable */
          if(SH[0] < MinEntropyCutoff) {
            SH[0] = MinEntropy;
            SH[1] = 0.0;
          }
          EntropyDPT(i,j) = SH[0];
          EnthalpyDPT(i,j) = SH[1];
        }
      }
    }
  }
  free(SH);
}


static void
calc_hairpin(
    int i,
    int j,
    double* EntropyEnthalpy,
    int traceback,
    unsigned char* numSeq1,
    int len1,
    unsigned char* numSeq2,
    int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT,
    double dplx_init_S,
    double dplx_init_H,
    double RC
) {
  int loopSize = j - i - 1;
  double G1, G2;
  G1 = G2 = -_INFINITY;
  double* SH;
  SH = (double*) safe_malloc(2 * sizeof(double), 0);
  SH[0] = -1.0;
  SH[1] = _INFINITY;

  if(loopSize < MIN_HRPN_LOOP) {
    EntropyEnthalpy[0] = -1.0;
    EntropyEnthalpy[1] = _INFINITY;
    return;
  }
  if (i <= len1 && len2 < j) {
    EntropyEnthalpy[0] = -1.0;
    EntropyEnthalpy[1] = _INFINITY;
    return;
  } else if (i > len2) {
    i -= len1;
    j -= len2;
  }
  if(loopSize <= 30) {
    EntropyEnthalpy[1] = hairpinLoopEnthalpies[loopSize - 1];
    EntropyEnthalpy[0] = hairpinLoopEntropies[loopSize - 1];
  } else {
    EntropyEnthalpy[1] = hairpinLoopEnthalpies[29];
    EntropyEnthalpy[0] = hairpinLoopEntropies[29];
  }

  if (loopSize > 3) { /* for loops 4 bp and more in length, terminal mm are accounted */
    EntropyEnthalpy[1] += tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
    EntropyEnthalpy[0] += tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq1[j]][numSeq1[j - 1]];
  } else if (loopSize == 3) { /* for loops 3 bp in length at-penalty is considered */
    EntropyEnthalpy[1] += atPenaltyH(numSeq1[i], numSeq1[j]);
    EntropyEnthalpy[0] += atPenaltyS(numSeq1[i], numSeq1[j]);
  }

  if (loopSize == 3) {  /* closing AT-penalty (+), triloop bonus, hairpin of 3 (+) */
    triloop_t* loop;
    if (numTriloops) {
      if ((loop = (triloop_t*) bsearch(numSeq1 + i, triloopEnthalpies,
                            numTriloops, sizeof(triloop_t), comp3loop))) {
        EntropyEnthalpy[1] += loop->value;
      }
      if ((loop = (triloop_t*) bsearch(numSeq1 + i, triloopEntropies,
                             numTriloops, sizeof(triloop_t), comp3loop))) {
          EntropyEnthalpy[0] += loop->value;
      }
    }
  } else if (loopSize == 4) { /* terminal mismatch, tetraloop bonus, hairpin of 4 */
    tetraloop_t* loop;
    if (numTetraloops) {
      if ((loop = (tetraloop_t*) bsearch(numSeq1 + i, tetraloopEnthalpies, numTetraloops, sizeof(tetraloop_t), comp4loop))) {
        EntropyEnthalpy[1] += loop->value;
      }
      if ((loop = (tetraloop_t*) bsearch(numSeq1 + i, tetraloopEntropies, numTetraloops, sizeof(tetraloop_t), comp4loop))) {
        EntropyEnthalpy[0] += loop->value;
      }
    }
  }
  if (!isFinite(EntropyEnthalpy[1])) {
    EntropyEnthalpy[1] = _INFINITY;
    EntropyEnthalpy[0] = -1.0;
  }
  if (
    isPositive(EntropyEnthalpy[1]) &&
    isPositive(EntropyEnthalpy[0]) &&
    (!isPositive(EnthalpyDPT(i, j)) || !isPositive(EntropyDPT(i, j)))
  ) { /* if both, S and H are positive */
    EntropyEnthalpy[1] = _INFINITY;
    EntropyEnthalpy[0] = -1.0;
  }
  RSH(i, j, SH);
  G1 = EntropyEnthalpy[1] + SH[1] -TEMP_KELVIN * (EntropyEnthalpy[0] + SH[0]);
  G2 = EnthalpyDPT(i, j) + SH[1] -TEMP_KELVIN * (EntropyDPT(i, j) + SH[0]);
  if(G2 < G1 && traceback == 0) {
    EntropyEnthalpy[0] = EntropyDPT(i, j);
    EntropyEnthalpy[1] = EnthalpyDPT(i, j);
  }
  free(SH);
  return;
}


static void
calc_bulge_internal(
    int i,
    int j,
    int ii,
    int jj,
    double* EntropyEnthalpy,
    int traceback,
    int maxLoop,
    unsigned char* numSeq1,
    int len1,
    unsigned char* numSeq2,
    int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT,
    double dplx_init_S,
    double dplx_init_H,
    double RC
) {
  int loopSize1, loopSize2, loopSize;
  double S, H, G1, G2;
  int N, N_loop;
  double* SH;
  SH = (double*) safe_malloc(2 * sizeof(double), 0);
  SH[0] = -1.0;
  SH[1] = _INFINITY;
  S = -1.0;
  H = _INFINITY;
  loopSize1 = ii - i - 1;
  loopSize2 = jj - j - 1;
  if(ii < jj) {
    N = ((2 * i)/2);
    N_loop = N;
    if(loopSize1 > 2) {N_loop -= (loopSize1 - 2);}
    if(loopSize2 > 2) {N_loop -= (loopSize2 - 2);}
  } else {
    N = ((2 * j)/2);
    N_loop = 2 * jj;
    if(loopSize1 > 2) {N_loop -= (loopSize1 - 2);}
    if(loopSize2 > 2) {N_loop -= (loopSize2 - 2);}
    N_loop = (N_loop/2) - 1;
  }
#ifdef DEBUG
  if (ii <= i) {
    fputs("Error in calc_bulge_internal(): ii is not greater than i\n", stderr);
  }
  if (jj <= j) {
    fputs("Error in calc_bulge_internal(): jj is not greater than j\n", stderr);
  }
#endif

#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      free(SH);
      return;
   }
#endif
#ifdef DEBUG
  if (loopSize1 == 0 && loopSize2 == 0) {
    fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
    return;
  }
#endif
  loopSize = loopSize1 + loopSize2 - 1;
  /* only bulges have to be considered */
  if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) {
    /* bulge loop of size one is treated differently the intervening nn-pair must be added */
    if(loopSize2 == 1 || loopSize1 == 1) {
      if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
        H = bulgeLoopEnthalpies[loopSize] +
            stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
        S = bulgeLoopEntropies[loopSize] +
            stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]];
      }
      if(isPositive(H) || isPositive(S)){
         H = _INFINITY;
         S = -1.0;
      }
      H += EnthalpyDPT(i, j);
      S += EntropyDPT(i, j);
      if(!isFinite(H)) {
        H = _INFINITY;
        S = -1.0;
      }
      RSH(ii, jj, SH);
      G1 = H + SH[1] - TEMP_KELVIN*(S + SH[0]);
      G2 = EnthalpyDPT(ii, jj) + SH[1] - TEMP_KELVIN*((EntropyDPT(ii, jj) + SH[0]));
      if( (G1 < G2) || (traceback == 1) ) {
        EntropyEnthalpy[0] = S;
        EntropyEnthalpy[1] = H;
      }
    } else {
      /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum
      * bulgeloop size is 30
      */
      H = (
        bulgeLoopEnthalpies[loopSize] +
        atPenaltyH(numSeq1[i], numSeq2[j]) +
        atPenaltyH(numSeq1[ii], numSeq2[jj])
      );
      H += EnthalpyDPT(i, j);

      S = (
        bulgeLoopEntropies[loopSize] +
        atPenaltyS(numSeq1[i], numSeq2[j]) +
        atPenaltyS(numSeq1[ii], numSeq2[jj])
      );
      S += EntropyDPT(i, j);
      if(!isFinite(H)) {
        H = _INFINITY;
        S = -1.0;
      }
      if(isPositive(H) && isPositive(S)) {
        H = _INFINITY;
        S = -1.0;
      }
      RSH(ii, jj, SH);
      G1 = H+SH[1] -TEMP_KELVIN*(S+SH[0]);
      G2 = EnthalpyDPT(ii, jj) + SH[1] - TEMP_KELVIN * (EntropyDPT(ii, jj) + SH[0]);
      if( (G1 < G2) || (traceback == 1) ) {
        EntropyEnthalpy[0] = S;
        EntropyEnthalpy[1] = H;
      }
    }
  } else if (loopSize1 == 1 && loopSize2 == 1) {
    S = (
      stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
      stackint2Entropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]
    );
    S += EntropyDPT(i, j);

    H = (
      stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
      stackint2Enthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]]
    );
    H += EnthalpyDPT(i, j);
    if(!isFinite(H)) {
      H = _INFINITY;
      S = -1.0;
    }
    if( isPositive(H) && isPositive(S) ) {
      H = _INFINITY;
      S = -1.0;
    }

    RSH(ii, jj, SH);
    G1 = H + SH[1] -TEMP_KELVIN * (S + SH[0]);
    G2 = EnthalpyDPT(ii, jj) + SH[1] - TEMP_KELVIN * (EntropyDPT(ii, jj) + SH[0]);
    if( (G1< G2) || (traceback == 1)) {
        EntropyEnthalpy[0] = S;
        EntropyEnthalpy[1] = H;
    }
    free(SH);
    return;
  } else { /* only internal loops */
    H = (
      interiorLoopEnthalpies[loopSize] +
      tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
      tstackEnthalpies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] +
      (ILAH * abs(loopSize1 - loopSize2))
    );
    H += EnthalpyDPT(i, j);

    S = (
      interiorLoopEntropies[loopSize] +
      tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
      tstackEntropies[numSeq2[jj]][numSeq2[jj-1]][numSeq1[ii]][numSeq1[ii-1]] +
      (ILAS * abs(loopSize1 - loopSize2))
    );
    S += EntropyDPT(i, j);

    if(!isFinite(H)) {
      H = _INFINITY;
      S = -1.0;
    }
    if(isPositive(H) && isPositive(S)) {
      H = _INFINITY;
      S = -1.0;
    }

    RSH(ii, jj, SH);
    G1 = H + SH[1] - TEMP_KELVIN*(S + SH[0]);
    G2 = EnthalpyDPT(ii, jj) + SH[1] - TEMP_KELVIN * (EntropyDPT(ii, jj) + SH[0]);

    if( (G1 < G2) || (traceback == 1) ) {
      EntropyEnthalpy[0] = S;
      EntropyEnthalpy[1] = H;
    }
  }
  free(SH);
  return;
}

static void
calc_bulge_internal2(
    int i,
    int j,
    int ii,
    int jj,
    double* EntropyEnthalpy,
    int traceback,
    int maxLoop,
    unsigned char* numSeq1,
    int len1,
    unsigned char* numSeq2,
    int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT,
    double dplx_init_S,
    double dplx_init_H,
    double RC
) {
  int loopSize1, loopSize2, loopSize;
  double T1, T2;
  double S,H;
  /* int N, N_loop; Triinu, please review */
  T1 = T2 = -_INFINITY;
  S = MinEntropy;
  H = 0.0;
  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
  if (loopSize1 + loopSize2 > maxLoop) {
    EntropyEnthalpy[0] = -1.0;
    EntropyEnthalpy[1] = _INFINITY;
    return;
  }
  /* Triinu, please review the statements below. */
  // if(i < (len1 -j)) {
   /* N  = i; */
    /* N_loop = (i - 1); */
  // } else {
   /* N = len1-j;  */
    /* N_loop = len1 - j - 1; */
  // }
#ifdef DEBUG
  if (ii <= i)
    fputs("Error in calc_bulge_internal(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in calc_bulge_internal(): jj isn't less than j\n", stderr);
  if (ii >= jj)
    fputs("Error in calc_bulge_internal(): jj isn't greater than ii\n", stderr);

  if ((i <= len1 && len1 < ii) || (jj <= len2 && len2 < j))  {
    EntropyEnthalpy[0] = -1.0;
    EntropyEnthalpy[1] = _INFINITY;
    return;
  }
#endif

#ifdef DEBUG
  if (loopSize1 + loopSize2 > maxLoop) {
    fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
    return;
  }
#endif
#ifdef DEBUG
  if (loopSize1 == 0 && loopSize2 == 0) {
    fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
    return;
  }
#endif

#ifdef DEBUG
  if (i > len1)
    i -= len1;
  if (ii > len1)
    ii -= len1;
  if (j > len2)
    j -= len2;
  if (jj > len2)
    jj -= len2;
#endif
  loopSize = loopSize1 + loopSize2 - 1; /* for indx only */
  if (
    (loopSize1 == 0 && loopSize2 > 0) ||
    (loopSize2 == 0 && loopSize1 > 0)
  ) { /* only bulges have to be considered */
    if(loopSize2 == 1 || loopSize1 == 1) {
      /* bulge loop of size one is treated differently the intervening
      * nn-pair must be added
      */
      if (
        (loopSize2 == 1 && loopSize1 == 0) ||
        (loopSize2 == 0 && loopSize1 == 1)
      ) {
        H = (
          bulgeLoopEnthalpies[loopSize] +
          stackEnthalpies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]]
        );
        S = (
          bulgeLoopEntropies[loopSize] +
          stackEntropies[numSeq1[i]][numSeq1[ii]][numSeq2[j]][numSeq2[jj]]
        );
      }
      if(traceback!=1) {
        H += EnthalpyDPT(ii, jj); /* bulge koos otsaga, st bulge i,j-ni */
        S += EntropyDPT(ii, jj);
      }

      if(!isFinite(H)) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);

      if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
        EntropyEnthalpy[0] = S;
        EntropyEnthalpy[1] = H;
      }

    } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

      H = (
        bulgeLoopEnthalpies[loopSize] +
        atPenaltyH(numSeq1[i], numSeq2[j]) +
        atPenaltyH(numSeq1[ii], numSeq2[jj])
      );
      if (traceback !=1 ) {
        H += EnthalpyDPT(ii, jj);
      }
      S = (
        bulgeLoopEntropies[loopSize] +
        atPenaltyS(numSeq1[i], numSeq2[j]) +
        atPenaltyS(numSeq1[ii], numSeq2[jj])
      );
      if (traceback != 1) {
        S += EntropyDPT(ii, jj);
      }
      if (!isFinite(H)) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
      T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);

      if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
        EntropyEnthalpy[0] = S;
        EntropyEnthalpy[1] = H;
      }
    }
  } /* end of calculating bulges */
  else if (loopSize1 == 1 && loopSize2 == 1) {
    /* mismatch nearest neighbor parameters */

    S = (
      stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
      stackint2Entropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]]
    );
    if(traceback != 1) {
      S += EntropyDPT(ii, jj);
    }

    H = (
      stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
      stackint2Enthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]]
    );
    if(traceback != 1) {
      H += EnthalpyDPT(ii, jj);
    }
    if(!isFinite(H)) {
      H = _INFINITY;
      S = -1.0;
    }
    T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
    T2 = (EnthalpyDPT(i, j) + dplx_init_H) / (EntropyDPT(i, j) + dplx_init_S + RC);
    if((DBL_EQ(T1,T2) == 2) || traceback) {
      if((T1 > T2) || ((traceback && T1 >= T2) || traceback==1)) {
        EntropyEnthalpy[0] = S;
        EntropyEnthalpy[1] = H;
      }
    }
    return;
  } else { /* only internal loops */

    H = (
      interiorLoopEnthalpies[loopSize] +
      tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
      tstackEnthalpies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]] +
      (ILAH * abs(loopSize1 - loopSize2))
    );
    if(traceback != 1) {
      H += EnthalpyDPT(ii, jj);
    }

    S = (
      interiorLoopEntropies[loopSize] +
      tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]] +
      tstackEntropies[numSeq2[jj]][numSeq2[jj+1]][numSeq1[ii]][numSeq1[ii-1]] +
      (ILAS * abs(loopSize1 - loopSize2))
    );
    if(traceback != 1) {
      S += EntropyDPT(ii, jj);
    }
    if(!isFinite(H)) {
      H = _INFINITY;
      S = -1.0;
    }
    T1 = (H + dplx_init_H) / ((S + dplx_init_S) + RC);
    T2 = (EnthalpyDPT(i, j) + dplx_init_H) / ((EntropyDPT(i, j)) + dplx_init_S + RC);
    if((T1 > T2) || ((traceback && T1 >= T2) || (traceback==1))) {
      EntropyEnthalpy[0] = S;
      EntropyEnthalpy[1] = H;
    }
  }
  return;
}

static void
calc_terminal_bp(
    double temp,
    double* send5,
    double* hend5,
    double* entropyDPT,
    double* enthalpyDPT,
    unsigned char* numSeq1,
    int len1,
    unsigned char* numSeq2,
    int len2,
    int len3,
    double dplx_init_S,
    double dplx_init_H,
    double RC
) { /* compute exterior loop */
  int i;
  int max;
  double T1, T2, T3, T4, T5;
  double G;
  SEND5(0) = SEND5(1) = -1.0;
  HEND5(0) = HEND5(1) = _INFINITY;
  for(i = 2; i<=(len1); i++) {
    SEND5(i) = MinEntropy;
    HEND5(i) = 0;
  }

  thalflex_packet_t tfpacket = make_thalflex_packet(
    send5, hend5,
    entropyDPT, enthalpyDPT,
    numSeq1, len1,
    numSeq2, len2,
    len3,
    dplx_init_S, dplx_init_H, RC
  );

  T1 = T2 = T3 = T4 = T5 = -_INFINITY;

  /* adding terminal penalties to 3' end and to 5' end */
  for(i = 2; i <= len1; ++i) {
    max = 0;
    T1 = T2 = T3 = T4 = T5 = -_INFINITY;
    T1 = (HEND5(i - 1) + dplx_init_H) / (SEND5(i - 1) + dplx_init_S + RC);
    T2 = (
      (END5_1(i, 1, &tfpacket) + dplx_init_H) /
      (END5_1(i,2, &tfpacket) + dplx_init_S + RC)
    );
    T3 = (
      (END5_2(i, 1, &tfpacket) + dplx_init_H) /
      (END5_2(i, 2, &tfpacket) + dplx_init_S + RC)
    );
    T4 = (
      (END5_3(i,1, &tfpacket) + dplx_init_H) /
      (END5_3(i, 2, &tfpacket) + dplx_init_S + RC)
    );
    T5 = (
      (END5_4(i, 1, &tfpacket) + dplx_init_H) /
      (END5_4(i, 2, &tfpacket) + dplx_init_S + RC)
    );
    max = max5(T1, T2, T3, T4, T5);
    switch (max) {
      case 1:
        SEND5(i) = SEND5(i - 1);
        HEND5(i) = HEND5(i - 1);
        break;
      case 2:
        G = (
          END5_1(i, 1, &tfpacket) -
          (temp * END5_1(i, 2, &tfpacket))
        );
        if(G < G2) {
          SEND5(i) = END5_1(i, 2, &tfpacket);
          HEND5(i) = END5_1(i, 1, &tfpacket);
        } else {
          SEND5(i) = SEND5(i - 1);
          HEND5(i) = HEND5(i - 1);
        }
        break;
      case 3:
        G = END5_2(i, 1, &tfpacket) - (temp * (END5_2(i, 2, &tfpacket)));
        if(G < G2) {
          SEND5(i) = END5_2(i, 2, &tfpacket);
          HEND5(i) = END5_2(i, 1, &tfpacket);
        } else {
          SEND5(i) = SEND5(i - 1);
          HEND5(i) = HEND5(i - 1);
        }
        break;
      case 4:
        G = END5_3(i, 1, &tfpacket) - (temp * (END5_3(i, 2, &tfpacket)));
        if(G < G2) {
          SEND5(i) = END5_3(i, 2, &tfpacket);
          HEND5(i) = END5_3(i, 1, &tfpacket);
        } else {
          SEND5(i) = SEND5(i - 1);
          HEND5(i) = HEND5(i - 1);
        }
        break;
      case 5:
        G = END5_4(i, 1, &tfpacket) - (temp * (END5_4(i, 2, &tfpacket)));
        if(G < G2) {
          SEND5(i) = END5_4(i, 2, &tfpacket);
          HEND5(i) = END5_4(i, 1, &tfpacket);
        } else {
          SEND5(i) = SEND5(i - 1);
          HEND5(i) = HEND5(i - 1);
        }
        break;
      default:
#ifdef DEBUG
        printf ("WARNING: max5 returned character code %d ??\n", max);
#endif
        break;
    }
  }
}

/* Return 1 if string is symmetrical, 0 otherwise. */
static int
symmetry_thermo(const unsigned char* seq)
{
  register char s;
  register char e;
  const unsigned char *seq_end = seq;
  int i = 0;
  int seq_len = length_unsig_char(seq);
  int mp = seq_len / 2;
  if((seq_len % 2) == 1) {
    return 0;
  }
  seq_end+=seq_len;
  seq_end--;
  while(i < mp) {
    i++;
    s = toupper(*seq);
    e = toupper(*seq_end);
    if ((s == 'A' && e != 'T')
        || (s == 'T' && e != 'A')
        || (e == 'A' && s != 'T')
        || (e == 'T' && s != 'A')) {
      return 0;
    }
    if ((s == 'C' && e != 'G')
        || (s == 'G' && e != 'C')
        || (e == 'C' && s !='G')
        || (e == 'G' && s != 'C')) {
      return 0;
    }
    seq++;
    seq_end--;
  }
  return 1;
}

static int
length_unsig_char(const unsigned char * str)
{
  int i = 0;
  while(*(str++)) {
    i++;
    if(i == INT_MAX) {
      return -1;
    }
  }
  return i;
}

static void
tracebacku(
    int* bp,
    int maxLoop,
    thal_results* o,
    double* send5,
    double* hend5,
    double* entropyDPT,
    double* enthalpyDPT,
    unsigned char* numSeq1,
    int len1,
    unsigned char* numSeq2,
    int len2,
    int len3,
    double dplx_init_S,
    double dplx_init_H,
    double RC
) { /* traceback for unimolecular structure */
  int i = 0, j = 0;
  int ii, jj, k;
  struct tracer *top, *stack = NULL;
  double* SH1;
  double* SH2;
  double* EntropyEnthalpy;
  SH1 = (double*) safe_malloc(2 * sizeof(double), o);
  SH2 = (double*) safe_malloc(2 * sizeof(double), o);
  EntropyEnthalpy = (double*) safe_malloc(2 * sizeof(double), o);
  push(&stack,len1, 0, 1, o);

  thalflex_packet_t tfpacket = make_thalflex_packet(
    send5, hend5,
    entropyDPT, enthalpyDPT,
    numSeq1, len1,
    numSeq2, len2,
    len3,
    dplx_init_S, dplx_init_H, RC
  );

  while(stack) {
    top = stack;
    stack = stack->next;
    i = top->i;
    j = top->j;
    if(top->mtrx == 1) {
      /* If previous structure is the same as this one */
      while (equal(SEND5(i), SEND5(i - 1)) && equal(HEND5(i), HEND5(i - 1))) {
        --i;
      }
      if (i == 0) {
        continue;
      }
      if (
          equal(SEND5(i), END5_1(i, 2, &tfpacket)) &&
          equal(HEND5(i), END5_1(i, 1, &tfpacket))
      ) {
        for (k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k) {
          if (
              equal(
                SEND5(i),
                (
                  atPenaltyS(numSeq1[k + 1], numSeq1[i]) +
                  EntropyDPT(k + 1, i)
                )
              ) &&
              equal(
                HEND5(i),
                (
                  atPenaltyH(numSeq1[k + 1], numSeq1[i]) +
                  EnthalpyDPT(k + 1, i)
                )
              )
          ) {
            push(&stack, k + 1, i,0, o);
            break;
          } else if (
              equal(
                SEND5(i),
                (
                  SEND5(k) +
                  atPenaltyS(numSeq1[k + 1], numSeq1[i]) +
                  EntropyDPT(k + 1, i)
                )
              ) &&
              equal(
                HEND5(i),
                (
                  HEND5(k) +
                  atPenaltyH(numSeq1[k + 1], numSeq1[i]) +
                  EnthalpyDPT(k + 1, i)
                )
              )
          ) {
            push(&stack, k + 1, i, 0, o);
            push(&stack, k, 0, 1, o);
            break;
          }
        }
      } else if (
          equal(
            SEND5(i),
            END5_2(i, 2, &tfpacket)
          ) &&
          equal(
            HEND5(i),
            END5_2(i, 1, &tfpacket)
          )
      ) {
        for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
          if (
              equal(
                SEND5(i),
                (
                  atPenaltyS(numSeq1[k + 2], numSeq1[i]) +
                  Sd5(i, k + 2, numSeq1) +
                  EntropyDPT(k + 2, i)
                )
              ) &&
              equal(
                HEND5(i),
                (
                  atPenaltyH(numSeq1[k + 2], numSeq1[i]) +
                  Hd5(i, k + 2, numSeq1) +
                  EnthalpyDPT(k + 2, i)
                )
              )
          ) {
            push(&stack, k + 2, i, 0, o);
            break;
          } else if (
              equal(
                SEND5(i),
                (
                  SEND5(k) +
                  atPenaltyS(numSeq1[k + 2], numSeq1[i]) +
                  Sd5(i, k + 2, numSeq1) +
                  EntropyDPT(k + 2, i)
                )
              ) &&
              equal(
                HEND5(i),
                (
                  HEND5(k) +
                  atPenaltyH(numSeq1[k + 2], numSeq1[i]) +
                  Hd5(i, k + 2, numSeq1) +
                  EnthalpyDPT(k + 2, i)
                )
              )
          ) {
            push(&stack, k + 2, i, 0, o);
            push(&stack, k, 0, 1, o);
            break;
          }
        }
      } else if (
          equal(SEND5(i), END5_3(i, 2, &tfpacket)) &&
          equal(HEND5(i), END5_3(i, 1, &tfpacket))
      ) {
        for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
          if (
              equal(
                SEND5(i),
                (
                  atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) +
                  Sd3(i - 1, k + 1, numSeq1) +
                  EntropyDPT(k + 1, i - 1)
                )
              ) &&
              equal(
                HEND5(i),
                (
                  atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) +
                  Hd3(i - 1, k + 1, numSeq1) +
                  EnthalpyDPT(k + 1, i - 1)
                )
              )
          ) {
            push(&stack, k + 1, i - 1, 0, o);
            break;
          } else if (
              equal(
                SEND5(i),
                (
                  SEND5(k) +
                  atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) +
                  Sd3(i - 1, k + 1, numSeq1) +
                  EntropyDPT(k + 1, i - 1)
                )
              ) &&
              equal(
                HEND5(i),
                (
                  HEND5(k) +
                  atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) +
                  Hd3(i - 1, k + 1, numSeq1) +
                  EnthalpyDPT(k + 1, i - 1)
                )
              )
          ) {
            push(&stack, k + 1, i - 1, 0, o); /* matrix 0  */
            push(&stack, k, 0, 1, o); /* matrix 3 */
            break;
          }
        }
      } else if(
          equal(
            SEND5(i),
            END5_4(i, 2, &tfpacket)
          ) &&
          equal(
            HEND5(i),
            END5_4(i, 1, &tfpacket)
          )
      ) {
        for (k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k) {
          if (
              equal(
                SEND5(i),
                (
                  atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) +
                  Ststack(i - 1, k + 2, numSeq1) +
                  EntropyDPT(k + 2, i - 1)
                )
              ) &&
              equal(
                HEND5(i),
                (
                  atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) +
                  Htstack(i - 1, k + 2, numSeq1) +
                  EnthalpyDPT(k + 2, i - 1)
                )
              )
          ) {
            push(&stack, k + 2, i - 1, 0, o);
            break;
          } else if (
              equal(
                SEND5(i),
                (
                  SEND5(k) +
                  atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) +
                  Ststack(i - 1, k + 2, numSeq1) +
                  EntropyDPT(k + 2, i - 1)
                )
              ) &&
              equal(
                HEND5(i),
                (
                  HEND5(k) +
                  atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) +
                  Htstack(i - 1, k + 2, numSeq1) +
                  EnthalpyDPT(k + 2, i - 1)
                )
              )
          ) {
            push(&stack, k + 2, i - 1, 0, o);
            push(&stack, k, 0, 1, o);
            break;
          }
        }
      }
    } else if(top->mtrx == 0) {
      bp[i - 1] = j;
      bp[j - 1] = i;
      SH1[0] = -1.0;
      SH1[1] = _INFINITY;
      calc_hairpin(
        i, j, SH1,
        1,  /* 1 means that we use this method in traceback */
        numSeq1, len1,
        numSeq2, len2,
        len3,
        entropyDPT, enthalpyDPT,
        dplx_init_S, dplx_init_H, RC
      );
      SH2[0] = -1.0;
      SH2[1] = _INFINITY;
      CBI(i, j, SH2, 2, maxLoop, &tfpacket);
      if (
            equal(
              EntropyDPT(i, j),
              (Ss(i, j, 2, numSeq1, len1, numSeq2, len2) + EntropyDPT(i + 1, j - 1))
            ) &&
            equal(
              EnthalpyDPT(i, j),
              (Hs(i, j, 2, numSeq1, len1, numSeq2, len2) + EnthalpyDPT(i + 1, j - 1))
            )
      ) {
        push(&stack, i + 1, j - 1, 0, o);
      } else if (equal(EntropyDPT(i, j), SH1[0]) && equal(EnthalpyDPT(i,j), SH1[1]));
      else if (
        equal(EntropyDPT(i, j), SH2[0]) &&
        equal(EnthalpyDPT(i, j), SH2[1])
      ) {
        int d, done;
        for (
          done = 0, d = (j - i - 3);
          (d >= MIN_HRPN_LOOP + 1) && (d >= j - i - 2 - maxLoop) && !done;
          --d
        ) {
          for (ii = i + 1; ii < j - d; ++ii) {
            jj = d + ii;
            EntropyEnthalpy[0] = -1.0;
            EntropyEnthalpy[1] = _INFINITY;
            calc_bulge_internal2(
              i, j,
              ii, jj,
              EntropyEnthalpy, 1, maxLoop,
              numSeq1, len1,
              numSeq2, len2,
              len3,
              entropyDPT, enthalpyDPT,
              dplx_init_S, dplx_init_H, RC
            );
            if (
                equal(EntropyDPT(i, j), EntropyEnthalpy[0] + EntropyDPT(ii, jj)) &&
                equal(EnthalpyDPT(i, j), EntropyEnthalpy[1] + EnthalpyDPT(ii, jj))
            ) {
              push(&stack, ii, jj, 0, o);
              ++done;
              break;
            }
          }
        }
      } else {
      }
    }
    free(top);
  }
  free(SH1);
  free(SH2);
  free(EntropyEnthalpy);
}


static void
traceback(
    int i,
    int j,
    int* ps1,
    int* ps2,
    int maxLoop,
    thal_results* o,
    double* entropyDPT,
    double* enthalpyDPT,
    unsigned char* numSeq1,
    int len1,
    unsigned char* numSeq2,
    int len2,
    int len3,
    double dplx_init_S,
    double dplx_init_H,
    double RC
) {
  int d, ii, jj, done;
  double* SH;
  SH = (double*) safe_malloc(2 * sizeof(double), o);
  ps1[i - 1] = j;
  ps2[j - 1] = i;
  while(1) {
    SH[0] = -1.0;
    SH[1] = _INFINITY;
    LSH(i, j, SH);
    if(equal(EntropyDPT(i, j), SH[0]) && equal(EnthalpyDPT(i, j), SH[1])) {
      break;
    }
    done = 0;
    if (
      (i > 1) &&
      (j > 1) &&
      equal(
        EntropyDPT(i, j),
        (Ss(i - 1, j - 1, 1, numSeq1, len1, numSeq2, len2) + EntropyDPT(i - 1, j - 1))
      ) &&
      equal(
        EnthalpyDPT(i, j),
        (Hs(i - 1, j - 1, 1, numSeq1, len1, numSeq2, len2) + EnthalpyDPT(i - 1, j - 1))
      )
    ) {
      i = i - 1;
      j = j - 1;
      ps1[i - 1] = j;
      ps2[j - 1] = i;
      done = 1;
    }
    for (d = 3; !done && d <= maxLoop + 2; ++d) {
      ii = i - 1;
      jj = -ii - d + (j + i);
      if (jj < 1) {
        ii -= abs(jj-1);
        jj = 1;
      }
      for (; !done && ii > 0 && jj < j; --ii, ++jj) {
        SH[0] = -1.0;
        SH[1] = _INFINITY;
        calc_bulge_internal(
          ii,
          jj,
          i,
          j,
          SH,
          1,
          maxLoop,
          numSeq1,
          len1,
          numSeq2,
          len2,
          len3,
          entropyDPT,
          enthalpyDPT,
          dplx_init_S,
          dplx_init_H,
          RC
        );
        if (
          equal(EntropyDPT(i, j), SH[0]) &&
          equal(EnthalpyDPT(i, j), SH[1])
        ) {
          i = ii;
          j = jj;
          ps1[i - 1] = j;
          ps2[j - 1] = i;
          done = 1;
          break;
        }
      }
    }
  }
  free(SH);
}


/* primer3-py special function */
static void
calcDimer(
    int* ps1,
    int* ps2,
    double temp,
    double H,
    double S,
    int temponly,
    double t37,
    thal_results *o,
    int len1,
    int len2,
    double saltCorrection,
    double RC
) {
    int i, N;
    double G, t;
    t = G = 0;
    if (!isFinite(temp)){
        if(temponly == 0) {
        }
        /* lets use generalization here; this should rather be very negative value */
        o->temp = 0.0;

        o->no_structure = 1;
        return;
    } else {
        N = 0;
        for(i = 0; i < len1; i++){
            if(ps1[i] > 0) ++N;
        }
        for(i = 0; i < len2; i++) {
            if(ps2[i] > 0) ++N;
        }
        N = (N/2) -1;
        t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
        if(temponly == 0) {
          G = (H) - (t37 * (S + (N * saltCorrection)));
          S = S + (N * saltCorrection);
          o->temp = (double) t;
          o->ds = (double) S;
          o->dh = (double) H;
          o->dg = (double) G;
        } else {
            o->temp = (double) t;
        }
    }
    return;
}


/* primer3-py special function */
static void
drawDimer(
    int* ps1,
    int* ps2,
    double temp,
    double H,
    double S,
    int temponly,
    double t37,
    thal_results *o,
    unsigned char* oligo1,
    int len1,
    unsigned char* oligo2,
    int len2,
    double saltCorrection,
    double RC
) {
  int i, j, k, numSS1, numSS2, N;
  double G = 0;
  double t = 0;
  char *output_buf = o->sec_struct;

  if (!isFinite(temp)) {
    if (temponly == 0 && output_buf == NULL) {
      printf("No predicted secondary structures for given sequences\n");
    }

    /* lets use generalization here; this should rather be very negative value */
    o->temp = 0.0;

    return;
  } else {
    N = 0;
    for(i = 0;i < len1; i++){
      if(ps1[i]>0) ++N;
    }
    for(i = 0; i < len2; i++) {
      if(ps2[i] > 0) ++N;
    }
    N = (N/2) -1;
    t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
    if(temponly==0) {
      G = (H) - (t37 * (S + (N * saltCorrection)));
      S = S + (N * saltCorrection);
      o->temp = (double) t;
      o->ds = (double) S;
      o->dh = (double) H;
      o->dg = (double) G;
      /* maybe user does not need as precise as that */
      /* printf("Thermodynamical values:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\tN = %d, SaltC=%f, RC=%f\n",
        len1, (double) S, (double) H, (double) G, (double) t, (int) N, saltCorrection, RC); */
      if (output_buf == NULL) {
        printf(
          "Calculated thermodynamical parameters for dimer:\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
          (double) S, (double) H, (double) G, (double) t
        );
      }
    } else {
      o->temp = (double) t;
      return;
    }
  }

  size_t str_block = (len1 + len2 + 1);

  /* primer3-py | MSVC C99 compiler does not support variable length arrays */
  #ifdef OS_WIN
  char* duplex_buffer = (char*)malloc(str_block * 4);
  #else
  char duplex_buffer[str_block * 4];
  #endif

  char* duplex_0 = duplex_buffer;
  char* duplex_1 = duplex_0 + str_block;
  char* duplex_2 = duplex_1 + str_block;
  char* duplex_3 = duplex_2 + str_block;
  size_t dlen_0 = 0;
  size_t dlen_1 = 0;
  size_t dlen_2 = 0;
  size_t dlen_3 = 0;

  i = 0;
  numSS1 = 0;
  while (ps1[i++] == 0) { ++numSS1; }
  j = 0;
  numSS2 = 0;
  while (ps2[j++] == 0) { ++numSS2; }

  if (numSS1 >= numSS2) {
    for (i = 0; i < numSS1; ++i) {
      duplex_0[dlen_0++] = oligo1[i];
      duplex_1[dlen_1++] = ' ';
      duplex_2[dlen_2++] = ' ';
    }
    for (j = 0; j < numSS1 - numSS2; ++j) { duplex_3[dlen_3++] = ' '; }
    for (j = 0; j < numSS2; ++j) { duplex_3[dlen_3++] = oligo2[j]; }
  } else {
    for (j = 0; j < numSS2; ++j) {
      duplex_3[dlen_3++] = oligo2[j];
      duplex_1[dlen_1++] = ' ';
      duplex_2[dlen_2++] = ' ';
    }
    for (i = 0; i < numSS2 - numSS1; ++i) {
      duplex_0[dlen_0++] = ' ';
    }
    for (i = 0; i < numSS1; ++i) {
      duplex_0[dlen_0++] = oligo1[i];
    }
  }
  i = numSS1 + 1;
  j = numSS2 + 1;

  while (i <= len1) {
    while (
      (i <= len1) &&
      (ps1[i - 1] != 0) &&
      (j <= len2) &&
      (ps2[j - 1] != 0)
    ) {
      duplex_0[dlen_0++] = ' ';
      duplex_1[dlen_1++] = oligo1[i - 1];
      duplex_2[dlen_2++] = oligo2[j - 1];
      duplex_3[dlen_3++] = ' ';
      ++i;
      ++j;
    }
    numSS1 = 0;
    while (i <= len1 && ps1[i - 1] == 0) {
      duplex_0[dlen_0++] = oligo1[i - 1];
      duplex_1[dlen_1++] = ' ';
      ++numSS1;
      ++i;
    }
    numSS2 = 0;
    while (j <= len2 && ps2[j - 1] == 0) {
      duplex_2[dlen_2++] = ' ';
      duplex_3[dlen_3++] = oligo2[j - 1];
      ++numSS2;
      ++j;
    }
    if (numSS1 < numSS2) {
      for (k = 0; k < numSS2 - numSS1; ++k) {
        duplex_0[dlen_0++] = '-';
        duplex_1[dlen_1++] = ' ';
      }
    } else if (numSS1 > numSS2) {
      for (k = 0; k < numSS1 - numSS2; ++k) {
        duplex_2[dlen_2++] = ' ';
        duplex_3[dlen_3++] = '-';
      }
    }
  }

  /* Trim trailing whitespace to ensure test output matches up */
  duplex_0[dlen_0+1] = 0;
  duplex_1[dlen_1+1] = 0;
  duplex_2[dlen_2+1] = 0;
  duplex_3[dlen_3+1] = 0;
  trim_trailing_whitespace(duplex_0, dlen_0);
  trim_trailing_whitespace(duplex_1, dlen_1);
  trim_trailing_whitespace(duplex_2, dlen_2);
  trim_trailing_whitespace(duplex_3, dlen_3);

  if (output_buf == NULL) {
    if (duplex_0[0] != 0) {
      printf("SEQ\t%s\n", duplex_0);
    } else {
      printf("SEQ\n");
    }

    if (duplex_1[0] != 0) {
      printf("SEQ\t%s\n", duplex_1);
    } else {
      printf("SEQ\n");
    }

    if (duplex_2[0] != 0) {
      printf("STR\t%s\n", duplex_2);
    } else {
      printf("STR\n");
    }

    if (duplex_3[0] != 0) {
      printf("STR\t%s\n", duplex_3);
    } else {
      printf("STR\n");
    }
  } else {
    snprintf(
      output_buf,
      str_block * 4 + 32,
      "SEQ\t%s\nSEQ\t%s\nSTR\t%s\nSTR\t%s\n",
      duplex_0,
      duplex_1,
      duplex_2,
      duplex_3
    );
  }
  /* primer3-py | MSVC C99 compiler does not support variable length arrays */
  #ifdef OS_WIN
  free(duplex_buffer);
  #endif
  return;
}


/* primer3-py special function */
static void
calcHairpin(
    int* bp,
    double mh,
    double ms,
    int temponly,
    double temp,
    thal_results *o,
    int len1,
    double saltCorrection
) {
    int i, N = 0;
    double mg, t;
    if (!isFinite(ms) || !isFinite(mh)) {
        if(temponly == 0) {
            #ifdef DEBUG
                fputs("No temperature could be calculated\\n",stderr);
            #endif
        } else {
            o->temp = 0.0;
            o->no_structure = 1;
        }
    } else {
        if(temponly == 0) {
            for (i = 1; i < len1; ++i) {
                if(bp[i-1] > 0) N++;
            }
        } else {
            for (i = 1; i < len1; ++i) {
                if(bp[i-1] > 0) N++;
            }
        }
        t = (mh / (ms + (((N / 2) - 1) * saltCorrection))) - ABSOLUTE_ZERO;

        if(temponly == 0) {
            mg = mh - (temp * (ms + (((N / 2) - 1) * saltCorrection)));
            ms = ms + (((N / 2) - 1) * saltCorrection);
            o->temp = (double) t;
            o->ds = (double) ms;
            o->dh = (double) mh;
            o->dg = (double) mg;
        } else {
            o->temp = (double) t;
        }
    }
    return;
}


/* primer3-py special function */
static void
drawHairpin(
    int* bp,
    double mh,
    double ms,
    int temponly,
    double temp,
    thal_results* o,
    unsigned char* oligo1,
    int len1,
    double saltCorrection
) {
  /* Plain text */
  int i, N = 0;
  double mg, t;
  char *output_buf = o->sec_struct;
  if (!isFinite(ms) || !isFinite(mh)) {
    if (temponly == 0 && output_buf == NULL) {
      printf("0\tdS = %g\tdH = %g\tinf\tinf\n", (double) ms, (double) mh);
#ifdef DEBUG
      fputs("No temperature could be calculated\n",stderr);
#endif
    } else {
      o->temp = 0.0; /* lets use generalization here */
    }
  } else {
    if(temponly == 0) {
      for (i = 1; i < len1; ++i) {
        if(bp[i-1] > 0) N++;
      }
    } else {
      for (i = 1; i < len1; ++i) {
        if(bp[i-1] > 0) N++;
      }
    }
    t = (mh / (ms + (((N / 2) - 1) * saltCorrection))) - ABSOLUTE_ZERO;

    if(temponly == 0) {
      mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));
      ms = ms + (((N/2)-1) * saltCorrection);
      o->ds = (double) ms;
      o->dh = (double) mh;
      o->dg = (double) mg;
      o->temp = (double) t;
      if (output_buf == NULL) {
        printf(
          "Calculated thermodynamical parameters for dimer:\t%d\tdS = %g\tdH = %g\tdG = %g\tt = %g\n",
          len1, (double) ms, (double) mh, (double) mg, (double) t
        );
      }
    } else {
      o->temp = (double) t;
      return;
    }
  }
  /* primer3-py | MSVC C99 compiler does not support variable length arrays */
  #ifdef OS_WIN
  /* plain-text output */
  char* asciiRow = (char*)malloc(len1);
  #else
  /* plain-text output */
  char asciiRow[len1];
  #endif

  for(i = 1; i < len1+1; ++i) {
    if(bp[i-1] == 0) {
      asciiRow[(i-1)] = '-';
    } else {
      if(bp[i-1] > (i-1)) {
        asciiRow[(bp[i-1]-1)]='\\';
      } else {
        asciiRow[(bp[i-1]-1)]='/';
      }
    }
  }
  if (output_buf != NULL) {
    snprintf(output_buf, 5, "SEQ\t");
    output_buf += 4;
  } else{
    printf("SEQ\t");
  }
  for(i = 0; i < len1; ++i) {
    if (output_buf != NULL) {
      snprintf(output_buf, 2, "%c", asciiRow[i]);
      output_buf++;
    } else {
      printf("%c", asciiRow[i]);
    }
  }
  if (output_buf != NULL) {
    snprintf(output_buf, len1 + 12, "\nSTR\t%s\n", oligo1);
  } else {
    printf("\nSTR\t%s\n", oligo1);
  }

  /* primer3-py | MSVC C99 compiler does not support variable length arrays */
  #ifdef OS_WIN
  free(asciiRow);
  #endif

  return;
}

/* NOTE: Commented out because primer3-py does not use the drawDimer and drawHairpin
 * code of upstream primer3 2.6.1 due to the fact it produces different strings
 * under different conditions.  TODO: refactor upstream drawDimer and drawHairpin
 * code to unify the output
 */
// static void
// // save_append_string(char** ret, int *space, thal_results *o, const char *str) {
//   int xlen, slen;
//   if (str == NULL) {
//     return;
//   }
//   if (*ret == NULL) {
//     *ret = (char *) safe_malloc(sizeof(char)*500, o);
//     *ret[0] = '\0';
//     *space = 500;
//   }
//   xlen = strlen(*ret);
//   slen = strlen(str);
//   if (xlen + slen + 1 > *space) {
//     *space += 4 * (slen + 1);
//     *ret = (char *) safe_realloc(*ret, *space, o);
//   }
//   strcpy(*ret + xlen, str);
//   return;
// }

// static void
// save_append_char(char** ret, int *space, thal_results *o, const char str) {
//   char fix[3];
//   fix[0] = str;
//   fix[1] = '\0';
//   save_append_string(ret, space, o, fix);
// }


static int
equal(double a, double b)
{
#ifdef INTEGER
  return a == b;
#endif

  if (!isfinite(a) || !isfinite(b))
    return 0;
  return fabs(a - b) < 1e-5;

  // NOTE: code not reachable, commented out NC 2018.07.11
  // if (a == 0 && b == 0)
  //   return 1;
}


int
trim_trailing_whitespace(char *msg, size_t msg_len)
{
  char *end_ptr;
  if(*msg == 0) {
    return 0;
  }

  /* Work back from the end of the string to find last non-whitespace char */
  end_ptr = msg + msg_len - 1;
  while(end_ptr >= msg && isspace((unsigned char) *end_ptr)) { end_ptr--; }
  /* Write new termination */
  end_ptr[1] = '\0';
  return 0;
}
