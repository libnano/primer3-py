/* libprimer3/thalflexsignatures.h

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

Modifications to thal.c code to create thalflexsignatures.h.
Copyright (C) 2023. Ben Pruitt & Nick Conway;

See thalflex.c for more info on thalflexsignatures.h
*/
#ifndef _P3PY_THALFLEXSIG_H
#define _P3PY_THALFLEXSIG_H

#include <stdio.h>
#include "thal.h"

/*** BEGIN STRUCTs ***/

typedef struct triloop_t {
  char loop[5];
  double value;
}  triloop_t;

typedef struct tetraloop_t {
  char loop[6];
  double value;
} tetraloop_t;

struct tracer /* structure for tracebacku - unimolecular str */ {
  int i;
  int j;
  int mtrx; /* [0 1] EntropyDPT/EnthalpyDPT*/
  struct tracer* next;
};

/*** END STRUCTs ***/

/* returns length of unsigned char; to avoid warnings while compiling */
static int length_unsig_char(const unsigned char * str);

/* converts DNA sequence to int; 0-A, 1-C, 2-G, 3-T, 4-whatever */
static unsigned char str2int(char c);

/* part of calculating salt correction for Tm by SantaLucia et al */
static double saltCorrectS (double mv, double dv, double dntp);

/* file of thermodynamic params */
static char* readParamFile(const char* dirname, const char* fname, thal_results* o);

/* get thermodynamic tables */
static double readDouble(char **str, thal_results* o);

static void readLoop(char **str, double *v1, double *v2, double *v3, thal_results *o);

static int readTLoop(char **str, char *s, double *v, int triloop, thal_results *o);

static void getStack(
  double stackEntropies[5][5][5][5],
  double stackEnthalpies[5][5][5][5],
  const thal_parameters *tp, thal_results* o);

static void getStackint2(
  double stackEntropiesint2[5][5][5][5],
  double stackint2Enthalpies[5][5][5][5],
  const thal_parameters *tp, thal_results* o);

static void getDangle(
  double dangleEntropies3[5][5][5],
  double dangleEnthalpies3[5][5][5],
  double dangleEntropies5[5][5][5],
  double dangleEnthalpies5[5][5][5],
  const thal_parameters *tp, thal_results* o);

static void getTstack(
  double tstackEntropies[5][5][5][5],
  double tstackEnthalpies[5][5][5][5],
  const thal_parameters *tp, thal_results* o);

static void getTstack2(
  double tstack2Entropies[5][5][5][5],
  double tstack2Enthalpies[5][5][5][5],
  const thal_parameters *tp, thal_results* o);

static void getTriloop(
  triloop_t*,
  triloop_t*,
  int* num,
  const thal_parameters *tp,
  thal_results* o
);

static void getTetraloop(
  tetraloop_t*,
  tetraloop_t*,
  int* num,
  const thal_parameters *tp,
  thal_results* o
);

static void getLoop(
  double hairpinLoopEnntropies[30],
  double interiorLoopEntropies[30],
  double bulgeLoopEntropiess[30],
  double hairpinLoopEnthalpies[30],
  double interiorLoopEnthalpies[30],
  double bulgeLoopEnthalpies[30],
  const thal_parameters *tp,
  thal_results* o
);

/* creates table of entropy values for nucleotides to which AT-penlty must be applied */
static void tableStartATS(double atp_value, double atp[5][5]);

static void tableStartATH(double atp_value, double atp[5][5]);

/* checks if sequence consists of specific triloop */
static int comp3loop(const void*, const void*);

/* checks if sequence consists of specific tetraloop */
static int comp4loop(const void*, const void*);

/* initiates thermodynamic parameter tables of entropy and enthalpy for dimer */
static void initMatrix(
    unsigned char* numSeq1, int len1,
    unsigned char* numSeq2, int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT
);

/* initiates thermodynamic parameter tables of entropy and enthalpy for monomer */
static void initMatrix2(
    unsigned char* numSeq1, int len1,
    unsigned char* numSeq2, int len2,
    int len3,
    double* entropyDPT,
    double* enthalpyDPT
);

/* calc-s thermod values into dynamic progr table (dimer) */
static void fillMatrix(
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
);

/* calc-s thermod values into dynamic progr table (monomer) */
static void fillMatrix2(
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
);


/* calculates bulges and internal loops for dimer structures */
static void calc_bulge_internal(
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
);

/* calculates bulges and internal loops for monomer structures */
static void calc_bulge_internal2(
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
);


/* finds monomer structure that has maximum Tm */
static void calc_hairpin(
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
);

/* Returns stack entropy */
static double Ss(
  int i, int j, int k,
  unsigned char* numSeq1, int len1,
  unsigned char* numSeq2, int len2
);

/* Returns stack enthalpy */
static double Hs(
    int i, int j, int k,
    unsigned char* numSeq1, int len1,
    unsigned char* numSeq2, int len2
);

static void reverse(unsigned char *s);

static int max5(double, double, double, double, double);

/* Is sequence symmetrical */
static int symmetry_thermo(const unsigned char* seq);

/* traceback for dimers */
static void traceback(
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
);

/* traceback for hairpins */
static void tracebacku(
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
);

/* calculates but does not print dimer structure */
/* primer3-py special function */
static void calcDimer(
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
);

/* prints ascii output of dimer structure */
/* primer3-py special function */
static void drawDimer(
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
);

/* calculates but does not print hairpin structure */
/* primer3-py special function */
static void calcHairpin(
    int* bp,
    double mh,
    double ms,
    int temponly,
    double temp,
    thal_results *o,
    int len1,
    double saltCorrection
);

/* prints ascii output of hairpin structure */
/* primer3-py special function */
static void drawHairpin(
    int* bp,
    double mh,
    double ms,
    int temponly,
    double temp,
    thal_results* o,
    unsigned char* oligo1,
    int len1,
    double saltCorrection
);

int trim_trailing_whitespace(char *msg, size_t msg_len);

static int equal(double a, double b);

/* To add elements to struct */
static void push(struct tracer**, int, int, int, thal_results*);

/* terminal bp for monomer structure */
static void calc_terminal_bp(
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
);


/* memory stuff */
static void* safe_calloc(size_t, size_t, thal_results* o);
static void* safe_malloc(size_t, thal_results* o);
static void* safe_realloc(void*, size_t, thal_results* o);
static double* safe_recalloc(double* ptr, int m, int n, thal_results* o);

#endif
