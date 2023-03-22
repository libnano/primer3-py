/* libprimer3/thalflex.h

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

Modifications to thal.c code to create thalflex.h.
Copyright (C) 2023. Ben Pruitt & Nick Conway;

See thalflex.c for more info on thalflex.h
*/
#ifndef _P3PY_THALFLEX_H
#define _P3PY_THALFLEX_H

typedef struct thalflex_packet_t {
    double* send5;
    double* hend5;
    double* entropyDPT;
    double* enthalpyDPT;
    unsigned char* numSeq1;
    int len1;
    unsigned char* numSeq2;
    int len2;
    int len3;
    double dplx_init_S;
    double dplx_init_H;
    double RC;
} thalflex_packet_t;


thalflex_packet_t make_thalflex_packet(
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
) {
  thalflex_packet_t tfpacket = {
      .send5 = send5,
      .hend5 = hend5,
      .entropyDPT = entropyDPT,
      .enthalpyDPT = enthalpyDPT,
      .numSeq1 = numSeq1,
      .len1 = len1,
      .numSeq2 = numSeq2,
      .len2 = len2,
      .len3 = len3,
      .dplx_init_S = dplx_init_S,
      .dplx_init_H = dplx_init_H,
      .RC = RC
  };
  return tfpacket;
}

/* Handle unpacking boiler plate */
#define UNPACK_TFPACKET_CBI()                   \
  double* entropyDPT = tf_packet->entropyDPT;   \
  double* enthalpyDPT = tf_packet->enthalpyDPT; \
  unsigned char* numSeq1 = tf_packet->numSeq1;  \
  int len1 = tf_packet->len1;                   \
  unsigned char* numSeq2 = tf_packet->numSeq2;  \
  int len2 = tf_packet->len2;                   \
  int len3 = tf_packet->len3;                   \
  double dplx_init_S = tf_packet->dplx_init_S;  \
  double dplx_init_H = tf_packet->dplx_init_H;  \
  double RC = tf_packet->RC;


#define UNPACK_TFPACKET_ALL()       \
  double* send5 = tf_packet->send5; \
  double* hend5 = tf_packet->hend5; \
  UNPACK_TFPACKET_CBI()


#define UNPACK_TFPACKET_END5()                    \
  double* send5 = tf_packet->send5;               \
  double* hend5 = tf_packet->hend5;               \
  double* entropyDPT = tf_packet->entropyDPT;     \
  double* enthalpyDPT = tf_packet->enthalpyDPT;   \
  unsigned char* numSeq1 = tf_packet->numSeq1;    \
  int len3 = tf_packet->len3;                     \
  double dplx_init_S = tf_packet->dplx_init_S;    \
  double dplx_init_H = tf_packet->dplx_init_H;    \
  double RC = tf_packet->RC;


/* Returns thermodynamic value (S) for 5' dangling end */
#define Sd5(i, j, numSeq1) (dangleEntropies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]])

/* Returns thermodynamic value (H) for 5' dangling end */
#define Hd5(i, j, numSeq1) (dangleEnthalpies5[numSeq1[i]][numSeq1[j]][numSeq1[j - 1]])

/* Returns thermodynamic value (S) for 3' dangling end */
#define Sd3(i, j, numSeq1) (dangleEntropies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]])

/* Returns thermodynamic value (H) for 3' dangling end */
#define Hd3(i, j, numSeq1) (dangleEnthalpies3[numSeq1[i]][numSeq1[i+1]][numSeq1[j]])

/* Returns entropy value for terminal stack */
#define Ststack(i, j, numSeq1) (tstack2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]])

/* Returns enthalpy value for terminal stack */
#define Htstack(i, j, numSeq1) (tstack2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq1[j]][numSeq1[j-1]])


/* table where bp-s enthalpies, that retrieve to the most stable Tm, are saved */
#define EnthalpyDPT(i, j) enthalpyDPT[(j) + ((i-1)*len3) - (1)]

/* table where bp-s entropies, that retrieve to the most stable Tm, are saved */
#define EntropyDPT(i, j) entropyDPT[(j) + ((i-1)*len3) - (1)]

/* entropies of most stable hairpin terminal bp */
#define SEND5(i) send5[i]

/* enthalpies of most stable hairpin terminal bp */
#define HEND5(i) hend5[i]


#define bpIndx(a, b) BPI[a][b] /* for traceing matrix BPI */
#define atPenaltyS(a, b) atpS[a][b]
#define atPenaltyH(a, b) atpH[a][b]

/* NOTE: END5_{n}
* Executed in calc_terminal_bp; to find structure that corresponds to max Tm for
* terminal bp
* END5_1(X,1/2) - 1=Enthalpy, 2=Entropy
*/

// static double
// END5_1(
//     int i,
//     int hs,
//     double* send5,
//     double* hend5,
//     double* entropyDPT,
//     double* enthalpyDPT,
//     unsigned char* numSeq1,
//     int len3,
//     double dplx_init_S,
//     double dplx_init_H,
//     double RC
// ) {

static double
END5_1(
    int i,
    int hs,
    thalflex_packet_t* tf_packet
) {
  /* Unpack the tf_packet boiler plate variable definition */
  UNPACK_TFPACKET_END5()

  int k;
  double max_tm; /* energy min */
  double T1, T2;
  double H, S;
  double H_max, S_max;
  H_max = H = _INFINITY;
  S_max = S = -1.0;
  T1 = T2 = -_INFINITY;
  max_tm = -_INFINITY;
  for(k = 0; k <= i - MIN_HRPN_LOOP - 2; ++k) {
    T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
    T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
    if(T1 >= T2) {
      H = (
        HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i]) +
        EnthalpyDPT(k + 1, i)
      );
      S = (
        SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i]) +
        EntropyDPT(k + 1, i)
      );

      if(!isFinite(H) || H > 0 || S > 0) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
    } else {
      H = 0 + atPenaltyH(numSeq1[k + 1], numSeq1[i]) + EnthalpyDPT(k + 1, i);
      S = 0 + atPenaltyS(numSeq1[k + 1], numSeq1[i]) + EntropyDPT(k + 1, i);

      if(!isFinite(H) || H > 0 || S > 0) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
    }
    if(max_tm < T1) {
      if(S > MinEntropyCutoff) {
        H_max = H;
        S_max = S;
        max_tm = T1;
      }
    }
  }
  if (hs == 1) return H_max;
  return S_max;
}


static double
END5_2(
    int i,
    int hs,
    thalflex_packet_t* tf_packet
) {
  /* Unpack the tf_packet boiler plate variable definition */
  UNPACK_TFPACKET_END5()

  int k;
  double max_tm;
  double T1, T2;
  double H, S;
  double H_max, S_max;
  H_max = H = _INFINITY;
  T1 = T2 = max_tm = -_INFINITY;
  S_max = S = -1.0;
  for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
    T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
    T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
    if(T1 >= T2) {
      H = (
        HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i]) +
        Hd5(i, k + 2, numSeq1) + EnthalpyDPT(k + 2, i)
      );
      S = (
        SEND5(k) + atPenaltyS(numSeq1[k + 2], numSeq1[i]) +
        Sd5(i, k + 2, numSeq1) + EntropyDPT(k + 2, i)
      );
      if(!isFinite(H) || H > 0 || S > 0) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
    } else {
      H = (
        0 + atPenaltyH(numSeq1[k + 2], numSeq1[i]) +
        Hd5(i, k + 2, numSeq1) + EnthalpyDPT(k + 2, i)
      );
      S = (
        0 + atPenaltyS(numSeq1[k + 2], numSeq1[i]) +
        Sd5(i, k + 2, numSeq1) + EntropyDPT(k + 2, i)
      );
      if(!isFinite(H) || H > 0 || S > 0) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
    }
    if(max_tm < T1) {
      if(S > MinEntropyCutoff) {
        H_max = H;
        S_max = S;
        max_tm = T1;
      }
    }
  }
  if (hs == 1) return H_max;
  return S_max;
}


static double
END5_3(
    int i,
    int hs,
    thalflex_packet_t* tf_packet
) {
  /* Unpack the tf_packet boiler plate variable definition */
  UNPACK_TFPACKET_END5();

  int k;
  double max_tm;
  double T1, T2;
  double H, S;
  double H_max, S_max;
  H_max = H = _INFINITY;;
  T1 = T2 = max_tm = -_INFINITY;
  S_max = S = -1.0;
  for (k = 0; k <= i - MIN_HRPN_LOOP - 3; ++k) {
    T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
    T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
    if(T1 >= T2) {
      H = (
        HEND5(k) + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) +
        Hd3(i - 1, k + 1, numSeq1) + EnthalpyDPT(k + 1, i - 1)
      );
      S = (
        SEND5(k) + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) +
        Sd3(i - 1, k + 1, numSeq1) + EntropyDPT(k + 1, i - 1)
      );
      if(!isFinite(H) || (H > 0) || (S > 0)) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
    } else {
      H = (
        0 + atPenaltyH(numSeq1[k + 1], numSeq1[i - 1]) +
        Hd3(i - 1, k + 1, numSeq1) + EnthalpyDPT(k + 1, i - 1)
      );
      S = (
        0 + atPenaltyS(numSeq1[k + 1], numSeq1[i - 1]) +
        Sd3(i - 1, k + 1, numSeq1) + EntropyDPT(k + 1, i - 1)
      );
      if(!isFinite(H) || (H > 0) || (S > 0)) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
    }
    if(max_tm < T1) {
      if(S > MinEntropyCutoff) {
        H_max = H;
        S_max = S;
        max_tm = T1;
      }
    }
  }
  if (hs == 1) { return H_max; }
  return S_max;
}


static double
END5_4(
    int i,
    int hs,
    thalflex_packet_t* tf_packet
) {
  /* Unpack the tf_packet boiler plate variable definition */
  UNPACK_TFPACKET_END5()

  int k;
  double max_tm;
  double T1, T2;
  double H, S;
  double H_max, S_max;
  H_max = H = _INFINITY;
  T1 = T2 = max_tm = -_INFINITY;
  S_max = S = -1.0;
  for(k = 0; k <= i - MIN_HRPN_LOOP - 4; ++k) {
    T1 = (HEND5(k) + dplx_init_H) /(SEND5(k) + dplx_init_S + RC);
    T2 = (0 + dplx_init_H) /(0 + dplx_init_S + RC);
    if(T1 >= T2) {
      H = (
        HEND5(k) + atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) +
        Htstack(i - 1, k + 2, numSeq1) +
        EnthalpyDPT(k + 2, i - 1)
      );
      S = (
        SEND5(k) +
        atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) +
        Ststack(i - 1, k + 2, numSeq1) + EntropyDPT(k + 2, i - 1)
      );
      if(!isFinite(H) || H > 0 || S > 0) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) / (S + dplx_init_S + RC);
    } else {
      H = (
        0 +
        atPenaltyH(numSeq1[k + 2], numSeq1[i - 1]) +
        Htstack(i - 1, k + 2, numSeq1) +
        EnthalpyDPT(k + 2, i - 1)
      );
      S = (
        0 +
        atPenaltyS(numSeq1[k + 2], numSeq1[i - 1]) +
        Ststack(i - 1, k + 2, numSeq1) +
        EntropyDPT(k + 2, i - 1)
      );
      if(!isFinite(H) || H > 0 || S > 0) {
        H = _INFINITY;
        S = -1.0;
      }
      T1 = (H + dplx_init_H) /(S + dplx_init_S + RC);
    }
    if(max_tm < T1) {
      if(S > MinEntropyCutoff) {
        H_max = H;
        S_max = S;
        max_tm = T1;
      }
    }
  }
  if (hs == 1) {return H_max;}
  return S_max;
}


/* Calculate terminal entropy S and terminal enthalpy H starting reading from
* 5'end (Left hand/3' end - Right end)
*/
#define LSH(i, j, EntropyEnthalpy)                                              \
{                                                                               \
  double S1, H1, T1, G1;                                                        \
  double S2, H2, T2, G2;                                                        \
  S1 = S2 = -1.0;                                                               \
  H1 = H2 = -_INFINITY;                                                         \
  T1 = T2 = -_INFINITY;                                                         \
  if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {                                    \
    EntropyDPT(i, j) = -1.0;                                                    \
    EnthalpyDPT(i, j) = _INFINITY;                                              \
  }                                                                             \
  else {                                                                        \
    S1 = (                                                                      \
      atPenaltyS(numSeq1[i], numSeq2[j]) +                                      \
      tstack2Entropies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]]      \
    );                                                                          \
    H1 = (                                                                      \
      atPenaltyH(numSeq1[i], numSeq2[j]) +                                      \
      tstack2Enthalpies[numSeq2[j]][numSeq2[j-1]][numSeq1[i]][numSeq1[i-1]]     \
    );                                                                          \
    G1 = H1 - (TEMP_KELVIN * S1);                                               \
    if(!isFinite(H1) || (G1 > 0)) {                                             \
      H1 = _INFINITY;                                                           \
      S1 = -1.0;                                                                \
      G1 = 1.0;                                                                 \
    }                                                                           \
    /** If there is two dangling ends at the same end of duplex **/             \
    if (                                                                        \
        (bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1 ) &&                           \
        isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) &&  \
        isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])     \
    ) {                                                                         \
      S2 = (                                                                    \
        atPenaltyS(numSeq1[i], numSeq2[j]) +                                    \
        dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +              \
        dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]]                \
      );                                                                        \
      H2 = (                                                                    \
        atPenaltyH(numSeq1[i], numSeq2[j]) +                                    \
        dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +             \
        dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]]               \
      );                                                                        \
      G2 = H2 - TEMP_KELVIN*S2;                                                 \
      if(!isFinite(H2) || (G2 > 0)) {                                           \
        H2 = _INFINITY;                                                         \
        S2 = -1.0;                                                              \
        G2 = 1.0;                                                               \
      }                                                                         \
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);                        \
      if(isFinite(H1) && (G1 < 0)) {                                            \
        T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);                      \
        if( (T1 < T2) && (G2 < 0) ) {                                           \
          S1 = S2;                                                              \
          H1 = H2;                                                              \
          T1 = T2;                                                              \
        }                                                                       \
      } else if(G2 < 0) {                                                       \
        S1 = S2;                                                                \
        H1 = H2;                                                                \
        T1 = T2;                                                                \
      }                                                                         \
    } else if (                                                                 \
        (bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) &&                            \
        isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])     \
    ) {                                                                         \
      S2 = (                                                                    \
          atPenaltyS(numSeq1[i], numSeq2[j]) +                                  \
          dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]              \
      );                                                                        \
      H2 = (                                                                    \
        atPenaltyH(numSeq1[i], numSeq2[j]) +                                    \
        dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]               \
      );                                                                        \
      G2 = H2 - TEMP_KELVIN * S2;                                               \
      if(!isFinite(H2) || (G2 > 0)) {                                           \
        H2 = _INFINITY;                                                         \
        S2 = -1.0;                                                              \
        G2 = 1.0;                                                               \
      }                                                                         \
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);                        \
      if(isFinite(H1) && (G1 < 0)) {                                            \
        T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);                      \
        if((T1 < T2) && (G2 < 0)) {                                             \
          S1 = S2;                                                              \
          H1 = H2;                                                              \
          T1 = T2;                                                              \
        }                                                                       \
      } else if (G2 < 0){                                                       \
        S1 = S2;                                                                \
        H1 = H2;                                                                \
        T1 = T2;                                                                \
      }                                                                         \
    } else if (                                                                 \
        (bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) &&                            \
        isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])     \
    ) {                                                                         \
      S2 = (                                                                    \
        atPenaltyS(numSeq1[i], numSeq2[j]) +                                    \
        dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]]                \
      );                                                                        \
      H2 = (                                                                    \
        atPenaltyH(numSeq1[i], numSeq2[j]) +                                    \
        dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]]               \
      );                                                                        \
      G2 = H2 - TEMP_KELVIN * S2;                                               \
      if(!isFinite(H2) || (G2 > 0)) {                                           \
        H2 = _INFINITY;                                                         \
        S2 = -1.0;                                                              \
        G2 = 1.0;                                                               \
      }                                                                         \
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);                        \
                                                                                \
      if(isFinite(H1) && (G1 < 0)) {                                            \
        T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);                      \
                                                                                \
        if(T1 < T2  && (G2 < 0)) {                                              \
          S1 = S2;                                                              \
          H1 = H2;                                                              \
          T1 = T2;                                                              \
        }                                                                       \
                                                                                \
      } else if(G2 < 0) {                                                       \
        S1 = S2;                                                                \
        H1 = H2;                                                                \
        T1 = T2;                                                                \
      }                                                                         \
    }                                                                           \
    S2 = atPenaltyS(numSeq1[i], numSeq2[j]);                                    \
    H2 = atPenaltyH(numSeq1[i], numSeq2[j]);                                    \
    T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);                          \
    G1 = H1 -TEMP_KELVIN*S1;                                                    \
    G2 = H2 -TEMP_KELVIN*S2;                                                    \
    if (isFinite(H1)) {                                                         \
      if(T1 < T2) {                                                             \
        EntropyEnthalpy[0] = S2;                                                \
        EntropyEnthalpy[1] = H2;                                                \
      } else {                                                                  \
        EntropyEnthalpy[0] = S1;                                                \
        EntropyEnthalpy[1] = H1;                                                \
      }                                                                         \
    } else {                                                                    \
      EntropyEnthalpy[0] = S2;                                                  \
      EntropyEnthalpy[1] = H2;                                                  \
    }                                                                           \
  }                                                                             \
}


#define RSH(i, j, EntropyEnthalpy)                                              \
{                                                                               \
  double G1, G2;                                                                \
  double S1, S2;                                                                \
  double H1, H2;                                                                \
  double T1, T2;                                                                \
  S1 = S2 = -1.0;                                                               \
  H1 = H2 = _INFINITY;                                                          \
  T1 = T2 = -_INFINITY;                                                         \
  if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {                                    \
    EntropyEnthalpy[0] = -1.0;                                                  \
    EntropyEnthalpy[1] = _INFINITY;                                             \
  } else {                                                                      \
    S1 = (                                                                      \
      atPenaltyS(numSeq1[i], numSeq2[j]) +                                      \
      tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]]  \
    );                                                                          \
    H1 = (                                                                      \
      atPenaltyH(numSeq1[i], numSeq2[j]) +                                      \
      tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]] \
    );                                                                          \
    G1 = H1 - (TEMP_KELVIN * S1);                                               \
    if(!isFinite(H1) || (G1 > 0)) {                                             \
      H1 = _INFINITY;                                                           \
      S1 = -1.0;                                                                \
      G1 = 1.0;                                                                 \
    }                                                                           \
    if (                                                                        \
        (bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0) &&                            \
        isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) &&  \
        isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])     \
    ) {                                                                         \
      S2 = (                                                                    \
        atPenaltyS(numSeq1[i], numSeq2[j]) +                                    \
        dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +              \
        dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]]                \
      );                                                                        \
      H2 = (                                                                    \
        atPenaltyH(numSeq1[i], numSeq2[j]) +                                    \
        dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +             \
        dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]]               \
      );                                                                        \
      G2 = H2 - (TEMP_KELVIN * S2);                                             \
      if(!isFinite(H2) || G2>0) {                                               \
        H2 = _INFINITY;                                                         \
        S2 = -1.0;                                                              \
        G2 = 1.0;                                                               \
      }                                                                         \
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);                        \
      if(isFinite(H1) && (G1 < 0)) {                                            \
        T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);                      \
        if((T1 < T2) && (G2 < 0)) {                                             \
          S1 = S2;                                                              \
          H1 = H2;                                                              \
          T1 = T2;                                                              \
        }                                                                       \
      } else if (G2 < 0) {                                                      \
        S1 = S2;                                                                \
        H1 = H2;                                                                \
        T1 = T2;                                                                \
      }                                                                         \
    } else if (                                                                 \
        (bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0) &&                            \
        isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])     \
    ) {                                                                         \
      S2 = (                                                                    \
        atPenaltyS(numSeq1[i], numSeq2[j]) +                                    \
        dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]                \
      );                                                                        \
      H2 = (                                                                    \
        atPenaltyH(numSeq1[i], numSeq2[j]) +                                    \
        dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]               \
      );                                                                        \
      G2 = H2 - TEMP_KELVIN*S2;                                                 \
      if(!isFinite(H2) || (G2 > 0)) {                                           \
        H2 = _INFINITY;                                                         \
        S2 = -1.0;                                                              \
        G2 = 1.0;                                                               \
      }                                                                         \
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);                        \
      if(isFinite(H1) && (G1 < 0)) {                                            \
        T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);                      \
        if(T1 < T2 && (G2 < 0)) {                                               \
          S1 = S2;                                                              \
          H1 = H2;                                                              \
          T1 = T2;                                                              \
        }                                                                       \
      } else if (G2 < 0){                                                       \
        S1 = S2;                                                                \
        H1 = H2;                                                                \
        T1 = T2;                                                                \
      }                                                                         \
    } else if (                                                                 \
        (bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0) &&                            \
        isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])     \
    ) {                                                                         \
      S2 = (                                                                    \
        atPenaltyS(numSeq1[i], numSeq2[j]) +                                    \
        dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]]                \
      );                                                                        \
      H2 = (                                                                    \
        atPenaltyH(numSeq1[i], numSeq2[j]) +                                    \
        dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]]               \
      );                                                                        \
      G2 = H2 - TEMP_KELVIN*S2;                                                 \
      if (!isFinite(H2) || (G2 > 0)) {                                          \
        H2 = _INFINITY;                                                         \
        S2 = -1.0;                                                              \
        G2 = 1.0;                                                               \
      }                                                                         \
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);                        \
      if (isFinite(H1) && (G1 < 0)) {                                           \
        T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);                      \
        if ((T1 < T2) && (G2 < 0)) {                                            \
          S1 = S2;                                                              \
          H1 = H2;                                                              \
          T1 = T2;                                                              \
        }                                                                       \
      } else if (G2 < 0) {                                                      \
        S1 = S2;                                                                \
        H1 = H2;                                                                \
        T1 = T2;                                                                \
      }                                                                         \
    }                                                                           \
    S2 = atPenaltyS(numSeq1[i], numSeq2[j]);                                    \
    H2 = atPenaltyH(numSeq1[i], numSeq2[j]);                                    \
    T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);                          \
    G1 = H1 -TEMP_KELVIN*S1;                                                    \
    G2 =  H2 -TEMP_KELVIN*S2;                                                   \
    if(isFinite(H1)) {                                                          \
      if(T1 < T2) {                                                             \
        EntropyEnthalpy[0] = S2;                                                \
        EntropyEnthalpy[1] = H2;                                                \
      } else {                                                                  \
        EntropyEnthalpy[0] = S1;                                                \
        EntropyEnthalpy[1] = H1;                                                \
      }                                                                         \
    } else {                                                                    \
      EntropyEnthalpy[0] = S2;                                                  \
      EntropyEnthalpy[1] = H2;                                                  \
    }                                                                           \
  }                                                                             \
}


/* Finds max Tm while filling the dyn progr table using stacking S and stacking
* H (dimer)
*/
static void
maxTM(
  int i,
  int j,
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
  double T0, T1;
  double S0, S1;
  double H0, H1;
  double* SH;
  SH = (double*) safe_malloc(2 * sizeof(double), 0);
  T0 = T1 = -_INFINITY;
  S0 = EntropyDPT(i, j);
  H0 = EnthalpyDPT(i, j);
  RSH(i, j, SH);
  T0 = (H0 + dplx_init_H + SH[1]) /(S0 + dplx_init_S + SH[0] + RC);
  if (
      isFinite(EnthalpyDPT(i - 1, j - 1)) &&
      isFinite(Hs(i - 1, j - 1, 1, numSeq1, len1, numSeq2, len2))
  ) {
    S1 = (
      EntropyDPT(i - 1, j - 1) +
      Ss(i - 1, j - 1, 1, numSeq1, len1, numSeq2, len2)
    );
    H1 = (
      EnthalpyDPT(i - 1, j - 1) +
      Hs(i - 1, j - 1, 1, numSeq1, len1, numSeq2, len2)
    );
    T1 = (H1 + dplx_init_H + SH[1]) /(S1 + dplx_init_S + SH[0] + RC);
  } else {
    S1 = -1.0;
    H1 = _INFINITY;
    T1 = (H1 + dplx_init_H) /(S1 + dplx_init_S + RC);
  }
  if (S1 < MinEntropyCutoff) {
    /* to not give dH any value if dS is unreasonable */
    S1 = MinEntropy;
    H1 = 0.0;
  }
  if (S0 < MinEntropyCutoff) {
    /* to not give dH any value if dS is unreasonable */
    S0 = MinEntropy;
    H0 = 0.0;
  }
  if(T1 > T0) {
    EntropyDPT(i, j) = S1;
    EnthalpyDPT(i, j) = H1;
  } else if (T0 >= T1) {
    EntropyDPT(i, j) = S0;
    EnthalpyDPT(i, j) = H0;
  }
  free(SH);
}


/* Finds max Tm while filling the dyn progr table using stacking S and stacking
* H (monomer)
*/
static void
maxTM2(
  int i,
  int j,
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
  double T0, T1;
  double S0, S1;
  double H0, H1;
  T0 = T1 = -_INFINITY;
  S0 = EntropyDPT(i, j);
  H0 = EnthalpyDPT(i, j);
  T0 = (H0 + dplx_init_H) /(S0 + dplx_init_S + RC);
  if (isFinite(EnthalpyDPT(i, j))) {
    S1 = (
      EntropyDPT(i + 1, j - 1) +
      Ss(i, j, 2, numSeq1, len1, numSeq2, len2)
    );
    H1 = (
      EnthalpyDPT(i + 1, j - 1) +
      Hs(i, j, 2, numSeq1, len1, numSeq2, len2)
    );
  } else {
    S1 = -1.0;
    H1 = _INFINITY;
  }
  T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
  if (S1 < MinEntropyCutoff) {
    S1 = MinEntropy;
    H1 = 0.0;
  }
  if (S0 < MinEntropyCutoff) {
    S0 = MinEntropy;
    H0 = 0.0;
  }

  if (T1 > T0) {
    EntropyDPT(i, j) = S1;
    EnthalpyDPT(i, j) = H1;
  } else {
    EntropyDPT(i, j) = S0;
    EnthalpyDPT(i, j) = H0;
  }
}


static double
Ss(
  int i, int j, int k,
  unsigned char* numSeq1, int len1,
  unsigned char* numSeq2, int len2
) {
  if (k == 2) {
    if (i >= j) {
      return -1.0;
    }
    if ((i == len1) || (j == len2 + 1)) {
      return -1.0;
    }

    if (i > len1) {
      i -= len1;
    }
    if (j > len2) {
      j -= len2;
    }
    return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]];
  } else {
    return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
  }
}


static double
Hs(
    int i, int j, int k,
    unsigned char* numSeq1, int len1,
    unsigned char* numSeq2, int len2
) {
  if(k==2) {
    if (i >= j) {
      return _INFINITY;
    }
    if (i == len1 || j == len2 + 1) {
      return _INFINITY;
    }

    if (i > len1) {
      i -= len1;
    }
    if (j > len2) {
      j -= len2;
    }
    double check = stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j - 1]];
    if(isFinite(check)) {
      return check;
    } else {
      return _INFINITY;
    }
  } else {
    return stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
  }
}


/* Carries out Bulge and Internal loop and stack calculations to hairpin */
static void
CBI(
  int i, int j,
  double* EntropyEnthalpy,
  int traceback, int maxLoop,
  thalflex_packet_t* tf_packet
) {
  /* Unpack the tf_packet boiler plate variable definition */
  UNPACK_TFPACKET_CBI()

  int d, ii, jj;
  for (d = j - i - 3; d >= MIN_HRPN_LOOP + 1 && d >= j - i - 2 - maxLoop; --d) {
    for (ii = i + 1; ii < j - d && ii <= len1; ++ii) {
      jj = d + ii;
      if(traceback==0) {
        EntropyEnthalpy[0] = -1.0;
        EntropyEnthalpy[1] = _INFINITY;
      }

      if (isFinite(EnthalpyDPT(ii, jj)) && isFinite(EnthalpyDPT(i, j))) {
        calc_bulge_internal2(
            i, j, ii, jj,
            EntropyEnthalpy,
            traceback, maxLoop,
            numSeq1, len1,
            numSeq2, len2,
            len3,
            entropyDPT, enthalpyDPT,
            dplx_init_S, dplx_init_H, RC
        );
        if(isFinite(EntropyEnthalpy[1])) {
          if(EntropyEnthalpy[0] < MinEntropyCutoff) {
            EntropyEnthalpy[0] = MinEntropy;
            EntropyEnthalpy[1] = 0.0;
          }
          if(traceback == 0) {
            EnthalpyDPT(i, j) = EntropyEnthalpy[1];
            EntropyDPT(i, j) = EntropyEnthalpy[0];
          }
        }
      }
    }
  }
}


#endif
