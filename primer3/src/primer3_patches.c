/******************************************************************************
** primer3_patches.c
** 
** This file contains C code patches for the primer3-2.3.6 distribution that
** allow for better python C API integration. See buildutil.py for a more
** thorough explanation of the patching system.
******************************************************************************/

//#FILE#src/thal.h##
//#BLOCK#thal_results##
//#PRESIG#/\* Structure for receiving results from the thermodynamic alignment calculation \*/\s+typedef struct##
//#POSTSIG#^} thal_results;\s+##
//#STARTBLOCK##
/* Structure for receiving results from the thermo alignment calculation */
typedef struct {
   char msg[255];
   int no_structure;    // Added no structure (1 if no structure found)
   double temp;
   double ds;           // Added entropy value
   double dh;           // Added enthalpy value
   double dg;           // Added gibbs free energy value
   int align_end_1;
   int align_end_2;
} thal_results;
//#ENDBLOCK##

//#FILE#src/thal.c##
//#BLOCK#drawDimer##
//#PRESIG#static void\s*drawDimer\(int\* ps1,##
//#POSTSIG#^}##
//#STARTBLOCK##
static void
drawDimer(int* ps1, int* ps2, double temp, double H, double S, int temponly, double t37, thal_results *o)
{
    int i, N;
    // char* duplex[4];
    double G, t;
    t = G = 0;
    if (!isFinite(temp)){
        if(temponly==0) {
        }
        o->temp = 0.0; /* lets use generalization here; this should rather be very negative value */
        o->no_structure = 1;
        strcpy(o->msg, "No predicted sec struc for given seq");
        return;
    } else {
        N=0;
        for(i=0;i<len1;i++){
            if(ps1[i]>0) ++N;
        }
        for(i=0;i<len2;i++) {
            if(ps2[i]>0) ++N;
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
        } else {
            o->temp = (double) t;
        }
    }
    return;
}
//#ENDBLOCK##

//#FILE#src/thal.c##
//#BLOCK#drawHairpin##
//#PRESIG#static void\s*drawHairpin\(int\* bp,##
//#POSTSIG#^}##
//#STARTBLOCK##
static void
drawHairpin(int* bp, double mh, double ms, int temponly, double temp, thal_results *o)
{
    int i, N;
    N = 0;
    double mg, t;
    if (!isFinite(ms) || !isFinite(mh)) {
        if(temponly == 0) {
            #ifdef DEBUG
                fputs("No temperature could be calculated\\n",stderr);
            #endif
        } else {
            o->temp = 0.0;
            o->no_structure = 1;
            strcpy(o->msg, "No predicted sec struc for given seq\\n");
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
        t = (mh / (ms + (((N/2)-1) * saltCorrection))) - ABSOLUTE_ZERO;
        if(temponly == 0) {
            mg = mh - (temp * (ms + (((N/2)-1) * saltCorrection)));
            ms = ms + (((N/2)-1) * saltCorrection);
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
//#ENDBLOCK##

//#FILE#src/thal.c##
//#BLOCK#NoSS##
//#REPLACE#fputs\("No secondary structure could be calculated\\n",stderr\);##
//#STARTBLOCK##
o->no_structure = 1;
//#ENDBLOCK##


//#FILE#src/p3_seq_lib.h##
//#BLOCK#create_seq_lib##
//#REPLACE##endif##
//#STARTBLOCK##
seq_lib *create_empty_seq_lib();

#endif
//#ENDBLOCK##

//#FILE#src/libprimer3.cpp##
//#BLOCK#c++ extern##
//#REPLACE##include "dpal.h"\s+#include "thal.h"\s+#include "oligotm.h"\s+#include "libprimer3.h"##
//#STARTBLOCK##
extern "C" {
    #include "dpal.h"
    #include "thal.h"
    #include "oligotm.h"
    #include "libprimer3.h"
}
//#ENDBLOCK##
