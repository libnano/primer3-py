/******************************************************************************
** primer3_py.c
** 
** Python-callable C API bindings for Primer3 (these functions are exposed
** to Python via the _primer3 module).
******************************************************************************/

#include    <Python.h>

#include    <stdio.h>
#include    <string.h> /* for NULL pointers */

// Primer 3 includes (_mod includes are generated via patching)
#include    <thal_mod.h>
#include    <oligotm.h>
#include    <libprimer3_mod.h>

// Helper functions for parameter + output parsing
#include "primer3_py_helpers.h"

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* module doc string */
PyDoc_STRVAR(primer3__doc__, "Python C API bindings for Primer3\n");

/* function doc strings */
PyDoc_STRVAR(getThermoParams__doc__,
"Load the Primer3 thermodynamic parameters into memory.\n\n"
"Should only need to be called once on module import\n"
"path: path to the parameter directory"
);

PyDoc_STRVAR(calcThermo__doc__,
"Compute the best thermodynamic alignment between two DNA sequences\n\n"
"seq1: sequence of the first oligo\n"
"seq2: sequence of the second oligo\n"
"align_type: alignment type (1: any, 2: end1, 3: end3, 4: hairpin)\n"
"mv_conc: concentration of monovalent cations\n"
"dv_conc: concentration of divalent cations\n"
"dntp_conc: concentration of dntps\n"
"dna_conc: concentration of oligonucleotides\n"
"temp: temperature at which hairpin structures will be calculated\n"
"max_loop: maximum size of a loop to consider (must be <=30)\n"
"temponly: calculate only the temperature\n"
"debug: debug printing turned on if true\n"
);

PyDoc_STRVAR(calcTm__doc__,
"Compute the melting temperature of a DNA sequence\n\n"
"seq: dna sequence\n"
"mv_conc: concentration of monovalent cations\n"
"dv_conc: concentration of divalent cations\n"
"dntp_conc: concentration of dntps\n"
"dna_conc: concentration of oligonucleotides\n"
"dna_conc: concentration of oligonucleotides\n"
"max_nn_length: max oligo length for which NN calculations will be used\n"
"tm_method: the method used for tm calculation (see oligotm.h)\n"
"salt_correction_method: the method used for salt corrections (see oligotm.h)\n"
);

PyDoc_STRVAR(setGlobals__doc__,
"Set the Primer3 global args and add a mispriming and/or mishyb libary\n\n"
"global_args: dictionary of Primer3 args\n"
"misprime_lib: mispriming library dictionary\n"
"mishyb_lib: mishybridization library dictionary\n"
);

PyDoc_STRVAR(setSeqArgs__doc__,
"Set the Primer3 sequence args\n\n"
"seq_args: dictionary of Primer3 sequence args\n"
);

PyDoc_STRVAR(runDesign__doc__,
"Design primers using the internal Primer3 design engine\n\n"
);



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HELPER FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

static PyObject*
getThermoParams(PyObject *self, PyObject *args) {
    /* This loads the thermodynamic parameters from the parameter files
     * found at the provided path and should only need to be called once
     * prior to running thermodynamic calculations. Returns boolean indicating
     * success.
     */

    char            *param_path=NULL;
    thal_results    o;

    if (!PyArg_ParseTuple(args, "s", &param_path)) {
        return NULL;
    }

    if (get_thermodynamic_values(param_path, &o)){
        PyErr_SetString(PyExc_IOError, o.msg);
        return NULL;
    } else {
        Py_RETURN_TRUE;
    }
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOW-LEVEL BINDINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~ */

static PyObject*
calcThermo(PyObject *self, PyObject *args){
    /* Wraps the main thal function in thal.h/thal.c. Arguments are parsed
     * into a thal_args struct (see thal.h).
     */

    char                    *oligo1=NULL;
    char                    *oligo2=NULL;
    thal_args               thalargs;
    thal_results            thalres;

    thalargs.dimer = 1; // this param is used by thal_main.c as a placeholder
                          // when determining whether the user has provided
                          // one or two sequences via the command line. if the
                          // user has only provided 1 sequence, then thal is
                          // run with the same sequence for oligo1 and oligo2

    // Initialize some thalres values to zero
    thalres.no_structure = 0;
    thalres.ds = thalres.dh = thalres.dg = 0.0;
    thalres.align_end_1 = thalres.align_end_2 = 0;

    if (!PyArg_ParseTuple(args, "ssidddddiii",
                          &oligo1, &oligo2, &thalargs.type, &thalargs.mv,
                          &thalargs.dv, &thalargs.dntp, &thalargs.dna_conc,
                          &thalargs.temp,  &thalargs.maxLoop,
                          &thalargs.temponly, &thalargs.debug)) {
        return NULL;
    }

    thal((const unsigned char *)oligo1, (const unsigned char *)oligo2,
         (const thal_args *)&thalargs, &thalres);

    return Py_BuildValue("siddddii", thalres.msg, thalres.no_structure,
                              thalres.temp, thalres.ds, thalres.dh,
                              thalres.dg, thalres.align_end_1,
                              thalres.align_end_2);
}

static PyObject*
calcTm(PyObject *self, PyObject *args){
    /* Wraps the seq function in oligotm.c/oligotm.h, which is used to
     * calculate the tm of a short sequence.
     */

    char            *oligo=NULL;
    double          dna_conc, mv_conc, dv_conc, dntp_conc, tm;
    int             max_nn_length, tm_method, salt_correction_method;

    if (!PyArg_ParseTuple(args, "sddddiii",
                          &oligo, &mv_conc, &dv_conc, &dntp_conc, &dna_conc,
                          &max_nn_length, &tm_method,
                          &salt_correction_method)) {
        return NULL;
    }

    tm = seqtm((const char *)oligo, dna_conc, mv_conc, dv_conc,
                 dntp_conc, max_nn_length, (tm_method_type)tm_method,
                 (salt_correction_type)salt_correction_method);

    return Py_BuildValue("d", tm);
}

/* ~~~~~~~~~~~~~~~~~~~~~ PRIMER / OLIGO DESIGN BINDINGS ~~~~~~~~~~~~~~~~~~~~ */

p3_global_settings           *pa=NULL;
seq_args                     *sa=NULL;

static PyObject*
setGlobals(PyObject *self, PyObject *args){
    /* Sets the Primer3 global settings from a Python dictionary containing
    ** key: value pairs that correspond to the documented Primer3 
    ** global parameters. Also accepts a mispriming or mishybridization library
    ** organized as `seq_name`:`seq_value` key:value pairs. 
    */

    PyObject                *global_args=NULL, *misprime_lib=NULL;
    PyObject                *mishyb_lib=NULL;
    seq_lib                 *mp_lib, *mh_lib;

    if (pa != NULL) {
        // Free memory for previous global settings
        p3_destroy_global_settings(pa);
    }
    // Allocate memory for global settings
    if (!(pa = (p3_global_settings *) malloc(sizeof(*pa)))) {
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O!OO", &PyDict_Type, &global_args, 
                          &misprime_lib, &mishyb_lib)) {
        return NULL;
    }

    if (misprime_lib != NULL && misprime_lib != Py_None) {
        if ((mp_lib = createSeqLib(misprime_lib)) == NULL) {
            return NULL;
        }
        pa->p_args.repeat_lib = mp_lib;
    }
    if (mishyb_lib != NULL && mishyb_lib != Py_None) {
        if ((mh_lib = createSeqLib(mishyb_lib))==NULL) {
            return NULL;
        }
        pa->p_args.repeat_lib = mh_lib;
    }

    if ((pa = _setGlobals(global_args)) == NULL){
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject*
setSeqArgs(PyObject *self, PyObject *args){
    /* Sets the Primer3 sequence args from a Python dictionary containing
    ** key: value pairs that correspond to the documented Primer3 
    ** sequence parameters.
    */

    PyObject                *seq_args=NULL;

    if (pa == NULL) {
        PyErr_SetString(PyExc_IOError, "Primer3 global args must be \
            set prior to sequence args.");
        return NULL;
    }

    if (sa != NULL) {
        // Free memory for previous seq args
        destroy_seq_args(sa);
    }

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &seq_args)) {
        return NULL;
    }

    if ((sa = _setSeqArgs(seq_args, pa))==NULL) {
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject*
runDesign(PyObject *self, PyObject *args){
    /* Wraps the primer design functionality of Primer3. Should be called 
    ** after setting the global and sequence-specific Primer3 parameters
    ** (see setGlobals and setSeqArgs, above)
     */

    PyObject                *results=NULL;
    p3retval                *retval=NULL;

    if (pa == NULL || sa == NULL) {
        PyErr_SetString(PyExc_IOError, "Primer3 global args and sequence\
            args must be set prior to calling runDesign.");
        return NULL;        
    }

    retval = choose_primers(pa, sa);
    if ((results = p3OutputToDict(pa, sa, retval)) == NULL){
        if (retval != NULL) {
            destroy_p3retval(retval);
        }
        return NULL;
    }

    destroy_p3retval(retval);
    // Commented out for now (causes "malloc: *** error for object 0x101a03e20:
    // pointer being freed was not allocated") error
    // destroy_dpal_thal_arg_holder();

    return results;
}

void
cleanUp(){
    /* Free any remaining global Primer3 objects */
    if (pa != NULL) {
        // Free memory for previous global settings
        p3_destroy_global_settings(pa);
    }

    if (sa != NULL) {
        // Free memory for previous seq args
        destroy_seq_args(sa);
    }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

static PyMethodDef primer3_methods[] = {
    { "getThermoParams", getThermoParams, METH_VARARGS, getThermoParams__doc__ },
    { "calcThermo", calcThermo, METH_VARARGS, calcThermo__doc__ },
    { "calcTm", calcTm, METH_VARARGS, calcTm__doc__ },
    { "setGlobals", setGlobals,  METH_VARARGS, setGlobals__doc__},
    { "setSeqArgs", setSeqArgs,  METH_VARARGS, setSeqArgs__doc__},
    { "runDesign", runDesign,  METH_VARARGS, runDesign__doc__},
    { NULL, NULL }
};

MOD_INIT(_primer3){
/* standard numpy compatible initiation */
#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_primer3",           /* m_name */
        primer3__doc__,      /* m_doc */
        -1,                 /* m_size */
        primer3_methods,     /* m_methods */
        NULL,               /* m_reload */
        NULL,               /* m_traverse */
        NULL,               /* m_clear */
        NULL,               /* m_free */
    };
    PyObject* m = PyModule_Create(&moduledef);
    return m;
#else
    Py_InitModule3("_primer3", primer3_methods, primer3__doc__);
#endif

Py_AtExit(&cleanUp);
}
