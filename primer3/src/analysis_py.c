#include <Python.h>
#include "structmember.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for NULL pointers */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

// Primer 3 includes (_mod includes are generated via patching)
#include    <thal.h>
#include    <oligotm.h>

// Helper functions for parameter + output parsing
#include "primer3_py_helpers.h"
#include "thermoresult_py.h"

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif

#if PY_MAJOR_VERSION >= 3
#define PYSTR_GET_STR_AND_SIZE(obj, buf, len) buf = PyUnicode_AsUTF8AndSize(obj, (&len))
#else
#define PYSTR_GET_STR_AND_SIZE(obj, buf, len) buf = PyString_AsString(obj)
#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C API functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

PyDoc_STRVAR(analysis__doc__,
              "Analyze many sequences using fixed input parameters\n");

typedef struct {
    PyObject_HEAD
    
    thal_args               thalargs;
    // for melting temperature
    int max_nn_length;
    int tm_method;
    int salt_correction_method;

} ThermoAnalysis;

static void
ThermoAnalysis_dealloc(ThermoAnalysis* self) {
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
ThermoAnalysis_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    ThermoAnalysis *self;
    self = (ThermoAnalysis *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static int Analyis_setParameters(ThermoAnalysis *self, PyObject *args, PyObject *kwds) {

    static char *kwlist[] = {   "type",
                                "mv_conc",
                                "dv_conc",
                                "dntp_conc",
                                "dna_conc",
                                "temp_c",
                                "max_loop",
                                "temp_only",
                                "debug",
                                "max_nn_length",
                                "tm_method",
                                "salt_correction_method",
                                NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|idddddiiiiii", kwlist, 
                            &(self->thalargs.type),
                            &(self->thalargs.mv),
                            &(self->thalargs.dv),
                            &(self->thalargs.dntp),
                            &(self->thalargs.dna_conc),
                            &(self->thalargs.temp),
                            &(self->thalargs.maxLoop),
                            &(self->thalargs.temponly),
                            &(self->thalargs.debug),
                            &(self->max_nn_length),
                            &(self->tm_method),
                            &(self->salt_correction_method)
                            )) {
        return -1;
    } else {
        return 0;
    }
}

static int
ThermoAnalysis_init(ThermoAnalysis *self, PyObject *args, PyObject *kwds) {
    if (self == NULL) {
        return -1;
    }
    // set defaults
    
    /*
    1 THAL_ANY, (by default)
    2 THAL_END1,
    3 THAL_END2,
    4 THAL_HAIRPIN */
    self->thalargs.type = 1;

    self->thalargs.mv = 50;
    self->thalargs.dv = 0;
    self->thalargs.dntp = 0.8;
    self->thalargs.dna_conc = 50;
    self->thalargs.temp = 37 + 273.15;
    self->thalargs.maxLoop = 30;
    self->thalargs.temponly = 0;
    self->thalargs.debug = 0;

    self->max_nn_length = 60;

    /*
    'breslauer': 0,
    'santalucia': 1
    */
    self->tm_method = 1;

    /*
    'schildkraut': 0,
    'santalucia': 1,
    'owczarzy': 2
    */
    self->salt_correction_method = 1;

    if (Analyis_setParameters(self, args, kwds) < 0) {
        return -1;
    }
    return 0;
}

PyDoc_STRVAR(heterodimer__doc__,
              "heterodimer\n");

static PyObject*
ThermoAnalysis_heterodimer(ThermoAnalysis *self, PyObject *args) {
    /* Wraps the main thal function in thal.h/thal.c. Arguments are parsed
     * into a thal_args struct (see thal.h).
     */

    char                    *oligo1=NULL, *oligo2=NULL;
    int                     oligo1_len, oligo2_len;
    ThermoResult *tr_obj = NULL;
    tr_obj = PyObject_New(ThermoResult, &ThermoResultType);
    if (tr_obj == NULL) {
        return NULL;
    }

    self->thalargs.dimer = 1;   // this param is used by thal_main.c as a placeholder
                                // when determining whether the user has provided
                                // one or two sequences via the command line. if the
                                // user has only provided 1 sequence, then thal is
                                // run with the same sequence for oligo1 and oligo2
    self->thalargs.type = 1;    // ANY

    if (!PyArg_ParseTuple(args, "s#s#",
                          &oligo1, &oligo1_len, &oligo2, &oligo2_len)) {
        return NULL;
    }

    if (oligo1_len > 60 && oligo2_len > 60) { 
        return PyErr_Format(PyExc_ValueError, 
                            "Only one input sequence may have a length > 60 "
                            "(lengths are %d and %d, respectively)", oligo1_len,
                            oligo2_len);
    }

    thal((const unsigned char *) oligo1, (const unsigned char *) oligo2,
         (const thal_args *) &(self->thalargs), &(tr_obj->thalres), 0);

    return (PyObject*) tr_obj;
}


PyDoc_STRVAR(homodimer__doc__,
              "homodimer\n");

static PyObject*
ThermoAnalysis_homodimer(ThermoAnalysis *self, PyObject *args) {
    /* Wraps the main thal function in thal.h/thal.c. Arguments are parsed
     * into a thal_args struct (see thal.h).
     */

    char                    *oligo1=NULL;
    int                     oligo1_len;
    ThermoResult *tr_obj = NULL;
    tr_obj = PyObject_New(ThermoResult, &ThermoResultType);
    if (tr_obj == NULL) {
        return NULL;
    }

    self->thalargs.dimer = 1;   // this param is used by thal_main.c as a placeholder
                                // when determining whether the user has provided
                                // one or two sequences via the command line. if the
                                // user has only provided 1 sequence, then thal is
                                // run with the same sequence for oligo1 and oligo2
    self->thalargs.type = 1;    // ANY

    if (!PyArg_ParseTuple(args, "s#",
                          &oligo1, &oligo1_len)) {
        return NULL;
    }

    if (oligo1_len > 60) { 
        return PyErr_Format(PyExc_ValueError, 
                            "Input sequence may not have a length > 60 "
                            "(length is %d respectively)", oligo1_len);
    }

    thal((const unsigned char *) oligo1, (const unsigned char *) oligo1,
         (const thal_args *) &(self->thalargs), &(tr_obj->thalres), 0);

    return (PyObject*) tr_obj;
}

PyDoc_STRVAR(hairpin__doc__,
              "hairpin\n");

static PyObject*
ThermoAnalysis_hairpin(ThermoAnalysis *self, PyObject *args) {
    /* Wraps the main thal function in thal.h/thal.c. Arguments are parsed
     * into a thal_args struct (see thal.h).
     */

    char                    *oligo1=NULL;
    int                     oligo1_len;
    ThermoResult *tr_obj = NULL;
    tr_obj = PyObject_New(ThermoResult, &ThermoResultType);
    if (tr_obj == NULL) {
        return NULL;
    }

    self->thalargs.dimer = 1;   // this param is used by thal_main.c as a placeholder
                                // when determining whether the user has provided
                                // one or two sequences via the command line. if the
                                // user has only provided 1 sequence, then thal is
                                // run with the same sequence for oligo1 and oligo2
    self->thalargs.type = 4;    // THAL_HAIRPIN

    if (!PyArg_ParseTuple(args, "s#",
                          &oligo1, &oligo1_len)) {
        return NULL;
    }

    if (oligo1_len > 60) { 
        return PyErr_Format(PyExc_ValueError, 
                            "Input sequence may not have a length > 60 "
                            "(length is %d respectively)", oligo1_len);
    }

    thal((const unsigned char *) oligo1, (const unsigned char *) oligo1,
     (const thal_args *) &(self->thalargs), &(tr_obj->thalres), 0);

    return (PyObject*) tr_obj;
}

PyDoc_STRVAR(meltingTemp__doc__,
              "meltingTemp\n");

static PyObject*
ThermoAnalysis_meltingTemp(ThermoAnalysis *self, PyObject *args) {
    /* Wraps the seq function in oligotm.c/oligotm.h, which is used to
     * calculate the tm of a short sequence.
     */

    char            *oligo=NULL;
    double          tm;
    int             oligo_len;

    if (!PyArg_ParseTuple(args, "s#",
                          &oligo, &oligo_len)) {
        return NULL;
    }

    thal_args ta = self->thalargs;
    tm = seqtm((const char *)oligo, ta.dna_conc, ta.mv, ta.dv,
                 ta.dntp, self->max_nn_length, (tm_method_type)self->tm_method,
                 (salt_correction_type)self->salt_correction_method);

    return PyFloat_FromDouble(tm);
}

static PyMemberDef ThermoAnalysis_members[] = {
    { "debug", T_INT, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, debug),
         0, "debug"},
    { "type", T_INT, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, type),
        0, "alignment type"},
    { "max_loop", T_INT, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, maxLoop),
        0, "max loop"},
    { "mv_conc", T_DOUBLE, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, mv),
        0, "mv concentration"},
    { "dv_conc", T_DOUBLE, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, dv),
        0, "dv_ concentration"},
    { "dntp_conc", T_DOUBLE, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, dntp),
        0, "dntp concentration"},
    { "dna_conc", T_DOUBLE, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, dna_conc),
        0, "dna concentration"},
    { "temp", T_DOUBLE, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, temp),
        0, "temperature"},
    { "dimer", T_INT, 
        offsetof(ThermoAnalysis, thalargs) + offsetof(thal_args, dimer),
        0, "dimer"},
    { "max_nn_length", T_INT, 
        offsetof(ThermoAnalysis, max_nn_length),
        0, "max_nn_length"},
    { "tm_method", T_INT, 
        offsetof(ThermoAnalysis, tm_method),
        0, "tm_method"},
    { "salt_correction_method", T_INT, 
        offsetof(ThermoAnalysis, salt_correction_method),
        0, "salt_correction_method"},
    {NULL}  /* Sentinel */
};

static PyGetSetDef ThermoAnalysis_getsetters[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef ThermoAnalysis_methods[] = {
    { "heterodimer", (PyCFunction) ThermoAnalysis_heterodimer, METH_VARARGS, heterodimer__doc__ },
    { "homodimer", (PyCFunction) ThermoAnalysis_homodimer, METH_VARARGS, homodimer__doc__ },
    { "hairpin", (PyCFunction) ThermoAnalysis_hairpin, METH_VARARGS, hairpin__doc__ },
    { "meltingTemp", (PyCFunction) ThermoAnalysis_meltingTemp, METH_VARARGS, meltingTemp__doc__ },
    {NULL}  /* Sentinel */
};

static PyTypeObject ThermoAnalysisType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "analysis.ThermoAnalysis",            /*tp_name*/
    sizeof(ThermoAnalysis),               /*tp_basicsize*/
    0,                              /*tp_itemsize*/
    (destructor)ThermoAnalysis_dealloc,   /*tp_dealloc*/
    0,                              /*tp_print*/
    0,                              /*tp_getattr*/
    0,                              /*tp_setattr*/
    0,                              /*tp_compare*/
    0,                              /*tp_repr*/
    0,                              /*tp_as_number*/
    0,                              /*tp_as_sequence*/
    0,                              /*tp_as_mapping*/
    0,                              /*tp_hash */
    0,                              /*tp_call*/
    0,                              /*tp_str*/
    0,                              /*tp_getattro*/
    0,                              /*tp_setattro*/
    0,                              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "ThermoAnalysis objects",             /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    ThermoAnalysis_methods,               /* tp_methods */
    ThermoAnalysis_members,               /* tp_members */
    ThermoAnalysis_getsetters,            /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc)ThermoAnalysis_init,        /* tp_init */
    0,                              /* tp_alloc */
    ThermoAnalysis_new,                   /* tp_new */
};

static PyMethodDef analysis_mod_methods[] = {
    {NULL}
};


static PyObject *setUp(void) {
    /* This loads the thermodynamic parameters from the parameter files
     * found at the provided path and should only need to be called once
     * prior to running thermodynamic calculations. Returns boolean indicating
     * success.
     */

    char            *param_path=NULL;
    thal_results    thalres;

    char* p3path = getenv("PRIMER3HOME");
    if (p3path == NULL) {
        return NULL;
    }

    param_path = (char *) calloc( (strlen(p3path)+20), sizeof(char));
    if (param_path == NULL) {
        return NULL;
    }
    sprintf(param_path, "%s/%s/", p3path, "primer3_config");

    if (get_thermodynamic_values(param_path, &thalres)){
        PyErr_SetString(PyExc_IOError, thalres.msg);
        free(param_path);
        return NULL;
    } else {
        free(param_path);
        Py_RETURN_TRUE;
    }
}

// void
// cleanUp(void){
//     /* Free any remaining global Primer3 objects */
//     if (pa != NULL) {
//         // Free memory for previous global settings
//         p3_destroy_global_settings(pa);
//     }

//     if (sa != NULL) {
//         // Free memory for previous seq args
//         destroy_seq_args(sa);
//     }
// }

MOD_INIT(analysis) {
    if (PyType_Ready(&ThermoAnalysisType) < 0) {
    #if PY_MAJOR_VERSION >= 3
        return NULL;
    #else
        return;
    #endif
    }

    if (PyType_Ready(&ThermoResultType) < 0) {
    #if PY_MAJOR_VERSION >= 3
        return NULL;
    #else
        return;
    #endif
    }

    if (setUp() == NULL) {
    #if PY_MAJOR_VERSION >= 3
        return NULL;
    #else
        return;
    #endif
    }
    #if PY_MAJOR_VERSION >= 3
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "analysis",            /* m_name */
            analysis__doc__,       /* m_doc */
            -1,                     /* m_size */
            analysis_mod_methods,  /* m_methods */
            NULL,                   /* m_reload */
            NULL,                   /* m_traverse */
            NULL,                   /* m_clear */
            NULL,                   /* m_free */
        };
        PyObject* m = PyModule_Create(&moduledef);
        import_array();
        if (m == NULL) { return NULL; }

        Py_INCREF(&ThermoAnalysisType);
        PyModule_AddObject(m, "ThermoAnalysis", (PyObject *)&ThermoAnalysisType);

        Py_INCREF(&ThermoResultType);
        PyModule_AddObject(m, "ThermoResult", (PyObject *)&ThermoResultType);
        // Py_AtExit(&cleanUp);
        return m;
    #else
        // Py_AtExit(&cleanUp);
        PyObject* m = Py_InitModule3("analysis", analysis_mod_methods, analysis__doc__);
        if (m == NULL) { return; }
        import_array();
        
        Py_INCREF(&ThermoAnalysisType);
        PyModule_AddObject(m, "ThermoAnalysis", (PyObject *)&ThermoAnalysisType);

        Py_INCREF(&ThermoResultType);
        PyModule_AddObject(m, "ThermoResult", (PyObject *)&ThermoResultType);
    #endif
};






