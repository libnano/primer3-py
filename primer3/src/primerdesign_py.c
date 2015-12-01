// Copyright (C) 2014. Ben Pruitt & Nick Conway; Wyss Institute
// See LICENSE for full GPLv2 license.
// #
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// #
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// #
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

/*
primer3.primerdesign

Python-callable C API bindings for the Primer3 primer design engine.

These functions are wrapped in ``primer3.bindings`` to provide better
abstraction / function signatures, but may be used directly if the additional
function call overhead is of concern.

*/


#include    <Python.h>

#include    <stdio.h>
#include    <string.h> /* for NULL pointers */

#define PYTHON_BINDING

// Primer 3 includes (_mod includes are generated via patching)
#include    <thal.h>
#include    <oligotm.h>
#include    <libprimer3.h>

// Helper functions for parameter + output parsing
#include "primerdesign_helpers.h"

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* module doc string */
PyDoc_STRVAR(primerdesign__doc__, "Python C API bindings for the Primer3 "
                                  "design engine\n");

/* function doc strings */
PyDoc_STRVAR(loadThermoParams__doc__,
"Load the Primer3 thermodynamic parameters into memory.\n\n"
"Should only need to be called once on module import\n"
"path: path to the parameter directory"
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
loadThermoParams(PyObject *self, PyObject *args) {
    /* This loads the thermodynamic parameters from the parameter files
     * found at the provided path and should only need to be called once
     * prior to running thermodynamic calculations. Returns boolean indicating
     * success.
     */

    char            *param_path=NULL;
    thal_results    thalres;

    if (!PyArg_ParseTuple(args, "s", &param_path)) {
        return NULL;
    }

    if (get_thermodynamic_values(param_path, &thalres)){
        PyErr_SetString(PyExc_IOError, thalres.msg);
        return NULL;
    } else {
        Py_RETURN_TRUE;
    }
}


/* ~~~~~~~~~~~~~~~~~~~~~ PRIMER / OLIGO DESIGN BINDINGS ~~~~~~~~~~~~~~~~~~~~ */

p3_global_settings          *pa=NULL;
seq_args                    *sa=NULL;

static PyObject*
setGlobals(PyObject *self, PyObject *args){
    /* Sets the Primer3 global settings from a Python dictionary containing
    ** key: value pairs that correspond to the documented Primer3
    ** global parameters. Also accepts a mispriming or mishybridization library
    ** organized as `seq_name`:`seq_value` key:value pairs.
    */

    PyObject        *global_args=NULL, *misprime_lib=NULL;
    PyObject        *mishyb_lib=NULL;
    seq_lib         *mp_lib, *mh_lib;


    if (pa != NULL) {
        // Free memory for previous global settings
        p3_destroy_global_settings(pa);
        pa = NULL;
    }

    // Allocate memory for global settings
    if ((pa = p3_create_global_settings()) == NULL) {
        PyErr_SetString(PyExc_IOError,
                        "Could not allocate memory for p3 globals");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O!OO", &PyDict_Type, &global_args,
                          &misprime_lib, &mishyb_lib)) {
        goto err_set_global;
    }


    if ((pdh_setGlobals(pa, global_args)) != 0) {
        goto err_set_global;
    }

    if ((misprime_lib != NULL) && (misprime_lib != Py_None)) {
        if ((mp_lib = pdh_createSeqLib(misprime_lib)) == NULL) {
            goto err_set_global;
        }
        pa->p_args.repeat_lib = mp_lib;
    }
    if ((mishyb_lib != NULL) && (mishyb_lib != Py_None)) {
        if ((mh_lib = pdh_createSeqLib(mishyb_lib))==NULL) {
            goto err_set_global;
        }
        pa->o_args.repeat_lib = mh_lib;
    }

    Py_RETURN_NONE;

err_set_global:
    p3_destroy_global_settings(pa);
    pa = NULL;
    return NULL;
}

static PyObject*
setSeqArgs(PyObject *self, PyObject *args){
    /* Sets the Primer3 sequence args from a Python dictionary containing
    ** key: value pairs that correspond to the documented Primer3
    ** sequence parameters.
    */

    PyObject        *seq_args=NULL;

    if (pa == NULL) {
        PyErr_SetString(PyExc_IOError, "Primer3 global args must be \
            set prior to sequence args.");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &seq_args)) {
        return NULL;
    }


    if (sa != NULL) {
        // Free memory for previous seq args
        destroy_seq_args(sa);
        sa = NULL;
    }

    if ((sa = create_seq_arg()) == NULL) {
        PyErr_SetString(PyExc_IOError, "Could not allocate memory for seq_args");
        return NULL;
    }

    if (pdh_setSeqArgs(seq_args, sa) != 0) {
        destroy_seq_args(sa);
        sa = NULL;
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

    int             *debug=0;
    PyObject        *results=NULL;
    p3retval        *retval=NULL;

    if (pa == NULL || sa == NULL) {
        PyErr_SetString(PyExc_IOError, "Primer3 global args and sequence "
            "args must be set prior to calling runDesign.");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "i", &debug)) {
        return NULL;
    }

    if (debug) {
        p3_print_args(pa, sa);
    }

    retval = choose_primers(pa, sa);
    if ((results = pdh_outputToDict(pa, sa, retval)) == NULL){
        if (retval != NULL) {
            destroy_p3retval(retval);
        }
        return NULL;
    }

    destroy_p3retval(retval);
    retval = NULL;
    // Commented out for now (causes "malloc: *** error for object 0x101a03e20:
    // pointer being freed was not allocated") error
    destroy_dpal_thal_arg_holder();

    return results;
}

void
cleanUp(void){
    /* Free any remaining global Primer3 objects */
    if (pa != NULL) {
        // Free memory for previous global settings
        p3_destroy_global_settings(pa);
        pa = NULL;
    }

    if (sa != NULL) {
        // Free memory for previous seq args
        destroy_seq_args(sa);
        sa = NULL;
    }
    destroy_thal_structures();
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

static PyMethodDef primerdesign_methods[] = {
    { "loadThermoParams", loadThermoParams, METH_VARARGS,
     loadThermoParams__doc__ },
    { "setGlobals", setGlobals,  METH_VARARGS, setGlobals__doc__},
    { "setSeqArgs", setSeqArgs,  METH_VARARGS, setSeqArgs__doc__},
    { "runDesign", runDesign,  METH_VARARGS, runDesign__doc__},
    { NULL, NULL }
};

MOD_INIT(primerdesign){
#if PY_MAJOR_VERSION >= 3
    PyObject* m;
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "primerdesign",           /* m_name */
        primerdesign__doc__,      /* m_doc */
        -1,                 /* m_size */
        primerdesign_methods,     /* m_methods */
        NULL,               /* m_reload */
        NULL,               /* m_traverse */
        NULL,               /* m_clear */
        NULL,               /* m_free */
    };
    Py_AtExit(&cleanUp);
    m = PyModule_Create(&moduledef);
    return m;
#else
    Py_AtExit(&cleanUp);
    Py_InitModule3("primerdesign", primerdesign_methods, primerdesign__doc__);
#endif

}
