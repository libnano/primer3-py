typedef struct {
    PyObject_HEAD
    thal_results thalres;
} ThermoResult;

static void
ThermoResult_dealloc(ThermoResult* self) {
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
ThermoResult_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    ThermoResult *self;
    self = (ThermoResult *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static int
ThermoResult_init(ThermoResult *self, PyObject *args, PyObject *kwds) {
    if (self == NULL) {
        return -1;
    }
    self->thalres.no_structure = 0;
    self->thalres.ds = self->thalres.dh = self->thalres.dg = 0.0;
    self->thalres.align_end_1 = self->thalres.align_end_2 = 0;
    return 0;
}

static PyMemberDef ThermoResult_members[] = {
    { "temp", T_DOUBLE, 
        offsetof(ThermoResult, thalres) + offsetof(thal_results, temp),
         0, "temperature"},
    { "ds", T_DOUBLE, 
        offsetof(ThermoResult, thalres) + offsetof(thal_results, ds),
        0, "ds"},
    { "dh", T_DOUBLE, 
        offsetof(ThermoResult, thalres) + offsetof(thal_results, dh),
        0, "dh"},
    { "dg", T_DOUBLE, 
        offsetof(ThermoResult, thalres) + offsetof(thal_results, dg),
        0, "dg"},
    { "align_end_1", T_INT, 
        offsetof(ThermoResult, thalres) + offsetof(thal_results, align_end_1),
        0, "align_end_1"},
    { "align_end_2", T_INT, 
        offsetof(ThermoResult, thalres) + offsetof(thal_results, align_end_2),
        0, "align_end_2"},
    { "no_structure", T_INT, 
        offsetof(ThermoResult, thalres) + offsetof(thal_results, no_structure),
        0, "no structure"},
    {NULL}  /* Sentinel */
};

static PyGetSetDef ThermoResult_getsetters[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef ThermoResult_methods[] = {
    {NULL}  /* Sentinel */
};

static PyTypeObject ThermoResultType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "analysis.ThermoResult",                        /*tp_name*/
    sizeof(ThermoResult),                           /*tp_basicsize*/
    0,                                              /*tp_itemsize*/
    (destructor)ThermoResult_dealloc,               /*tp_dealloc*/
    0,                                              /*tp_print*/
    0,                                              /*tp_getattr*/
    0,                                              /*tp_setattr*/
    0,                                              /*tp_compare*/
    0,                                              /*tp_repr*/
    0,                                              /*tp_as_number*/
    0,                                              /*tp_as_sequence*/
    0,                                              /*tp_as_mapping*/
    0,                                              /*tp_hash */
    0,                                              /*tp_call*/
    0,                                              /*tp_str*/
    0,                                              /*tp_getattro*/
    0,                                              /*tp_setattro*/
    0,                                              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,       /*tp_flags*/
    "ThermoResult objects",                         /* tp_doc */
    0,                                              /* tp_traverse */
    0,                                              /* tp_clear */
    0,                                              /* tp_richcompare */
    0,                                              /* tp_weaklistoffset */
    0,                                              /* tp_iter */
    0,                                              /* tp_iternext */
    ThermoResult_methods,                           /* tp_methods */
    ThermoResult_members,                           /* tp_members */
    ThermoResult_getsetters,                        /* tp_getset */
    0,                                              /* tp_base */
    0,                                              /* tp_dict */
    0,                                              /* tp_descr_get */
    0,                                              /* tp_descr_set */
    0,                                              /* tp_dictoffset */
    (initproc)ThermoResult_init,                    /* tp_init */
    0,                                              /* tp_alloc */
    ThermoResult_new,                               /* tp_new */
};
