#include <Python.h>
#include "tools.hpp" // header file for libfoo.a

static PyObject* py_foo(PyObject* self, PyObject* args) {
    // Parse arguments
    int arg1;
    if (!PyArg_ParseTuple(args, "i", &arg1))
        return NULL;

    // Call libfoo.a function
    int result = foo(arg1);

    // Build return value
    PyObject* ret = Py_BuildValue("i", result);
    return ret;
}

static PyMethodDef methods[] = {
    {"foo", py_foo, METH_VARARGS, "Call the foo() function from libfoo.a"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mylib", // name of the module
    NULL,
    -1,
    methods // methods table
};

PyMODINIT_FUNC PyInit_mylib(void) {
    return PyModule_Create(&moduledef);
}
