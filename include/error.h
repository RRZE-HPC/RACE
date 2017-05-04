#ifndef NAME_ERROR_H
#define NAME_ERROR_H

#include "print.h"

enum NAME_error{
    NAME_SUCCESS,
    NAME_ERR_INVALID_ARG,
    NAME_ERR_MATRIX_SYMM,
    NAME_ERR_D2_COLOR,
    NAME_ERR_D1_COLOR,
    NAME_ERR_GRAPH_TRAVERSAL,
    NAME_ERR_INCOMPATIBILITY,
    NAME_ERR_HWLOC,
    NAME_ERR_NOT_IMPLEMENTED
};

char const* NAME_error_string(NAME_error e);

#define NAME_FN(call) {\
    NAME_error ret = NAME_SUCCESS;\
    ret = call;\
    if (ret != NAME_SUCCESS) {\
        PRINT(NAME_ERROR,ANSI_COLOR_RED,"%s",NAME_error_string((NAME_error)ret));\
    }\
}\


#endif
