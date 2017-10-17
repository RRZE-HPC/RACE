#ifndef RACE_ERROR_H
#define RACE_ERROR_H

#include "print.h"

enum RACE_error{
    RACE_SUCCESS,
    RACE_ERR_INVALID_ARG,
    RACE_ERR_MATRIX_SYMM,
    RACE_ERR_D2_COLOR,
    RACE_ERR_D1_COLOR,
    RACE_ERR_GRAPH_TRAVERSAL,
    RACE_ERR_INCOMPATIBILITY,
    RACE_ERR_HWLOC,
    RACE_ERR_NOT_IMPLEMENTED
};

char const* RACE_error_string(RACE_error e);

#define RACE_FN(call) {\
    RACE_error ret = RACE_SUCCESS;\
    ret = call;\
    if (ret != RACE_SUCCESS) {\
        PRINT(RACE_ERROR,ANSI_COLOR_RED,"%s",RACE_error_string((RACE_error)ret));\
    }\
}\


#endif
