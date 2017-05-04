#include "error.h"

char const * NAME_error_string(NAME_error e)
{
    char const *ret;
    switch (e) {
        case NAME_ERR_INVALID_ARG:
            ret = "Invalid argument";
            break;
        case NAME_ERR_MATRIX_SYMM:
            ret = "Matrix not symmetric";
            break;
        case NAME_ERR_D2_COLOR:
            ret = "Conflict in D2 coloring";
            break;
        case NAME_ERR_D1_COLOR:
            ret = "Conflict in D1 coloring";
            break;
        case NAME_ERR_INCOMPATIBILITY:
            ret = "INCOMPATIBILITY ERROR";
            break;
        case NAME_ERR_HWLOC:
            ret = "HWLOC ERROR";
            break;
        case NAME_ERR_NOT_IMPLEMENTED:
            ret = "NOT IMPLEMENTED";
            break;
        default:
            ret = "Invalid";
            break;
    }

   return ret;
}


