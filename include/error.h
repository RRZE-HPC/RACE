/*
 * =======================================================================================
 *
 *   RACE: Recursicve Algebraic Coloring Engine
 *   Copyright (C) 2019, RRZE, Friedrich-Alexander-Universität Erlangen-Nürnberg
 *   Author: Christie Alappat
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * =======================================================================================
 */

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
    RACE_ERR_PIN,
    RACE_ERR_INTERNAL,
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
