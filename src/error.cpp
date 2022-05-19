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

#include "error.h"

char const * RACE_error_string(RACE_error e)
{
    char const *ret;
    switch (e) {
        case RACE_ERR_INVALID_ARG:
            ret = "Invalid argument";
            break;
        case RACE_ERR_MATRIX_SYMM:
            ret = "Matrix not symmetric";
            break;
        case RACE_ERR_D2_COLOR:
            ret = "Conflict in D2 coloring";
            break;
        case RACE_ERR_D1_COLOR:
            ret = "Conflict in D1 coloring";
            break;
        case RACE_ERR_INCOMPATIBILITY:
            ret = "INCOMPATIBILITY ERROR";
            break;
        case RACE_ERR_HWLOC:
            ret = "HWLOC ERROR";
            break;
        case RACE_ERR_PIN:
            ret = "PIN ERROR";
            break;
        case RACE_ERR_NOT_IMPLEMENTED:
            ret = "NOT IMPLEMENTED";
            break;
        default:
            ret = "Invalid";
            break;
    }

   return ret;
}


