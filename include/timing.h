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

#ifndef RACE_TIMING_H
#define RACE_TIMING_H

#include <sys/time.h>
#include "print.h"

#define TIMING

#ifdef TIMING
#define START_TIME(region)\
timeval start##region, end##region;\
double start_tym##region, end_tym##region;\
gettimeofday(&start##region, NULL); \
start_tym##region = start##region.tv_sec + start##region.tv_usec*1e-6;\

#define STOP_TIME(region)\
gettimeofday(&end##region, NULL);\
end_tym##region = end##region.tv_sec + end##region.tv_usec*1e-6;\

#define PRINT_TIME(region)\
TIME_PRINT("%s : Time = %f", #region, end_tym##region - start_tym##region);\

#define GET_TIME(region)\
(end_tym##region - start_tym##region);\


#else
#define START_TIME(region)\

#define STOP_TIME(region)\

#define PRINT_TIME(region)\

#define GET_TIME(region)\

#endif

#endif
