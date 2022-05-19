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

//if GAP is there use optimised version
#include "config.h"

#ifdef RACE_USE_GAP
    #ifdef RACE_USE_SOA_GRAPH
        #include "graph_SoA.cpp"
        #pragma message ("SoA graph being used")
    #else
        #include "graph_AoS.cpp"
        #pragma message ("AoS graph being used")
    #endif
#else
    #ifdef RACE_USE_SOA_GRAPH
        #pragma message ("SoA graph not supported with serial BFS. Switch on RACE_USE_GAP. Now AoS graph being used")
    #else
        #pragma message ("AoS graph being used")
    #endif
    #include "graph_AoS.cpp"
#endif
