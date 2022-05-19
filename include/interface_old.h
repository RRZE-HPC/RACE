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

#ifndef RACE_INTERFACE_H
#define RACE_INTERFACE_H
#include "graph.h"
#include "traverse.h"
#include "lb.h"
#include "type.h"
#include "levelData.h"

class RACEInterface{
    private:
        int nrow;
        RACE::dist dist;
        int requestedThreads;
        int availableThreads;
        int *initPerm;
        int *initInvPerm;
        int *rowPtr;
        int *col;
        int *zonePtr;
        int zonePtrLen;
        int *perm;
        int permLen;
        int *invPerm;
        int invPermLen;

    public:
        RACEInterface(int nrow_, int nthreads_, RACE::dist dist_, int *rowPtr_, int *col_, int *initPerm_=NULL, int *initInvPerm_=NULL);
        void RACEColor();
        void getZonePtr(int **zonePtr_, int *len_);
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getNumThreads();
};

#endif
