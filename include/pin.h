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

#ifndef RACE_PIN_H
#define RACE_PIN_H

#include "zone_tree.h"
#include "error.h"
#include "type.h"
#include "thpool.h"
#include <omp.h>

class Pin{
    private:
        ZoneTree *zoneTree;
        std::vector<std::pair<int,int>> pinMap;
        void* machine; //this is of type Machine, but casting to void to avoid external library dependency
        int SMT;
        RACE::PinMethod method;
        void pinOrderRecursive(int parentIdx);
        void calcPinOrder();
        RACE_error createPuNodeMapping();
        void getNodeId(int cpuId);
        //For pthread
        //void createPinnedThreadPool(int parentIdx);
        //For OMP: currently broken
        RACE_error pinApplicationRecursive(int parentIdx);
    public:
        Pin(ZoneTree* zoneTree_, int SMT_, RACE::PinMethod method_);
        RACE_error pinInit();
        RACE_error pinThread(int pinOrder);
        RACE_error pinApplication();
        void pinPowerThread(int nodes);
        void resetMaster();
};

#endif
