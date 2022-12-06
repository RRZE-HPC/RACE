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

#ifndef RACE_FUNCTION_MANAGER_H
#define RACE_FUNCTION_MANAGER_H

#include <functional>
#include "zone_tree.h"
#include "level_recursion.h"
#include "omp.h"
#include "level_pool.h"
#include "timing.h"
#include "matrixPowerRecursive.h"

class FuncManager;

void recursiveCall(FuncManager* funMan, int parent);

void recursivePowerCall(FuncManager* funMan, int parent);

class FuncManager
{
    private:
        bool rev;
        typedef std::function< void(int,int,void *) > funcType;
        typedef std::function< void(int,int,int,int,int,void *) > powerFuncType;
        typedef std::function< void(void *) > commFuncType;
        bool power_fn;
        bool numaSplit;
        funcType func;
        powerFuncType powerFunc;
        commFuncType commFunc;
        void* args;
        void* commArgs;
        //totPower = power*subPower
        int totPower;
        int power;
        int subPower;
        int activethreads;
        int origthreads;
        int totalNodes;
        int threadPerNode;
        int CL_pad;
        ZoneTree* zoneTree;
        LevelPool* pool;
        mtxPowerRecursive* matPower;
        std::vector<int> serialPart;
        volatile int* nodeBarrier;
        volatile int* barrierCount;
        std::vector<volatile int*> lockCtr;
        std::vector<volatile bool*> lockTable;
        std::vector<volatile int*> lockTableCtr;
        //int* unlockRow;
        //int* dangerRow;
        //int* unlockCtr;
        friend void recursiveCall(FuncManager* funMan, int parent);
        friend void recursivePowerCall(FuncManager* funMan, int parent);
        std::function<void(int)> recursiveFun;

        void initFuncColor();
        void initFuncPower();
        void recursivePowerCallSerial(int parent);
        inline void powerCallGeneral(int startLevel, int endLevel, int boundaryStart, int boundaryEnd, int startSlope, int endSlope, const std::vector<int> *levelPtr, const std::vector<std::map<int, std::vector<std::vector<int>>>>* boundaryLevelPtr, const std::vector<int> *unlockRow, const std::vector<std::map<int, std::vector<std::vector<int>>>>* boundaryUnlockRow, const std::vector<int> *unlockCtr, const std::vector<int> *dangerRow, const std::vector<std::map<int, std::vector<std::vector<int>>>>* boundaryDangerRow, int numaLocalArg, int offset, int parent);
        inline void powerCallNodeReminder(int startSlope, int endSlope, const std::vector<int> *levelPtr, const std::vector<int> *nodePtr, int numaLocalArg, int offset);
        inline void powerCallHopelessRightReminder(int leftmostLevel, const std::vector<int> *levelPtr, const std::vector<std::map<int, std::vector<std::vector<int>>>>* boundaryLevelPtr, const std::vector<int> *unlockRow, const std::vector<int> *unlockCtr, const std::vector<int> *dangerRow, int numaLocalArg, int offset, int parent);
        inline void powerCallHopelessLeftReminder(int rightmostLevel, const std::vector<int> *levelPtr, const std::vector<std::map<int, std::vector<std::vector<int>>>>* boundaryLevelPtr, const std::vector<int> *unlockRow, const std::vector<int> *unlockCtr, const std::vector<int> *dangerRow, int numaLocalArg, int offset, int parent);
        inline void mpiPreComputation(const std::vector<int> *distFromRemotePtr, int numaLocalArg, int offset);
        inline void mpiPostComputation(const std::vector<int> *distFromRemotePtr, int numaLocalArg, int offset);
        inline void powerCallReleaseHopelessRegionLocks(int hopelessStartLevel, int parent);
        std::function<void(int)> recursivePowerFun;

        //double *a, *b, *c, *d;
        //int len;
    public:
        //FuncManager(void (*f_) (int,int,void *), void* args_, ZoneTree *zoneTree_, LevelPool* pool_, std::vector<int> serialPart);
        FuncManager(funcType f_, void* args_, ZoneTree *zoneTree_, LevelPool* pool_, std::vector<int> serialPart);
        //power is the actual power
        //subPower is the sub iterations within a power
        //Used for example when working with preconditioners
        //FuncManager(void (*f_) (int,int,int,int,int,void *), void* args_, int power, int subPower, mtxPowerRecursive *matPower, bool numaSplit_);
        FuncManager(powerFuncType f_, void* args_, int power, int subPower, mtxPowerRecursive *matPower, bool numaSplit_);
        FuncManager(const FuncManager &obj);
        ~FuncManager();

        void registerCommFunc(commFuncType commFunc_, void* commArgs_);
        void SerialPartRun();
        void initPowerRun(int parent);
        void NUMAInitPower();
        // void powerRun();
        void Run(bool rev_=false);
        void setPower(int power);
        void setSubPower(int subPower);
        int getPower();
        int getSubPower();
        void setSerial();
        void unsetSerial();
        bool isNumaSplit();
        bool isCommRegistered();
        //void RunOMP();
        //	double barrierTime;
};

#endif
