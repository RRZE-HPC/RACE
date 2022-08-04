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

/**
 * @file interface.h
 * @brief User interface to RACE library.
 * @author Christie Alappat <christie.alappat@fau.de>
 */

#ifndef RACE_INTERFACE_H
#define RACE_INTERFACE_H
#include "graph.h"
#include "type.h"
#include "levelData.h"
#include "zone_tree.h"
#include "level_recursion.h"
#include "functionManager.h"
#include "type.h"
#include "level_pool.h"
#include "config.h"
#include "matrixPowerRecursive.h"
#include "string.h"
/**
 * @brief RACE namespace.
 */
 namespace RACE
{
    class Interface;
}

/**
 * @brief Interface class for RACE.
 */
class RACE::Interface{
    private:
        RACE::Graph* graph;
        int nrow;
        RACE::dist distance;
        RACE::d2Method d2Type;
        RACE::LBTarget lbTarget;
        int requestedThreads;
        int availableThreads;
        int SMT;
        RACE::PinMethod pinMethod;
        LevelPool* pool;
        int *initPerm;
        int *initInvPerm;
        int *rowPtr;
        int *col;
/*        int *zonePtr;
        int zonePtrLen;*/
        int *perm;
//        int permLen;
        int *invPerm;
//        int invPermLen;
        ZoneTree* zoneTree;
        mtxPowerRecursive* powerCalculator;
        std::vector<FuncManager*> funMan;
        int highestPower;
        int highestSubPower;
        //	FuncManager *funMan;
        bool detectConflict(std::vector<int> range1, std::vector<int> range2);
        bool recursiveChecker(int parent);
        bool D2Checker();
        std::map<int, int> tunedPowMap;
    public:
        /**
         * @brief Constructor to Interface class.
         *
         * @param[in] nrow_ number of rows (vertices) in the matrix (graph).
         * @param[in] nthreads_ required number of parallelism.
         * @param[in] dist_ Select type of coloring, possible values are:
         * <TABLE>
         * <TR><TD>Type</TD><TD>Value</TD></TR>
         * <TR><TD>Distance-1</TD><TD>RACE::ONE</TD></TR>
         * <TR><TD>Distance-2</TD><TD>RACE::TWO</TD></TR>
         * <TR><TD>Matrix Power</TD><TD>RACE::POWER</TD></TR>
         * </TABLE>
         * @param[in] rowPtr_ Provide rowPtr of the sparse matrix in CRS format.
         * @param[in] col_ Provide column index of the sparse matrix in CRS format.
         * @param[in] SMT_ Number of SMT threads required (optional, default 1).
         * @param[in] pinMethod_ Select the type of pinning required, options
         * are:
         * <TABLE>
         * <TR><TD>Type</TD><TD>Value</TD></TR>
         * <TR><TD>Close/Compact</TD><TD>RACE::FILL</TD></TR>
         * <TR><TD>Scattered</TD><TD>RACE::SCATTER (default)</TD></TR>
         * </TABLE>
         * @param[in] initPerm_ Provide intital permutation vector of the original graph (optional, default=NULL).
         * @param[in] initInvPerm_ Provide intital inverse permutation vector of the original graph (optional, default=NULL).
         * @param[in] d2Type_ if distance-2 coloring select between two or
         * three colors.
         * <TABLE>
         * <TR><TD>Type</TD><TD>Value</TD></TR>
         * <TR><TD>Two color</TD><TD>RACE::TWO_BLOCK (default)</TD></TR>
         * <TR><TD>Three color</TD><TD>RACE::THREE_BLOCK</TD></TR>
         * </TABLE>
         * @param[in] lbTarget_ Load balancing type, w.r.t ROWS or NNZ
         * <TABLE>
         * <TR><TD>Type</TD><TD>Value</TD></TR>
         * <TR><TD>Rows</TD><TD>RACE::ROW </TD></TR>
         * <TR><TD>Non-zeros</TD><TD>RACE::NNZ (default)</TD></TR>
         */
        Interface(int nrow_, int nthreads_, RACE::dist dist_, int *rowPtr_, int *col_, int SMT_=1, RACE::PinMethod method_=RACE::SCATTER, int *initPerm_=NULL, int *initInvPerm_=NULL, RACE::d2Method d2Type_=RACE::TWO_BLOCK, RACE::LBTarget lbTarget_=RACE::NNZ);
        ~Interface();
        //Pre-processing
        RACE_error RACEColor(int highestPower, int numSharedCache, double cacheSize, double safetyFactor=2, std::string mtxType="N", int highestSubPower=1);
        RACE_error RACEColor();
        double getEfficiency();
        int getMaxStageDepth();
        void printZoneTree();
        //void getZoneTree(int **zoneTree_, int *len);
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getNumThreads();

        //sleep all threads
        void sleep();

        //Execution
        //register function for coloring
        //TODO: we could probably disable support to raw function pointers. it
        //shouldn't change any interface from user-side since even if user
        //provides raw function it will be converted to std::functio
        int registerFunction(void (*f) (int,int,void *), void* args);
        int registerFunction(std::function<void (int,int,void *)> f, void* args);

        //register function for matrix power
        int registerFunction(std::function<void(int,int,int,int,int,void *)> f, void* args, int power=1, int subPower=1, int numaSplit=false);
        int registerFunction(void (*f) (int,int,int,int,int,void *), void* args, int power=1, int subPower=1, int numaSplit=false);
        //simple register function for matrix power, with subPower=1
        int registerFunction(std::function<void(int,int,int,int,void *)> f, void* args, int power=1, int numaSplit=false);
        int registerFunction(void (*f) (int,int,int,int,void *), void* args, int power=1, int numaSplit=false);

        void executeFunction(int funcId, bool rev=false);
        void setPower(int funcId, int pow);
        void setSubPower(int funcId, int subPow);
        int getPower(int funcId);
        int getSubPower(int funcId);
        void setSerial(int funcId);
        void unsetSerial(int funcId);
        int tuneFunction(int funcId, bool rev=false);
        //RACE_error powerRun(int power, int *rowPtr, int *col, double *A, double *x);

        void resetTime();

        bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, double* val, bool diagFirst=false);

        bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, float* val, bool diagFirst=false);

        bool simdifyD1(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, double* val, bool dist2Compatible=false);

        bool simdifyD1(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, float* val, bool dist2Compatible=false);


        //function to pin threads similar to RACE
        //from external applications like OpenMP
        void pinThread(int threadId);

        //normal NUMA init
        void numaInitRowPtr(int *rowPtr);
        void numaInitMtxVec(int *rowPtr, int *col, double *val, double *x, int power=1);

        void numaInitRowPtr(int **rowPtr);
        void numaInitMtxVec(int **rowPtr, int **col, double **val, int power=1);


        //NUMA with actual splitting
        void getNumaSplitting(int **split, int *splitLen);
        int getHighestPower();
};

namespace RACE
{
    void wrappedPowerFunc(std::function<void(int,int,int,int,void *)>, int start, int end, int pow, int subPow, int numaLocal, void* args);
    void powerInitRowPtrFunc(int start, int end, int pow, int subPow, int numa_domain, void* arg);
    void powerInitMtxVecFunc(int start, int end, int pow, int subPow, int numa_domain, void* arg);
    void powerInitRowPtrNumaLocalFunc(int start, int end, int pow, int subPow, int numa_domain, void* arg);
    void powerInitMtxVecNumaLocalFunc(int start, int end, int pow, int subPow, int numa_domain, void* arg);
}

struct matValArg
{
    int nrow;
    int *rowPtr;
    int *col;
    double *val;
    double *x;
};

#define ENCODE_ARG(_nrow_, _rowPtr_, _col_, _val_, _x_)\
    matValArg *newArg_ = new matValArg;\
    newArg_->nrow = _nrow_;\
    newArg_->rowPtr = _rowPtr_;\
    newArg_->col = _col_;\
    newArg_->val = _val_;\
    newArg_->x = _x_;\
    void *voidArg = (void*) newArg_;

#define DECODE_ARG(_voidArg_)\
    matValArg *nonVoidArg_ = (matValArg*) _voidArg_;\
    int nrow = nonVoidArg_->nrow;\
    int* rowPtr = nonVoidArg_->rowPtr;\
    int* col = nonVoidArg_->col;\
    double* val = nonVoidArg_->val;\
    double* x = nonVoidArg_->x;\

struct matValNumaLocalArg
{
    int nrow;
    int **rowPtr;
    int **col;
    double **val;
    double *x;
};


#define ENCODE_NUMA_LOCAL_ARG(_nrow_, _rowPtr_, _col_, _val_, _x_)\
    matValNumaLocalArg *newArg_ = new matValNumaLocalArg;\
    newArg_->nrow = _nrow_;\
    newArg_->rowPtr = _rowPtr_;\
    newArg_->col = _col_;\
    newArg_->val = _val_;\
    newArg_->x = _x_;\
    void *voidArg = (void*) newArg_;

#define DECODE_NUMA_LOCAL_ARG(_voidArg_)\
    matValNumaLocalArg *nonVoidArg_ = (matValNumaLocalArg*) _voidArg_;\
    int nrow = nonVoidArg_->nrow;\
    int** rowPtr = nonVoidArg_->rowPtr;\
    int** col = nonVoidArg_->col;\
    double** val = nonVoidArg_->val;\
    double* x = nonVoidArg_->x;\



#endif
