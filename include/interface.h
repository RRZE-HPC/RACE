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
#include "matrixPower.h"

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
        Graph* graph;
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
        mtxPower* powerCalculator;
        std::vector<FuncManager*> funMan;
        //	FuncManager *funMan;
        bool detectConflict(std::vector<int> range1, std::vector<int> range2);
        bool recursiveChecker(int parent);
        bool D2Checker();
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
        RACE_error RACEColor(int highestPower, int numSharedCache, double cacheSize, double safetyFactor=2);
        //Pre-processing
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
        int registerFunction(void (*f) (int,int,void *), void* args);
        int registerFunction(void (*f) (int,int,int,void *), void* args, int power=1);
        void executeFunction(int funcId, bool rev=false);
        //RACE_error powerRun(int power, int *rowPtr, int *col, double *A, double *x);

        void resetTime();

        bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, double* val, bool diagFirst=false);

        bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, float* val, bool diagFirst=false);

        bool simdifyD1(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, double* val, bool dist2Compatible=false);

        bool simdifyD1(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, float* val, bool dist2Compatible=false);


        //function to pin threads similar to RACE
        //from external applications like OpenMP
        void pinThread(int threadId);
        void numaInitRowPtr(int *rowPtr);
        void numaInitMtxVec(int *rowPtr, int *col, double *val, double *x, int power=1);
};

void powerInitRowPtrFunc(int start, int end, int pow, void* arg);
void powerInitMtxVecFunc(int start, int end, int pow, void* arg);

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


#endif
