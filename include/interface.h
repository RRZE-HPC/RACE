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


class RACEInterface{
    private:
        Graph* graph;
        int nrow;
        dist_t dist;
        d2Method d2Type;
        int requestedThreads;
        int availableThreads;
        int SMT;
        PinMethod pinMethod;
        LevelPool* pool;
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
        ZoneTree* zoneTree;
        std::vector<FuncManager*> funMan;
        //	FuncManager *funMan;
        bool detectConflict(std::vector<int> range1, std::vector<int> range2);
        bool recursiveChecker(int parent);
        bool D2Checker();
    public:
        RACEInterface(int nrow_, int nthreads_, dist_t dist_, int *rowPtr_, int *col_, int SMT_=1, PinMethod method_=SCATTER, int *initPerm_=NULL, int *initInvPerm_=NULL, d2Method d2Type_=TWO_BLOCK);
        ~RACEInterface();
        //Pre-processing
        void RACEColor();
        void printZoneTree();
        //void getZoneTree(int **zoneTree_, int *len);
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getNumThreads();

        //sleep all threads
        void sleep();

        //Execution
        int registerFunction(void (*f) (int,int,void *), void* args);
        void executeFunction(int funcId);
        void resetTime();

        bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, double* val);

        bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, float* val);

        //function to pin threads similar to RACE
        //from external applications like OpenMP
        void pinThread(int threadId);
};

#endif
