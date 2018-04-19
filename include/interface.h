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

namespace RACE
{
    class Interface;
}

class RACE::Interface{
    private:
        Graph* graph;
        int nrow;
        RACE::dist distance;
        RACE::d2Method d2Type;
        int requestedThreads;
        int availableThreads;
        int SMT;
        RACE::PinMethod pinMethod;
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
        Interface(int nrow_, int nthreads_, RACE::dist dist_, int *rowPtr_, int *col_, int SMT_=1, RACE::PinMethod method_=RACE::SCATTER, int *initPerm_=NULL, int *initInvPerm_=NULL, RACE::d2Method d2Type_=RACE::TWO_BLOCK);
        ~Interface();
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
        void executeFunction(int funcId, bool rev=false);
        void resetTime();

        bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, double* val, bool diagFirst=false);

        bool simdify(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, float* val, bool diagFirst=false);

        bool simdifyD1(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, double* val, bool dist2Compatible=false);

        bool simdifyD1(int simdWidth, int C, int nrows, int*col, int* chunkStart, int* rl, int* clp, float* val, bool dist2Compatible=false);


        //function to pin threads similar to RACE
        //from external applications like OpenMP
        void pinThread(int threadId);
};

#endif
