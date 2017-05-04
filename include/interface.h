#ifndef NAME_INTERFACE_H
#define NAME_INTERFACE_H
#include "graph.h"
#include "type.h"
#include "levelData.h"
#include "zone_tree.h"
#include "level_recursion.h"
#include "functionManager.h"
#include "machine.h"
#include "pin.h"
#include "type.h"

class NAMEInterface{
    private:
        Graph* graph;
        int nrow;
        dist_t dist;
        int requestedThreads;
        int availableThreads;
        bool SMT;
        PinMethod pinMethod;
        Pin *pin;
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
        std::vector<FuncManager> funMan;
        bool detectConflict(std::vector<int> range1, std::vector<int> range2);
        bool recursiveChecker(int parent);
        bool D2Checker();

    public:
        NAMEInterface(int nrow_, int nthreads_, dist_t dist_, int *rowPtr_, int *col_, bool SMT_=false, PinMethod method_=SCATTER, int *initPerm_=NULL, int *initInvPerm_=NULL);
        ~NAMEInterface();

        //Pre-processing
        void NAMEColor();
        void printZoneTree();
        //void getZoneTree(int **zoneTree_, int *len);
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getNumThreads();

        //Execution
        int registerFunction(void (*f) (int,int,void *), void* args);
        void executeFunction(int funcId);
};

#endif
