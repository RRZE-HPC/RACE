#ifndef NAME_INTERFACE_H
#define NAME_INTERFACE_H
#include "graph.h"
#include "traverse.h"
#include "lb.h"
#include "type.h"
#include "levelData.h"

class NAMEInterface{
    private:
        int nrow;
        dist_t dist;
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
        NAMEInterface(int nrow_, int nthreads_, dist_t dist_, int *rowPtr_, int *col_, int *initPerm_=NULL, int *initInvPerm_=NULL);
        void NAMEColor();
        void getZonePtr(int **zonePtr_, int *len_);
        void getPerm(int **perm_, int *len_);
        void getInvPerm(int **invPerm_, int *len_);
        int getNumThreads();
};

#endif
