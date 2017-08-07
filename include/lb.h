#ifndef NAME_LB_H
#define NAME_LB_H

#include "type.h"
#include "levelData.h"
#include "error.h"

enum neighbour_t{
    acquire,
    give
};

class Stat{
    private:
        int *arr;
        int len;

    public:
        Stat(int *arr_, int len_, int numPartitions_);
        ~Stat();
        int numPartitions;
        void calculate();
        double *mean;
        double *var;
        double *acquireWeight;
        double *giveWeight;
        double *weight;
};

class LB{
    private:
        int* levelPtr;
        int* zonePtr;
        LevelData* levelData;
        dist_t dist;
        d2Method d2Type;
        //gap to be left between 2 levels
        int minGap;

        int maxThreads;
        int nThreads;
        int blockPerThread;
        int nBlocks;

        LB_t lbTarget;

        void calcChunkSum(int *arr, bool forceRow=false);
        void calcZonePtr(int base);
        int  findNeighbour(const Stat &stats, neighbour_t type);
        void moveOneStep(int toIdx, int fromIdx);
        void splitZones();

    public:
        LB(int nThreads_, LevelData* levelData_, dist_t dist_, d2Method d2Type, LB_t lbTarget = NNZ); //constructor
        ~LB(); //destructor

        int getMaxThreads();
        int getNumThreads();
        void getZonePtr(int **zonePtr_, int *len, int base=0);
        NAME_error balance();
};


#endif
