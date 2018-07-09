#ifndef RACE_FUNCTION_MANAGER_H
#define RACE_FUNCTION_MANAGER_H

#include <functional>
#include "zone_tree.h"
#include "level_recursion.h"
#include "omp.h"
#include "level_pool.h"
#include "timing.h"

class FuncManager;

void recursiveCall(FuncManager* funMan, int parent);

class FuncManager
{
    private:
        bool rev;
        typedef std::function< void(int,int,void *) > funcType;
        funcType func;
        void* args;
        ZoneTree* zoneTree;
        LevelPool* pool;
        std::vector<int> serialPart;
        friend void recursiveCall(FuncManager* funMan, int parent);
        std::function<void(int)> recursiveFun;
        //double *a, *b, *c, *d;
        //int len;
    public:
        FuncManager(void (*f_) (int,int,void *), void* args_, ZoneTree *zoneTree_, LevelPool* pool_, std::vector<int> serialPart);
        FuncManager(const FuncManager &obj);
        ~FuncManager();
        void SerialPartRun();
        void Run(bool rev_=false);
        //void RunOMP();
        //	double barrierTime;
};

#endif
