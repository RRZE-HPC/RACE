#ifndef NAME_FUNCTION_MANAGER_H
#define NAME_FUNCTION_MANAGER_H

#include <functional>
#include "zone_tree.h"
#include "level_recursion.h"
#include "omp.h"
#include "level_pool.h"
#include "timing.h"

class FuncManager
{
    private:
        typedef std::function< void(int,int,void *) > funcType;
        funcType func;
        void* args;
        ZoneTree* zoneTree;
        LevelPool* pool;
        void recursiveCall(int parent);
        double *a, *b, *c, *d;
        int len;
    public:
        FuncManager(void (*f_) (int,int,void *), void* args_, ZoneTree *zoneTree_, LevelPool* pool_);
        ~FuncManager();
        void Run();
};

#endif
