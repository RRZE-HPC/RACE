#ifndef NAME_FUNCTION_MANAGER_H
#define NAME_FUNCTION_MANAGER_H

#include <functional>
#include "zone_tree.h"
#include "level_recursion.h"
#include "omp.h"
#include "pin.h"

class FuncManager
{
    private:
        typedef std::function< void(int,int,void *) > funcType;
        funcType func;
        void* args;
        ZoneTree* zoneTree;
        Pin* pin;
        void recursiveCall(int parent);
    public:
        FuncManager(void (*f_) (int,int,void *), void* args_, ZoneTree *zoneTree_, Pin* pin_);
        ~FuncManager();
        void Run();
};

#endif
