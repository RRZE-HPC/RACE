#ifndef RACE_FUNCTION_MANAGER_H
#define RACE_FUNCTION_MANAGER_H

#include <functional>
#include "zone_tree.h"
#include "level_recursion.h"
#include "omp.h"
#include "level_pool.h"
#include "timing.h"
#include "matrixPower.h"

class FuncManager;

void recursiveCall(FuncManager* funMan, int parent);

class FuncManager
{
    private:
        bool rev;
        typedef std::function< void(int,int,void *) > funcType;
        typedef std::function< void(int,int,int,void *) > powerFuncType;
        bool power_fn;
        funcType func;
        powerFuncType powerFunc;
        void* args;
        int power;
        ZoneTree* zoneTree;
        LevelPool* pool;
        mtxPower* matPower;
        std::vector<int> serialPart;
        volatile int* lockCtr;
        volatile bool* lockTable;
        volatile int* lockTableCtr;
        int* unlockRow;
        int* dangerRow;
        int* unlockCtr;
        friend void recursiveCall(FuncManager* funMan, int parent);
        std::function<void(int)> recursiveFun;
        //double *a, *b, *c, *d;
        //int len;
    public:
        FuncManager(void (*f_) (int,int,void *), void* args_, ZoneTree *zoneTree_, LevelPool* pool_, std::vector<int> serialPart);
        FuncManager(void (*f_) (int,int,int,void *), void* args_, int power, mtxPower *matPower);
        FuncManager(const FuncManager &obj);
        ~FuncManager();
        void SerialPartRun();
        void initPowerRun();
        void NUMAInitPower();
        void powerRun();
        void Run(bool rev_=false);
        //void RunOMP();
        //	double barrierTime;
};

#endif
