#ifndef RACE_FUNCTION_MANAGER_H
#define RACE_FUNCTION_MANAGER_H

#include <functional>
#include "zone_tree.h"
#include "level_recursion.h"
#include "omp.h"
#include "level_pool.h"
#include "timing.h"
#include "matrixPowerRecursive.h"

class FuncManager;

void recursiveCall(FuncManager* funMan, int parent);

void recursivePowerCall(FuncManager* funMan, int parent);

void recursivePowerCallSerial(FuncManager* funMan, int parent);

class FuncManager
{
    private:
        bool rev;
        typedef std::function< void(int,int,void *) > funcType;
        typedef std::function< void(int,int,int,int,void *) > powerFuncType;
        bool power_fn;
        bool numaSplit;
        funcType func;
        powerFuncType powerFunc;
        void* args;
        int power;
        int activethreads;
        ZoneTree* zoneTree;
        LevelPool* pool;
        mtxPowerRecursive* matPower;
        std::vector<int> serialPart;
        std::vector<volatile int*> lockCtr;
        std::vector<volatile bool*> lockTable;
        std::vector<volatile int*> lockTableCtr;
        //int* unlockRow;
        //int* dangerRow;
        //int* unlockCtr;
        friend void recursiveCall(FuncManager* funMan, int parent);
        std::function<void(int)> recursiveFun;

        friend void recursivePowerCallSerial(FuncManager* funMan, int parent);

        void powerCallGeneral(int startLevel, int endLevel, int boundaryStart, int boundaryEnd, int startSlope, int endSlope, int threadPerNode, const std::vector<int> *levelPtr, const std::vector<int> *unlockRow, const std::vector<int> *unlockCtr, const std::vector<int> *dangerRow, int numaLocalArg, int offset, int parent);
        void powerCallNodeReminder(int startSlope, int endSlope, int threadPerNode, const std::vector<int> *levelPtr, const std::vector<int> *nodePtr, int numaLocalArg, int offset);
        void powerCallHopelessReminder(int level, int startPower, int endPower, int direction, int threadPerNode, const std::vector<int> *levelPtr, int numaLocalArg, int offset, int parent);
            std::function<void(int)> recursivePowerFun;

        //double *a, *b, *c, *d;
        //int len;
    public:
        FuncManager(void (*f_) (int,int,void *), void* args_, ZoneTree *zoneTree_, LevelPool* pool_, std::vector<int> serialPart);
        FuncManager(void (*f_) (int,int,int,int,void *), void* args_, int power, mtxPowerRecursive *matPower, bool numaSplit_);
        FuncManager(const FuncManager &obj);
        ~FuncManager();
        void SerialPartRun();
        void initPowerRun(int parent);
        void NUMAInitPower();
       // void powerRun();
        void Run(bool rev_=false);
        bool isNumaSplit();
        //void RunOMP();
        //	double barrierTime;
};

#endif
