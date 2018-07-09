#ifndef RACE_PIN_H
#define RACE_PIN_H

#include "zone_tree.h"
#include "error.h"
#include "type.h"
#include "thpool.h"
#include <omp.h>

class Pin{
    private:
        ZoneTree *zoneTree;
        std::vector<std::pair<int,int>> pinMap;
        void* machine; //this is of type Machine, but casting to void to avoid external library dependency
        int SMT;
        RACE::PinMethod method;
        void pinOrderRecursive(int parentIdx);
        void calcPinOrder();
        RACE_error createPuNodeMapping();
        void getNodeId(int cpuId);
        //For pthread
        //void createPinnedThreadPool(int parentIdx);
        //For OMP: currently broken
        RACE_error pinApplicationRecursive(int parentIdx);
    public:
        Pin(ZoneTree* zoneTree_, int SMT_, RACE::PinMethod method_);
        RACE_error pinInit();
        RACE_error pinThread(int pinOrder);
        RACE_error pinApplication();
        void resetMaster();
};

#endif
