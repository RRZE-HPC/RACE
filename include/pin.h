#ifndef NAME_PIN_H
#define NAME_PIN_H

#include "machine.h"
#include "zone_tree.h"
#include "error.h"
#include "type.h"
#include "thpool.h"
#include <omp.h>

class Pin{
    private:
        ZoneTree *zoneTree;
        std::vector<std::pair<int,int>> pinMap;
        Machine machine;
        int SMT;
        PinMethod method;
        void pinOrderRecursive(int parentIdx);
        void calcPinOrder();
        void createPuNodeMapping();
        void getNodeId(int cpuId);
        //For pthread
        //void createPinnedThreadPool(int parentIdx);
        //For OMP: currently broken
        void pinApplicationRecursive(int parentIdx);
    public:
        Pin(ZoneTree* zoneTree_, int SMT_, PinMethod method_);
        void pinInit();
        void pinThread(int pinOrder);
        void pinApplication();
        void resetMaster();
};

#endif
