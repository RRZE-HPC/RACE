#ifndef RACE_GRAPH_H
#define RACE_GRAPH_H
#include "config.h"

#ifdef RACE_USE_GAP
    #ifdef RACE_USE_SOA_GRAPH
        #include "graph_SoA.h"
    #else
        #include "graph_AoS.h"
    #endif
#else
    #include "graph_AoS.h"
#endif

#endif
