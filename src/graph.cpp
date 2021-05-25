//if GAP is there use optimised version
#include "config.h"

#ifdef RACE_USE_GAP
    #ifdef RACE_USE_SOA_GRAPH
        #include "graph_SoA.cpp"
        #pragma message ("SoA graph being used")
    #else
        #include "graph_AoS.cpp"
        #pragma message ("AoS graph being used")
    #endif
#else
    #ifdef RACE_USE_SOA_GRAPH
        #pragma message ("SoA graph not supported with serial BFS. Switch on RACE_USE_GAP. Now AoS graph being used")
    #else
        #pragma message ("AoS graph being used")
    #endif
    #include "graph_AoS.cpp"
#endif
