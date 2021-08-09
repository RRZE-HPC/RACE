//if GAP is there use optimised version
#include "config.h"

#ifdef RACE_USE_GAP
    #ifdef RACE_USE_SOA_GRAPH
        #include "traverse_GAP_graphSoA.cpp"
        #pragma message ("GAP BFS with SoA graph being using")
    #else
        #include "traverse_GAP_graphAoS.cpp"
        #pragma message ("GAP BFS with AoS graph being using")
    #endif
#else
    #include "traverse_serial.cpp"
    #pragma message ("Serial BFS with AoS graph being used")
#endif
