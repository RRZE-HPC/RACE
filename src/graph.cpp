//if GAP is there use optimised version
#include "config.h"

#ifdef RACE_USE_GAP
    #include "graph_serial.cpp"
#else
    #include "graph_serial.cpp"
#endif
