//if GAP is there use optimised version
#include "config.h"

#ifdef RACE_USE_GAP
    #include "traverse_GAP.cpp"
#else
    #include "traverse_serial.cpp"
#endif
