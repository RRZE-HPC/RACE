#ifndef RACE_TRAVERSE_H
#define RACE_TRAVERSE_H
#include "config.h"

#ifdef RACE_USE_GAP
    #include "traverse_GAP.h"
    #pragma message ( "GAP BFS being used" )
#else
    #include "traverse_serial.h"
    #pragma message ( "Serial BFS being used" )
#endif

#endif
