#ifndef RACE_TRAVERSE_H
#define RACE_TRAVERSE_H
#include "config.h"

#ifdef RACE_USE_GAP
    #include "traverse_GAP.h"
#else
    #include "traverse_serial.h"
#endif

#endif
