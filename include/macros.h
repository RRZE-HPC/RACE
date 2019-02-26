#ifndef RACE_MACRO_H
#define RACE_MACRO_H

#include "type.h"

#define getBlockPerThread(dist, d2Type)\
   ( (dist==RACE::POWER) ? 1 : (dist==RACE::ONE)?2:(d2Type==RACE::TWO_BLOCK)?2:3 )


#define getMinGap(dist, d2Type)\
   ( (dist==RACE::POWER) ? 1 : (dist==RACE::ONE)?1:(d2Type==RACE::TWO_BLOCK)?2:1 )

#define getPossibleThreads(totalLevel, dist, d2Type)\
    ( (dist==RACE::POWER) ? static_cast<int>(totalLevel) : (dist==RACE::ONE)?static_cast<int>( (totalLevel) / 2.0):(d2Type==RACE::TWO_BLOCK)?static_cast<int>( (totalLevel) / 4.0): static_cast<int>( (totalLevel) /3.0) )

#define UNUSED(x) (void)(x)


#endif
