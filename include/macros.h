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

#define SPLIT_LEVEL_PER_THREAD(_level_)\
    int _startRow_ = levelPtr->at(_level_);\
    int _endRow_ = levelPtr->at(_level_+1);\
    int _RowPerThread_ = (_endRow_ - _startRow_)/threadPerNode;\
    int startRow_tid = _startRow_ + localTid*_RowPerThread_;\
    startRow_tid -= offset;\
    int endRow_tid = (localTid == (threadPerNode-1)) ? _endRow_ : _startRow_ + (localTid+1)*_RowPerThread_;\
    endRow_tid -= offset;\

#define SPLIT_LEVEL_PER_THREAD_P2P(_level_)\
    int _startRow_ = levelPtr->at(_level_);\
    int _endRow_ = levelPtr->at(_level_+1);\
    int _RowPerThread_ = (_endRow_ - _startRow_)/threadPerNode;\
    int startRow_tid = _startRow_ + localTid*_RowPerThread_;\
    startRow_tid -= offset;\
    int endRow_tid = (localTid == (threadPerNode-1)) ? _endRow_ : _startRow_ + (localTid+1)*_RowPerThread_;\
    endRow_tid -= offset;\
    int currUnlockRow = unlockRow->at(_level_);\
    currUnlockRow -= offset;\
    int dangerRowStart = dangerRow->at(_level_);\
    dangerRowStart -= offset;\

#define SPLIT_LEVEL_PER_THREAD_P2P_NOREF(_level_)\
    int _startRow_ = levelPtr[_level_];\
    int _endRow_ = levelPtr[_level_+1];\
    int _RowPerThread_ = (_endRow_ - _startRow_)/threadPerNode;\
    int startRow_tid = _startRow_ + localTid*_RowPerThread_;\
    startRow_tid -= offset;\
    int endRow_tid = (localTid == (threadPerNode-1)) ? _endRow_ : _startRow_ + (localTid+1)*_RowPerThread_;\
    endRow_tid -= offset;\
    int currUnlockRow = unlockRow[_level_];\
    currUnlockRow -= offset;\
    int dangerRowStart = dangerRow[_level_];\
    dangerRowStart -= offset;\



#endif
