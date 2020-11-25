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

#define SPLIT_LEVEL_PER_THREAD_BOUNDARY(_level_, _boundary_)\
    int _negativeStartRow_ = levelPtrNegativeBoundary->at(_boundary_)[_level_];\
    int _negativeEndRow_ = levelPtrNegativeBoundary->at(_boundary_)[_level_+1];\
    int _negativeRowPerThread_ = (_negativeEndRow_- _negativeStartRow_)/threadPerNode;\
    int negative_startRow_tid = _negativeStartRow_ + localTid*_negativeRowPerThread_;\
    negative_startRow_tid -= offset;\
    int negative_endRow_tid = (localTid == (threadPerNode-1)) ? _negativeEndRow_ : _negativeStartRow_+(localTid+1)*_negativeRowPerThread_;\
    negative_endRow_tid -= offset;\
    int _positiveStartRow_ = levelPtrPositiveBoundary->at(_boundary_)[_level_];\
    int _positiveEndRow_ = levelPtrPositiveBoundary->at(_boundary_)[_level_+1];\
    int _positiveRowPerThread_ = (_positiveEndRow_- _positiveStartRow_)/threadPerNode;\
    int positive_startRow_tid = _positiveStartRow_ + localTid*_positiveRowPerThread_;\
    positive_startRow_tid -= offset;\
    int positive_endRow_tid = (localTid == (threadPerNode-1)) ? _positiveEndRow_ : _positiveStartRow_+(localTid+1)*_positiveRowPerThread_;\
    positive_endRow_tid -= offset;\


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
