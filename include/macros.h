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

#define SPLIT_LEVEL_PER_THREAD_BOUNDARY(_level_)\
    int _StartRow_ = _val_[_level_];\
    int _EndRow_ = _val_[_level_+1];\
    int _RowPerThread_ = (_EndRow_- _StartRow_)/threadPerNode;\
    int startRow_tid = _StartRow_ + localTid*_RowPerThread_;\
    startRow_tid -= offset;\
    int endRow_tid = (localTid == (threadPerNode-1)) ? _EndRow_ : _StartRow_+(localTid+1)*_RowPerThread_;\
    endRow_tid -= offset;\


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


//last argument is the statement that you want to init for each boundary ranges
#define INIT_BOUNDARY_STRUCTURE(_boundaryRange_, _newVar_, ...)\
{\
    int _numBoundaries_ = (int)_boundaryRange_.size();\
    _newVar_.resize(_numBoundaries_);\
    for(int _region_=0; _region_<_numBoundaries_; ++_region_)\
    {\
        int _wbl_=_boundaryRange_[_region_].size();\
        _newVar_[_region_].resize(_wbl_);\
        for(int _workingRadius_=0; _workingRadius_<_wbl_; ++_workingRadius_)\
        {\
            std::map<int, Range> *_curBoundaryRange_ = &(_boundaryRange_[_region_][_workingRadius_]);\
            for(auto _mapIter_ = _curBoundaryRange_->begin(); _mapIter_ != _curBoundaryRange_->end(); ++_mapIter_)\
            {\
                int _radius_ = _mapIter_->first;\
                _newVar_[_region_][_workingRadius_][_radius_] = __VA_ARGS__;\
            }\
        }\
    }\
}\

//last argument is the statement that you want to execute for each boundary ranges
#define EXEC_BOUNDARY_STRUCTURE(_var_, ...)\
{\
    int _numBoundaries_ = (int)_var_.size();\
    for(int _region_=0; _region_<_numBoundaries_; ++_region_)\
    {\
        int _wbl_=_var_[_region_].size();\
        for(int _workingRadius_=0; _workingRadius_<_wbl_; ++_workingRadius_)\
        {\
            for(auto _mapIter_ = _var_[_region_][_workingRadius_].begin(); _mapIter_ != _var_[_region_][_workingRadius_].end(); ++_mapIter_)\
            {\
                int _radius_ = _mapIter_->first;\
                auto _val_ = _mapIter_->second;\
                __VA_ARGS__;\
            }\
        }\
    }\
}\


//last argument is the statement that you want to execute for each boundary ranges
//only difference with previous macro is: _radius_ is not defined in inner-loop,
//so avoid annoying compiler remarks
#define EXEC_BOUNDARY_STRUCTURE_wo_radius(_var_, ...)\
{\
    int _numBoundaries_ = (int)_var_.size();\
    for(int _region_=0; _region_<_numBoundaries_; ++_region_)\
    {\
        int _wbl_=_var_[_region_].size();\
        for(int _workingRadius_=0; _workingRadius_<_wbl_; ++_workingRadius_)\
        {\
            for(auto _mapIter_ = _var_[_region_][_workingRadius_].begin(); _mapIter_ != _var_[_region_][_workingRadius_].end(); ++_mapIter_)\
            {\
                auto _val_ = _mapIter_->second;\
                __VA_ARGS__;\
            }\
        }\
    }\
}\

#endif
