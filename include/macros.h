/*
 * =======================================================================================
 *
 *   RACE: Recursicve Algebraic Coloring Engine
 *   Copyright (C) 2019, RRZE, Friedrich-Alexander-Universität Erlangen-Nürnberg
 *   Author: Christie Alappat
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * =======================================================================================
 */

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


#define SPLIT_LEVEL_PER_THREAD_BOUNDARY(_level_)\
    int _StartRowBoundary_ = _val_->at(_level_);\
    int _EndRowBoundary_ = _val_->at(_level_+1);\
    int _RowPerThreadBoundary_ = (_EndRowBoundary_- _StartRowBoundary_)/threadPerNode;\
    int startRow_tid_b = _StartRowBoundary_ + localTid*_RowPerThreadBoundary_;\
    startRow_tid_b -= offset;\
    int endRow_tid_b = (localTid == (threadPerNode-1)) ? _EndRowBoundary_ : _StartRowBoundary_+(localTid+1)*_RowPerThreadBoundary_;\
    endRow_tid_b -= offset;\

#define SPLIT_LEVEL_PER_THREAD_BOUNDARY_NOREF(_level_)\
    int _StartRowBoundary_ = _val_.at(_level_);\
    int _EndRowBoundary_ = _val_.at(_level_+1);\
    int _RowPerThreadBoundary_ = (_EndRowBoundary_- _StartRowBoundary_)/threadPerNode;\
    int startRow_tid_b = _StartRowBoundary_ + localTid*_RowPerThreadBoundary_;\
    startRow_tid_b -= offset;\
    int endRow_tid_b = (localTid == (threadPerNode-1)) ? _EndRowBoundary_ : _StartRowBoundary_+(localTid+1)*_RowPerThreadBoundary_;\
    endRow_tid_b -= offset;\

#define SPLIT_LEVEL_PER_THREAD_BOUNDARY_w_UNLOCK_DANGER(_level_)\
    int _StartRowBoundary_ = _val_->at(_level_);\
    int _EndRowBoundary_ = _val_->at(_level_+1);\
    int _RowPerThreadBoundary_ = (_EndRowBoundary_- _StartRowBoundary_)/threadPerNode;\
    int startRow_tid_b = _StartRowBoundary_ + localTid*_RowPerThreadBoundary_;\
    startRow_tid_b -= offset;\
    int endRow_tid_b = (localTid == (threadPerNode-1)) ? _EndRowBoundary_ : _StartRowBoundary_+(localTid+1)*_RowPerThreadBoundary_;\
    endRow_tid_b -= offset;\
    int currUnlockRow_b = ((boundaryUnlockRow->at(_workingRadius_)).at(_radius_))[_region_][_level_];\
    currUnlockRow_b -= offset;\
    int dangerRowStart_b = ((boundaryDangerRow->at(_workingRadius_)).at(_radius_))[_region_][_level_];\
    dangerRowStart_b -= offset;

#define SPLIT_LEVEL_PER_THREAD_BOUNDARY_w_UNLOCK_DANGER_NOREF(_level_)\
    int _StartRowBoundary_ = _val_.at(_level_);\
    int _EndRowBoundary_ = _val_.at(_level_+1);\
    int _RowPerThreadBoundary_ = (_EndRowBoundary_- _StartRowBoundary_)/threadPerNode;\
    int startRow_tid_b = _StartRowBoundary_ + localTid*_RowPerThreadBoundary_;\
    startRow_tid_b -= offset;\
    int endRow_tid_b = (localTid == (threadPerNode-1)) ? _EndRowBoundary_ : _StartRowBoundary_+(localTid+1)*_RowPerThreadBoundary_;\
    endRow_tid_b -= offset;\
    int currUnlockRow_b = boundaryUnlockRow[_workingRadius_][_radius_][_region_][_level_];\
    currUnlockRow_b -= offset;\
    int dangerRowStart_b = boundaryDangerRow[_workingRadius_][_radius_][_region_][_level_];\
    dangerRowStart_b -= offset;


//last argument is the statement that you want to init for each boundary ranges
#define INIT_BOUNDARY_STRUCTURE(_boundaryRange_, _newVar_, ...)\
{\
    int _wbl_ = (int)_boundaryRange_.size();\
    _newVar_.resize(_wbl_);\
    for(int _workingRadius_=0; _workingRadius_<_wbl_; ++_workingRadius_)\
    {\
        std::map<int, std::vector<Range>> *_curBoundaryRange_ = &(_boundaryRange_[_workingRadius_]);\
        for(auto _mapIter_ = _curBoundaryRange_->begin(); _mapIter_ != _curBoundaryRange_->end(); ++_mapIter_)\
        {\
            int _radius_ = _mapIter_->first;\
            std::vector<Range>* _ranges_ = &(_mapIter_->second);\
            int _num_regions_ = (int)_ranges_->size();\
            _newVar_[_workingRadius_][_radius_].resize(_num_regions_);\
            for(int _region_=0; _region_<_num_regions_; ++_region_)\
            {\
                _newVar_[_workingRadius_][_radius_][_region_] = __VA_ARGS__;\
            }\
        }\
    }\
}

//last argument is the statement that you want to execute for each boundary ranges
#define EXEC_BOUNDARY_STRUCTURE(_var_, ...)\
{\
    int _wbl_ = (int)_var_.size();\
    for(int _workingRadius_=0; _workingRadius_<_wbl_; ++_workingRadius_)\
    {\
        for(auto _mapIter_ = _var_[_workingRadius_].begin(); _mapIter_ != _var_[_workingRadius_].end(); ++_mapIter_)\
        {\
            int _radius_ = _mapIter_->first;\
            auto* _entity_ = &(_mapIter_->second);\
            int _numBoundaries_ = (int)_entity_->size();\
            for(int _region_=0; _region_<_numBoundaries_; ++_region_)\
            {\
                auto _val_ = _entity_->at(_region_);\
                __VA_ARGS__;\
            }\
        }\
    }\
}

//last argument is the statement that you want to execute for each boundary ranges
//only difference with previous macro is: _radius_ is not defined in inner-loop,
//so avoid annoying compiler remarks
#define EXEC_BOUNDARY_STRUCTURE_wo_radius(_var_, ...)\
{\
    int _wbl_ = (int)_var_.size();\
    for(int _workingRadius_=0; _workingRadius_<_wbl_; ++_workingRadius_)\
    {\
        for(auto _mapIter_ = _var_[_workingRadius_].begin(); _mapIter_ != _var_[_workingRadius_].end(); ++_mapIter_)\
        {\
            auto _entity_ = &(_mapIter_->second);\
            int _numBoundaries_ = (int)_entity_->size();\
            for(int _region_=0; _region_<_numBoundaries_; ++_region_)\
            {\
                auto _val_ = _entity_->at(_region_);\
                __VA_ARGS__;\
            }\
        }\
    }\
}


#define EXEC_BOUNDARY_STRUCTURE_w_wave_shape_wo_radius(_var_, _pow_,...)\
{\
    int _wbl_ = (int)_var_.size();\
    for(int _workingRadius_=0; _workingRadius_<_wbl_; ++_workingRadius_)\
    {\
        if((_pow_ > _workingRadius_) && ( _pow_ < (power-(_workingRadius_+1)) ))\
        {\
            for(auto _mapIter_ = _var_[_workingRadius_].begin(); _mapIter_ != _var_[_workingRadius_].end(); ++_mapIter_)\
            {\
                auto* _entity_ = &(_mapIter_->second);\
                int _numBoundaries_ = (int)_entity_->size();\
                for(int _region_=0; _region_<_numBoundaries_; ++_region_)\
                {\
                    auto* _val_ = &(_entity_->at(_region_));\
                    __VA_ARGS__;\
                }\
            }\
        }\
    }\
}


#define EXEC_BOUNDARY_STRUCTURE_w_wave_shape(_var_, _pow_,...)\
{\
    int _wbl_ = (int)_var_.size();\
    for(int _workingRadius_=0; _workingRadius_<_wbl_; ++_workingRadius_)\
    {\
        if((_pow_ > _workingRadius_) && ( _pow_ < (power-(_workingRadius_+1)) ))\
        {\
            for(auto _mapIter_ = _var_[_workingRadius_].begin(); _mapIter_ != _var_[_workingRadius_].end(); ++_mapIter_)\
            {\
                int _radius_ = _mapIter_->first;\
                int _absRadius_ = std::abs(_radius_);\
                if( (_pow_ > (_absRadius_-1)) && (_pow_ < (power-(_absRadius_))) )\
                {\
                    auto* _entity_ = &(_mapIter_->second);\
                    int _numBoundaries_ = (int)_entity_->size();\
                    for(int _region_=0; _region_<_numBoundaries_; ++_region_)\
                    {\
                        auto* _val_ = &(_entity_->at(_region_));\
                        __VA_ARGS__;\
                    }\
                }\
            }\
        }\
    }\
}


#endif
