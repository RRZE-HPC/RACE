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

#include "level_pool.h"
#include "functional"

LevelPool::LevelPool(ZoneTree *zoneTree_, int SMT_, RACE::PinMethod pinMethod_):zoneTree(zoneTree_),pin(zoneTree_,SMT_, pinMethod_),tree(NULL)
{
    int totThreads = zoneTree->at(0).nthreadsZ;
    pool.init(totThreads);

    int ctr = 0;
    //create mappedIdx
    for(int i=0; i<zoneTree->size(); ++i)
    {
        mappedIdx.push_back(ctr);
        ctr += zoneTree->at(i).totalSubBlocks;
    }

    tree = new team<int> [ctr];
    pinInitSuccess = 1;
    RACE_error ret = pin.pinInit();

    if(ret!=RACE_SUCCESS)
    {
        pinInitSuccess = 0;
        ERROR_PRINT("Error in initing pin");
    }
}

LevelPool::~LevelPool()
{
    delete[] tree;
}

//Recursively spawn thread pools
void LevelPool::createPoolRecursive(int parentIdx)
{

    RECURSIVE_HELPER(LevelPool::createPoolRecursive,
            std::vector<int> gid(nThreads);
            for(unsigned i=0; i<gid.size(); ++i)
            {
                gid[i] = zoneTree->at(children->at(2*parentSubIdx)+blockPerThread*i).pinOrder;
            }
            tree[poolTreeIdx(parentIdx, parentSubIdx)].init(gid, &pool);
    );
}

RACE_error LevelPool::createPool()
{
    int root = 0;
    createPoolRecursive(root);
    RACE_error ret = pinPool();
    return ret;
}


//Recursively pin pools
RACE_error LevelPool::pinPoolRecursive(int parentIdx)
{
    RACE_error save_ret = RACE_SUCCESS;
    RECURSIVE_HELPER(LevelPool::pinPoolRecursive, 
            RACE_error ret = pin.pinThread(zoneTree->at(parentIdx).pinOrder);
            if(ret != RACE_SUCCESS)
            {
               save_ret = ret;
            }
        );

    return save_ret;
}


RACE_error LevelPool::pinPool()
{
    if(pinInitSuccess)
    {
        int root = 0;
        RACE_error ret = pinPoolRecursive(root);
        return ret;
    }
    else
    {
        return RACE_ERR_PIN;
    }
}

void LevelPool::resetMaster()
{
    pin.resetMaster();
}

void LevelPool::sleepPoolRecursive(int parentIdx)
{
    RECURSIVE_HELPER(LevelPool::sleepPoolRecursive,
           tree[poolTreeIdx(parentIdx,parentSubIdx)].sleep();
    );
}


void LevelPool::sleepPool()
{
    int root = 0;
    sleepPoolRecursive(root);
}


