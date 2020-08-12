#include "matrixPowerRecursive.h"
#include "utility.h"

#define UPDATE_INVPERM(_invPerm_, _perm_)\
{\
    for(int i=0; i<totalRows; ++i)\
    {\
        _invPerm_[_perm_[i]] = i;\
    }\
}\

mtxPowerRecursive::mtxPowerRecursive(Graph* graph_, int highestPower_, int numSharedCache_, double cacheSize_, double safetyFactor_):graph(graph_), levelPtr(NULL),nodePtr(NULL), unlockRow(NULL), dangerRow(NULL),unlockCtr(NULL), highestPower(highestPower_), numSharedCache(numSharedCache_), cacheSize(cacheSize_), safetyFactor(safetyFactor_), perm(NULL), invPerm(NULL)
{
    int default_cutoff = 50; //make 50 factor so not much cut-off happens by default
    cache_violation_cutoff.push_back(default_cutoff); //default
    getEnv("RACE_CACHE_VIOLATION_CUTOFF", cache_violation_cutoff);

    totalRows = graph->NROW;
    perm =  new int[totalRows];
    invPerm =  new int[totalRows];
    for(int i=0; i<totalRows; ++i)
    {
        perm[i] = i;
        invPerm[i] = i;
    }
}

mtxPowerRecursive::~mtxPowerRecursive()
{
    if(levelPtr)
    {
        delete[] levelPtr;
    }
    if(nodePtr)
    {
        delete[] nodePtr;
    }
    if(unlockRow)
    {
        delete[] unlockRow;
    }
    if(dangerRow)
    {
        delete[] dangerRow;
    }
    if(unlockCtr)
    {
        delete[] unlockCtr;
    }
    if(perm)
    {
        delete[] perm;
    }
    if(invPerm)
    {
        delete[] invPerm;
    }
}

int mtxPowerRecursive::get_cache_violation_cutoff(int stage)
{
    if(stage < (int)cache_violation_cutoff.size())
    {
        return cache_violation_cutoff[stage];
    }
    else
    {
        return cache_violation_cutoff[0];
    }
}

#define INIT_STAGE_DATA\
    std::vector<int> hopelessRegions_curStage;\
 std::vector<int> levelPtr_curStage;\
 std::vector<int> nodePtr_curStage;\
 std::vector<int> unlockRow_curStage;\
 std::vector<int> unlockCtr_curStage;\
 std::vector<int> dangerRow_curStage;\
 int *perm_curStage;\
 int nrows_curStage;\

#define READ_STAGE_DATA(_stage_mtxPower_)\
    hopelessRegions_curStage = _stage_mtxPower_.getHopelessRegions();\
    levelPtr_curStage = _stage_mtxPower_.getLevelPtr();\
    nodePtr_curStage = _stage_mtxPower_.getNodePtr();\
    unlockRow_curStage = _stage_mtxPower_.getUnlockRow();\
    unlockCtr_curStage = _stage_mtxPower_.getUnlockCtr();\
    dangerRow_curStage = _stage_mtxPower_.getDangerRow();\
    _stage_mtxPower_.getPerm(&perm_curStage, &nrows_curStage);\

#define COPY_TO_PLAIN_INT_PTR(_vec_, _plain_)\
    _plain_ = new int[_vec_.size()];\
    for(int i=0; i<(int)_vec_.size(); ++i)\
    {\
        _plain_[i] = _vec_[i];\
    }\


void mtxPowerRecursive::findPartition()
{
    INIT_STAGE_DATA;
    int  stage=1;
    //partition for first stage
    mtxPower curStage(graph, highestPower, numSharedCache, cacheSize, safetyFactor, get_cache_violation_cutoff(stage), 0, graph->NROW);
    curStage.findPartition();
    READ_STAGE_DATA(curStage);
    updatePerm(&perm, perm_curStage, totalRows, totalRows);
    delete[] perm_curStage;
    ++stage;

    while(!hopelessRegions_curStage.empty())
    {
        int num_hopeless = (int)hopelessRegions_curStage.size()/2;
        //do recursive treatment at hopeless regions
        for(int h=0; h<num_hopeless; ++h)
        {
            int start_idx = hopelessRegions_curStage[2*h];
            int end_idx = hopelessRegions_curStage[2*h+1];
            int start_row = levelPtr_curStage[start_idx];
            int end_row = levelPtr_curStage[end_idx];
            mtxPower curNewStage(graph, highestPower, numSharedCache, cacheSize, safetyFactor, get_cache_violation_cutoff(stage), start_row, end_row);
            curNewStage.findPartition();
            READ_STAGE_DATA(curNewStage);
            updatePerm(&perm, perm_curStage, totalRows, totalRows);
            delete[] perm_curStage;
        }
        ++stage;
    }

    UPDATE_INVPERM(invPerm, perm);
    totalLevel = (int)(levelPtr_curStage.size()-1);
    //if empty this is final
    COPY_TO_PLAIN_INT_PTR(levelPtr_curStage, levelPtr);
    COPY_TO_PLAIN_INT_PTR(nodePtr_curStage, nodePtr);
    COPY_TO_PLAIN_INT_PTR(unlockRow_curStage, unlockRow);
    COPY_TO_PLAIN_INT_PTR(unlockCtr_curStage, unlockCtr);
    COPY_TO_PLAIN_INT_PTR(dangerRow_curStage, dangerRow);
}


int mtxPowerRecursive::getTotalLevel()
{
    return totalLevel;
}

int mtxPowerRecursive::getTotalNodes()
{
    return numSharedCache;
}

int*  mtxPowerRecursive::getLevelPtrRef()
{
    return levelPtr;
}

int*  mtxPowerRecursive::getNodePtrRef()
{
    return nodePtr;
}

int*  mtxPowerRecursive::getUnlockRowRef()
{
    return unlockRow;
}


int*  mtxPowerRecursive::getDangerRowRef()
{
    return dangerRow;
}

int*  mtxPowerRecursive::getUnlockCtrRef()
{
    return unlockCtr;
}

void mtxPowerRecursive::getPerm(int **perm_, int *len)
{
    (*perm_) = perm;
    perm = NULL; //hand-over responsibility
    (*len) = totalRows;
}

void mtxPowerRecursive::getInvPerm(int **invPerm_, int *len)
{
    (*invPerm_) = invPerm;
    invPerm = NULL; //hand-over responsibility
    (*len) = totalRows;
}
