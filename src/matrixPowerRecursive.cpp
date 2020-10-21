#include "matrixPowerRecursive.h"
#include "utility.h"

#define UPDATE_INVPERM(_invPerm_, _perm_)\
{\
    for(int i=0; i<totalRows; ++i)\
    {\
        _invPerm_[_perm_[i]] = i;\
    }\
}\

mtxPowerRecursive::mtxPowerRecursive(Graph* graph_, int highestPower_, int numSharedCache_, double cacheSize_, double safetyFactor_):graph(graph_), highestPower(highestPower_), numSharedCache(numSharedCache_), cacheSize(cacheSize_), safetyFactor(safetyFactor_), perm(NULL), invPerm(NULL)
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
    int cache_cutoff = 50;
    if(stage < (int)(cache_violation_cutoff.size()-1))
    {
        cache_cutoff = cache_violation_cutoff[stage+1]; //Stage numbering starts with 0, therefore stage+1
    }
    else
    {
        cache_cutoff = cache_violation_cutoff[0];
    }

    return cache_cutoff;
}

#define READ_STAGE_DATA(_leaf_, _stage_mtxPower_)\
    _leaf_.hopelessRegions = _stage_mtxPower_.getHopelessRegions();\
    _leaf_.levelPtr = _stage_mtxPower_.getLevelPtr();\
    _leaf_.unlockRow = _stage_mtxPower_.getUnlockRow();\
    _leaf_.unlockCtr = _stage_mtxPower_.getUnlockCtr();\
    _leaf_.dangerRow = _stage_mtxPower_.getDangerRow();\
    _stage_mtxPower_.getPerm(&perm_curStage, &_leaf_.nrows);\
    std::vector<int> _nodePtr_ = _stage_mtxPower_.getNodePtr();\
    _leaf_.nodePtr = _nodePtr_;\
    /*create unitPtr*/\
    int _n_nodes_ = _nodePtr_.size()-1;\
    std::vector<int> _hopelessNodePtr_ = _stage_mtxPower_.getHopelessNodePtr();\
    _leaf_.unitNodePtr.push_back(0);\
    for(int n=0; n<_n_nodes_; ++n)\
    {\
        int levelNodeStart = _nodePtr_[n];\
        int levelNodeEnd = _nodePtr_[n+1];\
        _leaf_.unitPtr.push_back(levelNodeStart);\
        for(int h=_hopelessNodePtr_[n]; h<_hopelessNodePtr_[n+1]; ++h)\
        {\
            _leaf_.unitPtr.push_back(_leaf_.hopelessRegions[h]); \
            _leaf_.unitPtr.push_back(_leaf_.hopelessRegions[h]+1); \
        }\
        if(_leaf_.unitPtr.back() != levelNodeEnd) /*if there is no hopeless in last unit, make a dummy*/\
        {\
            _leaf_.unitPtr.push_back(levelNodeEnd);\
        }\
        _leaf_.unitNodePtr.push_back(_leaf_.unitPtr.size());\
        if(n==(_n_nodes_-1)) \
        {\
            _leaf_.unitPtr.push_back(levelNodeEnd);\
        }\
    }

#define COPY_TO_PLAIN_INT_PTR(_vec_, _plain_)\
    _plain_ = new int[_vec_.size()];\
    for(int i=0; i<(int)_vec_.size(); ++i)\
    {\
        _plain_[i] = _vec_[i];\
    }\

void mtxPowerRecursive::recursivePartition(int parentIdx)
{
    MPLeaf parentLeaf = tree[parentIdx];
    std::vector<int> parentHopelessRegions = parentLeaf.hopelessRegions;
    int n_hopeless = (int)(parentHopelessRegions.size());
    int parentStage = parentLeaf.stage;

    for(int h=0; h<n_hopeless; ++h)
    {
        MPLeaf curLeaf;
        curLeaf.parent = parentIdx;
        curLeaf.stage = parentStage + 1;
        int hopelessStart = parentHopelessRegions[h];
        printf("@@@ hopless = %d,%d\n", hopelessStart, hopelessStart+1);
        if(parentStage>0)
        {
            curLeaf.nodeId = parentLeaf.nodeId;
        }
        else
        {
            //search and find
            std::vector<int> parentNodePtr = parentLeaf.nodePtr;
            for(int n=0; n<(int)(parentNodePtr.size()-1); ++n)
            {
                if( (parentNodePtr[n] <= hopelessStart) && (parentNodePtr[n+1] > hopelessStart) )
                {
                    curLeaf.nodeId = n;
                    break;
                }
            }
        }
        curLeaf.range = std::vector<int>{parentLeaf.levelPtr[hopelessStart], parentLeaf.levelPtr[hopelessStart+1]};
        printf("@@@ range = %d,%d\n", curLeaf.range[0], curLeaf.range[1]);
        mtxPower curStage(graph, highestPower, 1, cacheSize, safetyFactor, get_cache_violation_cutoff(curLeaf.stage), curLeaf.range[0], curLeaf.range[1]);
        curStage.findPartition();
        int *perm_curStage;
        READ_STAGE_DATA(curLeaf, curStage);
        //permute
        updatePerm(&perm, perm_curStage, totalRows, totalRows);
        delete[] perm_curStage;
        //push leaf to tree
        tree.push_back(curLeaf);
        //make this a child of parent
        int curIdx = (int)tree.size()-1;
        tree[parentIdx].children.push_back(curIdx);

        if(!curLeaf.hopelessRegions.empty())
        {
            printf("recursively calling for parent = %d\n", curIdx);
            recursivePartition(curIdx);
        }
    }
}

void mtxPowerRecursive::findPartition()
{
    MPLeaf curLeaf;
    curLeaf.parent = -1;
    curLeaf.nodeId = -1;
    curLeaf.range = std::vector<int>{0, graph->NROW};
    curLeaf.stage = 0;
    //partition for first stage
    mtxPower curStage(graph, highestPower, numSharedCache, cacheSize, safetyFactor, get_cache_violation_cutoff(curLeaf.stage), curLeaf.range[0], curLeaf.range[1]);
    curStage.findPartition();
    int* perm_curStage;
    READ_STAGE_DATA(curLeaf, curStage);
    hopelessNodePtr = curStage.getHopelessNodePtr();
    updatePerm(&perm, perm_curStage, totalRows, totalRows);
    delete[] perm_curStage;
    tree.push_back(curLeaf);

    if(!curLeaf.hopelessRegions.empty())
    {
        recursivePartition(0); //my curId=0
    }

    UPDATE_INVPERM(invPerm, perm);
    //totalLevel = (int)(curLeaf.levelPtr.size()-1);
    //if empty this is final
    /*COPY_TO_PLAIN_INT_PTR(curLeaf.levelPtr, levelPtr);
    COPY_TO_PLAIN_INT_PTR(nodePtr_vec, nodePtr);
    COPY_TO_PLAIN_INT_PTR(curLeaf.unlockRow, unlockRow);
    COPY_TO_PLAIN_INT_PTR(curLeaf.unlockCtr, unlockCtr);
    COPY_TO_PLAIN_INT_PTR(curLeaf.dangerRow, dangerRow);
*/

    printTree();
}

void mtxPowerRecursive::printTree()
{
    for(int i=0; i<(int)tree.size(); ++i)
    {
        MPLeaf curLeaf = tree[i];
        printf("%d. Range:[%d, %d], ", i, curLeaf.range[0], curLeaf.range[1]);
        printf("Children:[");
        for(int j=0; j<(int)curLeaf.children.size(); ++j)
        {
            printf("%d, ", curLeaf.children[j]);
        }
        printf("], Parent: %d, node: %d, stage: %d, cache_cutoff: %d, hopelessRegions: [",
                curLeaf.parent, curLeaf.nodeId, curLeaf.stage, get_cache_violation_cutoff(curLeaf.stage));
        for(int j=0; j<(int)curLeaf.hopelessRegions.size(); ++j)
        {
            printf("%d:%d, ", curLeaf.levelPtr[curLeaf.hopelessRegions[j]], curLeaf.levelPtr[curLeaf.hopelessRegions[j]+1]);
        }
        printf("], unitPtr: [");
        for(int j=0; j<(int)curLeaf.unitPtr.size(); ++j)
        {
            printf("%d ", curLeaf.unitPtr[j]);
        }
        printf("], unitNodePtr: [");
        for(int j=0; j<(int)curLeaf.unitNodePtr.size(); ++j)
        {
            printf("%d ", curLeaf.unitNodePtr[j]);
        }
        printf("]\n");
    }
}

/*int mtxPowerRecursive::getTotalLevel()
{
    return totalLevel;
}
*/
/*
int mtxPowerRecursive::getTotalNodes()
{
    return numSharedCache;
}
*/
/*
int*  mtxPowerRecursive::getLevelPtrRef()
{
    return levelPtr;
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
*/
/*
std::vector<int> getNodePtr()
{
    return nodePtr;
}
*/
/*
std::vector<int> getHopelessNodePtr()
{
    return hopelessNodePtr;
}*/

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
