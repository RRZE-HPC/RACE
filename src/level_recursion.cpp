#include "level_recursion.h"
#include "utility.h"

LevelRecursion::LevelRecursion(Graph* graph_, int requestNThreads_, dist_t dist_):graph(graph_), dist(dist_), requestNThreads(requestNThreads_), perm(NULL), invPerm(NULL)
{
    zoneTree = new ZoneTree;
    perm = new int[graph->NROW];
    invPerm = new int[graph->NROW];

    for(int i=0; i<graph->NROW; ++i)
    {
        perm[i] = i;
        invPerm[i] = i;
    }
}

LevelRecursion::~LevelRecursion()
{
    if(perm)
    {
        delete[] perm;
    }

    if(invPerm)
    {
        delete[] invPerm;
    }

    if(zoneTree)
    {
        delete zoneTree;
    }
}

ZoneTree* LevelRecursion::getZoneTree()
{
    ZoneTree* toRet = zoneTree;
    zoneTree = NULL;
    return toRet;
}

bool LevelRecursion::recursivePartition(int parentIdx, int subRequestNThreads)
{
    int obtainedNThreads = zoneTree->at(parentIdx).nthreadsZ;
    bool glbBreak = false;

    while(obtainedNThreads < subRequestNThreads)
    {
        KeyChild keyChild;
        keyChild = zoneTree->findKeyChild(parentIdx);

        //TODO 3 block  method
        const int numBlockPerThread = 2;
        bool glbBreakBlock[numBlockPerThread];

        for(int i=0; i<numBlockPerThread; ++i)
        {
            glbBreakBlock[i] = false;

            int currIdx = keyChild.indices[i];
            //spawn from currIdx
            if(zoneTree->at(currIdx).reachedLimit == false)
            {
                std::vector<int> range = zoneTree->at(currIdx).valueZ;

                //Traverse
                Traverse traverse(graph, dist, range[0], range[1]);
                traverse.calculateDistance();
                int *levelPerm = NULL;
                int len;
                traverse.getPerm(&levelPerm, &len);
                updatePerm(&perm, levelPerm, len);
                delete[] levelPerm;

                LevelData* levelData = traverse.getLevelData();
                int currIdxNThreads = zoneTree->at(currIdx).nthreadsZ;
                bool locFlag = zoneTree->spawnChild(currIdx, currIdxNThreads+1, currIdxNThreads, levelData, dist);
                delete levelData;

                if(!locFlag){
                    int subAvailableNThreads = zoneTree->at(currIdx).nthreadsZ;
                    if(subAvailableNThreads == 1) {
                        zoneTree->at(currIdx).reachedLimit = true;
                        glbBreakBlock[i] = true;
                    }
                    if(subAvailableNThreads > 1) {
                        glbBreakBlock[i] = recursivePartition(currIdx, currIdxNThreads+1);
                    }
                }
            } else {
                glbBreakBlock[i] = true;
            }
        }
        int glbBreakCtr = 0;
        //check whether all glbBreakBlock is true
        for(int i=0; i<numBlockPerThread; ++i)
        {
            if(glbBreakBlock[i] == true)
            {
                ++glbBreakCtr;
            }
        }

        if(glbBreakCtr == numBlockPerThread)
        {
            glbBreak = true;
            break;
        }

        obtainedNThreads = zoneTree->at(0).nthreadsZ;
    }

    return glbBreak;
}

void LevelRecursion::levelBalancing()
{
    //Step 1: Construct Root
    ZoneLeaf root(0,graph->NROW,-1);
    zoneTree->push_back(root);
    int parentIdx = 0;

    //Step 2: Construct level 1
    //Traverse
    Traverse traverse(graph, dist);
    traverse.calculateDistance();
    int *levelPerm = NULL;
    int len;
    traverse.getPerm(&levelPerm, &len);
    updatePerm(&perm, levelPerm, len);
    delete[] levelPerm;
    LevelData* levelData = traverse.getLevelData();
    zoneTree->spawnChild(parentIdx, requestNThreads, 1, levelData, dist);
    delete levelData;

    //Step 3: Recursive partition until nthreads satisfied
    bool glbBreak = recursivePartition(parentIdx, requestNThreads);

    availableNThreads = zoneTree->at(0).nthreadsZ;

    if(glbBreak == true)
    {
        WARNING_PRINT("Could not spawn requested threads = %d. Threads limited to %d", requestNThreads, availableNThreads);
    }

    //update invPerm
    for(int i=0; i<graph->NROW; ++i)
    {
        invPerm[perm[i]] = i;
    }
}

//Getter functions
void LevelRecursion::getPerm(int **perm_, int *len)
{
    (*perm_) = perm;
    perm = NULL;
    (*len) = graph->NROW;
}

void LevelRecursion::getInvPerm(int **invPerm_, int *len)
{
    (*invPerm_) = invPerm;
    invPerm = NULL;
    (*len) = graph->NROW;
}

int LevelRecursion::getAvailableThreads()
{
    return availableNThreads;
}

