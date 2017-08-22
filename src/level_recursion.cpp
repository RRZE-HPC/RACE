#include "level_recursion.h"
#include "utility.h"

LevelRecursion::LevelRecursion(Graph* graph_, int requestNThreads_, dist_t dist_, d2Method d2Type_):graph(graph_), dist(dist_), d2Type(d2Type_), requestNThreads(requestNThreads_), perm(NULL), invPerm(NULL)
{
    zoneTree = new ZoneTree(dist, d2Type);
    perm = new int[graph->NROW];
    invPerm = new int[graph->NROW];

    for(int i=0; i<graph->NROW; ++i)
    {
        perm[i] = i;
        invPerm[i] = i;
    }

    double default_eff = 40;
    char *lvlEff = getenv("NAME_EFFICIENCY");
    char *lvlThreads = getenv("NAME_THREADS");

    if(lvlThreads != NULL)
    {
        char* token = strtok(lvlThreads, ",");
        while(token != NULL)
        {
            lvl_threads.push_back(atoi(token));
            token = strtok(NULL, ",");
        }
    }

    printf("NAME_THREADS = \n");
    for(unsigned i=0; i<lvl_threads.size(); ++i)
    {
        printf("%d\n",lvl_threads[i]);
    }


    if(lvlEff == NULL)
    {
        eff_vec.push_back(default_eff);
    } else {
        eff_vec.push_back(default_eff);
        char* token = strtok(lvlEff, ",");
        while(token != NULL)
        {
            eff_vec.push_back(atof(token));
            token = strtok(NULL, ",");
        }
    }

    printf("Efficiency Vector = \n");
    for(unsigned i=0; i<eff_vec.size(); ++i)
    {
        printf("%f\n",eff_vec[i]);
    }

    //manipulate efficiency vector so threads are explicitly created 
    for(unsigned lvl=1; lvl<lvl_threads.size()+1; ++lvl)
    {
        printf("Efficiency at lvl %d dropped due to THREAD specification\n", lvl);
        eff_vec[lvl] = 0;
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

//returns efficiency of each level
double LevelRecursion::efficiency(unsigned levelNum)
{
    if(levelNum < eff_vec.size())
    {
        return eff_vec[levelNum];
    } else {
        return eff_vec[0];
    }
}

//returns threads of each level
int LevelRecursion::lvlThreads(unsigned levelNum)
{
    if((levelNum-1) < lvl_threads.size())
    {
        return lvl_threads[levelNum-1];
    } else {
        return -1;
    }
}

ZoneTree* LevelRecursion::getZoneTree()
{
    ZoneTree* toRet = zoneTree;
    zoneTree = NULL;
    return toRet;
}

/*
void LevelRecursion::calculateIdealNthreads(int parentIdx, int currLvl)
{
    std::vector<int>* children = &(zoneTree->at(parentIdx).childrenZ);
    unsigned nBlocks = children->size()/2;
    if(lvlThreads(currLvl)==-1)
    {
        int parentThreads = zoneTree->at(parentIdx).idealNthreadsZ;
        int parentRow = zoneTree->at(parentIdx).valueZ[1] - zoneTree->at(parentIdx).valueZ[0];
        double* remainder = new double[nBlocks];
        int remainingThreads = parentThreads;
        //One could individually treat red and black, but this would lead to less locality
        //TODO : To try the method of treating red and black individually.
        for(unsigned i=0; i<nBlocks; ++i)
        {
            std::vector<int>* redRange = &(zoneTree->at(children->at(2*i)).valueZ);
            std::vector<int>* blackRange = &(zoneTree->at(children->at(2*i+1)).valueZ);
            int totalRow = blackRange->at(1)-redRange->at(0);
            double myWeight = totalRow/(static_cast<double>(parentRow));
            double myThreads_d = myWeight*parentThreads;
            //atleast 1 thread, but not greater than remainingThreads
            int myThreads = static_cast<int>(myThreads_d);
            zoneTree->at(children->at(2*i)).idealNthreadsZ = myThreads;
            zoneTree->at(children->at(2*i+1)).idealNthreadsZ = myThreads;
            //Give high priority for 0 assigned regions
            remainder[i] = (myThreads!=0)?(myThreads_d - myThreads):1;
            remainingThreads -= myThreads;
        }
        if(remainingThreads>0)
        {
            int* rankPerm = new int[nBlocks];
            for(unsigned i=0; i<nBlocks; ++i)
            {
                rankPerm[i] = i;
            }
            sortPerm(remainder, rankPerm, 0, nBlocks, true);
            int currRank = 0;
            while(remainingThreads != 0)
            {
                int target = rankPerm[currRank];
                zoneTree->at(children->at(2*target)).idealNthreadsZ += 1;
                zoneTree->at(children->at(2*target+1)).idealNthreadsZ += 1;
                remainingThreads -= 1;
                currRank++;
            }
        }
        else if(remainingThreads<0)
        {
            ERROR_PRINT("Thread assignment error");
        }
        delete[] remainder;
    }
    else
    {
        int myThreads = lvlThreads(currLvl);
        for(unsigned i=0; i<nBlocks; ++i)
        {
            zoneTree->at(children->at(2*i)).idealNthreadsZ = myThreads;
            zoneTree->at(children->at(2*i+1)).idealNthreadsZ = myThreads;
        }
    }
}
*/

void LevelRecursion::calculateIdealNthreads(int parentIdx, int currLvl)
{
    std::vector<int>* children = &(zoneTree->at(parentIdx).childrenZ);
    int blockPerThread =  getBlockPerThread(dist, d2Type);
    unsigned nThreads = children->size()/blockPerThread;

    if(lvlThreads(currLvl)==-1)
    {
        int parentThreads = zoneTree->at(parentIdx).idealNthreadsZ;
        int parentRow = zoneTree->at(parentIdx).valueZ[1] - zoneTree->at(parentIdx).valueZ[0];
        double* remainder = new double[nThreads];
        int remainingThreads = parentThreads;
        //One could individually treat red and black, but this would lead to less locality
        //TODO : To try the method of treating red and black individually.
        for(unsigned i=0; i<nThreads; ++i)
        {
            std::vector<int>* firstRange = &(zoneTree->at(children->at(blockPerThread*i)).valueZ);
            std::vector<int>* lastRange = &(zoneTree->at(children->at(blockPerThread*(i+1)-1)).valueZ);
            int totalRow = lastRange->at(1)-firstRange->at(0);
            double myWeight = totalRow/(static_cast<double>(parentRow));
            double myThreads_d = myWeight*parentThreads;
            //atleast 1 thread, but not greater than remainingThreads
            int myThreads = static_cast<int>(myThreads_d);
            for(int block=0; block<blockPerThread; ++block)
            {
                zoneTree->at(children->at(blockPerThread*i+block)).idealNthreadsZ = myThreads;
            }
            //Give high priority for 0 assigned regions
            remainder[i] = (myThreads!=0)?(myThreads_d - myThreads):1;
            remainingThreads -= myThreads;
        }
        if(remainingThreads>0)
        {
            int* rankPerm = new int[nThreads];
            for(unsigned i=0; i<nThreads; ++i)
            {
                rankPerm[i] = i;
            }
            sortPerm(remainder, rankPerm, 0, nThreads, true);
            int currRank = 0;
            while(remainingThreads != 0)
            {
                int target = rankPerm[currRank];
                for(int block=0; block<blockPerThread; ++block)
                {
                    zoneTree->at(children->at(blockPerThread*target+block)).idealNthreadsZ += 1;
                }
                remainingThreads -= 1;
                currRank++;
            }
        }
        else if(remainingThreads<0)
        {
            ERROR_PRINT("Thread assignment error");
        }
        delete[] remainder;
    }
    else
    {
        int myThreads = lvlThreads(currLvl);
        for(unsigned i=0; i<nThreads; ++i)
        {
            for(int block=0; block<blockPerThread; ++block)
            {
                zoneTree->at(children->at(blockPerThread*i+block)).idealNthreadsZ = myThreads;
            }
        }
    }
}

void LevelRecursion::recursivePartition(int parentIdx, int currLevel)
{
    calculateIdealNthreads(parentIdx, currLevel);
    std::vector<int> children = zoneTree->at(parentIdx).childrenZ;
    unsigned numChildren = children.size();
    for(unsigned i=0; i< numChildren; ++i)
    {
        int currIdx = children[i];
        int currIdxNThreads = zoneTree->at(currIdx).idealNthreadsZ;
        std::vector<int> range = zoneTree->at(currIdx).valueZ;

        //Traverse
        if(currIdxNThreads > 1)
        {
            Traverse traverse(graph, dist, range[0], range[1], currIdx);
            traverse.calculateDistance();
            int *levelPerm = NULL;
            int len;
            traverse.getPerm(&levelPerm, &len);
            updatePerm(&perm, levelPerm, len);
            delete[] levelPerm;

            LevelData* levelData = traverse.getLevelData();
            //Try to spawn all the threads required
            bool locFlag = zoneTree->spawnChild(currIdx, currIdxNThreads, levelData, efficiency(currLevel));
            delete levelData;

            //if failed to obtain threads, try going down the tree recursively
            if(!locFlag)
            {
                int subAvailableNThreads = zoneTree->at(currIdx).nthreadsZ;
                if(subAvailableNThreads == 1) {
                    zoneTree->at(currIdx).reachedLimit = true;
                }
                if(subAvailableNThreads > 1) {
                    recursivePartition(currIdx, currLevel+1);
                }
            }
        }
    }
}

void LevelRecursion::levelBalancing()
{
    //Step 1: Construct Root
    ZoneLeaf root(0,graph->NROW,-1);
    root.idealNthreadsZ = requestNThreads;
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
    int currLevel = 1;
    int lvlOneThreads = requestNThreads;
    if(lvlThreads(1)!=-1)
    {
        lvlOneThreads = lvlThreads(1);
    }
    zoneTree->spawnChild(parentIdx, lvlOneThreads, levelData, efficiency(currLevel));
    delete levelData;

    if( zoneTree->at(0).nthreadsZ != requestNThreads)
    {
        //Step 3: Recursive partition until nthreads satisfied
        recursivePartition(parentIdx, currLevel+1);
    }

    availableNThreads = zoneTree->at(0).nthreadsZ;

    if(availableNThreads != requestNThreads)
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

