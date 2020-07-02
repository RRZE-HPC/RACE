#include "level_recursion.h"
#include "utility.h"

LevelRecursion::LevelRecursion(Graph* graph_, int requestNThreads_, RACE::dist dist_, RACE::d2Method d2Type_, RACE::LBTarget lbTarget_):graph(graph_), dist(dist_), d2Type(d2Type_), lbTarget(lbTarget_), requestNThreads(requestNThreads_), perm(NULL), invPerm(NULL)
{
    zoneTree = new ZoneTree(dist, d2Type, lbTarget);
    int totalRows;
    graph->getInitialPerm(&perm, &totalRows);
    invPerm = new int[totalRows];

    for(int i=0; i<totalRows; ++i)
    {
        invPerm[perm[i]] = i;
    }

    double default_eff = 40;
/*    char *lvlEff = getenv("RACE_EFFICIENCY");
    char *lvlThreads = getenv("RACE_THREADS");

    if(lvlThreads != NULL)
    {
        char* token = strtok(lvlThreads, ",");
        while(token != NULL)
        {
            lvl_threads.push_back(atoi(token));
            token = strtok(NULL, ",");
        }
    }
*/
    getEnv("RACE_THREADS", lvl_threads);
    printf("RACE_THREADS = \n");
    for(unsigned i=0; i<lvl_threads.size(); ++i)
    {
        printf("%d\n",lvl_threads[i]);
    }

    eff_vec.push_back(default_eff);
    getEnv("RACE_EFFICIENCY", eff_vec);

    /*
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
    }*/

    printf("Efficiency Vector = \n");
    for(unsigned i=0; i<eff_vec.size(); ++i)
    {
        printf("%f\n",eff_vec[i]);
    }

    //manipulate efficiency vector so threads are explicitly created 
    for(unsigned lvl=1; lvl<lvl_threads.size()+1; ++lvl)
    {
        printf("Efficiency at lvl %d dropped due to THREAD specification\n", (int)lvl);
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

void LevelRecursion::calculateIdealNthreads(int parentIdx, int parentSubIdx, int currLvl)
{
    std::vector<int>* children = &(zoneTree->at(parentIdx).children);
    int blockPerThread =  getBlockPerThread(dist, d2Type);
    int nThreads = (zoneTree->at(parentIdx).children[2*parentSubIdx+1]-zoneTree->at(parentIdx).children[2*parentSubIdx])/blockPerThread;

    if(lvlThreads(currLvl)==-1)
    {
        int parentThreads = zoneTree->at(parentIdx).idealNthreadsZ;
        int startRow = zoneTree->at(parentIdx).valueZ[parentSubIdx];
        int endRow = zoneTree->at(parentIdx).valueZ[parentSubIdx+1];
        int parentRow = endRow - startRow;
        int parentNnz = 0;
        if(lbTarget == RACE::NNZ)
        {

            for(int row=startRow; row<endRow; ++row)
            {
                parentNnz += graph->at(row).children.size();
            }
        }

        double* remainder = new double[nThreads];
        int remainingThreads = parentThreads;
        //One could individually treat red and black, but this would lead to less locality
        //TODO : To try the method of treating red and black individually.
        for(int i=0; i<nThreads; ++i)
        {
            std::vector<int>* firstRange = &(zoneTree->at(children->at(2*parentSubIdx)+i*blockPerThread).valueZ);
            std::vector<int>* lastRange = &(zoneTree->at(children->at(2*parentSubIdx)+blockPerThread*(i+1)-1).valueZ);
            int startSubRow = firstRange->front();
            int endSubRow = lastRange->back();
            int totalRow = endSubRow - startSubRow;

            double myWeight = totalRow/(static_cast<double>(parentRow));
            if(lbTarget == RACE::NNZ)
            {
                int totalNnz = 0;
                for(int row=startSubRow; row<endSubRow; ++row)
                {
                    totalNnz += (int)(graph->at(row).children.size());
                }

                myWeight = totalNnz/(static_cast<double>(parentNnz));
            }

            double myThreads_d = myWeight*parentThreads;
            //atleast 1 thread, but not greater than remainingThreads
            int myThreads = static_cast<int>(myThreads_d);
            for(int block=0; block<blockPerThread; ++block)
            {
                zoneTree->at(children->at(2*parentSubIdx)+blockPerThread*i+block).idealNthreadsZ = myThreads;
            }
            //Give high priority for 0 assigned regions
            remainder[i] = (myThreads!=0)?(myThreads_d - myThreads):1;
            remainingThreads -= myThreads;
        }

        if(remainingThreads>0)
        {
            int* rankPerm = new int[nThreads];
            for(int i=0; i<nThreads; ++i)
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
                    zoneTree->at(children->at(2*parentSubIdx)+blockPerThread*target+block).idealNthreadsZ += 1;
                }
                remainingThreads -= 1;
                currRank++;
            }

            delete[] rankPerm;
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
        for(int i=0; i<nThreads; ++i)
        {
            for(int block=0; block<blockPerThread; ++block)
            {
                zoneTree->at(children->at(2*parentSubIdx)+blockPerThread*i+block).idealNthreadsZ = myThreads;
            }
        }
    }
}

void LevelRecursion::recursivePartition(int parentIdx, int parentSubIdx, int currLevel)
{
    calculateIdealNthreads(parentIdx, parentSubIdx, currLevel);
    std::vector<int> children = zoneTree->at(parentIdx).children;

    for(int currIdx=children[2*parentSubIdx]; currIdx<children[2*parentSubIdx+1]; ++currIdx)
    {
        int currIdxNThreads = zoneTree->at(currIdx).idealNthreadsZ;
        std::vector<int> range = zoneTree->at(currIdx).valueZ;

        //Traverse
        if(currIdxNThreads > 1)
        {
            for(unsigned j=0; j<range.size()-1; ++j)
            {
                Traverse traverse(graph, dist, range[j], range[j+1], currIdx);
                traverse.calculateDistance();
                int *levelPerm = NULL;
                int len;
                traverse.getPerm(&levelPerm, &len);
                updatePerm(&perm, levelPerm, len, len+graph->NROW_serial);
                delete[] levelPerm;

                LevelData* levelData = traverse.getLevelData();
                //Try to spawn all the threads required
                bool locFlag = zoneTree->spawnChild(currIdx, j, currIdxNThreads, levelData, efficiency(currLevel));

                delete levelData;

                //if failed to obtain threads, try going down the tree recursively
                if(!locFlag)
                {
                    //This might not be updated with blocking
                    //int subAvailableNThreads = zoneTree->at(currIdx).nthreadsZ;
                    int blockPerThread =  getBlockPerThread(dist, d2Type);
                    int subAvailableNThreads = (zoneTree->at(currIdx).children[2*j+1]-zoneTree->at(currIdx).children[2*j])/blockPerThread;
                    if(subAvailableNThreads == 1) {
                        zoneTree->at(currIdx).reachedLimit = true;
                    }
                    if(subAvailableNThreads > 1) {
                        recursivePartition(currIdx, j, currLevel+1);
                    }
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
    root.totalSubBlocks = 1;
    zoneTree->push_back(root);
    int parentIdx = 0;

    //Step 2: Construct level 1
    //Traverse
    Traverse traverse(graph, dist);
    traverse.calculateDistance();
    int *levelPerm = NULL;
    int len;
    traverse.getPerm(&levelPerm, &len);
    updatePerm(&perm, levelPerm, len, len+graph->NROW_serial);
    delete[] levelPerm;
    LevelData* levelData = traverse.getLevelData();
    int currLevel = 1;
    int lvlOneThreads = requestNThreads;
    if(lvlThreads(1)!=-1)
    {
        lvlOneThreads = lvlThreads(1);
    }
    zoneTree->spawnChild(parentIdx, 0, lvlOneThreads, levelData, efficiency(currLevel));
    delete levelData;

    if( zoneTree->at(0).nthreadsZ != requestNThreads)
    {
        //Step 3: Recursive partition until nthreads satisfied
        recursivePartition(parentIdx, 0, currLevel+1);
    }

    availableNThreads = zoneTree->at(0).nthreadsZ;

    if(availableNThreads != requestNThreads)
    {
        WARNING_PRINT("Could not spawn requested threads = %d. Threads limited to %d", requestNThreads, availableNThreads);
    }

    //update invPerm
    for(int i=0; i<graph->NROW+graph->NROW_serial; ++i)
    {
        invPerm[perm[i]] = i;
    }
}

//Getter functions
void LevelRecursion::getPerm(int **perm_, int *len)
{
    (*perm_) = perm;
    perm = NULL;
    (*len) = graph->NROW + graph->NROW_serial;
}

void LevelRecursion::getInvPerm(int **invPerm_, int *len)
{
    (*invPerm_) = invPerm;
    invPerm = NULL;
    (*len) = graph->NROW + graph->NROW_serial;
}

int LevelRecursion::getAvailableThreads()
{
    return availableNThreads;
}

