#include "traverse.h"
#include "utility.h"
#include <set>
#include <omp.h>

std::map<int, LevelData> Traverse::cachedData;
Traverse::Traverse(Graph *graph_, RACE::dist dist_, int rangeLo_, int rangeHi_, int parentIdx_, int numRoots_, std::vector<int> negativeBoundary_, std::vector<int> positiveBoundary_):graph(graph_),dist(dist_), rangeLo(rangeLo_),rangeHi(rangeHi_),parentIdx(parentIdx_), numRoots(numRoots_), graphSize(graph_->graphData.size()),distFromRoot(NULL),perm(NULL),invPerm(NULL), negativeBoundary(negativeBoundary_), positiveBoundary(positiveBoundary_), levelData(NULL)
{
    if(rangeHi == -1)
    {
        rangeHi = graph->NROW;
    }

    distFromRoot = new int[graphSize];
    for(int i=0; i<graphSize; ++i) {
        distFromRoot[i] = -1;
    }

    perm = new int[graph->NROW];
    for(int i=0; i<graph->NROW; ++i) {
        perm[i] = i;
    }

    invPerm = new int[graph->NROW];

    int numNegativeBoundary = negativeBoundary.empty()?0:negativeBoundary.size()-1;
    int numPositiveBoundary = positiveBoundary.empty()?0:positiveBoundary.size()-1;
    printf("Num negative boundary = %d, positive boundary = %d\n", numNegativeBoundary, numPositiveBoundary);
    levelDataNegativeBoundary.resize(numNegativeBoundary, NULL);
    levelDataPositiveBoundary.resize(numPositiveBoundary, NULL);
    //for main (target) region
    levelData = new LevelData;
    UNUSED(parentIdx);

}

Traverse::~Traverse()
{
    if(distFromRoot) {
        delete[] distFromRoot;
    }

    if(perm) {
        delete[] perm;
    }

    if(invPerm) {
        delete[] invPerm;
    }

    if(levelData) {
        delete levelData;
    }

    int numNegativeBoundary = negativeBoundary.empty()?0:negativeBoundary.size()-1;
    int numPositiveBoundary = positiveBoundary.empty()?0:positiveBoundary.size()-1;

    for(int b=0; b<numNegativeBoundary; ++b)
    {
        if(levelDataNegativeBoundary[b]) {
            delete levelDataNegativeBoundary[b];
        }
    }
    for(int b=0; b<numPositiveBoundary; ++b)
    {
        if(levelDataPositiveBoundary[b]) {
            delete levelDataPositiveBoundary[b];
        }
    }
}

//@brief This function finds children and assigns proper distance to them
//@in parent_id
//@out marked children and grandchildren id's
std::vector<int> Traverse::markChildren(int currChild, int currLvl)
{
    std::vector<int> nxtChildren;//next children
    int targetRangeLo = rangeLo;
    int targetRangeHi = rangeHi;

    if(dist==RACE::POWER) {
        targetRangeLo = negativeBoundary.empty()?rangeLo:negativeBoundary.back();
        targetRangeHi = positiveBoundary.empty()?rangeHi:positiveBoundary.back();
    }
    if( (targetRangeLo<=currChild) && (targetRangeHi>currChild) ) {
        if(distFromRoot[currChild] == -1) {
            Counter::add();
            distFromRoot[currChild] = currLvl;
            std::vector<int> *grandChildren = &(graph->graphData[currChild].children);
            //Push even if its not within Range; since dist 2 coloring 
            //requires them
            for(unsigned childIdx=0; childIdx < grandChildren->size(); ++childIdx) {
                int grandChild = grandChildren->at(childIdx);
                if(distFromRoot[grandChild] == -1) {
                    nxtChildren.push_back(grandChildren->at(childIdx));
                }
            }

        }
    } else if(dist==RACE::TWO) {
        std::vector<int> *grandChildren = &(graph->graphData[currChild].children);
        for(unsigned childIdx=0; childIdx < grandChildren->size(); ++childIdx) {
            int grandChild = grandChildren->at(childIdx);
            if( (distFromRoot[grandChild] == -1) && ((rangeLo <= grandChild) && (rangeHi > grandChild)) ) {
                nxtChildren.push_back(grandChild);
            }
        }
    }

    return nxtChildren;
}

int Counter::val = 0;

void Counter::add()
{
#pragma omp atomic update
    val++;
}

void Counter::reset()
{
    val = 0;
}

void Traverse::calculateDistance()
{
    //traverse only if level has not been cached
    /*    if(cachedData.find(parentIdx) != cachedData.end())	
          {
          printf("Retrieving from cache\n");
          (*levelData) = cachedData[parentIdx];
          }
          else*/
    {
        bool marked_all = false;
        int root = rangeLo;
        int targetRangeLo = rangeLo;
        int targetRangeHi = rangeHi;

        if(dist==RACE::POWER) {
            targetRangeLo = negativeBoundary.empty()?rangeLo:negativeBoundary.back();
            targetRangeHi = positiveBoundary.empty()?rangeHi:positiveBoundary.back();
        }

        std::vector<int> currChildren;
        //currChildren.push_back(root); //The root
        for(int i=0; i<numRoots; ++i)
        {
            currChildren.push_back(root+i);
        }
        int currLvl = 0;
        Counter::reset();

        while(!marked_all)
        {
            marked_all = true;
            std::set<int> children;
            //in parallel
            //	#pragma omp parallel for  shared(marked_all)
            for(unsigned i=0; i<currChildren.size(); ++i) {
                std::vector<int> nxtChildren = markChildren(currChildren[i],currLvl);
                //atomic needed

                if(!(nxtChildren.empty()))
                {
                    //		#pragma omp critical 
                    {
                        marked_all = false;
                        for(unsigned j=0; j<nxtChildren.size(); ++j) {
                            children.insert(nxtChildren[j]);
                        }
                    }
                }

                //now assign current children as parent
            }

            currChildren.assign(children.begin(),children.end());
            currLvl += 1;

            if( marked_all==true && (Counter::val != (targetRangeHi-targetRangeLo)) ) {
                //now process islands
                for(int i=targetRangeLo; i<targetRangeHi; ++i) {
                    if(distFromRoot[i] == -1) {
                        //Found him, mark him as distance 2 apart
                        currLvl += 1;
                        currChildren.push_back(i);
                        marked_all = false;
                        break;
                    }
                }
            }
        }

        levelData->totalLevel = currLvl;
        printf("Total Level = %d\n",levelData->totalLevel);

        createLevelData();
        printf("created Level Data\n");
        permuteGraph();
        printf("permuted graph\n");
    }
}

RACE_error Traverse::findLevelData(int lower_nrows, int upper_nrows, int totalLevel, LevelData* curLevelData)
{
    curLevelData->totalLevel = totalLevel;
    int* levelRow_ = new int[totalLevel];
    int* levelNnz_ = new int[totalLevel];
    for(int i=0; i<totalLevel; ++i) {
        levelRow_[i] = 0;
        levelNnz_[i] = 0;
    }

    for(int i=lower_nrows; i<upper_nrows; ++i)
    {
        int curr_dist = distFromRoot[i];
        if(curr_dist == -1)
        {
            ERROR_PRINT("There are orphan nodes; I thought this wouldn't happen");
            return RACE_ERR_GRAPH_TRAVERSAL;

        }

        levelRow_[curr_dist]+=1;
        levelNnz_[curr_dist] += graph->graphData[i].children.size();
        //levelNnz_[curr_dist] += graph->graphData[i].upperNnz;

    }

    curLevelData->levelRow = levelRow_;
    curLevelData->levelNnz = levelNnz_;

    curLevelData->nrow=0;
    curLevelData->nnz=0;
    for(int i=0; i<totalLevel; ++i)
    {
        curLevelData->nrow += levelRow_[i];
        curLevelData->nnz += levelNnz_[i];
    }
    return RACE_SUCCESS;
}

RACE_error Traverse::createLevelData()
{
    RACE_error err_flag = RACE_SUCCESS;
    //levelData for main body (region)
    err_flag = findLevelData(rangeLo, rangeHi, levelData->totalLevel, levelData);
    if(err_flag != RACE_SUCCESS)
    {
        ERROR_PRINT("Something went wrong in levelData calculation");
        return err_flag;
    }
    //handle bondary levels in case its power calculation
    if(dist == RACE::POWER)
    {
        int numNegativeBoundaries = negativeBoundary.empty()?0:(int)(negativeBoundary.size()-1);
        int numPositiveBoundaries = positiveBoundary.empty()?0:(int)(positiveBoundary.size()-1);

        //negative-boundary
        for(int b=0; b<numNegativeBoundaries; ++b)
        {
            int lower_nrows = negativeBoundary[b+1];
            int upper_nrows = negativeBoundary[b];

            LevelData* curLevelData = new LevelData;
            err_flag = RACE_SUCCESS;
            err_flag = findLevelData(lower_nrows, upper_nrows, levelData->totalLevel, curLevelData);
            if(err_flag != RACE_SUCCESS)
            {
                ERROR_PRINT("Something went wrong in levelData calculation for boundaries");
                return err_flag;
            }

            levelDataNegativeBoundary[b] = curLevelData;
        }
        //positive-boundary
        for(int b=0; b<numPositiveBoundaries; ++b)
        {
            int lower_nrows = positiveBoundary[b];
            int upper_nrows = positiveBoundary[b+1];

            LevelData* curLevelData = new LevelData;
            err_flag = RACE_SUCCESS;
            err_flag = findLevelData(lower_nrows, upper_nrows, levelData->totalLevel, curLevelData);
            if(err_flag != RACE_SUCCESS)
            {
                ERROR_PRINT("Something went wrong in levelData calculation for boundaries");
                return err_flag;
            }

            levelDataPositiveBoundary[b] = curLevelData;
        }
    }
    //cache this data for later use
    //Don't cache this won't happen in 
    //the current strategy
    //cachedData[parentIdx] = (*levelData);
    return err_flag;
}

void Traverse::permuteGraph()
{
    int numRegions = 1;
    std::vector<int> regionRange;
    regionRange.push_back(rangeLo);
    regionRange.push_back(rangeHi);
    int minRange = rangeLo;
    int maxRange = rangeHi;

    if(dist==RACE::POWER)
    {
        int numNegativeBoundaries = negativeBoundary.empty()?0:(int)(negativeBoundary.size()-1);
        int numPositiveBoundaries = positiveBoundary.empty()?0:(int)(positiveBoundary.size()-1);
        numRegions += numNegativeBoundaries + numPositiveBoundaries;

        //negative-boundary
        for(int b=0; b<numNegativeBoundaries; ++b)
        {
            int lower_nrows = negativeBoundary[b+1];
            int upper_nrows = negativeBoundary[b];
            regionRange.push_back(lower_nrows);
            regionRange.push_back(upper_nrows);
        }

        //positive-boundary
        for(int b=0; b<numPositiveBoundaries; ++b)
        {
            int lower_nrows = positiveBoundary[b];
            int upper_nrows = positiveBoundary[b+1];
            regionRange.push_back(lower_nrows);
            regionRange.push_back(upper_nrows);
        }

        minRange = negativeBoundary.empty()?rangeLo:negativeBoundary.back();
        maxRange = positiveBoundary.empty()?rangeHi:positiveBoundary.back();
    }

    for(int regionIdx=0; regionIdx<numRegions; ++regionIdx)
    {
        int targetRangeLo = regionRange[2*regionIdx];
        int targetRangeHi = regionRange[2*regionIdx+1];

        //create permutation vector First
        sortPerm(distFromRoot, perm, targetRangeLo, targetRangeHi);

        //create invPerm
        for(int i=targetRangeLo; i<targetRangeHi; ++i) {
            invPerm[perm[i]] = i;
        }
    }

    Graph permutedGraph(*(graph));


    //Permute rows
    for(int i=minRange; i<maxRange; ++i)
    {
        permutedGraph.graphData[i].children = graph->graphData[perm[i]].children;
    }

    //Permute columns
    //TODO only neighbors
    for(int i=0/*rangeLo*/; i<graph->NROW/*rangeHi*/; ++i)
    {
        std::vector<int> *children = &(permutedGraph.graphData[i].children);
        for(unsigned j=0; j<children->size(); ++j)
        {
            int child = children->at(j);
            if( (child>=minRange) && (child<maxRange) ) {
                children->at(j) = invPerm[child];
            }
        }
    }

    permutedGraph.swap(*(graph));
}

//Getter functions
void Traverse::getPerm(int **perm_, int *len)
{
    (*perm_) = perm;
    perm = NULL;
    (*len) = graph->NROW;
}

void Traverse::getInvPerm(int **invPerm_, int *len)
{
    (*invPerm_) = invPerm;
    invPerm = NULL;
    (*len) = graph->NROW;
}

LevelData* Traverse::getLevelData()
{
    LevelData* levelData_ = levelData;
    levelData = NULL;
    return levelData_;
}

std::vector<LevelData*> Traverse::getLevelDataNegativeBoundary()
{
     std::vector<LevelData*> retLevelDataNegativeBoundary = levelDataNegativeBoundary;
     for(int b=0; b<(int)levelDataNegativeBoundary.size(); ++b)
     {
         levelDataNegativeBoundary[b] = NULL;
     }
    return retLevelDataNegativeBoundary;
}

std::vector<LevelData*> Traverse::getLevelDataPositiveBoundary()
{
     std::vector<LevelData*> retLevelDataPositiveBoundary = levelDataPositiveBoundary;
     for(int b=0; b<(int)levelDataPositiveBoundary.size(); ++b)
     {
         levelDataPositiveBoundary[b] = NULL;
     }
    return retLevelDataPositiveBoundary;
}
