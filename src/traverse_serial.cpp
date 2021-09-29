#include "traverse.h"
#include "utility.h"
#include <set>
#include <omp.h>

std::map<int, LevelData> Traverse::cachedData;
Traverse::Traverse(Graph *graph_, RACE::dist dist_, int rangeLo_, int rangeHi_, int parentIdx_, int numRoots_, std::vector<std::map<int, std::vector<Range>>> boundaryRange_, std::string mtxType_):graph(graph_),dist(dist_), rangeLo(rangeLo_),rangeHi(rangeHi_),parentIdx(parentIdx_), numRoots(numRoots_), graphSize(graph_->graphData.size()),distFromRoot(NULL),perm(NULL),invPerm(NULL), boundaryRange(boundaryRange_), levelData(NULL), mtxType(mtxType_)
{
    if( (mtxType != "N") && ( (mtxType != "L" && mtxType != "U") ) )
    {

        ERROR_PRINT("Matrix type %s does not exist. Available options are: N, L, or U", mtxType.c_str());
        return;
    }

    if(rangeHi == -1)
    {
        rangeHi = graph->NROW;
    }

    colRangeLo = rangeLo;
    colRangeHi = rangeHi;

    distFromRoot = new int[graphSize];
    for(int i=0; i<graphSize; ++i) {
        distFromRoot[i] = -1;
    }

    perm = new int[graph->NROW];
    for(int i=0; i<graph->NROW; ++i) {
        perm[i] = i;
    }

    invPerm = new int[graph->NROW];

    //totalNodesIncBoundary=(rangeHi-rangeLo);//the main nodes
    if(dist == RACE::POWER)
    {
        INIT_BOUNDARY_STRUCTURE(boundaryRange, boundaryLevelData, NULL);
        //EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange, totalNodesIncBoundary += (_val_.hi - _val_.lo); );
    }

    //printf("Num boundary regions = %d\n", (int)boundaryRange.size());
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

    if(dist == RACE::POWER)
    {
        EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryLevelData,
                if(_val_!=NULL)
                {
                    delete _val_;
                }
        );
    }
}

//@brief This function finds children and assigns proper distance to them
//@in parent_id
//@out marked children and grandchildren id's
std::vector<int> Traverse::markChildren(int currChild, int currLvl)
{
    std::vector<int> nxtChildren;//next children

    bool flag_child = false;

    if( (rangeLo<=currChild) && (rangeHi>currChild) )
    {
        flag_child = true;
        if(distFromRoot[currChild] == -1) {
            Counter::add(); //add counter only if in main region
        }
    }
    else if(dist==RACE::POWER)
    {
        EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange,
                if((_val_.lo <= currChild) && (_val_.hi > currChild))
                {
                flag_child=true;
                _numBoundaries_ = -1; //break from regions
                _wbl_ = -1; //break from working radius
                _mapIter_ = boundaryRange[_workingRadius_].end();//break from radius
                --_mapIter_;
                break;
                }
                );
    }

    if( flag_child) {
        if(distFromRoot[currChild] == -1) {
           // Counter::add();
            distFromRoot[currChild] = currLvl;
            std::vector<int> *grandChildren = &(graph->graphData[currChild].children);
            //Push even if its not within Range; since dist 2 coloring 
            //requires them
            for(unsigned childIdx=0; childIdx < grandChildren->size(); ++childIdx) {
                int grandChild = grandChildren->at(childIdx);
                if(distFromRoot[grandChild] == -1) {
                    nxtChildren.push_back(grandChildren->at(childIdx));
                }
                colRangeLo = std::min(colRangeLo, grandChild);
                colRangeHi = std::max(colRangeHi, grandChild);
            }

        }
    } else if(dist==RACE::TWO) {
        std::vector<int> *grandChildren = &(graph->graphData[currChild].children);
        for(unsigned childIdx=0; childIdx < grandChildren->size(); ++childIdx) {
            int grandChild = grandChildren->at(childIdx);
            if( (distFromRoot[grandChild] == -1) && ((rangeLo <= grandChild) && (rangeHi > grandChild)) ) {
                nxtChildren.push_back(grandChild);
            }
            colRangeLo = std::min(colRangeLo, grandChild);
            colRangeHi = std::max(colRangeHi, grandChild);
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

    if(mtxType == "N")
    {
        //traverse only if level has not been cached
        /*    if(cachedData.find(parentIdx) != cachedData.end())	
              {
              printf("Retrieving from cache\n");
              (*levelData) = cachedData[parentIdx];
              }
              else*/
        bool marked_all = false;
        int root = rangeLo;

        std::vector<int> currChildren;
        //currChildren.push_back(root); //The root
        for(int i=0; i<numRoots; ++i)
        {
            currChildren.push_back(root+i);
        }
        int currLvl = 0;
        Counter::reset();

        int prevIslandStop = rangeLo;
        while(!marked_all)
        {
            //printf("in while loop counter::val = %d\n", Counter::val);
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

            /*if(dist==RACE::POWER)
              {
              if( marked_all==true && (Counter::val != totalNodesIncBoundary) ) {
              printf("Not yet full nodes: counter=%d, totalNodesIncBoundary = %d\n", Counter::val, totalNodesIncBoundary);

              for(int i=rangeLo; i<rangeHi; ++i) {
              if(distFromRoot[i] == -1) {
            //Found him, mark him as distance 2 apart
            currLvl += 1;
            currChildren.push_back(i);
            marked_all = false;
            break;
            }
            }

            //if still missing
            if(marked_all == true)
            {
            EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange,
            for(int i=_val_.lo; i<_val_.hi; ++i) {
            if(distFromRoot[i] == -1) {
            printf("I am at island stuff for region [%d, %d]\n", _val_.lo, _val_.hi);
            currLvl += 1;
            currChildren.push_back(i);
            marked_all = false;
            _numBoundaries_ = -1; //break from region loop
            _wbl_ = -1; //break from workingRadius loop
            break; //this breaks from current radius loop
            }
            }
            );
            }
            }
            }*/
            // else
            {
                if( marked_all==true && (Counter::val != (rangeHi-rangeLo)) ) {
                    //now process islands
                    for(int i=prevIslandStop; i<rangeHi; ++i) {
                        if(distFromRoot[i] == -1) {
                            //Found him, mark him as distance 2 apart
                            if(dist != RACE::POWER) //do this only for coloring, else it might make levels with 0 elements for power kernel
                            {
                                currLvl += 1;
                            }
                            currChildren.push_back(i);
                            marked_all = false;
                            //printf("Island stuff found one\n");
                            prevIslandStop=i;
                            break;
                        }
                    }
                    if(marked_all)
                    {
                        ERROR_PRINT("Some nodes went missing");
                        exit(-1);
                    }
                }
            }
        }

        printf("Range = [%d, %d], ColRange = [%d, %d]\n", rangeLo, rangeHi, colRangeLo, colRangeHi);
        levelData->totalLevel = currLvl;
    }
    else
    {
        if(mtxType == "L")
        {
            //forward
            for(int i=0; i<graph->NROW; ++i)
            {
                distFromRoot[i] = i;
            }
        }
        else if(mtxType == "U")
        {
            //backward
            for(int i=0; i<graph->NROW; ++i)
            {
                distFromRoot[i] = graph->NROW-i;
            }
        }
        levelData->totalLevel = graph->NROW;
    }



    printf("Total Level = %d\n",levelData->totalLevel);

    createLevelData();
    printf("created Level Data\n");
    permuteGraph();
    printf("permuted graph\n");
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
    bool untouchedBoundaryNodes = false;
    //handle bondary levels in case its power calculation
    if(dist == RACE::POWER)
    {
        //set distFromRoot at active boundaries to totalLevel+1, so if it is not
        //touched from main region, it becomes the last level
        EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange,
                for(int node=_val_.lo; node<_val_.hi; ++node)
                {
                    if(distFromRoot[node] == -1)
                    {
                        distFromRoot[node] = levelData->totalLevel;
                        untouchedBoundaryNodes = true;
                    }
                }
            );

        //increase total levels if there are untouched boundary nodes
        if(untouchedBoundaryNodes)
        {
            levelData->totalLevel += 1;
        }
        EXEC_BOUNDARY_STRUCTURE(boundaryLevelData,
                _val_ = new LevelData;
                err_flag = RACE_SUCCESS;
                Range curRange = boundaryRange[_workingRadius_][_radius_][_region_];
                err_flag = findLevelData(curRange.lo, curRange.hi, levelData->totalLevel, _val_);
                if(err_flag != RACE_SUCCESS)
                {
                    ERROR_PRINT("Something went wrong in levelData calculation for boundaries");
                    return err_flag;
                }
                boundaryLevelData[_workingRadius_][_radius_][_region_] = _val_;
            );
    }

    //levelData for main body (region)
    err_flag = findLevelData(rangeLo, rangeHi, levelData->totalLevel, levelData);
    if(err_flag != RACE_SUCCESS)
    {
        ERROR_PRINT("Something went wrong in levelData calculation");
        return err_flag;
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

    Graph permutedGraph(*(graph));
    if(dist==RACE::POWER)
    {
        EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange, regionRange.push_back(_val_.lo); regionRange.push_back(_val_.hi); numRegions++;);
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

        //Permute rows
        for(int i=targetRangeLo; i<targetRangeHi; ++i)
        {
            permutedGraph.graphData[i].children = graph->graphData[perm[i]].children;
        }

    }


    //Permute rows : moved inside previous loop, I dont see any problem now.
    //TODO: verify
    /*for(int i=minRange; i<maxRange; ++i)
    {
        permutedGraph.graphData[i].children = graph->graphData[perm[i]].children;
    }*/

    //Permute columns
    //TODO only neighbors
    for(int i=0/*rangeLo*/; i<graph->NROW/*rangeHi*/; ++i)
    {
        std::vector<int> *children = &(permutedGraph.graphData[i].children);
        for(unsigned j=0; j<children->size(); ++j)
        {
            int child = children->at(j);
            bool inNodesIncBoundaries = false;
            if( (child>=rangeLo) && (child<rangeHi) ) {
                inNodesIncBoundaries = true;
            }
            if((RACE::POWER) && (!inNodesIncBoundaries))
            {
                EXEC_BOUNDARY_STRUCTURE_wo_radius(boundaryRange,
                        if(( (child>=_val_.lo) && (child<_val_.hi) ))
                        {
                            inNodesIncBoundaries = true;
                            _numBoundaries_ = -1; //break from regions
                            _wbl_ = -1; //break from working radius
                            _mapIter_ = boundaryRange[_workingRadius_].end();//break from radius
                            --_mapIter_;
                            break;
                        }
                    );
            }
            if(inNodesIncBoundaries)
            {
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

std::vector<std::map<int, std::vector<LevelData*>>> Traverse::getBoundaryLevelData()
{
    std::vector<std::map<int, std::vector<LevelData*>>> retBoundaryLevelData = boundaryLevelData;
    EXEC_BOUNDARY_STRUCTURE(boundaryLevelData,
            _val_ = NULL;
            boundaryLevelData[_workingRadius_][_radius_][_region_] = _val_;
            );
    return retBoundaryLevelData;
}
