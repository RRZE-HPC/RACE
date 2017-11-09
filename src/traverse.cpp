#include "traverse.h"
#include "utility.h"
#include <set>
#include <omp.h>

std::map<int, LevelData> Traverse::cachedData;
Traverse::Traverse(Graph *graph_, RACE::dist dist_, int rangeLo_, int rangeHi_, int parentIdx_):graph(graph_),dist(dist_), rangeLo(rangeLo_),rangeHi(rangeHi_),parentIdx(parentIdx_),graphSize(graph_->graphData.size()),distFromRoot(NULL),perm(NULL),invPerm(NULL),levelData(NULL)
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


    levelData = new LevelData;
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
}



//@brief This function finds children and assigns proper distance to them
//@in parent_id
//@out marked children and grandchildren id's
std::vector<int> Traverse::markChildren(int currChild, int currLvl)
{
    std::vector<int> nxtChildren;//next children
    if( (rangeLo<=currChild) && (rangeHi>currChild) ) {
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
        std::vector<int> currChildren;
        currChildren.push_back(root); //The root

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

            if( marked_all==true && (Counter::val != (rangeHi-rangeLo)) ) {
                printf("We have islands in range [%d - %d]\n",rangeLo,rangeHi);
                //now process islands
                for(int i=rangeLo; i<rangeHi; ++i) {
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
        permuteGraph();
    }
}

RACE_error Traverse::createLevelData()
{
    int totalLevel = levelData->totalLevel;
    int* levelRow_ = new int[totalLevel];
    int* levelNnz_ = new int[totalLevel];
    for(int i=0; i<totalLevel; ++i) {
        levelRow_[i] = 0;
        levelNnz_[i] = 0;
    }

    for(int i=rangeLo; i<rangeHi; ++i)
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

    levelData->levelRow = levelRow_;
    levelData->levelNnz = levelNnz_;

    levelData->nrow=0;
    levelData->nnz=0;
    for(int i=0; i<totalLevel; ++i)
    {
        levelData->nrow += levelRow_[i];
        levelData->nnz += levelNnz_[i];
    }

    //cache this data for later use
    //Don't cache this won't happen in 
    //the current strategy
    //cachedData[parentIdx] = (*levelData);
    return RACE_SUCCESS;
}

void Traverse::permuteGraph()
{
    //create permutation vector First
    sortPerm(distFromRoot, perm, rangeLo, rangeHi);

    //create invPerm
    for(int i=rangeLo; i<rangeHi; ++i) {
        invPerm[perm[i]] = i;
    }

    Graph permutedGraph(*(graph));


    //Permute rows
    for(int i=rangeLo; i<rangeHi; ++i)
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
            if( (child>=rangeLo) && (child<rangeHi) ) {
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
