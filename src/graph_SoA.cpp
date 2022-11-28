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

//#include "graph_SoA.h"
#include "graph.h"
#include "error.h"
#include <cmath>
#include <set>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <limits>

// These two needed for collection of boundary nodes
#include "traverse_GAP.h"
#include "levelData.h"

#include "utility.h"
#include "config.h"
#include "type.h"

//TODO no need to put only diagonal elements in Graph
RACE_error RACE::Graph::createGraphFromCRS(int *rowPtr, int *col, int *initPerm, int *initInvPerm)
{
    // printf("I'm in Graph::createGraphFromCRS\n");
    /*if(NROW != NCOL) {
        ERROR_PRINT("NROW!=NCOL : Currently Graph BMC supports only undirected Graph");
        return RACE_ERR_MATRIX_SYMM;
    }*/


    std::vector<bool> pure_diag_flag(NROW,false);
    graphData = col;
    childrenStart = new int[NROW];
    childrenSize = new int[NROW];
    //upperNnz = new int[NROW];

    tmpGraphData = NULL;

    totalPerm = new int[NROW];
    totalInvPerm = new int[NROW];
#pragma omp parallel for schedule(runtime)
    for(int row=0; row<NROW; ++row)
    {
        if(initPerm)
        {
            totalPerm[row] = initPerm[row];
            totalInvPerm[row] = initInvPerm[row];
        }
        else
        {
            totalPerm[row] = row;
            totalInvPerm[row] = row;
        }
    }

#ifdef RACE_PERMUTE_ON_FLY
    //prohibit permuting now, since it is done on fly
    initPerm = NULL;
    initInvPerm = NULL;
#pragma omp parallel for schedule(runtime)
    for(int row=0; row<NROW; ++row)
    {
        //upperNnz[row]=0;
        int column = graphData[rowPtr[row]];
        int numChildren = rowPtr[row+1] - rowPtr[row];
        childrenSize[row] = numChildren;
        childrenStart[row] = rowPtr[row];
        if((numChildren == 1) && column == row)
        {
#pragma omp critical
            {
                pureDiag.push_back(row);
            }
        }
    }

#else
    tmpGraphData = new int[NNZ];
    int perm_nnz = 0;

    for(int row=0; row<NROW; ++row)
    {
        int permRow = row;
        if(initPerm)
        {
            permRow = initPerm[row];
        }
        int permCol = 0;
        permCol = graphData[rowPtr[permRow]];
        if(initInvPerm)
        {
            permCol = initInvPerm[permCol];
        }

        int numChildren = (rowPtr[permRow+1]-rowPtr[permRow]);
        if((numChildren==1) && (permCol == row))
        {
//#pragma omp critical
            {
                pureDiag.push_back(row);
            }
        }
        childrenSize[row] = numChildren;
        childrenStart[row] = perm_nnz;
        for(int idx=rowPtr[permRow]; idx<rowPtr[permRow+1]; ++idx) {
            permCol = graphData[idx];
            if(initInvPerm)
            {
                permCol = initInvPerm[graphData[idx]];
            }
            tmpGraphData[perm_nnz] = permCol;
            //                if(permCol >= row) // I think this is correct
            /*if(permCol >= permRow)
            {
                ++graphData[row].upperNnz;
            }*/
            ++perm_nnz;
        }
    }

    graphData = tmpGraphData;
    if(NNZ != perm_nnz)
    {
        ERROR_PRINT("NNZ of permuted matrix do not match");
        exit(-1);
    }

#endif

    if(!pureDiag.empty())
    {
        WARNING_PRINT("Graph is not connected or irreducible. For coloring this is necessary, however for matrix power kernels this will work");
    //    return RACE_ERR_NOT_IMPLEMENTED;
    }

    bool serialPermuted = getStatistics();
    if(serialPermuted)
    {
        int ctr=0;
        int ctrSerial=0;
        serialPerm.resize(NROW);
        //bring serial part to bottom and construct a graph out of rest
        for(int row=0; row<NROW; ++row)
        {
            if(std::find(serialPartRow.begin(), serialPartRow.end(), row)==serialPartRow.end())
            {
                serialPerm[ctr] = row;
                ++ctr;
            }
            else
            {
                serialPerm[(NROW) - NROW_serial + ctrSerial] = row;
                ++ctrSerial;
            }
        }
        //permute the graph
        permuteAndRemoveSerialPart();
        serialPart.resize(2);
        //Now delete the serialRows from graph
        NROW = NROW-NROW_serial;
        NNZ =  NNZ-NNZ_serial;
        serialPart[0] = NROW;
        serialPart[1] = NROW+NROW_serial;
    }
    else
    {
        serialPart.resize(2);
        serialPart[0] = NROW;
        serialPart[1] = NROW;
    }
    /*
    totalInvPerm = new int[NROW+NROW_serial];
#pragma omp parallel for schedule(static)
    for(int i=0; i<NROW+NROW_serial; ++i) {
        totalInvPerm[totalPerm[i]] = i;
    }*/

    return RACE_SUCCESS;
}

std::vector<int> RACE::Graph::collectBoundaryNodes(int powerMax){

    // These are assumes to already be "compressed" column indices
    // int localColLimitLo = 0; // lowest row, therefore, left wall of local cols
    // int localColLimitHi = NROW; // highest row, therefore, right wall of local cols

    int currentCol;

    // int *sizePtr[POWER] = {0};
    // bool isRemote; //, colInBoundaryNodes, rowInBoundaryNodes;

    // Only collects distance-1 nodes from halo elements
    for(int row = 0; row < NROW; ++row){
        for(int nzIdx=childrenStart[row]; nzIdx<childrenStart[row+1]; ++nzIdx) {
            currentCol = graphData[nzIdx]; // Extract column of this particular element 

            // isRemote = ((currentCol < localColLimitLo) || (currentCol >= localColLimitHi)); // Possibly leave in as sanity check
            // colInBoundaryNodes = std::find(collectedBoundaryNodes.begin(), collectedBoundaryNodes.end(), currentCol) != collectedBoundaryNodes.end();
            // rowInBoundaryNodes = std::find(collectedBoundaryNodes.begin(), collectedBoundaryNodes.end(), row) != collectedBoundaryNodes.end();
            // if(( isRemote || colInBoundaryNodes) && !rowInBoundaryNodes ){

            if(currentCol >= NROW) {
                boundaryNodes.push_back(row); //TODO: change to insert when unordered set
                break; // move to the next row
            }
        }
        
    }

    int startRow = 0;
    int endRow = NROW;
    int parentIdx=0;
    std::vector<std::map<int, std::vector<Range>>> boundaryRange = {};
    std::string mtxType="N";

    RACE::Traverse *traverser = new RACE::Traverse(this, RACE::POWER, startRow, endRow, parentIdx, boundaryNodes, boundaryRange, mtxType);
    traverser->calculateDistance(powerMax - 2, true);  
    LevelData* curLevelData = traverser->getLevelData();

    int totalLevel = powerMax;

    std::vector<int> distFromRemotePtr(totalLevel+1);
    distFromRemotePtr[0] = 0; //TODO: verify?

    // Takes the size of the "chunks" of how matrix is partitioned, and cumsums them
    // levelRow : # nodes in each corresponding ring
    // distFromRemotePtr -> cumsum of levelRows
    for(int level=0; level<totalLevel; ++level)
    {
        distFromRemotePtr[level+1] = distFromRemotePtr[level] + curLevelData->levelRow[level];
    }

    return distFromRemotePtr;
}

//only for debugging puposes
void RACE::Graph::writePattern(char* name)
{
#ifdef RACE_PERMUTE_ON_FLY
    WARNING_PRINT("Writing unpermuted original matrix. Writing permuted matrix not yet implemented with RACE_PERMUTE_ON_FLY");
#endif
    std::string base_name(name);
    std::stringstream fileName;
    fileName<<base_name <<".mtx";
    std::fstream file;
    file.open( fileName.str().c_str() ,std::fstream::out);
    file<<"%";
    file<<"%MatrixMarket matrix coordinate real general\n";
    file<<"%MatrixMarketm file generated by RACE\n";

    file<<NROW<<" "<<NCOL<<" "<<NNZ<<"\n";

    //TODO write permuted file too
    for(int row=0; row<NROW; ++row) {
        int* children = &(graphData[childrenStart[row]]);
        for (int j=0; j<childrenSize[row]; ++j) {
            file<<row+1<<" "<<children[j]+1<<" "<<10<<"\n";
        }
    }
}

bool RACE::Graph::getStatistics()
{

    double maxDense = std::numeric_limits<double>::max();//default 8 times mean nnzr
    char *maxDenseEnv = getenv("RACE_MAX_DENSE");

    bool permuted = false;

    if(maxDenseEnv != NULL)
    {
        maxDense = atof(maxDenseEnv);
    }

    std::vector<int> blackList;
    std::map<int,int> rowBucket;
    int MEAN = NNZ/NROW;
    //TODO to get permuted value if initPerm is there
    for(int row=0; row<NROW; ++row)
    {
        int rowLen = childrenSize[row];
        if(rowLen > maxDense*MEAN)
        {
            blackList.push_back(row);
        }
        rowBucket[rowLen] += 1;
    }
/*
    for(auto it=rowBucket.begin(); it!=rowBucket.end(); ++it)
    {
        printf("%d -> %d\n",it->first, it->second);
    }
*/
    //print black-listed rows
    if(!blackList.empty())
    {
        INFO_PRINT("Total of %lu rows separated due to high density", blackList.size());

        int nnzCtr = 0;
        //contribution of black-listed rows
        for(auto it = blackList.begin(); it!=blackList.end(); ++it)
        {
            int rowLen = childrenSize[*it];
            nnzCtr += rowLen;
        }

        NNZ_serial = nnzCtr;
        NROW_serial = blackList.size();
        INFO_PRINT("%5.2f %% is serialized", nnzCtr*100/((double)NNZ));
        permuted = true;
    }
    else
    {
        NNZ_serial = 0;
        NROW_serial = 0;
        permuted = false;
    }

    serialPartRow = blackList;

    return permuted;
}


void RACE::Graph::permuteAndRemoveSerialPart()
{
    // printf("I'm in Graph::permuteAndRemoveSerialPart\n");

    updatePerm(&totalPerm, serialPerm, NROW, NROW);
#pragma omp parallel for schedule(static)
    for(int i=0; i<NROW; ++i) {
        totalInvPerm[totalPerm[i]] = i;
    }
    WARNING_PRINT("DENSE_ROW removal will not work correctly with RACE_PERMUTE_ON_FLY switched on. This is not yet implemented");
#ifndef RACE_PERMUTE_ON_FLY
    //permute only if on fly permuting is disabled

    int *invPerm = new int[NROW];
    //create invPerm
#pragma omp parallel for schedule(static)
    for(int i=0; i<NROW; ++i) {
        invPerm[serialPerm[i]] = i;
    }


    RACE::Graph permutedGraph(*(this));

    //Permute rows
#pragma omp parallel for schedule(static)
    for(int i=0; i<NROW; ++i)
    {
        permutedGraph.childrenSize[i] = childrenSize[serialPerm[i]];
    }

    int nnz=0;
    //Permute columns
    for(int i=0; i<NROW; ++i)
    {
        int *children = &(graphData[serialPerm[i]]);
        permutedGraph.childrenStart[i] = nnz;
        for(int j=0; j<permutedGraph.childrenSize[i]; ++j)
        {
            int child = children[j];
            int invChild = invPerm[child];

            if(invChild < (NROW-NROW_serial))
            {
                permutedGraph.graphData[nnz] = invPerm[child];
                ++nnz;
            }
            //remove children in serialPartRow
        }
    }

    permutedGraph.swap(*(this));
    delete[] invPerm;
#endif
}

RACE::Graph::Graph(int nrow, int ncol, int *rowPtr, int *col, RACE::dist distance, bool symm_hint, int *initPerm, int *initInvPerm):graphData(NULL), totalPerm(NULL), totalInvPerm(NULL), NROW(nrow), NCOL(ncol)
{
    NNZ = rowPtr[NROW];
    int *outRowPtr=NULL;
    int *outCol=NULL;

    // printf("I'm in Graph::Graph, constructor\n");

    //check if a rec stage is there when performing MPK
    //if not then we can avoid creating a symmetric matrix
    //as there is only a forward dependencies.
    //TODO: might need to tweak for MPI, depending on the implementation
    std::vector<int> maxRecStagesVec;
    getEnv("RACE_MAX_RECURSION_STAGES", maxRecStagesVec);
    int maxRecStages = 10; //some value
    if(distance == RACE::POWER)
    {
        if(!maxRecStagesVec.empty())
        {
            maxRecStages = maxRecStagesVec[0];
        }
    }
    if((maxRecStages > 0) && (!symm_hint))
    {
        RACE::makeSymmetricGraph(nrow, ncol, rowPtr, col, &outRowPtr, &outCol);
    }

    if(outRowPtr) //made the graph symmetric, new structure ==> read the new graph
    {
        RACE_FN(createGraphFromCRS(outRowPtr, outCol, initPerm, initInvPerm));
    }
    else
    {
        RACE_FN(createGraphFromCRS(rowPtr, col, initPerm, initInvPerm));
    }

    manageGraphData = false;
    if(outRowPtr)
    {
        delete[] outRowPtr;
        manageGraphData = true;
    }
}

RACE::Graph::~Graph()
{
    if(tmpGraphData)
    {
        delete[] tmpGraphData;
    }
    if(childrenStart)
    {
        delete[] childrenStart;
    }
    if(childrenSize)
    {
        delete[] childrenSize;
    }
    if(totalPerm)
    {
        delete[] totalPerm;
    }
    if(totalInvPerm)
    {
        delete[] totalInvPerm;
    }
    if(manageGraphData)
    {
        delete[] graphData;
    }
}

/*RACE::Graph::Graph(const Graph& srcGraph):graphData(srcGraph.graphData),pureDiag(srcGraph.pureDiag),serialPartRow(srcGraph.serialPartRow),perm(srcGraph.perm),NROW(srcGraph.NROW),NCOL(srcGraph.NCOL), NROW_serial(srcGraph.NROW_serial), NNZ_serial(srcGraph.NNZ_serial)
{
}*/

RACE::Graph::Graph(const Graph& srcGraph)
{

    pureDiag = srcGraph.pureDiag;
    serialPartRow = srcGraph.serialPartRow;
    serialPerm = srcGraph.serialPerm;
    NROW = srcGraph.NROW;
    NCOL = srcGraph.NCOL;
    NNZ = srcGraph.NNZ;
    NROW_serial = srcGraph.NROW_serial;
    NNZ_serial = srcGraph.NNZ_serial;
    totalPerm = new int[NROW+NROW_serial];
    totalInvPerm = new int[NROW+NROW_serial];
    childrenStart = new int[NROW+NROW_serial];
    childrenSize = new int[NROW+NROW_serial];
#pragma omp parallel for schedule(static)
    for(int i=0; i<NROW+NROW_serial; ++i)
    {
        childrenSize[i] = srcGraph.childrenSize[i];
        childrenStart[i] = srcGraph.childrenStart[i];
    }

    if(srcGraph.tmpGraphData)
    {
        tmpGraphData = new int[NNZ+NNZ_serial];
        graphData = tmpGraphData;
        int nnz=0;
#pragma omp parallel for schedule(static) reduction(+:nnz)
        for(int j=0; j<NNZ+NNZ_serial; ++j)
        {
            tmpGraphData[j] = srcGraph.tmpGraphData[j];
            ++nnz;
        }
        if(nnz != NNZ)
        {
            ERROR_PRINT("Error in copying");
        }
    }
    else
    {
        graphData = srcGraph.graphData;
    }

#pragma omp parallel for schedule(static)
    for(int i=0; i<NROW+NROW_serial; ++i)
    {
        totalPerm[i] = srcGraph.totalPerm[i];
        totalInvPerm[i] = srcGraph.totalInvPerm[i];
    }
}


RACE_error RACE::Graph::swap(Graph& other)
{
    if( ((NROW != other.NROW) || (NCOL != other.NCOL)) || (NNZ != other.NNZ) ) {
        ERROR_PRINT("Incompatible Graphs");
        return RACE_ERR_INCOMPATIBILITY;
    }

    //swap graphData
    int *GraphData_tmp = graphData;
    graphData = other.graphData;
    other.graphData = GraphData_tmp;

    int *tmpGraphData_tmp = tmpGraphData;
    tmpGraphData = other.tmpGraphData;
    other.tmpGraphData = tmpGraphData_tmp;

    int* tmpChildrenStart = childrenStart;
    childrenStart = other.childrenStart;
    other.childrenStart = tmpChildrenStart;

    int* tmpChildrenSize = childrenSize;
    childrenSize = other.childrenSize;
    other.childrenSize = tmpChildrenSize;

    int* tmpTotalPerm = totalPerm;
    totalPerm = other.totalPerm;
    other.totalPerm = tmpTotalPerm;

    int* tmpInvTotalPerm = totalInvPerm;
    totalInvPerm = other.totalInvPerm;
    other.totalInvPerm = tmpInvTotalPerm;

    pureDiag.swap(other.pureDiag);
    serialPartRow.swap(other.serialPartRow);
    serialPerm.swap(other.serialPerm);

    return RACE_SUCCESS;
}

//this is permutation by removing serial part, does not
//include the initPerm passed to graph
void RACE::Graph::getSerialPerm(int **perm_, int *len_)
{
    (*perm_) = new int[NROW+NROW_serial];

    for(int i=0; i<NROW+NROW_serial; ++i)
    {
        (*perm_)[i] = serialPerm[i];
    }

    (*len_) = NROW + NROW_serial;
}

void RACE::Graph::getPerm(int **perm_, int *len_)
{
    (*perm_) = totalPerm;
    (*len_) = NROW + NROW_serial;
    totalPerm = NULL;
}

void RACE::Graph::getInvPerm(int **invPerm_, int *len_)
{
    (*invPerm_) = totalInvPerm;
    (*len_) = NROW + NROW_serial;
    totalInvPerm = NULL;
}

int RACE::Graph::getChildrenSize(int row)
{
    return childrenSize[row];
}

int* RACE::Graph::getChildren(int row)
{
    return &(graphData[childrenStart[row]]);
}
