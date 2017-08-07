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

//TODO no need to put only diagonal elements in Graph
NAME_error Graph::createGraphFromCRS(int *rowPtr, int *col, int *initPerm, int *initInvPerm)
{
    if(NROW != NCOL) {
        ERROR_PRINT("NROW!=NCOL : Currently Graph BMC supports only undirected Graph");
        return NAME_ERR_MATRIX_SYMM;
    }

    int nodeWithChildren = 0;

    std::vector<bool> pure_diag_flag(NROW,false);

    int nnz = 0;
#pragma omp parallel for schedule(runtime) reduction(+:nodeWithChildren) reduction(+:nnz)
    for(int row=0; row<NROW; ++row) {
        graphData[row].upperNnz=0;
        int permRow = row;
        if(initPerm)
        {
            permRow = initPerm[row];
        }
        int permCol = 0;
        if( (rowPtr[permRow+1] - rowPtr[permRow]) > 1) {
            for(int idx=rowPtr[permRow]; idx<rowPtr[permRow+1]; ++idx) {
                permCol = col[idx];
                if(initInvPerm)
                {
                    permCol = initInvPerm[col[idx]];
                }
                graphData[row].children.push_back(permCol);
                if(permCol >= row)
                {
                    ++graphData[row].upperNnz;
                }
                ++nnz;
            }
            nodeWithChildren++;
        } else {
#pragma omp critical
            {
                pureDiag.push_back(row);
            }
        }
    }

    NNZ = nnz;
    //resize Graph to the small size
    if(nodeWithChildren != NROW) {
        ERROR_PRINT("Currently Graph has to be connected and irreducible");
        return NAME_ERR_NOT_IMPLEMENTED;
        int ctr = 0;
        for(intIter iter=pureDiag.begin(); iter!=pureDiag.end(); ++iter) {
            graphData.erase(graphData.begin()+(*iter)-ctr);
            ++ctr;
        }
    }

    if(!pureDiag.empty())
    {
        ERROR_PRINT("Currently Graph has to be connected and irreducible");
        return NAME_ERR_NOT_IMPLEMENTED;
    }
    return NAME_SUCCESS;
}

//only for debugging puposes
void Graph::writePattern(char* name)
{
    std::string base_name(name);
    std::stringstream fileName;
    fileName<<base_name <<".mtx";
    std::fstream file;
    file.open( fileName.str().c_str() ,std::fstream::out);
    file<<"%";
    file<<"%MatrixMarket matrix coordinate real general\n";
    file<<"%MatrixMarketm file generated by NAME\n";

    file<<NROW<<" "<<NCOL<<" "<<NNZ<<"\n";

    for(int row=0; row<NROW; ++row) {
        for (unsigned j=0; j<graphData[row].children.size(); ++j) {
            file<<row+1<<" "<<graphData[row].children[j]+1<<" "<<10<<"\n";
        }
    }
}


Graph::Graph(int nrow, int ncol, int *rowPtr, int *col, int *initPerm, int *initInvPerm):graphData(nrow),NROW(nrow),NCOL(ncol)
{
    NAME_FN(createGraphFromCRS(rowPtr, col, initPerm, initInvPerm));
}

Graph::Graph(const Graph& srcGraph):graphData(srcGraph.graphData),pureDiag(srcGraph.pureDiag),NROW(srcGraph.NROW),NCOL(srcGraph.NCOL)
{
}

NAME_error Graph::swap(Graph& other)
{
    if( (NROW != other.NROW) || (NCOL != other.NCOL)) {
        ERROR_PRINT("Incompatible Graphs");
        return NAME_ERR_INCOMPATIBILITY;
    }

    //swap graphData
    graphData.swap(other.graphData);
    pureDiag.swap(other.pureDiag);
    return NAME_SUCCESS;
}

Node& Graph::at(unsigned Idx)
{
    return graphData[Idx];
}


