#ifndef RACE_GRAPH_H
#define RACE_GRAPH_H

#include "print.h"
#include "error.h"
#include <vector>

typedef std::vector<int>::iterator intIter;

struct Node{
    std::vector<int> children;
    //number of elements in upper part
    int upperNnz;
};


class Graph{
    private:
       /**
         * @brief Graph of the matrix.
         */
        std::vector<Node> graphData;
        /**
         * @brief Diagonal Elements of the matrix (may be removed in future).
         */
        std::vector<int> pureDiag;
        /**
         * @brief Serially procedssed row list
         */
        std::vector<int> serialPartRow;
        /**
         * @brief Create Graph from sparsematrix in CRS fromat.
         *
         * @param[in] rowPtr rowPtr of the matrix.
         * @param[in] col column index of the matrix.
         */
        RACE_error createGraphFromCRS(int *rowPtr, int *col, int *initPerm=NULL, int *initInvPerm=NULL);

        std::vector<int> perm;
        void permuteAndRemoveSerialPart();
    public:
        /**
         * @brief Number of Rows in the matrix.
         */
        int NROW;
        int NCOL;
        int NROW_serial;
        /**
         * @brief Number of Columns in the matrix.
         */
        int NNZ;
        int NNZ_serial;
        std::vector<int> serialPart;

        Graph(int nrow, int ncol, int *row_ptr, int *col, int *initPerm=NULL, int *initInvPerm=NULL);//constructor
        Graph(const Graph &srcGraph);//copy constructor
        void writePattern(char *name);
        //returns a list of nodes having load more than RACE_MAX_LOAD
        void getStatistics();
        void getInitialPerm(int **perm_, int *len_);
        RACE_error swap(Graph& other);
        Node& at(unsigned idx);
        /**
         * @brief Traverse class is responsible for traversing through graph
         * and building levelPtr
         */
        friend class Traverse;
};
#endif
