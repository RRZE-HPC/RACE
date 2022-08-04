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

#ifndef RACE_GRAPH_AOS_H
#define RACE_GRAPH_AOS_H

#include "print.h"
#include "error.h"
#include <vector>

typedef std::vector<int>::iterator intIter;

struct Node{
    std::vector<int> children;
    //number of elements in upper part
    int upperNnz;
};


/**
 * @brief RACE namespace.
 */
 namespace RACE
{
    class Graph;
}

class RACE::Graph{
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

        std::vector<int> serialPerm;
        void permuteAndRemoveSerialPart();
    public:
        int* totalPerm;
        int* totalInvPerm;
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
        ~Graph();
        void writePattern(char *name);
        //returns a list of nodes having load more than RACE_MAX_LOAD
        bool getStatistics();
        //returns only permutation for serial part
        void getSerialPerm(int **perm_, int *len_);
        //if used control will be passed to calling class
        void getPerm(int **perm_, int *len);
        void getInvPerm(int **invPerm_, int *len);

        RACE_error swap(Graph& other);
        Node& at(unsigned idx);
        /**
         * @brief Traverse class is responsible for traversing through graph
         * and building levelPtr
         */
        friend class Traverse;
};
#endif
