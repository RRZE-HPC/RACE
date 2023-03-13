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

#include "interface.h"
#include "simdify.h"
#include <algorithm>
#include <iostream>
#include "utility.h"
#include "timing.h"
#include <algorithm>
#include "config.h"

RACE::Interface::Interface(int nrow_,int nthreads_, RACE::dist dist_, int *rowPtr_, int *col_, bool symm_hint, int SMT_, RACE::PinMethod pinMethod_, int *initPerm_, int *initInvPerm_, RACE::d2Method d2Type_, RACE::LBTarget lbTarget_):graph(NULL),nrow(nrow_),distance(dist_),d2Type(d2Type_),lbTarget(lbTarget_),requestedThreads(nthreads_),availableThreads(-1),SMT(SMT_),pinMethod(pinMethod_),pool(NULL),initPerm(initPerm_),initInvPerm(initInvPerm_),rowPtr(rowPtr_),col(col_),zoneTree(NULL),powerCalculator(NULL), highestPower(1), highestSubPower(1)
{
    graph = new Graph(nrow, nrow, rowPtr, col, distance, symm_hint, initPerm, initInvPerm);
}

RACE::Interface::~Interface()
{
    if(zoneTree) {
        delete zoneTree;
    }

    if(graph) {
        delete graph;
    }

    if(pool) {
        delete pool;
    }

    if(powerCalculator) {
        delete powerCalculator;
    }

    for(unsigned i=0; i<funMan.size(); ++i)
    {
        delete funMan[i];
    }
}

// TODO: needed when?
void RACE::Interface::passGlobalRows(int globalStartRow, int globalEndRow){
    // TODO: Affix these to Graph object 
    // graph->globalStartRow = globalStartRow;
    // graph->globalEndRow = globalEndRow;
}

// TODO: needed when?
void RACE::Interface::compressColIdx(){
    // template <typename VT, typename IT>
// void collect_local_needed_heri(
//     std::vector<std::vector<IT>> *communication_recv_idxs,
//     std::vector<IT> *recv_counts_cumsum,
//     ScsData<VT, IT> *local_scs,
//     const IT *work_sharing_arr,
//     int my_rank,
//     int comm_size)
// {
    // int from_proc, to_proc, elem_col;
    // int needed_heri_count = 0;
    // int amnt_lhs_halo_elems = 0;

    // // To remember which columns have already been accounted for
    // std::unordered_set<int> remote_elem_col_bk;

    // std::vector<int> remote_elem_idxs;
    // std::vector<int> original_col_idxs(local_scs->col_idxs.data(), local_scs->col_idxs.data() + local_scs->n_elements);

    // // TODO: these sizes are too large. Should be my_rank and comm_size - my_rank
    // int lhs_needed_heri_counts[comm_size] = {0};
    // int rhs_needed_heri_counts[comm_size] = {0};
    // // int lhs_needed_heri_counts[my_rank] = {0};
    // // int rhs_needed_heri_counts[comm_size - my_rank] = {0};

    // // COUNTING LOOP / PROC
    // for (IT i = 0; i < local_scs->n_elements; ++i)
    // {
    //     // If true, this is a remote element, and needs to be added to vector
    //     elem_col = local_scs->col_idxs[i];

    //     // TODO: still need to verify what is happening with padded elements
    //     // if this column corresponds to a padded element, continue to next nnz
    //     // if(elem_col == 0 && local_scs->values[i] == 0) {continue;}

    //     if (elem_col < graph->globalStartRow)
    //     {
    //         remote_elem_idxs.push_back(i);
    //         if(remote_elem_col_bk.find(elem_col) == remote_elem_col_bk.end()){
    //             // if this column has not yet been seen
    //             for (IT j = 0; j < comm_size; ++j) //TODO: change to only go until my_rank
    //             {
    //                 if (elem_col >= work_sharing_arr[j] && elem_col < work_sharing_arr[j + 1])
    //                 {
    //                     // Remember column corresponding to remote element
    //                     // remote_elem_col_bk.push_back(remote_elem_col);
    //                     remote_elem_col_bk.insert(elem_col);

    //                     ++lhs_needed_heri_counts[j]; // This array describes how many remote elements

    //                     break;
    //                 }
    //             }
    //             //if nothing found, error
    //         }
    //     }
    //     else if (elem_col > graph->globalEndRow - 1)
    //     { // i.e. if RHS remote element
    //         remote_elem_idxs.push_back(i);
    //         if(remote_elem_col_bk.find(elem_col) == remote_elem_col_bk.end()){
    //             for (IT j = 0; j < comm_size; ++j) // TODO: change to start at my_rank
    //             {
    //                 if (elem_col >= work_sharing_arr[j] && elem_col < work_sharing_arr[j + 1])
    //                 {
    //                     // Remember column corresponding to remote element
    //                     remote_elem_col_bk.insert(elem_col);

    //                     ++rhs_needed_heri_counts[j];
    //                     break;
    //                 }
    //             }

    //         }
    //     }
    //     else
    //     { // i.e. local element
    //         local_scs->col_idxs[i] -= work_sharing_arr[my_rank];
    //     }
    // }

    // IT local_elem_offset = work_sharing_arr[my_rank + 1] - graph->globalStartRow;

    // int lhs_heri_ctr[comm_size] = {0}; //TODO: change to only have size my_rank
    // int rhs_heri_ctr[comm_size] = {0}; //TODO: change to have size comm_size - my_rank

    // std::map<int, int> remote_cols;

    // // ASSIGNMENT LOOP
    // for (auto remote_elem_idx : remote_elem_idxs)
    // {
    //     elem_col = original_col_idxs[remote_elem_idx];

    //     // if this column corresponds to a padded element, continue to next nnz
    //     // if(elem_col == 0 && local_scs->values[i] == 0) {continue;}

    //     if (elem_col < graph->globalStartRow)
    //     { 
    //         for (IT j = 0; j < comm_size; ++j)//TODO: change to only go until my_rank
    //         {
    //             if (elem_col >= work_sharing_arr[j] && elem_col < work_sharing_arr[j + 1])
    //             {
    //                 // So, on the current process, I will know from which process this 
    //                 // new col_idx for lhs halo element
                    
    //                 if(remote_cols.find(elem_col) == remote_cols.end()){
    //                     remote_cols[elem_col] = local_elem_offset + lhs_cumsum_heri_counts[j] + lhs_heri_ctr[j];
    //                     ++lhs_heri_ctr[j];
    //                 }

    //                 local_scs->col_idxs[remote_elem_idx] = remote_cols[elem_col];
    //             }
    //         }

    //     }
    //     else if (elem_col > graph->globalEndRow - 1)
    //     { // i.e. if RHS remote element
    //         // The rank of where this needed element resides is deduced from the work sharing array.
    //         for (IT j = 0; j < comm_size; ++j)// TODO: change to start at my_rank
    //         {
    //             if (elem_col >= work_sharing_arr[j] && elem_col < work_sharing_arr[j + 1])
    //             {
    //                 if(remote_cols.find(elem_col) == remote_cols.end()){
    //                     remote_cols[elem_col] = local_elem_offset + lhs_cumsum_heri_counts[my_rank] + rhs_cumsum_heri_counts[j] + rhs_heri_ctr[j];
    //                     ++rhs_heri_ctr[j];
    //                 }
    //                 local_scs->col_idxs[remote_elem_idx] = remote_cols[elem_col];
    //             }
    //         }
    //     }
    // }

}

// For matrix power kernel
RACE_error RACE::Interface::RACEColor(int highestPower_, double cacheSize, double safetyFactor, std::string mtxType, int highestSubPower_)
{
    int numSharedCache = 1; //always 1, greater than 1 NUMA domain deprecated for OMP support
    if(numSharedCache != 1)
    {
        ERROR_PRINT("RACE currently does not support greater than one shared cache or NUMA domain. Please configure 'RACE_ENABLE_MPI_MPK' to 'ON' and use one MPI per NUMA domain configuration");
        return RACE_ERR_INVALID_ARG;
    }
    highestPower = highestPower_;
    highestSubPower = highestSubPower_;
    if(highestSubPower < 1)
    {
        ERROR_PRINT("Highest sub power is less than one. Expected a value greater than one");
        return RACE_ERR_INVALID_ARG;
    }
    if(distance != RACE::POWER)
    {
        ERROR_PRINT("If you need to calculate power set distance to RACE::POWER");
        return RACE_ERR_INVALID_ARG;
    }
    else
    {
        //std::vector<int> distFromRemotePtr;
        bool useMPI = true;
#ifndef RACE_ENABLE_MPI_MPK
        useMPI = false;
#endif
        // TODO: define this is config.h.in
        // #ifdef MPI_INCLUDED
        //     printf("MPI IS DEFINED\n");
        // #endif
        if(useMPI == true){
            // TODO: misnomer, this function does more than collect only mpi boundary nodes. It finds all the rings. Change name
            graph->collectBoundaryNodes(highestPower*highestSubPower);
        }
        powerCalculator = new mtxPowerRecursive(graph, highestPower*highestSubPower, numSharedCache, cacheSize, safetyFactor, mtxType);
        //sanity check
        if(requestedThreads%numSharedCache)
        {
            ERROR_PRINT("Threads (=%d) not a multiple of requested nodes (=%d)\n", requestedThreads, numSharedCache);
            exit(-1);
        }
        powerCalculator->findPartition();
        int len;
        powerCalculator->getPerm(&perm, &len);
        powerCalculator->getInvPerm(&invPerm, &len);
        //Pin pin(NULL, 1, RACE::FILL);
        //pin.pinPowerThread(numSharedCache);

        return RACE_SUCCESS;

    }
}


RACE_error RACE::Interface::RACEColor()
{

    if(distance == RACE::POWER)
    {
        ERROR_PRINT("If you need to color the matrix specify distance other that RACE::POWER");
        return RACE_ERR_INVALID_ARG;
    }
    else
    {
        //1. Construct Graph
        /*if(graph == NULL)
          {
          graph = new Graph(nrow, nrow, rowPtr, col, initPerm, initInvPerm);
          }*/
        //2. Call level_recursion
        LevelRecursion lr(graph, requestedThreads, distance, d2Type, lbTarget);
        lr.levelBalancing();
        availableThreads = lr.getAvailableThreads();
        int len;
        lr.getPerm(&perm, &len);
        lr.getInvPerm(&invPerm, &len);
        zoneTree = lr.getZoneTree();

        pool = new LevelPool(zoneTree, SMT, pinMethod);

        printZoneTree();

#ifdef RACE_KERNEL_THREAD_OMP
        return pool->pin.pinApplication();
#else
        return pool->createPool();//creates pinned thread pools
#endif
        /*    printf("Checking Coloring\n");
              if(D2Checker())
              {
              ERROR_PRINT("Conflict in coloring\n");
              }
              printf("Checking Finished\n");
              */
    }
}


double RACE::Interface::getEfficiency()
{
    if(distance == RACE::POWER){
        // MPK 
        // TODO: Expose/use all ring data for relevant statistics
        int powerMax = graph->distFromRemotePtr.size() - 1;
        int mainRows = graph->distFromRemotePtr[powerMax] - graph->distFromRemotePtr[powerMax - 1];
        // printf("%d - %d main rows\n", graph->distFromRemotePtr[powerMax], graph->distFromRemotePtr[powerMax - 1]);
        double eff = (double)mainRows / (double)nrow;
        return eff;
    }
    else{
        // Coloring
        double effThreads = ((double)nrow/(double)zoneTree->at(0).effRowZ);
        double eff = effThreads/requestedThreads;
        return eff;
    }
}

int RACE::Interface::getMaxStageDepth()
{
    return zoneTree->maxStages();
}

void RACE::Interface::printZoneTree()
{
   zoneTree->printTree();
}

void RACE::Interface::getPerm(int **perm_, int *len_)
{
#if 0 //now included
    if(initPerm)
    {
        int *totPerm = new int [nrow];
#pragma omp parallel for schedule(runtime)
        for(int i=0; i<nrow; ++i)
        {
            totPerm[i] = initPerm[perm[i]];
        }

        (*perm_) = totPerm;
        delete[] perm;
    }
    else
    {
        (*perm_) = perm;
    }
#else
    //initPerm included
    (*perm_) = perm;
#endif
    (*len_) = nrow;
}

void RACE::Interface::getInvPerm(int **invPerm_, int *len_)
{
#if 0 //now included
    if(initInvPerm)
    {
        int *totInvPerm = new int [nrow];
#pragma omp parallel for schedule(runtime)
        for(int i=0; i<nrow; ++i)
        {
            totInvPerm[i] = invPerm[initInvPerm[i]];
        }


        (*invPerm_) = totInvPerm;
        delete[] invPerm;
    }
    else
    {
        (*invPerm_) = invPerm;
    }
#else
    (*invPerm_) = invPerm;
#endif
    (*len_) = nrow;
}

int RACE::Interface::getNumThreads()
{
    return availableThreads;
}


int RACE::Interface::registerFunction(void (*f) (int,int,void *), void *args)
{
    std::function<void (int,int,void *)> f_func = f;
    return registerFunction(f_func, args);
}

int RACE::Interface::registerFunction(std::function<void (int,int,void *)> f, void *args)
{
    if(distance != RACE::POWER)
    {
        //    pool->pin.pinApplication();

        FuncManager* currFun = new FuncManager(f, args, zoneTree, pool, graph->serialPart);
        funMan.push_back(currFun);

        //     funMan = new FuncManager(f, args, zoneTree, pool);

        /* for(int i=0; i<10; ++i) {
           START_TIME(omp_fn);
           funMan->RunOMP();
           STOP_TIME(omp_fn);
           }*/
        //TODO just pin for the first time; but now OpenMP does not work with this
        //model
        //Pin threads
        //Pin pin(zoneTree, true);
        //pin.pinApplication();

        return (funMan.size()-1);
    }
    else
    {
        if(highestSubPower == 1)
        {
            ERROR_PRINT("To calcaulte power call back function prototype is void (*f) (int,int,int,int,void *)");
        }
        else
        {
            ERROR_PRINT("To calcaulte power call back function prototype is void (*f) (int,int,int,int,int,void *)");
        }
        return RACE_ERR_INVALID_ARG;
    }
}

int RACE::Interface::registerFunction(void (*f) (int,int,int,int,int,void *), void *args, int power, int subPower, int numaSplit)
{
    std::function<void (int,int,int,int,int,void *)> f_function = f;
    return registerFunction(f_function, args, power, subPower, numaSplit);
}

int RACE::Interface::registerFunction(std::function<void (int,int,int,int,int,void *)> f, void *args, int power, int subPower, int numaSplit)
{
    if(distance != RACE::POWER)
    {
        ERROR_PRINT("To parallel execution of dependent kernels callback  function prototype is void (*f) (int,int,void *)");
        return RACE_ERR_INVALID_ARG;
    }
    else
    {
        if(power*subPower > highestPower*highestSubPower)
        {
            ERROR_PRINT("Specified power*subPower (%d) greater than highestPower*highestSubPower passed (%d) during the pre-processing phase (call to RACEColor function). Performance might not be optimal and in some cases might cause segfaults.", power*subPower, highestPower*highestSubPower);
        }
        FuncManager* currFun = new FuncManager(f, args, power, subPower, powerCalculator, numaSplit);
        funMan.push_back(currFun);
        return (funMan.size()-1);
    }
}

void RACE::wrappedPowerFunc(std::function<void(int,int,int,int,void *)> f, int start, int end, int pow, int subPow, int numaLocal, void* args)
{
    f(start, end, pow, numaLocal, args);
    UNUSED(subPow);
}

int RACE::Interface::registerFunction(std::function<void (int,int,int,int,void *)> f, void *args, int power, int numaSplit)
{
    if(highestSubPower == 1)
    {
        using namespace std::placeholders;
        std::function<void (int,int,int,int,int,void *)> bindedFunc = std::bind(RACE::wrappedPowerFunc, f, _1, _2, _3, _4, _5, _6);
        return registerFunction(bindedFunc, args, power, 1, numaSplit);
    }
    else
    {
        ERROR_PRINT("Sub power value is higher than one. This means you have to supply a kernel that takes subPow into account. The funtion prototype is (*f) (int,int,int,int,void *)");
        return RACE_ERR_INCOMPATIBILITY;
    }
}

int RACE::Interface::registerFunction(void (*f) (int,int,int,int,void *), void *args, int power, int numaSplit)
{
    std::function<void (int,int,int,int,void *)> f_function = f;
    return registerFunction(f_function, args, power, numaSplit);
}

int RACE::Interface::attachCommunicationToFunction(int funcId, std::function<void (int, int, void *)> f, void* args)
{
    if(distance != RACE::POWER)
    {
        ERROR_PRINT("Communication currently supported only for RACE::POWER");
        return RACE_ERR_NOT_IMPLEMENTED;

    }
    else
    {
        funMan[funcId]->registerCommFunc(f, args);
        return RACE_SUCCESS;
    }
}

int RACE::Interface::attachCommunicationToFunction(int funcId, void (*f) (int, int, void*), void* args)
{
    std::function<void (int, int, void *)> comm_function = f;
    return attachCommunicationToFunction(funcId, comm_function, args);
}

void RACE::wrappedCommFunc(std::function<void (int, void *)> f, int pow, int subPow, void *args)
{
    f(pow, args);
    UNUSED(subPow);
}

int RACE::Interface::attachCommunicationToFunction(int funcId, std::function<void (int, void *)> f, void* args)
{
    if(highestSubPower == 1)
    {
        using namespace std::placeholders;
        std::function<void (int, int, void *)> bindedFunc = std::bind(RACE::wrappedCommFunc, f, _1, _2, _3);
        return attachCommunicationToFunction(funcId, bindedFunc, args);
    }
    else
    {
        ERROR_PRINT("Sub power value is higher than one. This means you have to supply a communication kernel that takes subPow into account. The funtion prototype is (*f) (int,int,void *)");
        return RACE_ERR_INCOMPATIBILITY;
    }
}

int RACE::Interface::attachCommunicationToFunction(int funcId, void (*f) (int, void*), void* args)
{
    std::function<void (int, void *)> comm_function = f;
    return attachCommunicationToFunction(funcId, comm_function, args);
}

inline void dummyComm(int power, void* args)
{
    UNUSED(power);
    UNUSED(args);
}

int RACE::Interface::executeFunction(int funcId, bool rev)
{
    RACE_error status = RACE_SUCCESS;
    if(distance == RACE::POWER)
    {
        //check if communication(if required)  is set properly
        bool haveMPI=false;
#ifdef RACE_ENABLE_MPI_MPK
        if(!graph->distFromRemotePtr.empty())
        {
            int totPower = highestPower*highestSubPower;
            int totRemoteElems = graph->distFromRemotePtr[totPower-1];

            if(totRemoteElems > 0)
            {
                haveMPI=true;
            }
        }
#endif
        if(haveMPI)
        {
            if(!funMan[funcId]->isCommRegistered())
            {
                attachCommunicationToFunction(funcId, dummyComm, NULL);
                WARNING_PRINT("You haven't registered any communication function to function with id=%d. Although it seems you have MPI partitioning. We are attaching a dummy communication function", funcId);
            }
        }
    }
    funMan[funcId]->Run(rev);
    //  funMan->Run();
    return status;
}

void RACE::Interface::setPower(int funcId, int pow)
{
    funMan[funcId]->setPower(pow);
}

void RACE::Interface::setSubPower(int funcId, int subPow)
{
    funMan[funcId]->setSubPower(subPow);
}

int RACE::Interface::getPower(int funcId)
{
    return funMan[funcId]->getPower();
}

int RACE::Interface::getSubPower(int funcId)
{
    return funMan[funcId]->getSubPower();
}

void RACE::Interface::setSerial(int funcId)
{
    funMan[funcId]->setSerial();
}

void RACE::Interface::unsetSerial(int funcId)
{
    funMan[funcId]->unsetSerial();
}


int RACE::Interface::tuneFunction(int funcId, bool rev)
{
    int bestPow;
    if(tunedPowMap.find(funcId) == tunedPowMap.end())
    {
        std::vector<double> time_vec(highestPower, 0);
        //int origPower = funMan[funcId]->getPower();
        for(int pow=1; pow<=highestPower; ++pow)
        {
            setPower(funcId, pow);

            timeval start, end;
            double start_tym, end_tym;
            gettimeofday(&start, NULL);
            start_tym = start.tv_sec + start.tv_usec*1e-6;
            double time = 0;
            int iter=0;
            while(time < 0.2)
            {
                funMan[funcId]->Run(rev);

                gettimeofday(&end, NULL);
                end_tym = end.tv_sec + end.tv_usec*1e-6;

                time = end_tym-start_tym;
                ++iter;
            }
            time_vec[pow-1] = time/static_cast<double>(iter*pow);
        }
        int minIndex = std::distance(time_vec.begin(), std::min_element(time_vec.begin(), time_vec.end()));

        bestPow = minIndex + 1;
        tunedPowMap[funcId] = bestPow;
        INFO_PRINT("Best tuned power for function with id %d is %d", funcId, bestPow);
    }
    else
    {
        bestPow = tunedPowMap[funcId];
    }

    setPower(funcId, bestPow);
    return bestPow;
}

void RACE::Interface::resetTime()
{
    zoneTree->resetTime();
}

//TODO: move to graph
bool RACE::Interface::detectConflict(std::vector<int> range1, std::vector<int> range2)
{
#ifdef RACE_PERMUTE_ON_FLY
    WARNING_PRINT("Detect conflict experimental with RACE_PERMUTE_ON_FLY");
#endif
    bool conflict = false;
    for(int i=range1[0]; i<range1[1]; ++i)
    {
        int perm_i = i;
#ifdef RACE_PERMUTE_ON_FLY
        perm_i = graph->totalPerm[i];
#endif
#ifdef RACE_USE_SOA_GRAPH
        //TODO:optimise
        int* children_data = graph->getChildren(perm_i);
        int childrenSize = graph->getChildrenSize(perm_i);
        std::vector<int> children_vec_i(childrenSize);
        for(int c=0; c<childrenSize; ++c)
        {
            int permCol = children_data[c];
#ifdef RACE_PERMUTE_ON_FLY
            permCol = graph->totalInvPerm[permCol];
#endif
            children_vec_i[c] = permCol;
        }
        std::vector<int>* children_i = &children_vec_i;
#else
#ifdef RACE_PERMUTE_ON_FLY
        std::vector<int> children_data = graph->at(perm_i).children;
        int childrenSize = children_data.size();
        std::vector<int> children_vec_i(childrenSize);
        for(int c=0; c<childrenSize; ++c)
        {
            int permCol = children_data[c];
            permCol = graph->totalInvPerm[permCol];
            children_vec_i[c] = permCol;
        }
        std::vector<int>* children_i = &children_vec_i;
#else
        std::vector<int>* children_i = &(graph->at(perm_i).children);
#endif

#endif
        for(int j=range2[0]; j<range2[1]; ++j)
        {
            int perm_j =  j;
#ifdef RACE_PERMUTE_ON_FLY
            perm_j = graph->totalPerm[j];
#endif
#ifdef RACE_USE_SOA_GRAPH
            //TODO:optimise
            children_data = graph->getChildren(perm_j);
            childrenSize = graph->getChildrenSize(perm_j);
            std::vector<int> children_vec_j(childrenSize);
            for(int c=0; c<childrenSize; ++c)
            {
                int permCol = children_data[c];
#ifdef RACE_PERMUTE_ON_FLY
                permCol = graph->totalInvPerm[permCol];
#endif
                children_vec_j[c] = permCol;
            }
            std::vector<int>* children_j = &children_vec_j;
#else
#ifdef RACE_PERMUTE_ON_FLY
        std::vector<int> children_data = graph->at(perm_j).children;
        int childrenSize = children_data.size();
        std::vector<int> children_vec_j(childrenSize);
        for(int c=0; c<childrenSize; ++c)
        {
            int permCol = children_data[c];
            permCol = graph->totalInvPerm[permCol];
            children_vec_j[c] = permCol;
        }
        std::vector<int>* children_j = &children_vec_j;
#else
        std::vector<int>* children_j = &(graph->at(perm_j).children);
#endif
#endif
            conflict = (std::find_first_of(children_i->begin(), children_i->end(), children_j->begin(), children_j->end()) != children_i->end());
            if(conflict)
            {
                ERROR_PRINT("Conflict at row: %d %d\nCheck i: [",i,j);
                for(unsigned child_i=0; child_i<children_i->size(); ++child_i)
                {
                    printf("%d, ",children_i->at(child_i));
                }
                printf("]\nCheck j: [");
                for(unsigned child_j=0; child_j<children_j->size(); ++child_j)
                {
                    printf("%d, ",children_j->at(child_j));
                }
                printf("]\n");

                auto result = std::find_first_of(children_i->begin(), children_i->end(), children_j->begin(), children_j->end());
                std::cout<<"conflicted values = " << std::distance(children_i->begin(),result)<< "and "<< std::distance(children_j->begin(),result)<<std::endl;

                break;
            }
        }
        if(conflict)
        {
            break;
        }
    }

    return conflict;
}

bool RACE::Interface::recursiveChecker(int parent)
{
    std::vector<int> *children = &(zoneTree->at(parent).children);
    int totalSubBlocks = zoneTree->at(parent).totalSubBlocks;
    int blockPerThread = getBlockPerThread(distance, d2Type);
    bool conflict = false;

    for(int subBlock=0; subBlock<totalSubBlocks; ++subBlock)
    {
        int currNThreads = (children->at(2*subBlock+1)-children->at(2*subBlock))/blockPerThread;
        printf("checking idx = %d\n",parent);
        if(currNThreads > 1)
        {
            for(int start=0; start<blockPerThread; ++start)
            {
                for(int i=start; i<currNThreads; i+=blockPerThread)
                {
                    std::vector<int> rangeOuter = zoneTree->at(children->at(2*subBlock)+i).valueZ;
                    for(int j=start; j<=currNThreads; j+=2)
                    {
                        if(i!=j)
                        {
                            std::vector<int> rangeInner = zoneTree->at(children->at(2*subBlock)+j).valueZ;
                            conflict = detectConflict(rangeOuter, rangeInner);
                            if(conflict)
                            {
                                break;
                            }
                        }
                    }
                    if(conflict)
                    {
                        break;
                    }
                }
            }
            for(int child=0; child< blockPerThread*currNThreads; ++child)
            {
                recursiveChecker(children->at(2*subBlock)+child);
            }
        }
    }
    return conflict;
}

bool RACE::Interface::D2Checker()
{
    bool conflict = false;
    conflict = recursiveChecker(5);

    return conflict;
}

void RACE::Interface::sleep()
{
#ifndef RACE_KERNEL_THREAD_OMP
    pool->sleepPool();
#endif
}

bool RACE::Interface::simdify(int simdWidth, int C, int nrows, int* col_new, int* chunkStart, int* rl, int* clp, double* val, bool diagFirst)
{
    return simdifyTemplate<double> (simdWidth, C, nrows, col_new, chunkStart, rl, clp, val, this, diagFirst);
}

bool RACE::Interface::simdify(int simdWidth, int C, int nrows, int* col_new, int* chunkStart, int* rl, int* clp, float* val, bool diagFirst)
{
    return simdifyTemplate<float> (simdWidth, C, nrows, col_new, chunkStart, rl, clp, val, this, diagFirst);
}

bool RACE::Interface::simdifyD1(int simdWidth, int C, int nrows, int* col_new, int* chunkStart, int* rl, int* cl, double* val, bool d2Compatible)
{
    return simdifyD1Template<double> (simdWidth, C, nrows, col_new, chunkStart, rl, cl, val, d2Compatible);
}

bool RACE::Interface::simdifyD1(int simdWidth, int C, int nrows, int* col_new, int* chunkStart, int* rl, int* cl, float* val, bool d2Compatible)
{
    return simdifyD1Template<float> (simdWidth, C, nrows, col_new, chunkStart, rl, cl, val, d2Compatible);
}


void RACE::Interface::pinThread(int threadId)
{
    pool->pin.pinThread(threadId);
}

void RACE::powerInitRowPtrFunc(int start, int end, int pow, int subPow, int numa_domain, void* arg)
{
     DECODE_ARG(arg);
     if((col != NULL && val != NULL) && x != NULL)
     {
         printf("Something went wrong, I shouldnn' t be here\n");
     }

     /*
     if( lock_memory( (char*)&(rowPtr[start]), (end-start)*sizeof(int) ) == -1 )
     {
         rowPtr[0] = 0*pow*nrow;
         ERROR_PRINT("Memory locking didn't work");
     }*/

//#pragma omp parallel for schedule(static)

     if((pow == 1) && (subPow == 1))
     {
         for(int row=start; row<end; ++row)
         {
             rowPtr[row] = 0*pow*nrow;
             rowPtr[row+1] = 0*pow*nrow;
         }
     }
     UNUSED(numa_domain);
}

void RACE::Interface::numaInitRowPtr(int *rowPtr_)
{
    if(distance == RACE::POWER)
    {
        ENCODE_ARG(nrow, rowPtr_, NULL, NULL, NULL);
        int fn_id = registerFunction(RACE::powerInitRowPtrFunc, voidArg, 1);
        executeFunction(fn_id);
    }
    else
    {
        ERROR_PRINT("NUMA not implemented for coloring")
    }
}


void RACE::powerInitRowPtrNumaLocalFunc(int start, int end, int pow, int subPow, int numa_domain, void* arg)
{
     DECODE_NUMA_LOCAL_ARG(arg);
     if((col != NULL && val != NULL) && x != NULL)
     {
         printf("Something went wrong, I shouldnn' t be here\n");
     }

     /*
     if( lock_memory( (char*)&(rowPtr[start]), (end-start)*sizeof(int) ) == -1 )
     {
         rowPtr[0] = 0*pow*nrow;
         ERROR_PRINT("Memory locking didn't work");
     }*/

//#pragma omp parallel for schedule(static)

     if((pow == 1) && (subPow == 1))
     {
         for(int row=start; row<end; ++row)
         {
             rowPtr[numa_domain][row] = 0*pow*nrow;
             //rowPtr[numa_domain][row+1] = 0*pow*nrow;
         }
     }
}


void RACE::Interface::numaInitRowPtr(int **rowPtr_)
{
    if(distance == RACE::POWER)
    {
        ENCODE_NUMA_LOCAL_ARG(nrow, rowPtr_, NULL, NULL, NULL);
        int fn_id = registerFunction(RACE::powerInitRowPtrNumaLocalFunc, voidArg, 1, true);
        executeFunction(fn_id);
    }
    else
    {
        ERROR_PRINT("NUMA not implemented for coloring")
    }
}

//without NUMA split to local parts
void RACE::powerInitMtxVecFunc(int start, int end, int pow, int subPow, int numa_domain, void* arg)
{
    DECODE_ARG(arg);

/*
    if(pow == 0)
    {
        if( lock_memory( (char*)&(val[rowPtr[start]]), (rowPtr[end]-rowPtr[start])*sizeof(double) ) == -1)
        {
            ERROR_PRINT("Memory locking didn't work");
        }

        if( lock_memory( (char*)&(col[rowPtr[start]]), (rowPtr[end]-rowPtr[start])*sizeof(int) ) == -1)
        {
            ERROR_PRINT("Memory locking didn't work");
        }
    }
    if(x != NULL)
    {
        if( lock_memory( (char*)&(x[pow*nrow+start]), (end-start)*sizeof(double) ) == -1)
        {
            ERROR_PRINT("Memory locking didn't work");
        }
    }
    */

//#pragma omp parallel for schedule(static)
    for(int row=start; row<end; ++row)
    {
        if((pow==1) && (subPow==1))
        {
            for(int idx=rowPtr[row]; idx<rowPtr[row+1]; ++idx)
            {
                val[idx] = 0;
                col[idx] = 0;
            }
        }
        if(x!=NULL)
        {
            x[(pow)*nrow+row] = 0;
            x[(pow-1)*nrow+row] = 0;
        }
    }

    UNUSED(numa_domain);
}

void RACE::Interface::numaInitMtxVec(int *rowPtr_, int *col_, double *val_, double *x_, int power_)
{
    if(distance == RACE::POWER)
    {
        ENCODE_ARG(nrow, rowPtr_, col_, val_, x_);
        int fn_id = registerFunction(RACE::powerInitMtxVecFunc, voidArg, 1);
        executeFunction(fn_id);
        //funMan[fn_id]->NUMAInitPower();

    }
    else
    {
        ERROR_PRINT("NUMA not implemented for coloring")
    }
    UNUSED(power_);
}


//NUMA local sparse matrix
void RACE::powerInitMtxVecNumaLocalFunc(int start, int end, int pow, int subPow, int numa_domain, void* arg)
{
    DECODE_NUMA_LOCAL_ARG(arg);

//#pragma omp parallel for schedule(static)
    for(int row=start; row<end; ++row)
    {
        if((pow==1) && (subPow==1))
        {
            for(int idx=rowPtr[numa_domain][row]; idx<rowPtr[numa_domain][row+1]; ++idx)
            {
                val[numa_domain][idx] = 0;
                col[numa_domain][idx] = 0;
            }
        }
    }
    UNUSED(nrow);
    UNUSED(x);
}


void RACE::Interface::numaInitMtxVec(int **rowPtr_, int **col_, double **val_, int power_)
{
    if(distance == RACE::POWER)
    {
        ENCODE_NUMA_LOCAL_ARG(nrow, rowPtr_, col_, val_, NULL);
        int fn_id = registerFunction(RACE::powerInitMtxVecNumaLocalFunc, voidArg, 1, true);
        executeFunction(fn_id);
        //funMan[fn_id]->NUMAInitPower();

    }
    else
    {
        ERROR_PRINT("NUMA not implemented for coloring")
    }
    UNUSED(power_);
}



void RACE::Interface::getNumaSplitting(int **split, int *splitLen)
{
    std::vector<int> levelGroupPtr = powerCalculator->tree[0].nodePtr;
    int NUMAdomains = levelGroupPtr.size()-1;
    std::vector<int> levelPtr = powerCalculator->tree[0].lp;
    (*split) = (int*) malloc(sizeof(int)*(NUMAdomains+1));
    for(int i=0; i<=NUMAdomains; ++i)
    {
        (*split)[i] = levelPtr[levelGroupPtr[i]];
    }
    (*splitLen) = NUMAdomains+1;
}

int RACE::Interface::getHighestPower()
{
    return highestPower;
}

inline void MAT_NUM_VEC_ACCESSES(int start, int end, int pow, int numa_domain, void* args)
{
    DECODE_ARG(args);

    for(int row=start; row<end; ++row)
    {
        x[row]++;
        if(x[row] != pow)
        {
            if(x[row] > pow)
            {
                ERROR_PRINT("Oh oh we have duplicate computations, error at pow=%d, for row=%d. Value I got at x is %f, expected %d. Level start =%d, Level end=%d", pow, row, x[row], pow, start, end);
            }
            else
            {
                ERROR_PRINT("Oh oh have some missing computations, error at pow=%d, for row=%d. Value I got at x is %f, expected %d. Level start =%d, Level end=%d", pow, row, x[row], pow, start, end);
            }
        }
    }
}

void RACE::Interface::checkPowerCoverage()
{
    double* x = new double[nrow];
    for(int row=0; row<nrow; ++row)
    {
        x[row] = 0;
    }
    ENCODE_ARG(nrow, NULL, NULL, NULL, x);
    int race_power_id = registerFunction(&MAT_NUM_VEC_ACCESSES, voidArg, getHighestPower());
    {
        executeFunction(race_power_id);
    }
    DELETE_ARG();
    delete[] x;
}
