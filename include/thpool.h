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

#ifndef RACE_THPOOL_H
#define RACE_THPOOL_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <errno.h>
#include <time.h>
#if defined(__linux__)
#include <sys/prctl.h>
#endif
#include "error.h"
#include "utility.h"
#include <functional>
#include <vector>

#include "signal.h"

/* Job */
template<typename arg_t> struct job{
    std::function<void(arg_t)> function;
    arg_t  arg;     /* function's argument       */
};

template<typename arg_t>
struct thpool;

template<typename arg_t>
struct team;


/* Thread */
template<typename arg_t>
struct thread{
    int gid;                         /* friendly id               */
    pthread_t *pthread;              /* pointer to actual thread  */
    job<arg_t> myJob;
    thpool<arg_t>* pool;
    team<arg_t>* workerTeam;      /*Current team in which the thread works*/
    Signal* signal;                 /* Signal for job wake      */

    volatile int  jobPresent;
    int RACE_BLOCKCTR;
    //constructor
    thread();
    ~thread();
    void kill();
    void init(int id, thpool<arg_t>* pool_);
    void addJob(std::function<void(arg_t)>, arg_t arg_, team<arg_t> *currTeam_);
    void finishJob();
};

template<typename arg_t>
void* run(void* thread);

/* Threadpool */
template<typename arg_t>
struct thpool{
    thread<arg_t>*   threads;            /*threads in the pool        */
    int num_threads;                     /* threads currently alive   */
    int num_slaves;                      /* number of slave threads   */
    volatile int num_threads_alive;      /* threads currently working */
    pthread_mutex_t*  thcount_lock;      /* used for thread count etc */
    volatile int interrupt;
    int RACE_BLOCKCTR;

    thpool();
    void init(int numThreads_);
    ~thpool();
    //TODO sleep and wake functionality
    bool initialized;
};

/* This is the team which executes a given job*/
template<typename arg_t>
struct team{
    thread<arg_t>** taskForce;
    int num_threads;
    int num_slaves;
    int RACE_BLOCKCTR;
    bool initialized;
    volatile int num_jobs;

    Signal* barrierSignal;
    double barrierTime;
    team();
    ~team();
    //tid are global id referring which thrads belong to the team
    void init(std::vector<int> tid, thpool<arg_t>* pool);
    void addJob(int tid, std::function<void(arg_t)>, arg_t arg_);
    void barrier();
    void sleep();
};

#include "thpool.tpp"

#endif
