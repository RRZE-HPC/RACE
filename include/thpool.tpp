/* vim: set filetype=cpp: */
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


#include <functional>
#include <sys/time.h>
#include "timing.h"
#include "test.h"
#include<malloc.h>
#include <unistd.h>

#define DEFAULT_RACE_BLOCKCTR 200

#define DEBUG 0

template< typename arg_t>
thread<arg_t>::thread():pthread(NULL),workerTeam(NULL),jobPresent(false)
{
    pthread = new pthread_t;
}

    template<typename arg_t>
thread<arg_t>::~thread()
{
    delete pthread;
    delete signal;
}

/* interrupt has to be called before; coming to this*/
    template<typename arg_t>
void thread<arg_t>::kill()
{
    jobPresent=1;
    //send signal to tell the thread to wake and destroy
    sendSignal(signal);
}

    template<typename arg_t>
void thread<arg_t>::init (int gid_, thpool<arg_t>* pool_)
{
    gid = gid_;
    pool = pool_;

    signal = new Signal(pool->RACE_BLOCKCTR);

    if(gid != 0)
    {
        pthread_create(pthread, NULL, (void *(*)(void *)) run<arg_t>, this);
        pthread_detach(*pthread);
    }
    else
    {
        //master thread
        (*pthread) = pthread_self();
    }
}

    template<typename arg_t>
void* run(void* thread_void)
{
    thread<arg_t>* curThread = static_cast<thread<arg_t> *>(thread_void);
    thpool<arg_t>* pool = curThread->pool;
    pthread_mutex_lock(pool->thcount_lock);
    pool->num_threads_alive++;
    pthread_mutex_unlock(pool->thcount_lock);

    while(!(__sync_fetch_and_add(&(pool->interrupt),0)))
    {
        //wait till jobPresent becomes true
        waitSignal(curThread->signal, curThread->jobPresent, true);
        //do RECV if master of the team is in different proc

        if(curThread->jobPresent && (!pool->interrupt))
        {
#if DEBUG
            printf("tid = %u doing job\n",(unsigned) (* (curThread->pthread)));
#endif
            //do current job
            //START_TIME(work);
            curThread->myJob.function(curThread->myJob.arg);
            //STOP_TIME(work);
            //finish job
            curThread->finishJob();
        }
    }

    pthread_mutex_lock(pool->thcount_lock);
    pool->num_threads_alive--;
    pthread_mutex_unlock(pool->thcount_lock);

    return NULL;
}

    template<typename arg_t>
void thread<arg_t>::addJob(std::function<void(arg_t)> function_, arg_t arg_, team<arg_t>* currTeam_)
{
    bool master = ( pthread_self() == (*pthread) );
    //START_TIME(add_job);

    if(jobPresent && !master)
    {
        ERROR_PRINT("thread = %d(%u), Previous job incomplete", gid, (unsigned)pthread_self());
    }
    else
    {
        myJob.function = function_;
        myJob.arg = arg_;

        /* pthread_mutex_lock(signal->lock);
           jobPresent = true;
           pthread_mutex_unlock(signal->lock);
           */
        //casted to void to avoid warning/error "value computed is not used"
        (void) __sync_lock_test_and_set(&jobPresent, 1);
        /*all the slaves belong to the current team,
         * master belongs to the parent team, this
         * is because master is a worker of previous
         * team and here he is the master; so don't
         * rewrite his team */
        if(!master)
        {
            workerTeam = currTeam_;
            sendSignal(signal);
            //if my proc id not equal to curr child; then do SEND
        }

    }
    //STOP_TIME(add_job);
}

   template<typename arg_t>
void thread<arg_t>::finishJob()
{
    bool master = ( pthread_self() == (*pthread) );

    if(!jobPresent && !master)
    {
        ERROR_PRINT("Nothing to finish; job is already complete");
    }
    else
    {
        //jobPresent = false;
        __sync_lock_release(&jobPresent);
#if DEBUG
        printf("tid = %d(%u) finished job\n", gid, (unsigned) (* (pthread)));
#endif

        if(__sync_add_and_fetch(&(workerTeam->num_jobs), 1) == workerTeam->num_slaves)
        {
#if DEBUG
            printf("tid = %d(%u) sending barrier signal\n", gid, (unsigned) (* (pthread)));
#endif
            sendSignal(workerTeam->barrierSignal);
            //if master not in same proc; do a MPI_SEND
#if DEBUG
            printf("tid = %d(%u) sent barrier signal\n", gid, (unsigned) (* (pthread)));
#endif
        }
    }
}

    template<typename arg_t>
thpool<arg_t>::thpool():initialized(false)
{
    char* RACE_BLOCKCTR_str = getenv("RACE_BLOCKCTR");
    if(RACE_BLOCKCTR_str == NULL)
    {
        RACE_BLOCKCTR = DEFAULT_RACE_BLOCKCTR;
    }
    else {
        RACE_BLOCKCTR = atoi(RACE_BLOCKCTR_str);
    }
}

    template<typename arg_t>
void thpool<arg_t>::init(int numThreads_)
{
    if(numThreads_ > 0)
    {
        initialized = true;
        num_threads  = std::max(numThreads_,1);
        num_threads_alive = 0;
        interrupt = 0;

        thcount_lock = new pthread_mutex_t;
        pthread_mutex_init(thcount_lock, NULL);

        //allocating threads
        threads = new thread<arg_t>[num_threads];

        for(int tid=0; tid<num_threads; ++tid)
        {
            threads[tid].init(tid, this);
        }

        while(__sync_fetch_and_add(&num_threads_alive,0) != (num_threads-1)){/*wait till all spawned threads are initialized*/ }
    }

}


    template<typename arg_t>
thpool<arg_t>::~thpool()
{
    //clean if pool is initialized before
    if(initialized)
    {
        interrupt = 1;

        //give Signal
        for(int i=1; i<num_threads; ++i)
        {
            threads[i].kill();
        }

        while(__sync_fetch_and_add(&num_threads_alive,0) != 0) {/*wait till all are destroyed*/ }
        delete thcount_lock;
        delete[] threads;
    }
}

    template<typename arg_t>
team<arg_t>::team():initialized(false)
{
}

    template<typename arg_t>
void team<arg_t>::init(std::vector<int> tid, thpool<arg_t>* pool)
{
    num_threads = tid.size();
    if(num_threads > 0)
    {
        initialized = true;

        num_slaves = num_threads-1;
        num_jobs = 0;
        //Construct task force
        taskForce = new thread<arg_t>*[num_threads];

        for(int i=0; i<num_threads; ++i)
        {
            taskForce[i] = &(pool->threads[tid[i]]);
        }

        RACE_BLOCKCTR = pool->RACE_BLOCKCTR;

        barrierSignal = new Signal(RACE_BLOCKCTR);
    }
}

    template<typename arg_t>
team<arg_t>::~team()
{
    if(initialized)
    {
        delete[] taskForce;
        delete barrierSignal;
    }
}

    template<typename arg_t>
void team<arg_t>::addJob(int tid, std::function<void(arg_t)> function_, arg_t arg_)
{
    if(tid >= 0 && tid < num_threads)
    {
#if DEBUG
        printf("tid = %d(%u), job added\n",tid, (unsigned)*(taskForce[tid]->pthread));
#endif
        taskForce[tid]->addJob(function_, arg_, this);
    }
    else
    {
        ERROR_PRINT("No job added; specified thread gid does not exists");
    }
}

/* Master does his job and waits
 * until all the submitted jobs have finished */
    template<typename arg_t>
void team<arg_t>::barrier()
{
    //current master does his job
#if DEBUG
    printf("master: tid = %u doing job\n",(unsigned) (pthread_self()));
#endif

    taskForce[0]->myJob.function(taskForce[0]->myJob.arg);
    taskForce[0]->jobPresent = 0;//false;

#if DEBUG
    printf("master: tid = %u finished job\n",(unsigned) (pthread_self()));
#endif

    /*timeval start, end;
      double startTime, endTime;

      gettimeofday(&start, NULL);		
      startTime = start.tv_sec + start.tv_usec*1e-6;*/
    //START_TIME(Barrier);

    waitSignal(barrierSignal, num_jobs, num_slaves);
    //If we have children in different procs do MPI_RECV; when counting check
    //only for children that do not belong to me

    __sync_lock_test_and_set(&num_jobs, 0);

    /*gettimeofday(&end, NULL);		
      endTime = end.tv_sec + end.tv_usec*1e-6;
      barrierTime = endTime - startTime;*/

#if DEBUG
    printf("master: tid = %u barrier done\n",(unsigned) (pthread_self()));
#endif


    // STOP_TIME(Barrier);
}

/*Make all threads in team to switch from
 * active to idle waiting*/
    template<typename arg_t>
void team<arg_t>::sleep()
{
    for(int tid=1; tid<num_threads; ++tid)
    {
        //added with spin, since spin value is 0
        if(__sync_fetch_and_add(&((taskForce[tid]->signal->mode)),spin)!=idle)
        {
            //wait till it starts spinning;
            //i.e its job is finished
            while(__sync_fetch_and_add(&(taskForce[tid]->signal->mode),spin)==released)
            {
            }
            sleepSignal(taskForce[tid]->signal);
            //wait  till thread is idle
            while(__sync_fetch_and_add(&(taskForce[tid]->signal->mode),spin)!=idle)
            {
            }
        }

    }

    //reset BLOCKCTR value
    taskForce[0]->pool->RACE_BLOCKCTR = DEFAULT_RACE_BLOCKCTR;
//    __sync_synchronize();
}
