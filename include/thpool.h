#ifndef NAME_THPOOL_H
#define NAME_THPOOL_H

#include <unistd.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <errno.h>
#include <time.h>
#if defined(__linux__)
#include <sys/prctl.h>
#endif
#include "error.h"
#include <functional>
#include <vector>

#define PTHREAD_BARRIER

//#define FAST_BARRIER

#ifdef FAST_BARRIER
#include "fast_barrier.h"
#define pthread_barrier_t fast_barrier_t
#define pthread_barrier_init(B, A, N) (fast_barrier_init((B),(A), (N)))
#define pthread_barrier_destroy(B) (fast_barrier_destroy(B))
#define pthread_barrier_wait fast_barrier_wait
#endif


/* Job */
template<typename arg_t> struct job{
    std::function<void(arg_t)> function;
    arg_t  arg;     /* function's argument       */
};

template<typename arg_t>
struct thpool;

/* Thread */
template<typename arg_t>
struct thread{
    int id;                              /* friendly id               */
    pthread_t *pthread;                  /* pointer to actual thread  */
    job<arg_t> myJob;
    thpool<arg_t>* pool;
    pthread_mutex_t*  jobLock;         /* used for thread count etc */
    pthread_cond_t*   jobArrived;       /* signal to wake           */
    volatile bool jobPresent;
    //constructor
    thread();
    ~thread();
    void init(int id, thpool<arg_t>* pool_);
    void addJob(std::function<void(arg_t)>, arg_t arg_);
    void finishJob();
};

template<typename arg_t>
void run(thread<arg_t>* thread);


/* Threadpool */
template<typename arg_t>
struct thpool{
    thread<arg_t>*   threads;            /* pointer to threads        */
    int num_threads;                     /* threads currently alive   */
    int num_slaves;                      /* number of slave threads   */
    volatile int num_threads_alive;      /* threads currently working */
    volatile int num_jobs;               /* jobs  currently done      */
    pthread_mutex_t*  thcount_lock;      /* used for thread count etc */
    pthread_cond_t*  threads_all_idle;   /* signal to wait            */
    pthread_mutex_t*  jobLock;        /* used for thread count etc */
    pthread_cond_t*   jobArrived;        /* signal to wake           */
    pthread_barrier_t* pool_barrier;
    volatile bool interrupt;
    void init(int numThreads_);
    bool initialized;

    thpool();
    ~thpool();
    void addJob(int tid, std::function<void(arg_t)>, arg_t arg_);
    void doJob();
    void barrier();
};


#include "thpool.tpp"

#endif
