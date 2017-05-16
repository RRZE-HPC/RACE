/* vim: set filetype=cpp: */

#include <functional>
#include <sys/time.h>
#include "timing.h"
#include "test.h"
#include<malloc.h>

#define DEBUG 0
//#define SIGNAL_JOB

    template< typename arg_t>
thread<arg_t>::thread():pthread(NULL),jobPresent(false)//,myJob(NULL)
{
    pthread = new pthread_t;
}

    template<typename arg_t>
thread<arg_t>::~thread()
{
    delete pthread;
    /*if(myJob)
    {
        delete myJob;
    }*/

    delete jobLock;
    delete jobArrived;
}

    template<typename arg_t>
void thread<arg_t>::init (int id_, thpool<arg_t>* pool_)
{
    id = id_;
    pool = pool_;

    jobLock = new pthread_mutex_t;
    pthread_mutex_init(jobLock, NULL);
    jobArrived = new pthread_cond_t;
    pthread_cond_init(jobArrived, NULL);

    if(id != 0)
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
void run(thread<arg_t> * thread)
{
    thpool<arg_t>* pool = thread->pool;
    pthread_mutex_lock(pool->thcount_lock);
    pool->num_threads_alive++;
    pthread_mutex_unlock(pool->thcount_lock);

    while(!(pool->interrupt))
    {
        /*timeval signal_tym;
          gettimeofday(&signal_tym, NULL);
          printf("tid = %d(%u) waits at %f\n",thread->id, (unsigned) pthread_self(), signal_tym.tv_sec + signal_tym.tv_usec*1e-6);
          */
         /*gettimeofday(&signal_tym, NULL);
           printf("tid = %d(%u) recieves signal at %f\n",thread->id, (unsigned) pthread_self(), signal_tym.tv_sec + signal_tym.tv_usec*1e-6);
         */

#ifdef SIGNAL_JOB
        pthread_mutex_lock(thread->jobLock);
        while((!thread->jobPresent) && (!pool->interrupt))
        {
            //wait till job arrives
#if DEBUG
            printf("tid = %d(%u) arrives to wait\n",thread->id, (unsigned) pthread_self());
#endif
            //   timeval signal_tym;
            //   gettimeofday(&signal_tym, NULL);
            //    printf("tid = %d(%u) waits at %f\n",thread->id, (unsigned) pthread_self(), signal_tym.tv_sec + signal_tym.tv_usec*1e-6);

            pthread_cond_wait(thread->jobArrived, thread->jobLock);
#if DEBUG
            gettimeofday(&signal_tym, NULL);
            printf("tid = %d(%u) recieves signal at %f\n",thread->id, (unsigned) pthread_self(), signal_tym.tv_sec + signal_tym.tv_usec*1e-6);
#endif
        }
        pthread_mutex_unlock(thread->jobLock);
#else
        pthread_barrier_wait(pool->pool_barrier);
#endif

        if(thread->jobPresent)
        {
#if DEBUG
            printf("tid = %u doing job\n",(unsigned) (* (thread->pthread)));
#endif
            //do current job
            //START_TIME(work);
            thread->myJob.function(thread->myJob.arg);
            //STOP_TIME(work);
            //finish job
            thread->finishJob();
#if DEBUG
            printf("tid = %u finished job\n",(unsigned) (* (thread->pthread)));
#endif


#ifdef PTHREAD_BARRIER
            pthread_barrier_wait(pool->pool_barrier);
#else
            pthread_mutex_lock(pool->thcount_lock);
            pool->num_jobs++;
            if(pool->num_jobs == pool->num_slaves) {
                pthread_cond_signal(pool->threads_all_idle);
            }
            pthread_mutex_unlock(pool->thcount_lock);
#endif
        }
    }

    pthread_mutex_lock(pool->thcount_lock);
    pool->num_threads_alive--;
    pthread_mutex_unlock(pool->thcount_lock);
}

    template<typename arg_t>
void thread<arg_t>::addJob(std::function<void(arg_t)> function_, arg_t arg_)
{
    //START_TIME(add_job);
   /* if(myJob)
    {
        ERROR_PRINT("thread = %u, Previous job incomplete", (unsigned)pthread_self());
    }
    else*/
    {
        //        pthread_mutex_lock(pool->jobLock);
    //    myJob = new job<arg_t>;
       // myJob = (job<arg_t>*)memalign(128, sizeof(job<arg_t>));
        myJob.function = function_;
        myJob.arg = arg_;
        jobPresent = true;

        //master thread does not require this signal
        /*  if(id > 0)
            {
        //START_TIME(signal_job_arrival);
        //signal arrival of Job
#if DEBUG
printf("tid = %u sending signal to %u\n",(unsigned) pthread_self(), (unsigned) (*pthread));
#endif
pthread_cond_signal(jobArrived);
#if DEBUG
printf("tid = %u sent signal to %u\n",(unsigned) pthread_self(), (unsigned) (*pthread));
#endif
        //STOP_TIME(signal_job_arrival)
        }*/

        //START_TIME(mutex_unlock);
        //        pthread_mutex_unlock(pool->jobLock);
        //STOP_TIME(mutex_unlock);
#if DEBUG
        printf("tid = %u mutex unlocked\n",(unsigned) pthread_self());
#endif

        //STOP_TIME(add_job)
    }
}

    template<typename arg_t>
void thread<arg_t>::finishJob()
{
    /*if(myJob==NULL)
    {
        ERROR_PRINT("Nothing to finish; job is already complete");
    }
    else*/
    {
        //pthread_mutex_init(jobMutex, NULL);
        //pthread_cond_init(jobArrived, NULL);
     //   delete myJob;
       // free(myJob);
       // myJob = NULL;
         jobPresent = false;
    }
}

    template<typename arg_t>
thpool<arg_t>::thpool():initialized(false)
{
}

    template<typename arg_t>
void thpool<arg_t>::init(int numThreads_)
{
    if(numThreads_ > 0)
    {
        initialized = true;
        num_threads  = std::max(numThreads_,1);
        num_threads_alive = 0;
        num_jobs = 0;
        interrupt = false;
        num_slaves = num_threads-1; //except master

        thcount_lock = new pthread_mutex_t;
        threads_all_idle = new pthread_cond_t;
        jobArrived = new pthread_cond_t;
        jobLock = new pthread_mutex_t;

        pool_barrier = new pthread_barrier_t;
        pthread_barrier_init(pool_barrier, NULL, num_threads);
        pthread_mutex_init(thcount_lock, NULL);
        pthread_cond_init(threads_all_idle, NULL);
        pthread_mutex_init(jobLock, NULL);
        pthread_cond_init(jobArrived, NULL);
        //allocating threads
        threads = new thread<arg_t>[num_threads];

        for(int tid=0; tid<num_threads; ++tid)
        {
            threads[tid].init(tid, this);
        }

        while(num_threads_alive != (num_slaves)){/*wait till all spawned threads are initialized*/ }
    }
}

    template<typename arg_t>
thpool<arg_t>::~thpool()
{
    //clean if pool is initialized before
    if(initialized)
    {
        interrupt = true;

#ifdef SIGNAL_JOB
        /*pthread_mutex_lock(jobLock);
          pthread_cond_broadcast(jobArrived);
          pthread_mutex_unlock(jobLock);*/

        for(int i=1; i<num_threads; ++i)
        {
            pthread_mutex_lock(threads[i].jobLock);
            pthread_cond_signal(threads[i].jobArrived);
            pthread_mutex_unlock(threads[i].jobLock);
        }
#else
        pthread_barrier_wait(pool_barrier);
#endif

        while(num_threads_alive != 0) { /*wait till all are destroyed*/ }
        delete thcount_lock;
        delete threads_all_idle;
        delete jobLock;
        delete jobArrived;
        pthread_barrier_destroy(pool_barrier);
        delete pool_barrier;
        delete[] threads;
    }
}

    template<typename arg_t>
void thpool<arg_t>::addJob(int tid, std::function<void(arg_t)> function_, arg_t arg_)
{
    if(tid >= 0 && tid < num_threads)
    {
        threads[tid].addJob(function_, arg_);
    }
    else
    {
        ERROR_PRINT("No job added; specified thread id does not exists");
    }
}

    template<typename arg_t>
void thpool<arg_t>::doJob()
{
#ifdef SIGNAL_JOB
    for(int i=1; i<num_threads; ++i)
    {
        pthread_mutex_lock(threads[i].jobLock);
        pthread_cond_signal(threads[i].jobArrived);
        pthread_mutex_unlock(threads[i].jobLock);
    }
    /* timeval signal_tym;
       gettimeofday(&signal_tym, NULL);
       printf("tid = %d(%u) waits at %f\n",0, (unsigned) pthread_self(), signal_tym.tv_sec + signal_tym.tv_usec*1e-6);
       */
#else
    pthread_barrier_wait(pool_barrier);
#endif
    /*
       gettimeofday(&signal_tym, NULL);
       printf("tid = %d(%u) recieves signal at %f\n",0, (unsigned) pthread_self(), signal_tym.tv_sec + signal_tym.tv_usec*1e-6);
       */

    //    if(threads[0].myJob)
    {
        //do master job
        //START_TIME(work);
#if DEBUG
        printf("tid = %u doing job\n",(unsigned)pthread_self());
#endif
        threads[0].myJob.function(threads[0].myJob.arg);
#if DEBUG
        printf("tid = %u finished job\n",(unsigned)pthread_self());
#endif
        //STOP_TIME(work);

        threads[0].finishJob();
    }

}

/* Wait until all the submitted jobs have finished */

    template<typename arg_t>
void thpool<arg_t>::barrier()
{
    doJob();

    START_TIME(Barrier);

#ifdef PTHREAD_BARRIER
    pthread_barrier_wait(pool_barrier);
#else
    pthread_mutex_lock(thcount_lock);
    while(num_jobs != num_slaves) {
        pthread_cond_wait(threads_all_idle, thcount_lock);
    }
    num_jobs = 0;
    pthread_cond_init(threads_all_idle, NULL);
    pthread_mutex_unlock(thcount_lock);
#endif

    STOP_TIME(Barrier);
}

