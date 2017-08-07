/* vim: set filetype=cpp: */

#include <functional>
#include <sys/time.h>
#include "timing.h"
#include "test.h"
#include<malloc.h>

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
    jobPresent=true;
    //send signal to tell the thread to wake and destroy
    sendSignal(signal);
}

    template<typename arg_t>
void thread<arg_t>::init (int gid_, thpool<arg_t>* pool_)
{
    gid = gid_;
    pool = pool_;

    signal = new Signal(pool->NAME_BLOCKCTR);

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
void run(thread<arg_t> * thread)
{
    thpool<arg_t>* pool = thread->pool;
    pthread_mutex_lock(pool->thcount_lock);
    pool->num_threads_alive++;
    pthread_mutex_unlock(pool->thcount_lock);

    while(!(__sync_fetch_and_add(&(pool->interrupt),0)))
    {
        //wait till jobPresent becomes true
        waitSignal(thread->signal, thread->jobPresent, true);

        if(thread->jobPresent && (!pool->interrupt))
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
        }
    }

    pthread_mutex_lock(pool->thcount_lock);
    pool->num_threads_alive--;
    pthread_mutex_unlock(pool->thcount_lock);
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
        jobPresent = false;
#if DEBUG
        printf("tid = %d(%u) finished job\n", gid, (unsigned) (* (pthread)));
#endif

        if(__sync_add_and_fetch(&(workerTeam->num_jobs), 1) == workerTeam->num_slaves)
        {
#if DEBUG
            printf("tid = %d(%u) sending barrier signal\n", gid, (unsigned) (* (pthread)));
#endif
            sendSignal(workerTeam->barrierSignal);
#if DEBUG
            printf("tid = %d(%u) sent barrier signal\n", gid, (unsigned) (* (pthread)));
#endif
        }
    }
}

    template<typename arg_t>
thpool<arg_t>::thpool():initialized(false)
{
    char* NAME_BLOCKCTR_str = getenv("NAME_BLOCKCTR");
    if(NAME_BLOCKCTR_str == NULL)
    {
        NAME_BLOCKCTR = 200000;
    }
    else {
        NAME_BLOCKCTR = atoi(NAME_BLOCKCTR_str);
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
        interrupt = false;

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
        interrupt = true;

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

        NAME_BLOCKCTR = pool->NAME_BLOCKCTR;

        barrierSignal = new Signal(NAME_BLOCKCTR);
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
    taskForce[0]->jobPresent = false;

#if DEBUG
    printf("master: tid = %u finished job\n",(unsigned) (pthread_self()));
#endif

    /*timeval start, end;
      double startTime, endTime;

      gettimeofday(&start, NULL);		
      startTime = start.tv_sec + start.tv_usec*1e-6;*/
    //START_TIME(Barrier);

    waitSignal(barrierSignal, num_jobs, num_slaves);
    __sync_lock_test_and_set(&num_jobs, 0);

    /*gettimeofday(&end, NULL);		
      endTime = end.tv_sec + end.tv_usec*1e-6;
      barrierTime = endTime - startTime;*/

#if DEBUG
    printf("master: tid = %u barrier done\n",(unsigned) (pthread_self()));
#endif


    // STOP_TIME(Barrier);
}

