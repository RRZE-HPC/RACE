#ifndef NAME_SPIN_BARRIER
#define NAME_SPIN_BARRIER

#include "utility.h"
#include <algorithm>
#include <sys/time.h>

//#define DEBUG_SPIN

struct spin_cond_t
{
    volatile unsigned int spinner;
    unsigned int ctrLimit;
};

inline void spin_cond_init(spin_cond_t *obj, int maxCtr)
{
    obj->spinner = 1;
    obj->ctrLimit = maxCtr;
}

//ends the current spinning
inline void spin_cond_signal(spin_cond_t *obj)
{
    __sync_lock_release(&(obj->spinner));

#ifdef DEBUG_SPIN
    timeval signal_time;
    gettimeofday(&signal_time, NULL);
    printf("%u cpu = %d signalled at %f, signal value = %u\n",(unsigned) pthread_self(), sched_getcpu(), signal_time.tv_sec+signal_time.tv_usec*1e-6, obj->spinner);
#endif
    /*    __sync_synchronize();//fetch_and_add(&(obj->spinner),0);
    //reset only if partner/s is still waiting
    printf("%u spinner value = %u\n",(unsigned)pthread_self(), (unsigned)obj->spinner);
    if(obj->spinner != 1)
    {
        timeval signal_time;
        gettimeofday(&signal_time, NULL);
        printf("%u really signalled at %f, pre-signal value = %u\n",(unsigned) pthread_self(), signal_time.tv_sec+signal_time.tv_usec*1e-6, obj->spinner);
        obj->spinner = 0;
    }*/
}

inline void spin_cond_wait(spin_cond_t *obj)
{
    //reset value in any case; don't consider signals issued before
    __sync_lock_test_and_set(&(obj->spinner),1);

#ifdef DEBUG_SPIN
    printf("%u cpu=%d started spinning\n", (unsigned) pthread_self(), sched_getcpu());
#endif

    while(__sync_fetch_and_add(&(obj->spinner), 1))
    {
        DUMMY(obj->spinner, true);
        if(obj->spinner > obj->ctrLimit)
        {
#ifdef DEBUG_SPIN
            printf("%u cpu=%d broke from loop\n", (unsigned) pthread_self(), sched_getcpu());
#endif
            break;
        }
    }
#ifdef DEBUG_SPIN
    printf("%u cpu=%d stopped spinning\n", (unsigned) pthread_self(), sched_getcpu());
#endif


#ifdef DEBUG_SPIN
    timeval signal_time;
    gettimeofday(&signal_time, NULL);
    printf("%u cpu=%d recieved at %f\n",(unsigned) pthread_self(), sched_getcpu(), signal_time.tv_sec+signal_time.tv_usec*1e-6);
#endif
}


#endif
