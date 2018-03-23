#ifndef RACE_SIGNAL_H
#define RACE_SIGNAL_H

#include <pthread.h>
#include  "spin_cond.h"

enum SignalMode { spin=0, idle=1, released=2 };

struct Signal{
    pthread_mutex_t*  lock;
    pthread_cond_t*   signal;
    spin_cond_t* preSignal;
    double signalTime;
    Signal(int RACE_BLOCKCTR);
    ~Signal();
    SignalMode mode;
};


/*wait till predicate_lhs != predicate_rhs
 * last argument is for resetting any data*/
#define waitSignal(signalObj, predicate_lhs, predicate_rhs)\
    /*if(__sync_fetch_and_add(&predicate_lhs,0) != predicate_rhs) {*/\
        __sync_lock_test_and_set(&(signalObj->mode), spin);\
        spin_cond_wait(signalObj->preSignal);\
    /* }*/\
    if(__sync_fetch_and_add(&predicate_lhs,0) != predicate_rhs)\
    {\
        pthread_mutex_lock(signalObj->lock);\
        __sync_lock_test_and_set(&(signalObj->mode), idle);\
        while(__sync_fetch_and_add(&predicate_lhs,0) != predicate_rhs) {\
            pthread_cond_wait(signalObj->signal, signalObj->lock);\
        }\
        pthread_mutex_unlock(signalObj->lock);\
    }\
    /*reset value for next iteration*/\
    __sync_lock_test_and_set(&(signalObj->mode), released);\
    __sync_lock_test_and_set(&(signalObj->preSignal->spinner),1);\


#define sendSignal(signalObj)\
    spin_cond_signal(signalObj->preSignal);\
    pthread_mutex_lock(signalObj->lock);\
    pthread_cond_signal(signalObj->signal);\
    pthread_mutex_unlock(signalObj->lock);\

#define sleepSignal(signalObj)\
    spin_cond_signal(signalObj->preSignal);\



#endif
