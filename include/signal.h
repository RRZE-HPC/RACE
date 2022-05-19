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
        signalObj->preSignal->spinner*=2;\
        __sync_fetch_and_add(&signalObj->preSignal->spinner,0);\
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
