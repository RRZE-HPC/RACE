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

#ifndef RACE_SPIN_BARRIER
#define RACE_SPIN_BARRIER

#include "utility.h"
#include <algorithm>
#include <sys/time.h>

//#define DEBUG_SPIN

struct spin_cond_t
{
    volatile unsigned int spinner;
    unsigned int maxCtr;
    unsigned int ctrLimit;
};

inline void spin_cond_init(spin_cond_t *obj, int maxCtr)
{
    obj->maxCtr = maxCtr;
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
//    if(__sync_lock_test_and_set(&(obj->spinner),1))
//	{

#ifdef DEBUG_SPIN
    printf("%u cpu=%d started spinning\n", (unsigned) pthread_self(), sched_getcpu());
#endif
    while(__sync_fetch_and_add(&(obj->spinner), 1))
    {
        DUMMY_spin(&(obj->spinner));
        //DUMMY(obj->spinner, true);
        if(obj->spinner > obj->ctrLimit)
        {
#ifdef DEBUG_SPIN
            printf("%u cpu=%d broke from loop\n", (unsigned) pthread_self(), sched_getcpu());
#endif
            //increase counter if it breaks
            break;
        }
    }
//}
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
