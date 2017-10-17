#ifndef RACE_LOCK_LESS_H
#define RACE_LOCK_LESS_H

#include <pthread.h>
#include <linux/futex.h>
#include <stdint.h>
#include <limits.h>
#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define _GNU_SOURCE
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/time.h>

//Modified from https://locklessinc.com/articles/barriers/

/* Compile Barrier */
//#define barrier() asm volatile("": : :"memory")

/* Atomic add, returning the new value after the addition */
#define atomic_add(P, V) __sync_add_and_fetch((P), (V))

/* Atomic add, returning the value before the addition */
#define atomic_xadd(P, V) __sync_fetch_and_add((P), (V))

/* Atomic add, returning the new value after the addition */
#define atomic_sub(P, V) __sync_sub_and_fetch((P), (V))

/* Atomic add, returning the value before the addition */
#define atomic_xsub(P, V) __sync_fetch_and_sub((P), (V))


/* Atomic or */
#define atomic_or(P, V) __sync_or_and_fetch((P), (V))

/* Force a read of the variable */
#define atomic_read(V) (*(volatile typeof(V) *)&(V))

/* Atomic 32 bit exchange */
static inline unsigned xchg_32(void *ptr, unsigned x)
{
        __asm__ __volatile__("xchgl %0,%1"
                                :"=r" (x)
                                                :"m" (*(volatile unsigned *)ptr), "0" (x)
                                                                :"memory");

            return x;
}

/* Atomic 64 bit exchange */
static inline unsigned long long xchg_64(void *ptr, unsigned long long x)
{
        __asm__ __volatile__("xchgq %0,%1"
                                :"=r" (x)
                                                :"m" (*(volatile unsigned long long *)ptr), "0" (x)
                                                                :"memory");

            return x;
}


static long sys_futex(void *addr1, int op, int val1, struct timespec *timeout, void *addr2, int val3)
{
        return syscall(SYS_futex, addr1, op, val1, timeout, addr2, val3);
}

#endif
