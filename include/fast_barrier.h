#ifndef RACE_FAST_BARRIER
#define RACE_FAST_BARRIER

#include "lock_less.h"

//Taken from https://locklessinc.com/articles/barriers/



typedef struct fast_barrier_t fast_barrier_t;
struct fast_barrier_t
{
    union
    {
        struct
        {
            unsigned seq;
            unsigned count;
        };
        unsigned long long reset;
    };
    unsigned refcount;
    unsigned total;
    int spins;
    unsigned flags;
};

inline void fast_barrier_init(fast_barrier_t *b, pthread_barrierattr_t *a, unsigned count)
{
    b->seq = 0;
    b->count = 0;
    b->refcount = 1;
    /* Total of waiting threads */
    b->total = count - 1;

    if (count < sysconf(_SC_NPROCESSORS_ONLN))
    {
        b->spins = 1000;
    }
    else
    {
        b->spins = 1;
    }

    /* Default to process private */
    b->flags = FUTEX_PRIVATE_FLAG;
    if (a)
    {
        /* Check for a shared
         * barrier */
        int shared = 0;
        pthread_barrierattr_getpshared(a, &shared);
        if (shared) b->flags = 0;
    }
}

inline void fast_barrier_destroy(fast_barrier_t *b)
{
    /* Trigger futex wake */
    atomic_sub(&b->refcount, 1);

    /* Wait until all references to the barrier are dropped */
    while (1)
    {
        unsigned ref = atomic_read(b->refcount);

        if (!ref) return;

        sys_futex(&b->refcount, FUTEX_WAIT | b->flags, ref, NULL, NULL, 0);
    }
}

inline int fast_barrier_wait(fast_barrier_t *b)
{
    int ret;

    atomic_add(&b->refcount, 1);

    while (1)
    {
        unsigned seq = atomic_read(b->seq);
        unsigned count = atomic_xadd(&b->count, 1);

        if (count < b->total)
        {
            int i;
            seq |= 1;

            for (i = 0; i < b->spins; i++)
            {
                if ((atomic_read(b->seq) | 1) != seq) break;
            }

            /* Can
             * we
             * proceed?
             * */
            while ((atomic_read(b->seq) | 1) == seq)
            {
                /* Hack
                 * - set
                 *   a
                 *   flag
                 *   that
                 *   says
                 *   we
                 *   are
                 *   sleeping
                 *   */
                *(volatile char *) &b->seq = 1;

                /* Sleep
                 * on
                 * it
                 * instead
                 * */
                sys_futex(&b->seq, FUTEX_WAIT | b->flags, seq, NULL, NULL, 0);
            }

            ret = 0;
            break;
        }

        if (count == b->total)
        {
            /* Simultaneously
             * clear
             * count,
             * increment
             * sequence
             * number,
             * and
             * clear
             * wait
             * flag
             * */
            seq = atomic_read(b->seq);

            if (xchg_64(&b->reset, (seq | 1) + 255) & 1)
            {
                /* Wake
                 * up
                 * sleeping
                 * threads
                 * */
                sys_futex(&b->seq, FUTEX_WAKE | b->flags, INT_MAX, NULL, NULL, 0);
            }

            ret = PTHREAD_BARRIER_SERIAL_THREAD;
            break;
        }

        seq |= 1;
        /* Hack
         * - set
         *   a
         *   flag
         *   that
         *   says
         *   we
         *   are
         *   sleeping
         *   */
        *(volatile char *) &b->seq = 1;

        /* We
         * were
         * too
         * slow...
         * wait
         * for
         * the
         * barrier
         * to
         * be
         * released
         * */
        sys_futex(&b->seq, FUTEX_WAIT | b->flags, seq, NULL, NULL, 0);
    }

    /* Are we the last to wake up? */
    if (atomic_xsub(&b->refcount, 1) == 1)
    {
        /* Wake destroying thread */
        sys_futex(&b->refcount, FUTEX_WAKE | b->flags, 1, NULL, NULL, 0);
    }
    return ret;
}

#endif
