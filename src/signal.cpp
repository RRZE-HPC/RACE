#include "signal.h"

Signal::Signal(int NAME_BLOCKCTR)
{
    lock = new pthread_mutex_t;
    pthread_mutex_init(lock, NULL);
    signal = new pthread_cond_t;
    pthread_cond_init(signal, NULL);
    preSignal = new spin_cond_t;
    spin_cond_init(preSignal, NAME_BLOCKCTR);
}

Signal::~Signal()
{
    delete lock;
    delete signal;
    delete preSignal;
}
