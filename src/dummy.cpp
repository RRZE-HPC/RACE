#include "utility.h"
#include "macros.h"

void DUMMY(int ctr, bool flag)
{
    //This is just a dummy fn to avoid compiler optimisation
    UNUSED(ctr);
    UNUSED(flag);
}

void DUMMY_spin(volatile unsigned int *a)
{
    UNUSED(a);
}
