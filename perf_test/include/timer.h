#ifndef RACE_TIMER_H
#define RACE_TIMER_H

#include <sys/time.h>

#ifdef TIMING
#define INIT_TIMER(region)\
timeval start##region, end##region;\
double start_tym##region, end_tym##region;\

#define START_TIMER(region)\
gettimeofday(&start##region, NULL); \
start_tym##region = start##region.tv_sec + start##region.tv_usec*1e-6;\

#define STOP_TIMER(region)\
gettimeofday(&end##region, NULL);\
end_tym##region = end##region.tv_sec + end##region.tv_usec*1e-6;\


#define PRINT_TIMER(region)\
printf("%s : Time = %\nf", #region, GET_TIMER(region));\

#define GET_TIMER(region)\
(end_tym##region-start_tym##region)\


#else
#define START_TIMER(region)\

#define STOP_TIMER(region)\

#define PRINT_TIMER(region)\

#define GET_TIMER(region)\

#endif

#endif
