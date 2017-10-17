#ifndef NAME_TIMING_H
#define NAME_TIMING_H

#include <sys/time.h>


#define TIMING

#ifdef TIMING
#define INIT_TIME(region)\
timeval start##region, end##region;\
double start_tym##region, end_tym##region;\

#define START_TIME(region)\
gettimeofday(&start##region, NULL); \
start_tym##region = start##region.tv_sec + start##region.tv_usec*1e-6;\

#define STOP_TIME(region)\
gettimeofday(&end##region, NULL);\
end_tym##region = end##region.tv_sec + end##region.tv_usec*1e-6;\


#define PRINT_TIME(region)\
printf("%s : Time = %\nf", #region, GET_TIME(region));\

#define GET_TIME(region)\
(end_tym##region-start_tym##region)\


#else
#define START_TIME(region)\

#define STOP_TIME(region)\

#define PRINT_TIME(region)\

#define GET_TIME(region)\

#endif

#endif
