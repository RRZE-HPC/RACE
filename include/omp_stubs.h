// DL 27.10.2025 To build without OpenMP

#ifndef RACE_OMP_STUBS_H
#define RACE_OMP_STUBS_H

/*
 * OpenMP compatibility stubs
 * These provide single-threaded, no-op replacements for OpenMP functions
 * when OpenMP is not enabled (_OPENMP not defined).
 */

#ifndef _OPENMP  /* Only define these if OpenMP is disabled */

#include <stddef.h>  /* for NULL */
#ifdef __cplusplus
extern "C" {
#endif

/* --- Environment control --- */
static inline void omp_set_num_threads(int x)          { (void)x; }
static inline int  omp_get_num_threads(void)           { return 1; }
static inline int  omp_get_max_threads(void)           { return 1; }
static inline int  omp_get_thread_num(void)            { return 0; }
static inline int  omp_get_num_procs(void)             { return 1; }
static inline int  omp_in_parallel(void)               { return 0; }

/* --- Dynamic threads --- */
static inline void omp_set_dynamic(int x)              { (void)x; }
static inline int  omp_get_dynamic(void)               { return 0; }

/* --- Nested parallelism --- */
static inline void omp_set_nested(int x)               { (void)x; }
static inline int  omp_get_nested(void)                { return 0; }
static inline void omp_set_max_active_levels(int x)    { (void)x; }
static inline int  omp_get_max_active_levels(void)     { return 1; }
static inline int  omp_get_level(void)                 { return 0; }
static inline int  omp_get_active_level(void)          { return 0; }
static inline int  omp_get_ancestor_thread_num(int l)  { (void)l; return 0; }
static inline int  omp_get_team_size(int l)            { (void)l; return 1; }

/* --- Scheduling --- */
static inline void omp_set_schedule(int kind, int chunk)
{ (void)kind; (void)chunk; }
static inline void omp_get_schedule(int *kind, int *chunk)
{ if (kind) *kind = 0; if (chunk) *chunk = 0; }
static inline int  omp_get_thread_limit(void)          { return 1; }

/* --- Timing --- */
static inline double omp_get_wtime(void)               { return 0.0; }
static inline double omp_get_wtick(void)               { return 1.0; }

/* --- Locks --- */
typedef void* omp_lock_t;
typedef void* omp_nest_lock_t;

static inline void omp_init_lock(omp_lock_t *l)        { (void)l; }
static inline void omp_destroy_lock(omp_lock_t *l)     { (void)l; }
static inline void omp_set_lock(omp_lock_t *l)         { (void)l; }
static inline void omp_unset_lock(omp_lock_t *l)       { (void)l; }
static inline int  omp_test_lock(omp_lock_t *l)        { (void)l; return 1; }

static inline void omp_init_nest_lock(omp_nest_lock_t *l)    { (void)l; }
static inline void omp_destroy_nest_lock(omp_nest_lock_t *l) { (void)l; }
static inline void omp_set_nest_lock(omp_nest_lock_t *l)     { (void)l; }
static inline void omp_unset_nest_lock(omp_nest_lock_t *l)   { (void)l; }
static inline int  omp_test_nest_lock(omp_nest_lock_t *l)    { (void)l; return 1; }

/* --- Places, affinity, proc bindings --- */
static inline int  omp_get_num_places(void)            { return 1; }
static inline int  omp_get_place_num_procs(int x)      { (void)x; return 1; }
static inline void omp_get_place_proc_ids(int x, int *y)
{ (void)x; (void)y; }
static inline int  omp_get_place_num(void)             { return 0; }
static inline int  omp_get_partition_num_places(void)  { return 1; }
static inline void omp_get_partition_place_nums(int *x){ (void)x; }
static inline void omp_set_affinity_format(const char *fmt) { (void)fmt; }
static inline void omp_get_affinity_format(char *buf, size_t size)
{ (void)buf; (void)size; }
static inline void omp_display_affinity(const char *fmt) { (void)fmt; }
static inline void omp_capture_affinity(char *buf, size_t size, const char *fmt)
{ (void)buf; (void)size; (void)fmt; }

/* --- Parallel region state --- */
static inline int  omp_in_final(void)                  { return 0; }
static inline int  omp_get_cancellation(void)          { return 0; }
static inline int  omp_get_proc_bind(void)             { return 0; }

/* --- Tasks --- */
static inline int  omp_get_num_teams(void)             { return 1; }
static inline int  omp_get_team_num(void)              { return 0; }
static inline int  omp_get_default_device(void)        { return 0; }
static inline void omp_set_default_device(int x)       { (void)x; }
static inline int  omp_get_initial_device(void)        { return 0; }

/* --- Memory allocators / target offload --- */
static inline void* omp_target_alloc(size_t x, int y)  { (void)x; (void)y; return NULL; }
static inline void  omp_target_free(void *x, int y)    { (void)x; (void)y; }
static inline int   omp_target_is_present(const void *x, int y)
{ (void)x; (void)y; return 0; }
static inline void  omp_target_memcpy(void *dst, const void *src,
                                      size_t n, size_t off1, size_t off2,
                                      int dev1, int dev2)
{ (void)dst; (void)src; (void)n; (void)off1; (void)off2; (void)dev1; (void)dev2; }
static inline void  omp_target_memcpy_rect(void)       { }
static inline void  omp_target_associate_ptr(void)     { }
static inline void  omp_target_disassociate_ptr(void *x, int y)
{ (void)x; (void)y; }

/* --- Cancellation / errors --- */
static inline int  omp_get_max_task_priority(void)     { return 0; }

/* --- Nothing parallel actually happens --- */
static inline void omp_parallel_for(void)              { /* no-op */ }

/* Intel OpenMP runtime stub (non-standard) */
static inline void kmp_set_warnings_off(void)          { }

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* !_OPENMP */

#endif /* RACE_OMP_STUBS_H */
