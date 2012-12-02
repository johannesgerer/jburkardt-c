#ifndef _OMP_H_DEF
#define _OMP_H_DEF
/*
  Define the lock data types
*/
typedef void *omp_lock_t;
typedef void *omp_nest_lock_t;
/*
  Define the schedule kinds.
  Add vendor specific schedule constants at end.
*/
typedef enum omp_sched_t
{
  omp_sched_static = 1,
  omp_sched_dynamic = 2,
  omp_sched_guided = 3,
  omp_sched_auto = 4
} omp_sched_t;
/*
  Exported OpenMP functions
*/
#ifdef __cplusplus
extern "C"
{
#endif

extern void   omp_destroy_lock ( omp_lock_t *lock );
extern void   omp_destroy_nest_lock ( omp_nest_lock_t *lock );
extern int    omp_get_active_level ( void );
extern int    omp_get_ancestor_thread_num ( int level );
extern int    omp_get_dynamic ( void );
extern int    omp_get_level ( void );
extern int    omp_get_max_active_levels ( void );
extern int    omp_get_max_threads ( void );
extern int    omp_get_nested ( void );
extern int    omp_get_num_procs ( void );
extern int    omp_get_num_threads ( void );
extern void   omp_get_schedule ( omp_sched_t *kind, int *modifier );
extern int    omp_get_team_size ( int level );
extern int    omp_get_thread_limit ( void );
extern int    omp_get_thread_num ( void );
extern double omp_get_wtick ( void );
extern double omp_get_wtime ( void );
extern int    omp_in_final ( void );
extern int    omp_in_parallel ( void );
extern void   omp_init_lock ( omp_lock_t *lock );
extern void   omp_init_nest_lock ( omp_nest_lock_t *lock );
extern void   omp_set_dynamic ( int dynamic_threads );
extern void   omp_set_lock ( omp_lock_t *lock );
extern void   omp_set_max_active_levels ( int max_active_levels );
extern void   omp_set_nest_lock ( omp_nest_lock_t *lock );
extern void   omp_set_nested ( int nested );
extern void   omp_set_num_threads ( int num_threads );
extern void   omp_set_schedule ( omp_sched_t kind, int modifier );
extern int    omp_test_lock ( omp_lock_t *lock );
extern int    omp_test_nest_lock ( omp_nest_lock_t *lock );
extern void   omp_unset_lock ( omp_lock_t *lock );
extern void   omp_unset_nest_lock ( omp_nest_lock_t *lock );

#ifdef __cplusplus
}
#endif
#endif