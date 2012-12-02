# include <stdio.h>
# include <stdlib.h>
# include <time.h>

# include "omp.h"

struct __omp_lock
{
  int lock;
};

enum { UNLOCKED = -1, INIT, LOCKED };

struct __omp_nest_lock
{
  short owner;
  short count;
};

enum { NOOWNER = -1, MASTER = 0 };

/******************************************************************************/

void omp_destroy_lock ( omp_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_DESTROY_LOCK destroys a simple lock.

  Discussion:

    The routine is intended to return the state of the lock to the 
    uninitialized state.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, omp_lock_t *ARG, the simple lock.
*/
{
  struct __omp_lock *lock = ( struct __omp_lock * ) arg;
  lock->lock = INIT;

  return;
}
/******************************************************************************/

void omp_destroy_nest_lock ( omp_nest_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_DESTROY_NEST_LOCK destroys a nestable lock.

  Discussion:

    The routine is intended to return the state of the lock to the 
    uninitialized state.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, omp_nest_lock_t *ARG, the nestable lock.
*/
{
  struct __omp_nest_lock *nlock = ( struct __omp_nest_lock * ) arg;
  nlock->owner = NOOWNER;
  nlock->count = UNLOCKED;
  return;
}
/******************************************************************************/

int omp_get_active_level ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_ACTIVE_LEVEL returns the number of nested active parallel regions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_ACTIVE_LEVEL, the number of nested active parallel 
    regions enclosing the task that contains this call.
*/
{
  return 0;
}
/******************************************************************************/

int omp_get_ancestor_thread_num ( int level )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_ANCESTOR_THREAD_NUM returns the thread number of the ancestor of this thread.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input, int LEVEL, the nested level.

    Output, int OMP_GET_ANCESTOR_THREAD_NUM, the thread number of the 
    ancestor of this thread.
*/
{
  if ( level == 0 )
  {
    return 0;
  }
  else
  {
    return -1;
  }
}
/******************************************************************************/

int omp_get_dynamic ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_DYNAMIC reports if dynamic adjustment of thread number is allowed.

  Discussion:

    The user can request dynamic thread adjustment by calling OMP_SET_DYNAMIC.

    For this stub library, the value FALSE is always returned.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_DYNAMIC, is TRUE (1) if dynamic adjustment of thread
    number has been enable, by default or by a user call, and FALSE (0)
    otherwise.
*/
{
  return 0;
}
/******************************************************************************/

int omp_get_level ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_LEVEL returns the number of nested parallel regions enclosing this task.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_LEVEL, the number of nested parallel regions 
    enclosing this task.
*/
{
  return 0;
}
/******************************************************************************/

int omp_get_max_active_levels ( void )

/******************************************************************************/
/*
  Purpose:

     OMP_GET_MAX_ACTIVE_LEVELS gets the maximum number of nested active parallel regions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_MAX_ACTIVE_LEVELS gets the maximum number of 
    nested active parallel regions.
*/
{
  return 0;
}
/******************************************************************************/

int omp_get_max_threads ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_MAX_THREADS returns the default number of threads.

  Discussion:

    If a parallel region is reached, and no number of threads has been
    specified explicitly, there is a default number of threads that will
    be used to form the new team.  That value is returned by this function.

    For this stub library, the value 1 is always returned.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_MAX_THREADS, the default number of threads.
*/
{
  return 1;
}
/******************************************************************************/

int omp_get_nested ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_NESTED reports if nested parallelism has been enabled.

  Discussion:

    The user can request nested parallelism by calling OMP_SET_NESTED.

    For this stub library, the value FALSE is always returned.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_NESTED, is TRUE (1) if nested parallelism has been
    enable by default or by a user call, and FALSE (0) otherwise.
*/
{
  return 0;
}
/******************************************************************************/

int omp_get_num_procs ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_NUM_PROCS returns the number of processors available to the program.

  Discussion:

    For this stub library, the value 1 is always returned.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_NUM_PROCS, the number of processors available.
*/
{
  return 1;
}
/******************************************************************************/

int omp_get_num_threads ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_NUM_THREADS returns the number of threads in the current team.

  Discussion:

    For this stub library, the value 1 is always returned.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_NUM_THREADS, the number of threads in the 
    current team.
*/
{
  return 1;
}
/******************************************************************************/

void omp_get_schedule ( omp_sched_t *kind, int *modifier )

/******************************************************************************/
/*
  Purpose:

     OMP_GET_SCHEDULE returns information about the "runtime" schedule.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, omp_sched_t *KIND, may be
    1, omp_sched_static,
    2, omp_sched_dynamic,
    3, omp_sched_guided,
    4, omp_sched_auto.

    Output, int *MODIFIER; this contains the "chunk_size" information for
    static, dynamic, or guided schedules, and is ignored for the auto schedule.
    If the chunk_size is less than 1, then the default value is used instead.
*/
{
  *kind = omp_sched_static;
  *modifier = 0;
  return;
}
/******************************************************************************/

int omp_get_team_size ( int level )

/******************************************************************************/
/*
  Purpose:

     OMP_GET_TEAM_SIZE returns the size of the thread team for a given level.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input, int LEVEL, the nested level.

    Output, int OMP_GET_TEAM_SIZE, the size of the thread team for 
    this level.
*/
{
  if ( level == 0 )
  {
    return 1;
  }
  else
  {
    return -1;
  }
}
/******************************************************************************/

int omp_get_thread_limit ( void )

/******************************************************************************/
/*
  Purpose:

     OMP_GET_THREAD_LIMIT returns the maximum number of OpenMP threads available.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_THREAD_LIMIT, the maximum number of OpenMP
    threads available.
*/
{
  return 1;
}
/******************************************************************************/

int omp_get_thread_num ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_THREAD_NUM returns the thread number of a thread in a team.

  Discussion:

    Thread numbers start at 0.

    If this function is not called from a parallel region, then only one
    thread is executing, so the value returned is 0.

    For this stub library, the value 0 is always returned.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_GET_THREAD_NUM, the thread number.
*/
{
  return 0;
}
/******************************************************************************/

double omp_get_wtick ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_WTICK returns the precision of the timer used by OMP_GET_WTIME.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, double OMP_GET_WTICK, the number of seconds between
    successive "ticks" of the wall clock timer.
*/
{
  double value;

  value = 1.0 / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

double omp_get_wtime ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_GET_WTIME returns elapsed wall clock time in seconds.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, double OMP_GET_WTIME, the current reading of the
    wall clock timer.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

int omp_in_final ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_IN_FINAL is true if the routine is executed in a final task region.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_IN_FINAL, is true if the routine is executed in a
    final task region.
*/
{
  return 1;
}
/******************************************************************************/

int omp_in_parallel ( void )

/******************************************************************************/
/*
  Purpose:

    OMP_IN_PARALLEL returns TRUE if the call is made from a parallel region.

  Discussion:

    For this stub library, the value FALSE is always returned.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, int OMP_IN_PARALLEL, is "TRUE" (1) if the routine was called
    from a parallel region and "FALSE" (0) otherwise.
*/
{
  return 0;
}
/******************************************************************************/

void omp_init_lock ( omp_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_INIT_LOCK initializes a simple lock.

  Discussion:

    This routine is intended to initialize the lock to the unlocked state.

    For this stub library, the lock points to -1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, omp_lock_t *ARG, the simple lock.
*/
{
  struct __omp_lock *lock = ( struct __omp_lock * ) arg;
  lock->lock = UNLOCKED;
  return;
}
/******************************************************************************/

void omp_init_nest_lock ( omp_nest_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_INIT_NEST_LOCK initializes a nestable lock.

  Discussion:

    This routine is intended to initialize the lock to the unlocked state.

    For this stub library, the lock points to -1.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Output, omp_nest_lock_t *ARG, the lock.
*/
{
  struct __omp_nest_lock *nlock = ( struct __omp_nest_lock * ) arg;
  nlock->owner = NOOWNER;
  nlock->count = 0;
  return;
}
/******************************************************************************/

void omp_set_dynamic ( int dynamic_threads )

/******************************************************************************/
/*
  Purpose:

    OMP_SET_DYNAMIC enables dynamic adjustment of the number of threads.

  Discussion:

    If DYNAMIC_THREADS is TRUE, then the number of threads available for
    execution in a parallel region may be dynamically adjusted.

    For this stub library, the input value is ignored.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input, int DYNAMIC_THREADS, is TRUE if the user wishes to allow
    dynamic adjustment of the number of threads available for execution
    in any parallel region.
*/
{
  return;
}
/******************************************************************************/

void omp_set_lock ( omp_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_SET_LOCK sets a simple lock.

  Discussion:

    The lock must already have been initialized.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input/output, omp_lock_t *ARG, the simple lock.
*/
{
  struct __omp_lock *lock = ( struct __omp_lock * ) arg;
  if ( lock->lock == UNLOCKED )
  {
    lock->lock = LOCKED;
  }
  else if ( lock->lock == LOCKED )
  {
    fprintf ( stderr, "error: deadlock in using lock variable\n" );
    exit ( 1 );
  }
  else
  {
    fprintf ( stderr, "error: lock not initialized\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void omp_set_max_active_levels ( int max_active_levels )

/******************************************************************************/
/*
  Purpose:

    OMP_SET_MAX_ACTIVE_LEVELS limits the number of nested active parallel regions.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input, int MAX_LEVELS, the maximum number of nested active parallel
    regions.
*/
{
  return;
}
/******************************************************************************/

void omp_set_nest_lock ( omp_nest_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_SET_NEST_LOCK sets a nestable lock.

  Discussion:

    The lock must already have been initialized.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input/output, omp_nest_lock_t *ARG, the nestable lock.
*/
{
  struct __omp_nest_lock *nlock = ( struct __omp_nest_lock * ) arg;
  if ( nlock->owner == MASTER && nlock->count >= 1 )
  {
    nlock->count++;
  }
  else if ( nlock->owner == NOOWNER && nlock->count == 0 )
  {
    nlock->owner = MASTER;
    nlock->count = 1;
  }
  else
  {
    fprintf ( stderr, "error: lock corrupted or not initialized\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void omp_set_nested ( int nested )

/******************************************************************************/
/*
  Purpose:

    OMP_SET_NESTED controls the use of nested parallelism.

  Discussion:

    For this stub library, the input value is ignored.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input, int NESTED, is TRUE (1) if nested parallelism is to be enabled.
*/
{
  return;
}
/******************************************************************************/

void omp_set_num_threads ( int num_threads )

/******************************************************************************/
/*
  Purpose:

    OMP_SET_NUM_THREADS sets the number of threads.

  Discussion:

    This routine sets the number of threads to be used in all subsequent
    parallel regions for which an explicit number of threads is not
    specified.

    For this stub library, the input value is ignored.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input, int NUM_THREADS, the number of threads to be used in all
    subsequent parallel regions for which an explicit number of threads
    is not specified.  0 < NUM_THREADS.
*/
{
  return;
}
/******************************************************************************/

void omp_set_schedule ( omp_sched_t kind, int modifier )

/******************************************************************************/
/*
  Purpose:

    OMP_SET_SCHEDULE chooses the schedule when "runtime" is the schedule kind.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input, omp_sched_t KIND, may be
    1, omp_sched_static,
    2, omp_sched_dynamic,
    3, omp_sched_guided,
    4, omp_sched_auto.

    Input, int MODIFIER; this contains the "chunk_size" information for
    static, dynamic, or guided schedules, and is ignored for the auto schedule.
    If the chunk_size is less than 1, then the default value is used instead.
*/
{
  return;
}
/******************************************************************************/

int omp_test_lock ( omp_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_TEST_LOCK tests a simple lock.

  Discussion:

    Calling this routine with an uninitialized lock causes an error.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input/output, omp_lock_t *ARG, the simple lock.
    If the lock was initialized but not set, it is set by this call.

    Output, int OMP_TEST_LOCK:
    TRUE (1), on input, the lock was initialized and not set;
    FALSE (0), on input the lock was initialized, and set.
*/
{
  struct __omp_lock *lock = ( struct __omp_lock * ) arg;
  if ( lock->lock == UNLOCKED )
  {
    lock->lock = LOCKED;
    return 1;
  }
  else if ( lock->lock == LOCKED )
  {
    return 0;
  }
  else
  {
    fprintf ( stderr, "error: lock not initialized\n" );
    exit ( 1 );
  }
}
/******************************************************************************/

int omp_test_nest_lock ( omp_nest_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_TEST_NEST_LOCK tests a nestable lock.

  Discussion:

    Calling this routine with an uninitialized lock causes an error.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input/output, omp_nest_lock_t *ARG, the nestable lock.
    If the lock was initialized but not set, it is set by this call.

    Output, int OMP_TEST_NEST_LOCK, returns the new nesting count,
    if the call was successful.  Otherwise, the value 0 is returned.
*/
{
  struct __omp_nest_lock *nlock = ( struct __omp_nest_lock * ) arg;
  omp_set_nest_lock ( arg) ;
  return nlock->count;
}
/******************************************************************************/

void omp_unset_lock ( omp_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_UNSET_LOCK unsets a simple lock.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input/output, omp_lock_t *ARG, the simple lock.
*/
{
  struct __omp_lock *lock = ( struct __omp_lock * ) arg;
  if ( lock->lock == LOCKED )
  {
    lock->lock = UNLOCKED;
  }
  else if ( lock->lock == UNLOCKED )
  {
    fprintf ( stderr, "error: lock not set\n" );
    exit ( 1 );
  }
  else
  {
    fprintf ( stderr, "error: lock not initialized\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void omp_unset_nest_lock ( omp_nest_lock_t *arg )

/******************************************************************************/
/*
  Purpose:

    OMP_UNSET_NEST_LOCK unsets a nestable lock.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 May 2012

  Author:

    John Burkardt

  Reference:

    OpenMP Application Program Interface,
    Version 3.1,
    July 2011.

  Parameters:

    Input/output, omp_nest_lock_t *ARG, the nestable lock.
*/
{
  struct __omp_nest_lock *nlock = ( struct __omp_nest_lock * ) arg;
  if ( nlock->owner == MASTER && nlock->count >= 1 )
  {
    nlock->count--;
    if ( nlock->count == 0 )
    {
      nlock->owner = NOOWNER;
    }
  }
  else if ( nlock->owner == NOOWNER && nlock->count == 0 )
  {
    fprintf ( stderr, "error: lock not set\n" );
    exit ( 1 );
  }
  else
  {
    fprintf ( stderr, "error: lock corrupted or not initialized\n" );
    exit ( 1 );
  }
  return;
}
