# include <stdio.h>
# include <stdlib.h>
# include <pthread.h>

pthread_mutex_t count_mutex     = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t condition_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  condition_cond  = PTHREAD_COND_INITIALIZER;

int main ( void );
void *functionCount1 ( void );
void *functionCount2 ( void );

int count = 0;

# define COUNT_DONE  10
# define COUNT_HALT1  3
# define COUNT_HALT2  6

/******************************************************************************/

int main ( void )

/******************************************************************************/
{
  pthread_t thread1;
  pthread_t thread2;

  pthread_create( &thread1, NULL, &functionCount1, NULL);
  pthread_create( &thread2, NULL, &functionCount2, NULL);
  pthread_join( thread1, NULL);
  pthread_join( thread2, NULL);

  return 0;
}
/******************************************************************************/

void *functionCount1 ( void )

/******************************************************************************/
{
  for ( ; ; )
  {
    pthread_mutex_lock( &condition_mutex );
    while ( count >= COUNT_HALT1 && count <= COUNT_HALT2 )
    {
      pthread_cond_wait( &condition_cond, &condition_mutex );
    }
    pthread_mutex_unlock( &condition_mutex );

    pthread_mutex_lock( &count_mutex );
    count++;
    printf("Counter value functionCount1: %d\n",count);
    pthread_mutex_unlock( &count_mutex );

    if(count >= COUNT_DONE) 
    {
      return(NULL);
    }
  }
}
/******************************************************************************/

void *functionCount2 ( void )

/******************************************************************************/
{
  for ( ; ; )
  {
    pthread_mutex_lock( &condition_mutex );
    if ( count < COUNT_HALT1 || count > COUNT_HALT2 )
    {
      pthread_cond_signal( &condition_cond );
    }
    pthread_mutex_unlock( &condition_mutex );

    pthread_mutex_lock( &count_mutex );
    count++;
    printf("Counter value functionCount2: %d\n",count);
    pthread_mutex_unlock( &count_mutex );

    if(count >= COUNT_DONE) 
    {
      return(NULL);
    }
  }
}
