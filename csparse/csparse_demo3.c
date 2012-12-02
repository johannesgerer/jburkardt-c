# include <stdlib.h>
# include <limits.h>
# include <math.h>
# include <stdio.h>

# include "csparse.h"

typedef struct problem_struct
{
  cs *A ;
  cs *C ;
  int sym ;
  double *x ;
  double *b ;
  double *r ;
} problem ;

problem *get_problem (FILE *f, double tol) ;
int demo3 (problem *Prob) ;
problem *free_problem (problem *Prob) ;


/* cs_demo3: read a matrix and test Cholesky update/downdate */
int main (void)
{
    problem *Prob = get_problem (stdin, 0) ;
    demo3 ( Prob ) ;
    free_problem (Prob) ;
    return (0) ;
}
