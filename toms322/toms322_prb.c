/*ftest.c - driver program for fisher.c.
  Copyright (c) 1995 by Matthew Belmonte.  Permission for use and distribution
  is hereby granted, subject to the restrictions that this copyright notice be
  included in its entirety, and that any and all changes made to the program be
  clearly noted in the program text.

  This software is provided 'as is', with no warranty, express or implied,
  including but not limited to warranties of merchantability or fitness for a
  particular purpose.  The user of this software assumes liability for any and
  all damages, whether direct or consequential, arising from its use.  The
  author of this software will not be liable for any such damages.
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>

# include "toms322.h"

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TOMS322_PRB.

  Discussion:

    TOMS322_PRB calls the TOMS322 tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2013

  Author:

    John Burkardt
*/
{
  int df;
  int df_num;
  double f;
  double t;
  char test_type;

  printf ( "\n" );
  printf ( "TOMS322_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the TOMS322 library.\n" );
  printf ( "\n" );
  printf("Enter test type, degrees of freedom, and value, or 'q' to quit -\nfor example 'f 1 10 10.04' or 't 10 3.169'\n");
  test_type = getchar();
  while(test_type != 'q')
    {
    if((test_type == 't') || (test_type == 'f'))
      {
      if(test_type != 't')
	scanf("%d", &df_num);
      else
	df_num = 1;
      scanf("%d", &df);
      if(test_type == 't')
	{
	scanf("%lf", &t);
	f = t*t;
	}
      else
	{
	scanf("%lf", &f);
	if(df_num == 1)
	  t = sqrt(f);
	}
      printf("\tF(%d,%d)=%f, ", df_num, df, f);
      if(df_num == 1)
	printf("t(%d)=%f, ", df, t);
      printf("p=%f\n", 1.0-fisher(df_num, df, f));
      }
    test_type = getchar();
    }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TOMS322_PRB:\n" );
  printf ( "  Normal end of execution\n" );

  return 0;
}
