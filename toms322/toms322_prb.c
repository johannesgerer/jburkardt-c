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

int main ( void )
{
  int df;
  int df_num;
  double f;
  double t;
  char test_type;

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
  return 0;
}
