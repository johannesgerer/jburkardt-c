# include <stdlib.h>
# include <stdio.h>

# include "asa172.h"

int main ( void );
void test01 ( void );
void test02 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ASA172_PRB.

  Discussion:

    ASA172_PRB tests the ASA172 library.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 July 2008

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "ASA172_PRB:\n" );
  printf ( "  C version\n" );
  printf ( "  Test the ASA172 library.\n" );

  test01 ( );
  test02 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ASA172_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01 compares indices computed by a triple loop.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 July 2008

  Author:

    John Burkardt
*/
{
# define KDIM 3

  int i;
  int i1;
  int i2;
  int i3;
  int ifault;
  int iprod[KDIM];
  int ivec[KDIM];
  int j;
  int jsub;
  int kdim = KDIM;
  int n;
  int nr[KDIM] = { 3, 2, 4 };
  int qfor;
  int qind;

  printf ( "\n" );
  printf ( "TEST01:\n" );
  printf ( "  SIMDO can convert between compressed and\n" );
  printf ( "  vector indices representing a nested loop.\n" );
  printf ( "\n" );
  printf ( "  Here, we set QFOR = FALSE, meaning we do\n" );
  printf ( "  NOT want to convert from FORTRAN ordering\n" );
  printf ( "  to lexical ordering.\n" );
  printf ( "\n" );
  printf ( "  Here, we actually carry out a triple loop\n" );
  printf ( "  list the indices, and then compare.\n" );

  qfor = 0;
/*
  If QFOR is FALSE, then the definition of IPROD is reversed...
*/
  iprod[0] = nr[kdim-1];
  for ( i = 1; i < kdim; i++ )
  {
    iprod[i] = iprod[i-1] * nr[kdim-1-i];
  }

  n = iprod[kdim-1];
/*
  Carry out the nested loops, and use JSUB to count each iteration.
  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
*/
  jsub = 0;

  printf ( "\n" );
  printf ( "  #1: Generate JSUB by counting as we DO the loops:\n" );
  printf ( "\n" );
  printf ( "  DO I1 = 1, N1\n" );
  printf ( "    DO I2 = 1, N2\n" );
  printf ( "      DO I3 = 1, N3\n" );
  printf ( "\n" );
  printf ( "      JSUB            I1        I2        I3\n" );
  printf ( "\n" );
  for ( i1 = 1; i1 <= nr[0]; i1++ )
  {
    ivec[0] = i1;
    for ( i2 = 1; i2 <= nr[1]; i2++ )
    {
      ivec[1] = i2;
      for ( i3 = 1; i3 <= nr[2]; i3++ )
      {
        ivec[2] = i3;
        jsub = jsub + 1;
        printf ( "  %8d      %8d  %8d  %8d\n", jsub, i1, i2, i3 );
      }
    }
  }
/*
  Now for each value of JSUB, retrieve the corresponding index subscript.
  In order to use the QFOR = .FALSE. switch, I apparently have to reverse
  the sense of the NR vector.
*/
  qind = 1;

  printf ( "\n" );
  printf ( "  #2: Loop on JSUB, retrieve loop indices\n" );
  printf ( "      QIND = TRUE J ->I(J)\n" );
  printf ( "      QFOR = FALSE\n" );
  printf ( "\n" );
  printf ( "      JSUB            I1        I2        I3\n" );
  printf ( "\n" );

  for ( j = 1; j <= n; j++ )
  {
    jsub = j;
    ifault= simdo ( qind, qfor, iprod, kdim, &jsub, ivec );
    printf ( "  %8d      %8d  %8d  %8d\n", jsub, ivec[0], ivec[1], ivec[2] );
  }
/*
  Carry out the nested loops, and DO NOT compute JSUB.
  Have SIMDO determine JSUB.
*/
  qind = 0;

  printf ( "\n" );
  printf ( "  #3: For any set of loop indices, retrieve JSUB\n" );
  printf ( "      QIND = FALSE I(J) -> J\n" );
  printf ( "      QFOR = FALSE\n" );
  printf ( "\n" );
  printf ( "      JSUB            I1        I2        I3\n" );
  printf ( "\n" );
  for ( i1 = 1; i1 <= nr[0]; i1++ )
  {
    ivec[0] = i1;
    for ( i2 = 1; i2 <= nr[1]; i2++ )
    {
      ivec[1] = i2;
      for ( i3 = 1; i3 <= nr[2]; i3++ )
      {
        ivec[2] = i3;
        ifault= simdo ( qind, qfor, iprod, kdim, &jsub, ivec );
        printf ( "  %8d      %8d  %8d  %8d\n", jsub, i1, i2, i3 );
      }
    }
  }
  return;
# undef KDIM
}
/******************************************************************************/

void test02 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST02 compares indices computed by a triple loop.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    27 July 2008

  Author:

    John Burkardt
*/
{
# define KDIM 3

  int i;
  int i1;
  int i2;
  int i3;
  int ifault;
  int iprod[KDIM];
  int ivec[KDIM];
  int j;
  int jsub;
  int kdim = KDIM;
  int n;
  int nr[KDIM] = { 3, 2, 4 };
  int qfor;
  int qind;

  printf ( "\n" );
  printf ( "TEST02:\n" );
  printf ( "  SIMDO can convert between compressed and\n" );
  printf ( "  vector indices representing a nested loop.\n" );
  printf ( "\n" );
  printf ( "  Here, we set QFOR = TRUE, meaning we DO\n" );
  printf ( "  want to convert from the FORTRAN \n" );
  printf ( "  ordering to lexical convention.\n" );
  printf ( "\n" );
  printf ( "  Here, we actually carry out a triple loop\n" );
  printf ( "  list the indices, and then compare.\n" );

  qfor = 1;

  iprod[0] = nr[0];
  for ( i = 1; i < kdim; i++ )
  {
    iprod[i] = iprod[i-1] * nr[i];
  }

  n = iprod[kdim-1];
/*
  Carry out the nested loops, and use JSUB to count each iteration.
  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
*/
  jsub = 0;

  printf ( "\n" );
  printf ( "  #1: Generate JSUB by counting as we do the loops.\n" );
  printf ( "\n" );
  printf ( "  DO I3 = 1, N3\n" );
  printf ( "    DO I2 = 1, N2\n" );
  printf ( "      DO I1 = 1, N1\n" );
  printf ( "\n" );
  printf ( "      JSUB            I1        I2        I3\n" );
  printf ( "\n" );
  for ( i3 = 1; i3 <= nr[2]; i3++ )
  {
    ivec[2] = i3;
    for ( i2 = 1; i2 <= nr[1]; i2++ )
    {
      ivec[1] = i2;
      for ( i1 = 1; i1 <= nr[0]; i1++ )
      {
        ivec[0] = i1;
        jsub = jsub + 1;
        printf ( "  %8d      %8d  %8d  %8d\n", jsub, i1, i2, i3 );
      }
    }
  }
/*
  Reverse the order, so that the loop indices are generated in lexical order.
*/
  qind = 1;

  printf ( "\n" );
  printf ( "  #2: Setting QFOR false means loop indices\n" );
  printf ( "  are generated in lexical order.\n" );
  printf ( "      QIND = TRUE J -> I(J)\n" );
  printf ( "      QFOR = TRUE\n" );
  printf ( "\n" );
  printf ( "      JSUB            I1        I2        I3\n" );
  printf ( "\n" );

  for ( j = 1; j <= n; j++ )
  {
    jsub = j;
    ifault= simdo ( qind, qfor, iprod, kdim, &jsub, ivec );
    printf ( "  %8d      %8d  %8d  %8d\n", jsub, ivec[0], ivec[1], ivec[2] );
  }
/*
  Carry out the nested loops, and DO NOT compute JSUB.
  Have SIMDO determine JSUB.
*/
  qind = 0;

  printf ( "\n" );
  printf ( "  #3: For any set of loop indices, retrieve JSUB\n" );
  printf ( "      QIND = FALSE I(J) -> J\n" );
  printf ( "      QFOR = TRUE\n" );
  printf ( "\n" );
  printf ( "      JSUB            I1        I2        I3\n" );
  printf ( "\n" );
  for ( i3 = 1; i3 <= nr[2]; i3++ )
  {
    ivec[2] = i3;
    for ( i2 = 1; i2 <= nr[1]; i2++ )
    {
      ivec[1] = i2;
      for ( i1 = 1; i1 <= nr[0]; i1++ )
      {
        ivec[0] = i1;
        ifault= simdo ( qind, qfor, iprod, kdim, &jsub, ivec );
        printf ( "  %8d      %8d  %8d  %8d\n", jsub, i1, i2, i3 );
      }
    }
  }
  return;
}
