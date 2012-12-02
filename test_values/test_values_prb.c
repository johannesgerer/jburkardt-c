# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex.h>

# include "test_values.h"

int main ( void );

void test001 ( void );
void test002 ( void );
void test003 ( void );
void test0035 ( void );
void test004 ( void );
void test005 ( void );
void test006 ( void );
void test007 ( void );
void test008 ( void );
void test009 ( void );
void test0093 ( void );
void test0095 ( void );

void test010 ( void );
void test011 ( void );
void test0114 ( void );
void test01145 ( void );
void test0115 ( void );
void test01155 ( void );
void test0116 ( void );
void test012 ( void );
void test0123 ( void );
void test0127 ( void );
void test0128 ( void );
void test013 ( void );
void test0134 ( void );
void test0135 ( void );
void test014 ( void );
void test015 ( void );
void test016 ( void );
void test017 ( void );
void test018 ( void );
void test0185 ( void );
void test019 ( void );
void test0195 ( void );

void test020 ( void );
void test0205 ( void );
void test021 ( void );
void test022 ( void );
void test023 ( void );
void test024 ( void );
void test025 ( void );
void test026 ( void );
void test0265 ( void );
void test027 ( void );
void test028 ( void );
void test029 ( void );

void test030 ( void );
void test0305 ( void );
void test031 ( void );
void test032 ( void );
void test033 ( void );
void test034 ( void );
void test035 ( void );
void test036 ( void );
void test0365 ( void );
void test037 ( void );
void test038 ( void );
void test039 ( void );
void test0395 ( void );

void test040 ( void );
void test041 ( void );
void test042 ( void );
void test0425 ( void );
void test043 ( void );
void test044 ( void );
void test0445 ( void );
void test045 ( void );
void test046 ( void );
void test0465 ( void );
void test047 ( void );
void test048 ( void );
void test049 ( void );

void test050 ( void );
void test051 ( void );
void test05125 ( void );
void test0515 ( void );
void test0517 ( void );
void test0519 ( void );
void test052 ( void );
void test053 ( void );
void test054 ( void );
void test055 ( void );
void test056 ( void );
void test057 ( void );
void test058 ( void );
void test059 ( void );

void test060 ( void );
void test061 ( void );
void test062 ( void );
void test063 ( void );
void test064 ( void );
void test065 ( void );
void test066 ( void );
void test0665 ( void );
void test067 ( void );
void test068 ( void );
void test0685 ( void );
void test069 ( void );

void test070 ( void );
void test071 ( void );
void test072 ( void );
void test073 ( void );
void test074 ( void );
void test075 ( void );
void test0755 ( void );
void test0756 ( void );
void test076 ( void );
void test077 ( void );
void test078 ( void );
void test079 ( void );

void test080 ( void );
void test081 ( void );
void test082 ( void );
void test083 ( void );
void test0835 ( void );
void test084 ( void );
void test0843 ( void );
void test0845 ( void );
void test085 ( void );
void test0855 ( void );
void test086 ( void );
void test087 ( void );
void test088 ( void );
void test089 ( void );

void test090 ( void );
void test091 ( void );
void test092 ( void );
void test093 ( void );
void test094 ( void );
void test0945 ( void );
void test095 ( void );
void test096 ( void );
void test097 ( void );
void test0972 ( void );
void test0973 ( void );
void test0974 ( void );
void test0975 ( void );
void test098 ( void );
void test099 ( void );
void test0995 ( void );

void test100 ( void );
void test101 ( void );
void test1015 ( void );
void test1016 ( void );
void test102 ( void );
void test103 ( void );
void test1035 ( void );
void test1037 ( void );
void test104 ( void );
void test105 ( void );
void test106 ( void );
void test107 ( void );
void test108 ( void );
void test10875 ( void );
void test109 ( void );

void test110 ( void );
void test1105 ( void );
void test111 ( void );
void test112 ( void );
void test113 ( void );
void test1135 ( void );
void test114 ( void );
void test115 ( void );
void test116 ( void );
void test117 ( void );
void test118 ( void );
void test1185 ( void );
void test119 ( void );

void test120 ( void );
void test121 ( void );
void test122 ( void );
void test123 ( void );
void test124 ( void );
void test125 ( void );
void test1255 ( void );
void test126 ( void );
void test127 ( void );
void test1275 ( void );
void test128 ( void );
void test1283 ( void );
void test1285 ( void );
void test129 ( void );

void test131 ( void );
void test132 ( void );
void test1325 ( void );
void test130 ( void );
void test133 ( void );
void test134 ( void );
void test135 ( void );
void test136 ( void );
void test137 ( void );
void test138 ( void );
void test139 ( void );

void test140 ( void );
void test141 ( void );
void test1415 ( void );
void test142 ( void );
void test143 ( void );
void test144 ( void );
void test1445 ( void );
void test1447 ( void );
void test145 ( void );
void test146 ( void );
void test1465 ( void );
void test147 ( void );
void test148 ( void );
void test149 ( void );

void test150 ( void );
void test151 ( void );
void test152 ( void );
void test153 ( void );
void test154 ( void );
void test1545 ( void );
void test155 ( void );
void test156 ( void );
void test157 ( void );
void test1575 ( void );
void test158 ( void );
void test159 ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for TEST_VALUES_PRB.

  Discussion:

    TEST_VALUES_PRB calls the TEST_VALUE routines.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "TEST_VALUES_PRB:\n" );
  printf ( "  C version,\n" );
  printf ( "  Test the TEST_VALUES library.\n" );

  test001 ( );
  test002 ( );
  test003 ( );
  test0035 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test008 ( );
  test009 ( );
  test0093 ( );
  test0095 ( );

  test010 ( );
  test011 ( );
  test0114 ( );
  test01145 ( );
  test0115 ( );
  test01155 ( );
  test0116 ( );
  test012 ( );
  test0123 ( );
  test0127 ( );
  test0128 ( );
  test013 ( );
  test0134 ( );
  test0135 ( );
  test014 ( );
  test015 ( );
  test016 ( );
  test017 ( );
  test018 ( );
  test0185 ( );
  test019 ( );
  test0195 ( );

  test020 ( );
  test0205 ( );
  test021 ( );
  test022 ( );
  test023 ( );
  test024 ( );
  test025 ( );
  test026 ( );
  test0265 ( );
  test027 ( );
  test028 ( );
  test029 ( );

  test030 ( );
  test0305 ( );
  test031 ( );
  test032 ( );
  test033 ( );
  test034 ( );
  test035 ( );
  test036 ( );
  test0365 ( );
  test037 ( );
  test038 ( );
  test039 ( );
  test0395 ( );

  test040 ( );
  test041 ( );
  test042 ( );
  test0425 ( );
  test043 ( );
  test044 ( );
  test0445 ( );
  test045 ( );
  test046 ( );
  test0465 ( );
  test047 ( );
  test048 ( );
  test049 ( );

  test050 ( );
  test051 ( );
  test05125 ( );
  test0515 ( );
  test0517 ( );
  test0519 ( );
  test052 ( );
  test053 ( );
  test054 ( );
  test055 ( );
  test056 ( );
  test057 ( );
  test058 ( );
  test059 ( );

  test060 ( );
  test061 ( );
  test062 ( );
  test063 ( );
  test064 ( );
  test065 ( );
  test066 ( );
  test0665 ( );
  test067 ( );
  test068 ( );
  test0685 ( );
  test069 ( );

  test070 ( );
  test071 ( );
  test072 ( );
  test073 ( );
  test074 ( );
  test075 ( );
  test0755 ( );
  test0756 ( );
  test076 ( );
  test077 ( );
  test078 ( );
  test079 ( );

  test080 ( );
  test081 ( );
  test082 ( );
  test083 ( );
  test0835 ( );
  test084 ( );
  test0843 ( );
  test0845 ( );
  test085 ( );
  test0855 ( );
  test086 ( );
  test087 ( );
  test088 ( );
  test089 ( );

  test090 ( );
  test091 ( );
  test092 ( );
  test093 ( );
  test094 ( );
  test0945 ( );
  test095 ( );
  test096 ( );
  test097 ( );
  test0972 ( );
  test0973 ( );
  test0974 ( );
  test0975 ( );
  test098 ( );
  test099 ( );
  test0995 ( );

  test100 ( );
  test101 ( );
  test1015 ( );
  test1016 ( );
  test102 ( );
  test103 ( );
  test1035 ( );
  test104 ( );
  test1037 ( );
  test105 ( );
  test106 ( );
  test107 ( );
  test108 ( );
  test10875 ( );
  test109 ( );

  test110 ( );
  test1105 ( );
  test111 ( );
  test112 ( );
  test113 ( );
  test1135 ( );
  test114 ( );
  test115 ( );
  test116 ( );
  test117 ( );
  test118 ( );
  test1185 ( );
  test119 ( );

  test120 ( );
  test121 ( );
  test122 ( );
  test123 ( );
  test124 ( );
  test125 ( );
  test1255 ( );
  test126 ( );
  test127 ( );
  test1275 ( );
  test128 ( );
  test1283 ( );
  test1285 ( );
  test129 ( );

  test131 ( );
  test132 ( );
  test1325 ( );
  test130 ( );
  test133 ( );
  test134 ( );
  test135 ( );
  test136 ( );
  test137 ( );
  test138 ( );
  test139 ( );

  test140 ( );
  test141 ( );
  test1415 ( );
  test142 ( );
  test143 ( );
  test144 ( );
  test1445 ( );
  test1447 ( );
  test145 ( );
  test146 ( );
  test1465 ( );
  test147 ( );
  test148 ( );
  test149 ( );

  test150 ( );
  test151 ( );
  test152 ( );
  test153 ( );
  test154 ( );
  test1545 ( );
  test155 ( );
  test156 ( );
  test157 ( );
  test1575 ( );
  test158 ( );
  test159 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "TEST_VALUES_PRB:\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test001 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST001 tests ABRAM0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST001:\n" );
  printf ( "  ABRAM0_VALUES stores values of \n" );
  printf ( "  the Abramowitz function of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                   ABRAM0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    abram0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test002 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST002 tests ABRAM1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST002:\n" );
  printf ( "  ABRAM1_VALUES stores values of \n" );
  printf ( "  the Abramowitz function of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                   ABRAM1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    abram1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test003 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST003 tests ABRAM2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST003:\n" );
  printf ( "  ABRAM2_VALUES stores values of \n" );
  printf ( "  the Abramowitz function of order 2.\n" );
  printf ( "\n" );
  printf ( "                X                   ABRAM3(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    abram2_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0035 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0035 tests AGM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 February 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;

  printf ( "\n" );
  printf ( "TEST0035:\n" );
  printf ( "  AGM_VALUES stores values of \n" );
  printf ( "  the arithmetic geometric mean function.\n" );
  printf ( "\n" );
  printf ( "           A          B              AGM(A,B)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    agm_values ( &n_data, &a, &b, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %14.6f  %14.6f  %24.16e\n", a, b, fx );
  }
  return;
}
/******************************************************************************/

void test004 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST004 tests AIRY_AI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double ai;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST004:\n" );
  printf ( "  AIRY_AI_VALUES stores values of \n" );
  printf ( "  the Airy functions Ai(X).\n" );
  printf ( "\n" );
  printf ( "                X                     Ai(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_values ( &n_data, &x, &ai );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, ai );
  }
  return;
}
/******************************************************************************/

void test005 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST005 tests AIRY_AI_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST005:\n" );
  printf ( "  AIRY_AI_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Airy Ai function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test006 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST006 tests AIRY_AI_PRIME_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double aip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST006:\n" );
  printf ( "  AIRY_AI_PRIME_VALUES stores values of \n" );
  printf ( "  the derivative of the Airy function Ai'(X).\n" );
  printf ( "\n" );
  printf ( "                X                    Ai'\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_ai_prime_values ( &n_data, &x, &aip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, aip );
  }
  return;
}
/******************************************************************************/

void test007 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST007 tests AIRY_BI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double bi;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST007:\n" );
  printf ( "  AIRY_BI_VALUES stores values of \n" );
  printf ( "  the Airy function Bi.\n" );
  printf ( "\n" );
  printf ( "                X                     Bi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_values ( &n_data, &x, &bi );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bi );
  }
  return;
}
/******************************************************************************/

void test008 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST008 tests AIRY_BI_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST008:\n" );
  printf ( "  AIRY_BI_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Airy Bi function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test009 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST009 tests AIRY_BI_PRIME_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST009:\n" );
  printf ( "  AIRY_BI_PRIME_VALUES stores values of \n" );
  printf ( "  the derivative of Airy function Bi'(X).\n" );
  printf ( "\n" );
  printf ( "                X                     Bi'\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_bi_prime_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void test0093 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0093 tests AIRY_CAI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double complex cai;
  int n_data;
  double complex x;

  printf ( "\n" );
  printf ( "TEST0093:\n" );
  printf ( "  AIRY_CAI_VALUES stores values of \n" );
  printf ( "  the Airy functions Ai(X) for complex argument.\n" );
  printf ( "\n" );
  printf ( "                X                     Ai\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_cai_values ( &n_data, &x, &cai );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %14.6f  %14.6f  %14.6f  %14.6f\n", x, cai );
  }
  return;
}
/******************************************************************************/

void test0095 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0095 tests AIRY_CBI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double complex cbi;
  int n_data;
  double complex x;

  printf ( "\n" );
  printf ( "TEST0095:\n" );
  printf ( "  AIRY_CBI_VALUES stores values of \n" );
  printf ( "  the Airy functions Bi(X) for complex argument.\n" );
  printf ( "\n" );
  printf ( "                X                     Bi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_cbi_values ( &n_data, &x, &cbi );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %14.6f  %14.6f  %14.6f  %14.6f\n", x, cbi );
  }
  return;
}
/******************************************************************************/

void test010 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST010 tests AIRY_GI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST010:\n" );
  printf ( "  AIRY_GI_VALUES stores values of \n" );
  printf ( "  the modified Airy function Gi(X).\n" );
  printf ( "\n" );
  printf ( "                X                     Gi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_gi_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void test011 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST011 tests AIRY_HI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST011:\n" );
  printf ( "  AIRY_HI_VALUES stores values of \n" );
  printf ( "  the modified Airy function Hi(X).\n" );
  printf ( "\n" );
  printf ( "                X                     Hi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    airy_hi_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void test0114 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0114 tests ARCCOS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0114:\n" );
  printf ( "  ARCCOS_VALUES stores values of the arc cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arccos_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test01145 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01145 tests ARCCOSH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST01145:\n" );
  printf ( "  ARCCOSH_VALUES stores values of the hyperbolic arc cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arccosh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0115 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0115 tests ARCSIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0115:\n" );
  printf ( "  ARCSIN_VALUES stores values of the arc sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arcsin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test01155 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01155 tests ARCSINH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST01155:\n" );
  printf ( "  ARCSINH_VALUES stores values of the hyperbolic arc sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arcsinh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0116 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0116 tests ARCTAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0116:\n" );
  printf ( "  ARCTAN_VALUES stores values of the arc tangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctan_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test012 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST012 tests ARCTAN_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST012:\n" );
  printf ( "  ARCTAN_INT_VALUES stores values of \n" );
  printf ( "  the arctangent integral.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctan_int_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void test0123 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST01235 tests ARCTANH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0123:\n" );
  printf ( "  ARCTANH_VALUES stores values of the hyperbolic arc tangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    arctanh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0127 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0127 tests BEI0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0127:\n" );
  printf ( "  BEI0_VALUES stores values of \n" );
  printf ( "  the Kelvin function BEI of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bei0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0128 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0128 tests BEI1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0128:\n" );
  printf ( "  BEI1_VALUES stores values of \n" );
  printf ( "  the Kelvin function BEI of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bei1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test013 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST013 tests BELL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST013:\n" );
  printf ( "  BELL_VALUES returns values of \n" );
  printf ( "  the Bell numbers.\n" );
  printf ( "\n" );
  printf ( "     N        BELL(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %10d\n", n, c );
  }
  return;
}
/******************************************************************************/

void test0134 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0134 tests BER0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0134:\n" );
  printf ( "  BER0_VALUES stores values of \n" );
  printf ( "  the Kelvin function BER of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ber0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0135 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0135 tests BER1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0135:\n" );
  printf ( "  BER1_VALUES stores values of \n" );
  printf ( "  the Kelvin function BER of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ber1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test014 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST014 tests BERNOULLI_NUMBER_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST014:\n" );
  printf ( "  BERNOULLI_NUMBER_VALUES returns values of \n" );
  printf ( "  the Bernoulli numbers.\n" );
  printf ( "\n" );
  printf ( "     N              B(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12g\n", n, c );
  }
  return;
}
/******************************************************************************/

void test015 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST015 tests BERNOULLI_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double b;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST015:\n" );
  printf ( "  BERNOULLI_POLY_VALUES returns values of \n" );
  printf ( "  the Bernoulli Polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     X      BERNOULLI(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_poly_values ( &n_data, &n, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, b );
  }
  return;
}
/******************************************************************************/

void test016 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST016 tests BERNSTEIN_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double b;
  int k;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST016:\n" );
  printf ( "  BERNSTEIN_POLY_VALUES returns values of \n" );
  printf ( "  the Bernstein Polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     K       X      BERNSTEIN(N,K)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bernstein_poly_values ( &n_data, &n, &k, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12e  %12e\n", n, k, x, b );
  }
  return;
}
/******************************************************************************/

void test017 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST017 tests BESSEL_I0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST017:\n" );
  printf ( "  BESSEL_I0_VALUES stores values of \n" );
  printf ( "  the Bessel I0 function.\n" );
  printf ( "\n" );
  printf ( "      X         I0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test018 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST018 tests BESSEL_I0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST018:\n" );
  printf ( "  BESSEL_I0_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Bessel I0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_int_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void test0185 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0185 tests BESSEL_I0_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0185:\n" );
  printf ( "  BESSEL_I0_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel i0 function.\n" );
  printf ( "\n" );
  printf ( "      X            i0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_i0_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test019 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST019 tests BESSEL_I1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST019:\n" );
  printf ( "  BESSEL_I1_VALUES stores values of \n" );
  printf ( "  the Bessel I1 function.\n" );
  printf ( "\n" );
  printf ( "      X         I1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_i1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0195 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0195 tests BESSEL_I1_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    06 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0195:\n" );
  printf ( "  BESSEL_I1_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel i1 function.\n" );
  printf ( "\n" );
  printf ( "      X            i1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_i1_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test020 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST020 tests BESSEL_IN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST020:\n" );
  printf ( "  BESSEL_IN_VALUES stores values of \n" );
  printf ( "  the Bessel In function.\n" );
  printf ( "\n" );
  printf ( "      N     X         IN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_in_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0205 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0205 tests BESSEL_IX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double nu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0205:\n" );
  printf ( "  BESSEL_IX_VALUES stores values of \n" );
  printf ( "  the Bessel In function for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      NU      X         IN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_ix_values ( &n_data, &nu, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12g  %12e  %12e\n", nu, x, fx );
  }
  return;
}
/******************************************************************************/

void test021 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST021 tests BESSEL_J0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST021:\n" );
  printf ( "  BESSEL_J0_VALUES stores values of \n" );
  printf ( "  the Bessel J0 function.\n" );
  printf ( "\n" );
  printf ( "      X         J0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test022 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST022 tests BESSEL_J0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST022:\n" );
  printf ( "  BESSEL_J0_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Bessel J0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j0_int_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void test023 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST023 tests BESSEL_J0_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST023:\n" );
  printf ( "  BESSEL_J0_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel j0 function.\n" );
  printf ( "\n" );
  printf ( "      X            j0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_j0_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test024 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST024 tests BESSEL_J1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST024:\n" );
  printf ( "  BESSEL_J1_VALUES stores values of \n" );
  printf ( "  the Bessel J1 function.\n" );
  printf ( "\n" );
  printf ( "      X         J1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_j1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test025 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST025 tests BESSEL_J1_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST025:\n" );
  printf ( "  BESSEL_J1_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel j1 function.\n" );
  printf ( "\n" );
  printf ( "      X            j1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_j1_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test026 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST026 tests BESSEL_JN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST026:\n" );
  printf ( "  BESSEL_JN_VALUES stores values of \n" );
  printf ( "  the Bessel Jn function.\n" );
  printf ( "\n" );
  printf ( "      N     X         JN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_jn_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0265 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0265 tests BESSEL_JX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  double fx;
  double nu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0265:\n" );
  printf ( "  BESSEL_JX_VALUES stores values of \n" );
  printf ( "  the Bessel Jn function for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      NU        X         JN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_jx_values ( &n_data, &nu, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12g  %12e  %12e\n", nu, x, fx );
  }
  return;
}
/******************************************************************************/

void test027 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST027 tests BESSEL_K0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST027:\n" );
  printf ( "  BESSEL_K0_VALUES stores values of \n" );
  printf ( "  the Bessel K0 function.\n" );
  printf ( "\n" );
  printf ( "      X         K0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_k0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test028 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST028 tests BESSEL_K0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double bip;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST028:\n" );
  printf ( "  BESSEL_K0_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Bessel K0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_k0_int_values ( &n_data, &x, &bip );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, bip );
  }
  return;
}
/******************************************************************************/

void test029 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST029 tests BESSEL_K1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST029:\n" );
  printf ( "  BESSEL_K1_VALUES stores values of \n" );
  printf ( "  the Bessel K1 function.\n" );
  printf ( "\n" );
  printf ( "      X         K1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_k1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test030 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST030 tests BESSEL_KN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST030:\n" );
  printf ( "  BESSEL_KN_VALUES stores values of \n" );
  printf ( "  the Bessel Kn function.\n" );
  printf ( "\n" );
  printf ( "      N     X         KN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_kn_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0305 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0305 tests BESSEL_KX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double nu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0305:\n" );
  printf ( "  BESSEL_KX_VALUES stores values of \n" );
  printf ( "  the Bessel Kn function for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      NU      X         KN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_kx_values ( &n_data, &nu, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12g  %12e  %12e\n", nu, x, fx );
  }
  return;
}
/******************************************************************************/

void test031 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST031 tests BESSEL_Y0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST031:\n" );
  printf ( "  BESSEL_Y0_VALUES stores values of \n" );
  printf ( "  the Bessel Y0 function.\n" );
  printf ( "\n" );
  printf ( "      X         Y0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_y0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test032 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST032 tests BESSEL_Y0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST032:\n" );
  printf ( "  BESSEL_Y0_INT_VALUES stores values of \n" );
  printf ( "  the integral of the Bessel Y0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_y0_int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test033 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST033 tests BESSEL_Y0_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST033:\n" );
  printf ( "  BESSEL_Y0_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel y0 function.\n" );
  printf ( "\n" );
  printf ( "                X                      y0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_y0_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test034 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST034 tests BESSEL_Y1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST034:\n" );
  printf ( "  BESSEL_Y1_VALUES stores values of \n" );
  printf ( "  the Bessel Y1 function.\n" );
  printf ( "\n" );
  printf ( "                X                   Y1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_y1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test035 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST035 tests BESSEL_Y1_SPHERICAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST035:\n" );
  printf ( "  BESSEL_Y1_SPHERICAL_VALUES stores values of\n" );
  printf ( "  the spherical Bessel y1 function.\n" );
  printf ( "\n" );
  printf ( "                X                      y1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   bessel_y1_spherical_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test036 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST036 tests BESSEL_YN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST036:\n" );
  printf ( "  BESSEL_YN_VALUES stores values of \n" );
  printf ( "  the Bessel Yn function.\n" );
  printf ( "\n" );
  printf ( "      N     X         YN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_yn_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0365 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0365 tests BESSEL_YX_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 April 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double nu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0365:\n" );
  printf ( "  BESSEL_YX_VALUES stores values of \n" );
  printf ( "  the Bessel Yn function for NONINTEGER order.\n" );
  printf ( "\n" );
  printf ( "      NU    X         YN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bessel_yx_values ( &n_data, &nu, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e  %12e\n", nu, x, fx );
  }
  return;
}
/******************************************************************************/

void test037 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST037 tests BETA_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST037:\n" );
  printf ( "  BETA_CDF_VALUES stores values of\n" );
  printf ( "  the Beta CDF.\n" );
  printf ( "\n" );
  printf ( "      A            B            X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void test038 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST038 tests BETA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST038:\n" );
  printf ( "  BETA_INC_VALUES stores values of\n" );
  printf ( "  the incomplete Beta function.\n" );
  printf ( "\n" );
  printf ( "      A            B            X            BETA_INC(A,B)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void test039 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST039 tests BETA_LOG_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fxy;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST039:\n" );
  printf ( "  BETA_LOG_VALUES stores values of\n" );
  printf ( "  the logarithm of the Beta function.\n" );
  printf ( "\n" );
  printf ( "      X              Y         BETA_LOG(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_log_values ( &n_data, &x, &y, &fxy );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", x, y, fxy );
  }
  return;
}
/******************************************************************************/

void test0395 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0395 tests BETA_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2008

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0395:\n" );
  printf ( "  BETA_NONCENTRAL_CDF_VALUES stores values of\n" );
  printf ( "  the noncentral Beta CDF.\n" );
  printf ( "\n" );
  printf ( "      A            B       LAMBDA             X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %12f  %24.16e\n", a, b, lambda, x, fx );
  }
  return;
}
/******************************************************************************/

void test040 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST040 tests BETA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fxy;
  int n_data;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST040:\n" );
  printf ( "  BETA_VALUES stores values of\n" );
  printf ( "  the Beta function.\n" );
  printf ( "\n" );
  printf ( "      X              Y         BETA(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( &n_data, &x, &y, &fxy );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", x,  y, fxy );
  }
  return;
}
/******************************************************************************/

void test041 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST041 tests BINOMIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  int c;
  int n_data;

  printf ( "\n" );
  printf ( "TEST041:\n" );
  printf ( "  BINOMIAL_VALUES returns values of\n" );
  printf ( "  the binomial numbers.\n" );
  printf ( "\n" );
  printf ( "     A     B        C(A,B)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    binomial_values ( &n_data, &a, &b, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", a, b, c );
  }
  return;
}
/******************************************************************************/

void test042 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST042 tests BINOMIAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  int a;
  double b;
  double fx;
  int n_data;
  int x;

  printf ( "\n" );
  printf ( "TEST042:\n" );
  printf ( "  BINOMIAL_CDF_VALUES returns values of \n" );
  printf ( "  the Binomial Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     A      B        X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8f  %4d  %24.16e\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void test0425 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0425 tests BETA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 April 2012

  Author:

    John Burkardt
*/
{
  double fxy;
  int n_data;
  double r;
  double x;
  double y;

  printf ( "\n" );
  printf ( "TEST0425:\n" );
  printf ( "  BIVARIATE_NORMAL_CDF_VALUES stores values of\n" );
  printf ( "  the bivariate normal CDF.\n" );
  printf ( "\n" );
  printf ( "      X            Y            R            F(R)(X,Y)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    bivariate_normal_cdf_values ( &n_data, &x, &y, &r, &fxy );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", x, y, r, fxy );
  }
  return;
}
/******************************************************************************/

void test043 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST043 tests CATALAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST043:\n" );
  printf ( "  CATALAN_VALUES returns values of \n" );
  printf ( "  the Catalan numbers.\n" );
  printf ( "\n" );
  printf ( "     N        C(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %10d\n", n, c );
  }
  return;
}
/******************************************************************************/

void test044 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST044 tests CAUCHY_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TEST044:\n" );
  printf ( "  CAUCHY_CDF_VALUES returns values of \n" );
  printf ( "  the Cauchy Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Mu      Sigma        X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cauchy_cdf_values ( &n_data, &mu, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", mu, sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void test0445 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0445 tests CBRT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0445:\n" );
  printf ( "  CBRT_VALUES stores values of the cube root function.\n" );
  printf ( "\n" );
  printf ( "      X            CBRT(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cbrt_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test045 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST045 tests CHEBY_T_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 March 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST045:\n" );
  printf ( "  CHEBY_T_POLY_VALUES returns values of\n" );
  printf ( "  the Chebyshev T polynomials.\n" );
  printf ( "\n" );
  printf ( "     N       X      T(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_t_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test046 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST046 tests CHEBY_U_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST046:\n" );
  printf ( "  CHEBY_U_POLY_VALUES returns values of\n" );
  printf ( "  the Chebyshev U polynomials.\n" );
  printf ( "\n" );
  printf ( "     N       X      U(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cheby_u_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0465 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0465 tests CHI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0465:\n" );
  printf ( "  CHI_VALUES stores values of\n" );
  printf ( "  the Hyperbolic Cosine Integral function CHI(X).\n" );
  printf ( "\n" );
  printf ( "      X            CHI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test047 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST047 tests CHI_SQUARE_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST047:\n" );
  printf ( "  CHI_SQUARE_CDF_VALUES returns values of \n" );
  printf ( "  the Chi-Squared Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     N       X    CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8f  %12e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void test048 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST048 tests CHI_SQUARE_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int df;
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST048:\n" );
  printf ( "  CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of\n" );
  printf ( "  the noncentral Chi-Squared Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "      X      LAMBDA     DF     CDF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    chi_square_noncentral_cdf_values ( &n_data, &df, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %10f  %8f  %4d  %12e\n", x, lambda, df, fx );
  }
  return;
}
/******************************************************************************/

void test049 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST049 tests CI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST049:\n" );
  printf ( "  CI_VALUES stores values of\n" );
  printf ( "  the Cosine Integral function CI(X).\n" );
  printf ( "\n" );
  printf ( "      X            CI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ci_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test050 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST050 tests CIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST050:\n" );
  printf ( "  CIN_VALUES stores values of\n" );
  printf ( "  the Cosine Integral function CIN(X).\n" );
  printf ( "\n" );
  printf ( "      X            CIN(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test051 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST051 tests BESSEL_Y0_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST051:\n" );
  printf ( "  CLAUSEN_VALUES stores values of \n" );
  printf ( "  Clausen's integral function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    clausen_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test05125 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST05125 tests CLEBSCH_GORDAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double j1;
  double j2;
  double j3;
  double m1;
  double m2;
  double m3;
  int n_data;

  printf ( "\n" );
  printf ( "TEST05125:\n" );
  printf ( "  CLEBSCH_GORDAN_VALUES returns values of\n" );
  printf ( "  the Clebsch Gordan coefficient.\n" );
  printf ( "\n" );
  printf ( "      J1      J2      J3      M1      M2      M3        CG\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    clebsch_gordan_values ( &n_data, &j1, &j2, &j3, &m1, &m2, &m3, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f  %24.16e\n",
    j1, j2, j3, m1, m2, m3, fx );
  }
  return;
}
/******************************************************************************/

void test0515 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0515 tests COLLATZ_COUNT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 March 2006

  Author:

    John Burkardt
*/
{
  int count;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST0515:\n" );
  printf ( "  COLLATZ_COUNT_VALUES returns values of\n" );
  printf ( "  the length of the Collatz sequence that\n" );
  printf ( "  starts at N.\n" );
  printf ( "\n" );
  printf ( "         N      COLLATZ_COUNT(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    collatz_count_values ( &n_data, &n, &count );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %12d\n", n, count );
  }

  return;
}
/******************************************************************************/

void test0517 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0517 tests COS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0517:\n" );
  printf ( "   COS_VALUES stores values of the cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cos_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0519 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0519 tests COSH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0519:\n" );
  printf ( "   COSH_VALUES stores values of the hyperbolic cosine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cosh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test052 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST052 tests CP_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double cp;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "TEST052:\n" );
  printf ( "  CP_VALUES stores values of\n" );
  printf ( "  the specific heat CP\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            CP(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    cp_values ( &n_data, &tc, &p, &cp );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, cp );
  }
  return;
}
/******************************************************************************/

void test053 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST053 tests DAWSON_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST053:\n" );
  printf ( "  DAWSON_VALUES stores values of\n" );
  printf ( "  Dawson's integral function.\n" );
  printf ( "\n" );
  printf ( "      X          DAWSON(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    dawson_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test054 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST054 tests DEBYE1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST054:\n" );
  printf ( "  DEBYE1_VALUES stores values of \n" );
  printf ( "  the Debye function of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    debye1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test055 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST055 tests DEBYE2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST055:\n" );
  printf ( "  DEBYE2_VALUES stores values of \n" );
  printf ( "  the Debye function of order 2.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    debye2_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test056 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST056 tests DEBYE3_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST056:\n" );
  printf ( "  DEBYE3_VALUES stores values of \n" );
  printf ( "  the Debye function of order 3.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    debye3_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test057 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST057 tests DEBYE4_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST057:\n" );
  printf ( "  DEBYE4_VALUES stores values of \n" );
  printf ( "  the Debye function of order 4.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    debye4_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test058 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST058 tests DIELECTRIC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double eps;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "TEST058:\n" );
  printf ( "  DIELECTRIC_VALUES stores values of\n" );
  printf ( "  the dielectric function.\n" );
  printf ( "\n" );
  printf ( "      T           P            EPS(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    dielectric_values ( &n_data, &tc, &p, &eps );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, eps );
  }
  return;
}
/******************************************************************************/

void test059 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST059 tests DILOGARITHM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST059:\n" );
  printf ( "  DILOGARITHM_VALUES stores values of\n" );
  printf ( "  the dilogarithm function.\n" );
  printf ( "\n" );
  printf ( "      X          DILOGARITHM(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    dilogarithm_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test060 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST060 tests E1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST060:\n" );
  printf ( "  E1_VALUES stores values of\n" );
  printf ( "  the exponential integral function E1(X).\n" );
  printf ( "\n" );
  printf ( "      X          E1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    e1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test061 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST061 tests EI_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST061:\n" );
  printf ( "  EI_VALUES stores values of\n" );
  printf ( "  the exponential integral function EI(X).\n" );
  printf ( "\n" );
  printf ( "      X          EI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ei_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test062 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST062 tests ELLIPTIC_EA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST062:\n" );
  printf ( "  ELLIPTIC_EA_VALUES stores values of\n" );
  printf ( "  the complete elliptic integral of the second\n" );
  printf ( "  kind, with parameter angle ALPHA in degrees.\n" );
  printf ( "\n" );
  printf ( "    ALPHA        EA(ALPHA)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    elliptic_ea_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test063 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST063 tests ELLIPTIC_EM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST063:\n" );
  printf ( "  ELLIPTIC_EM_VALUES stores values of\n" );
  printf ( "  the complete elliptic integral of the second\n" );
  printf ( "  kind, with parameter modulus M.\n" );
  printf ( "\n" );
  printf ( "      M            EM(M)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    elliptic_em_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test064 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST064 tests ELLIPTIC_KA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST064:\n" );
  printf ( "  ELLIPTIC_KA_VALUES stores values of\n" );
  printf ( "  the complete elliptic integral of the first\n" );
  printf ( "  kind, with parameter angle ALPHA in degrees.\n" );
  printf ( "\n" );
  printf ( "    ALPHA        KA(ALPHA)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    elliptic_ka_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test065 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST065 tests ELLIPTIC_KM_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST065:\n" );
  printf ( "  ELLIPTIC_KM_VALUES stores values of\n" );
  printf ( "  the complete elliptic integral of the first\n" );
  printf ( "  kind, with parameter modulus M.\n" );
  printf ( "\n" );
  printf ( "      M            KM(M)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    elliptic_km_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test066 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST066 tests ERF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST066:\n" );
  printf ( "  ERF_VALUES stores values of\n" );
  printf ( "  the error function ERF(X).\n" );
  printf ( "\n" );
  printf ( "      X          ERF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0665 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0665 tests ERFC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0665:\n" );
  printf ( "  ERFC_VALUES stores values of\n" );
  printf ( "  the complementary error function ERFC(X).\n" );
  printf ( "\n" );
  printf ( "      X          ERFC(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    erfc_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %12e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test067 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST067 tests EULER_NUMBER_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int c;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST067:\n" );
  printf ( "  EULER_NUMBER_VALUES returns values of\n" );
  printf ( "  the Euler numbers.\n" );
  printf ( "\n" );
  printf ( "     N        EULER_NUMBER(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %10d\n", n, c );
  }
  return;
}
/******************************************************************************/

void test068 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST068 tests EULER_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST068:\n" );
  printf ( "  EULER_POLY_VALUES returns values of\n" );
  printf ( "  the Euler numbers.\n" );
  printf ( "\n" );
  printf ( "     N     X       EULER_POLY(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    euler_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %12e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0685 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0685 tests EXP_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0685:\n" );
  printf ( "   EXP_VALUES stores values of the exponential function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    exp_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test069 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST069 tests EXP3_INT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST069:\n" );
  printf ( "  EXP3_INT_VALUES stores values of \n" );
  printf ( "  the exponential integral function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    exp3_int_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test070 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST070 tests EXPONENTIAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST070:\n" );
  printf ( "  EXPONENTIAL_CDF_VALUES stores values of \n" );
  printf ( "  the exponential CDF.\n" );
  printf ( "\n" );
  printf ( "      LAMBDA          X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    exponential_cdf_values ( &n_data, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.8f  %24.8f  %24.16e\n", lambda, x, fx );
  }
  return;
}
/******************************************************************************/

void test071 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST071 tests EXTREME_VALUES_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST071:\n" );
  printf ( "  EXTREME_VALUES_CDF_VALUES stores values of \n" );
  printf ( "  the extreme values CDF.\n" );
  printf ( "\n" );
  printf ( "        Alpha    Beta        X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    extreme_values_cdf_values ( &n_data, &alpha, &beta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", alpha, beta, x, fx );
  }
  return;
}
/******************************************************************************/

void test072 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST072 tests F_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST072:\n" );
  printf ( "   F_CDF_VALUES stores values of\n" );
  printf ( "   the F cumulative density function.\n" );
  printf ( "\n" );
  printf ( "     A       B            X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12f  %12f\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void test073 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST073 tests F_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int a;
  int b;
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST073:\n" );
  printf ( "   F_NONCENTRAL_CDF_VALUES stores values of\n" );
  printf ( "   the F cumulative density function.\n" );
  printf ( "\n" );
  printf ( "     A       B            LAMBDA    X            CDF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    f_noncentral_cdf_values ( &n_data, &a, &b, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %8f  %12f  %12e\n", a, b, lambda, x, fx );
  }
  return;
}
/******************************************************************************/

void test074 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST074 tests FRESNEL_COS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST074:\n" );
  printf ( "   FRESNEL_COS_VALUES stores values of\n" );
  printf ( "   the Fresnel cosine integral C(X).\n" );
  printf ( "\n" );
  printf ( "      X           C(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    fresnel_cos_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test075 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST075 tests FRESNEL_SIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST075:\n" );
  printf ( "   FRESNEL_SIN_VALUES stores values of\n" );
  printf ( "   the Fresnel sine integral S(X).\n" );
  printf ( "\n" );
  printf ( "      X           S(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    fresnel_sin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0755 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0755 tests FROBENIUS_NUMBER_ORDER2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  int c1;
  int c2;
  int f;
  int n_data;

  printf ( "\n" );
  printf ( "TEST0755:\n" );
  printf ( "  FROBENIUS_NUMBER_ORDER2_VALUES returns values of \n" );
  printf ( "  the Frobenius number of order 2.\n" );
  printf ( "\n" );
  printf ( "         C1        C2          F(C1,C2)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    frobenius_number_order2_values ( &n_data, &c1, &c2, &f );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %8d  %8d  %8d\n", c1, c2, f );
  }

  return;
}
/******************************************************************************/

void test0756 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0756 tests FROBENIUS_NUMBER_ORDER_VALUES, FROBENIUS_NUMBER_DATA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 November 2007

  Author:

    John Burkardt
*/
{
  int *c;
  int f;
  int i;
  int n_data;
  int order;

  printf ( "\n" );
  printf ( "TEST0756:\n" );
  printf ( "  FROBENIUS_NUMBER_ORDER_VALUES returns the order for\n" );
  printf ( "  a Frobenius problem;\n" );
  printf ( "  FROBENIUS_NUMBER_DATA_VALUES returns the corresponding\n" );
  printf ( "  coin denominations.\n" );

  n_data = 0;

  for ( ; ; )
  {
    frobenius_number_order_values ( &n_data, &order );

    if ( n_data == 0 )
    {
      break;
    }

    c = malloc ( order * sizeof ( int ) );

    frobenius_number_data_values ( &n_data, order, c, &f );

    printf ( "\n" );
    printf ( "  Order = %d\n", order );
    for ( i = 0; i < order; i++ )
    {
      printf ( "  %8d", c[i] );
    }
    printf ( "\n" );
    printf ( "  Frobenius number = %d\n", f );

    free ( c );
  }
  return;
}
/******************************************************************************/

void test076 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST076 tests GAMMA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST076:\n" );
  printf ( "   GAMMA_VALUES stores values of the Gamma function.\n" );
  printf ( "\n" );
  printf ( "      X            GAMMA(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test077 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST077 tests GAMMA_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( " TEST077:\n" );
  printf ( "   GAMMA_CDF_VALUES stores values of\n" );
  printf ( "   the Gamma CDF.\n" );
  printf ( "\n" );
  printf ( "      M    Sigma      X            CDF((X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_cdf_values ( &n_data, &mu, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", mu, sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void test078 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST078 tests GAMMA_INC_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST078:\n" );
  printf ( "   GAMMA_INC_VALUES stores values of\n" );
  printf ( "   the incomplete Gamma function.\n" );
  printf ( "\n" );
  printf ( "      A            X            GAMMA_INC(A)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}

/******************************************************************************/

void test079 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST079 tests GAMMA_LOG_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST079:\n" );
  printf ( "   GAMMA_LOG_VALUES stores values of\n" );
  printf ( "   the logarithm of the Gamma function.\n" );
  printf ( "\n" );
  printf ( "      X            GAMMA_LOG(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test080 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST080 tests GEGENBAUER_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST080:\n" );
  printf ( "  GEGENBAUER_POLY_VALUES returns values of\n" );
  printf ( "  the Gegenbauer polynomials.\n" );
  printf ( "\n" );
  printf ( "       N       A       X       G(N,A)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gegenbauer_poly_values ( &n_data, &n, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12f  %12f  %24.16e\n", n, a, x, fx );
  }

  return;
}
/******************************************************************************/

void test081 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST081 tests GEOMETRIC_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double cdf;
  int n_data;
  double p;
  int x;

  printf ( "\n" );
  printf ( " TEST081:\n" );
  printf ( "   GEOMETRIC_CDF_VALUES stores values of\n" );
  printf ( "   the Geometric Probability Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "      X      P       CDF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    geometric_cdf_values ( &n_data, &x, &p, &cdf );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12f  %24.16e\n", x, p, cdf );
  }
  return;
}
/******************************************************************************/

void test082 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST082 tests GOODWIN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST082:\n" );
  printf ( "  GOODWIN_VALUES stores values of \n" );
  printf ( "  the Goodwin function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    goodwin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test083 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST083 tests GUD_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST083:\n" );
  printf ( "   GUD_VALUES stores values of\n" );
  printf ( "   the Gudermannian function.\n" );
  printf ( "\n" );
  printf ( "      X            GUD(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    gud_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0835 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0835 tests HERMITE_FUNCTION_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST0835\n" );
  printf ( "   HERMITE_FUNCTION_VALUES stores values of\n" );
  printf ( "   the Hermite function.\n" );
  printf ( "\n" );
  printf ( "     N      X            Hf(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hermite_function_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test084 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST084 tests HERMITE_POLY_PHYS_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST084\n" );
  printf ( "   HERMITE_POLY_PHYS_VALUES stores values of\n" );
  printf ( "   the physicist's Hermite polynomials.\n" );
  printf ( "\n" );
  printf ( "     N      X            H(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_phys_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0843 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0843 tests HERMITE_POLY_PROB_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2012

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST0843\n" );
  printf ( "   HERMITE_POLY_PROB_VALUES stores values of\n" );
  printf ( "   the probabilist's Hermite polynomials.\n" );
  printf ( "\n" );
  printf ( "     N      X            He(N,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_prob_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0845 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0845 tests HYPER_2F1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 September 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double c;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( " TEST0845:\n" );
  printf ( "   HYPER_2F1_VALUES stores values of\n" );
  printf ( "   the hypergeometric function 2F1.\n" );
  printf ( "\n" );
  printf ( "      A      B     C      X   Hyper_2F1(A,B,C,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hyper_2f1_values ( &n_data, &a, &b, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %8f  %24.16e\n", a, b, c, x, fx ); 
  }
  return;
}
/******************************************************************************/

void test085 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST085 tests HYPERGEOMETRIC_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 January 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  int pop;
  int sam;
  int succ;
  int x;

  printf ( "\n" );
  printf ( " TEST085:\n" );
  printf ( "   HYPERGEOMETRIC_CDF_VALUES stores values of\n" );
  printf ( "   the Hypergeometric CDF.\n" );
  printf ( "\n" );
  printf ( "     SAM    SUC   POP     X   HyperCDF(S,S,P)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_cdf_values ( &n_data, &sam, &succ, &pop, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %8d  %8d  %8d  %24.16e\n", sam, succ, pop, x, fx );
  }
  return;
}
/******************************************************************************/

void test0855 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0855 tests HYPERGEOMETRIC_PDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 January 2008

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  int pop;
  int sam;
  int succ;
  int x;

  printf ( "\n" );
  printf ( " TEST0855:\n" );
  printf ( "   HYPERGEOMETRIC_PDF_VALUES stores values of\n" );
  printf ( "   the Hypergeometric PDF.\n" );
  printf ( "\n" );
  printf ( "     SAM    SUC   POP     X   HyperPDF(S,S,P)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    hypergeometric_pdf_values ( &n_data, &sam, &succ, &pop, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %8d  %8d  %8d  %24.16e\n", sam, succ, pop, x, fx );
  }
  return;
}
/******************************************************************************/

void test086 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST086 tests FACTORIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 March 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( " TEST086:\n" );
  printf ( "   FACTORIAL_VALUES return;s values of\n" );
  printf ( "   the factorial function.\n" );
  printf ( "\n" );
  printf ( "      N         Factorial(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test087 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST087 tests FACTORIAL2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( " TEST087:\n" );
  printf ( "   FACTORIAL2_VALUES return;s values of\n" );
  printf ( "   the double factorial function.\n" );
  printf ( "\n" );
  printf ( "      N         DoubleFactorial(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    factorial2_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test088 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST088 tests FACTORIAL_RISING_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 January 2012

  Author:

    John Burkardt
*/
{
  int fmn;
  int m;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST088:\n" );
  printf ( "  FACTORIAL_RISING_VALUES returns some exact values\n" );
  printf ( "  of the Pochhammer symbol:\n" );
  printf ( "\n" );
  printf ( "     M     N      Factorial_rising(M,N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    factorial_rising_values ( &n_data, &m, &n, &fmn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", m, n, fmn );
  }
  return;
}
/******************************************************************************/

void test089 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST089 tests I0ML0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST089:\n" );
  printf ( "  I0ML0_VALUES stores values of \n" );
  printf ( "  the I0-L0 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i0ml0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test090 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST090 tests I1ML1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST090:\n" );
  printf ( "  I1ML1_VALUES stores values of \n" );
  printf ( "  the I1-L1 function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    i1ml1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test091 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST091 tests JACOBI_CN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST091:\n" );
  printf ( "  JACOBI_CN_VALUES returns values of \n" );
  printf ( "  the Jacobi elliptic CN function.\n" );
  printf ( "\n" );
  printf ( "      A         X       CN(A,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_cn_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void test092 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST092 tests JACOBI_DN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST092:\n" );
  printf ( "  JACOBI_DN_VALUES returns values of \n" );
  printf ( "  the Jacobi elliptic DN function.\n" );
  printf ( "\n" );
  printf ( "      A         X       DN(A,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_dn_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void test093 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST093 tests JACOBI_POLY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 April 2012

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST093:\n" );
  printf ( "  JACOBI_POLY_VALUES returns values of\n" );
  printf ( "  the Jacobi polynomial.\n" );
  printf ( "\n" );
  printf ( "       N         A         B      X       J(N,A,B)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_poly_values ( &n_data, &n, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8f  %8f  %12f  %24.16e\n", n, a, b, x, fx );
  }

  return;
}
/******************************************************************************/

void test094 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST094 tests JACOBI_SN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST094:\n" );
  printf ( "  JACOBI_SN_VALUES returns values of \n" );
  printf ( "  the Jacobi elliptic SN function.\n" );
  printf ( "\n" );
  printf ( "      A         X       SN(A,X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jacobi_sn_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void test0945 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0945 tests JED_CE_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 May 2006

  Author:

    John Burkardt
*/
{
  int d;
  double f;
  double jed;
  int n_data;
  int m;
  int y;

  printf ( "\n" );
  printf ( "TEST0945:\n" );
  printf ( "  JED_CE_VALUES returns:\n" );
  printf ( "  JED, a Julian Ephemeris Date, and\n" );
  printf ( "  YMDF, the corresponding year, month, day, fraction.\n" );
  printf ( "\n" );
  printf ( "        JED          Y   M   D    F\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jed_ce_values ( &n_data, &jed, &y, &m, &d, &f );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %6d  %2d  %2d  %6f\n", jed, y, m, d, f );
  }
  return;
}
/******************************************************************************/

void test095 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST095 tests JED_MJD_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double jed;
  int n_data;
  double mjd;

  printf ( "\n" );
  printf ( "TEST095:\n" );
  printf ( "  JED_MJD_VALUES returns:\n" );
  printf ( "  JED, a Julian Ephemeris Date, and\n" );
  printf ( "  MJD, the corresponding Modified Julian Day count.\n" );
  printf ( "\n" );
  printf ( "   JED      MJD\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jed_mjd_values ( &n_data, &jed, &mjd );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", jed, mjd );
  }
  return;
}
/******************************************************************************/

void test096 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST096 tests JED_RD_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double jed;
  int n_data;
  double rd;

  printf ( "\n" );
  printf ( "TEST096:\n" );
  printf ( "  JED_RD_VALUES returns:\n" );
  printf ( "  JED, a Julian Ephemeris Date, and\n" );
  printf ( "  RD, the corresponding Reingold Dershowitz Day count.\n" );
  printf ( "\n" );
  printf ( "   JED      RD\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jed_rd_values ( &n_data, &jed, &rd );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", jed, rd );
  }
  return;
}
/******************************************************************************/

void test097 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST097 tests JED_WEEKDAY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double jed;
  int n_data;
  int weekday;

  printf ( "\n" );
  printf ( "TEST097:\n" );
  printf ( "  JED_WEEKDAY_VALUES returns Julian Ephemeris Dates \n" );
  printf ( "  (JED) and the corresponding weekday\n" );
  printf ( "\n" );
  printf ( "   JED      #  Weekday\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    jed_weekday_values ( &n_data, &jed, &weekday );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %d\n", jed, weekday );
  }
  return;
}
/******************************************************************************/

void test0972 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0972 tests KEI0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0972:\n" );
  printf ( "  KEI0_VALUES stores values of \n" );
  printf ( "  the Kelvin function KEI of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    kei0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0973 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0973 tests KEI1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0973:\n" );
  printf ( "  KEI1_VALUES stores values of \n" );
  printf ( "  the Kelvin function KEI of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    kei1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0974 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0974 tests KER0_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0974:\n" );
  printf ( "  KER0_VALUES stores values of \n" );
  printf ( "  the Kelvin function KER of order 0.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ker0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test0975 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0975 tests KER1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2006

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0975:\n" );
  printf ( "  KER1_VALUES stores values of \n" );
  printf ( "  the Kelvin function KER of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    ker1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test098 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST098 tests LAGUERRE_ASSOCIATED_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST098:\n" );
  printf ( "  LAGUERRE_ASSOCIATED_VALUES stores values of\n" );
  printf ( "  the associated Laguerre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     M    X             L(N,M)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laguerre_associated_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12f  %24.16e\n", n, m, x, fx );
  }
  return;
}
/******************************************************************************/

void test099 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST099 tests LAGUERRE_POLYNOMIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST099:\n" );
  printf ( "  LAGUERRE_POLYNOMIAL_VALUES stores values of \n" );
  printf ( "  the Laguerre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     X            L(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laguerre_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24.16e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test0995 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST0995 tests LAMBERT_W_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST0995:\n" );
  printf ( "  LAMBERT_W_VALUES stores values of \n" );
  printf ( "  the Lambert W function.\n" );
  printf ( "\n" );
  printf ( "                X                     W(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lambert_w_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test100 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST100 tests LAPLACE_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double beta;
  double fx;
  double mu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST100:\n" );
  printf ( "  LAPLACE_CDF_VALUES returns values of \n" );
  printf ( "  the Laplace Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Mu      Beta         X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    laplace_cdf_values ( &n_data, &mu, &beta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", mu, beta, x, fx );
  }
  return;
}
/******************************************************************************/

void test101 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST101 tests LEGENDRE_ASSOCIATED_VALUES.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST101:\n" );
  printf ( "  LEGENDRE_ASSOCIATED_VALUES stores values of\n" );
  printf ( "  the associated Legendre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     M    X             P(N,M)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12e  %24.16e\n", n, m, x, fx );
  }
  return;
}
/******************************************************************************/

void test1015 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1015 tests LEGENDRE_ASSOCIATED_NORMALIZED_VALUES.

  Modified:

    01 September 2010

  Author:

    John Burkardt
*/
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1015:\n" );
  printf ( "  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES stores values of\n" );
  printf ( "  the normalized associated Legendre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N     M    X             P(N,M)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12e  %24.16e\n", n, m, x, fx );
  }
  return;
}
/******************************************************************************/

void test1016 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1016 tests LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES.

  Modified:

    12 March 2012

  Author:

    John Burkardt
*/
{
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1016:\n" );
  printf ( "  LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES stores values of\n" );
  printf ( "  the associated Legendre polynomials normalized for the unit sphere.\n" );
  printf ( "\n" );
  printf ( "     N     M    X             P(N,M)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_sphere_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12e  %24.16e\n", n, m, x, fx );
  }
  return;
}
/******************************************************************************/

void test102 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST102 tests LEGENDRE_POLY_VALUES.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST102:\n" );
  printf ( "  LEGENDRE_POLY_VALUES stores values of \n" );
  printf ( "  the Legendre polynomials.\n" );
  printf ( "\n" );
  printf ( "     N    X             P(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24.16e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test103 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST103 tests LEGENDRE_FUNCTION_Q_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST103:\n" );
  printf ( "  LEGENDRE_FUNCTION_Q_VALUES stores values of\n" );
  printf ( "  the Legendre Q function.\n" );
  printf ( "\n" );
  printf ( "     N    X             Q(N)(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    legendre_function_q_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12e  %24.16e\n", n, x, fx );
  }
  return;
}
/******************************************************************************/

void test1035 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1035 tests LERCH_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  int s;
  double z;

  printf ( "\n" );
  printf ( "TEST1035:\n" );
  printf ( "  LERCH_VALUES returns values of\n" );
  printf ( "  the Lerch transcendent function.\n" );
  printf ( "\n" );
  printf ( "      Z      S      A      Fx\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lerch_values ( &n_data, &z, &s, &a, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %6d  %12f  %24.16e\n", z, s, a, fx );
  }
  return;
}
/******************************************************************************/

void test104 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST104 tests LOBACHEVSKY_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST104:\n" );
  printf ( "  LOBACHEVSKY_VALUES stores values of \n" );
  printf ( "  the Lobachevsky function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    lobachevsky_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test1037 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1037 tests LOG_VALUES.

  Modified:

    19 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1037:\n" );
  printf ( "   LOG_VALUES stores values of the natural logarithm function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test105 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST105 tests LOG_NORMAL_CDF_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TEST105:\n" );
  printf ( "  LOG_NORMAL_CDF_VALUES returns values of \n" );
  printf ( "  the Log Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Mu      Sigma        X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    log_normal_cdf_values ( &n_data, &mu, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", mu, sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void test106 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST106 tests LOG_SERIES_CDF_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double t;

  printf ( "\n" );
  printf ( "TEST106:\n" );
  printf ( "  LOG_SERIES_CDF_VALUES returns values of \n" );
  printf ( "  the Log Series Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     T      N   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    log_series_cdf_values ( &n_data, &t, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %6d  %24.16e\n", t, n, fx );
  }
  return;
}
/******************************************************************************/

void test107 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST107 tests LOGARITHMIC_INTEGRAL_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST107:\n" );
  printf ( "  LOGARITHMIC_INTEGAL_VALUES stores values of\n" );
  printf ( "  the logarithmic integral function.\n" );
  printf ( "\n" );
  printf ( "      X            LI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    logarithmic_integral_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test108 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST108 tests LOGISTIC_CDF_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double beta;
  double fx;
  double mu;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST108:\n" );
  printf ( "  LOGISTIC_CDF_VALUES returns values of \n" );
  printf ( "  the Logistic Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Mu      Beta         X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    logistic_cdf_values ( &n_data, &mu, &beta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", mu, beta, x, fx );
  }
  return;
}
/******************************************************************************/

void test10875 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST10875 tests MERTENS_VALUES.

  Modified:

    17 October 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST10875:\n" );
  printf ( "  MERTENS_VALUES returns values of\n" );
  printf ( "  the Mertens function.\n" );
  printf ( "\n" );
  printf ( "     N         MERTENS(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    mertens_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    } 
    printf ( "  %12d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test109 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST109 tests MOEBIUS_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST109:\n" );
  printf ( "  MOEBIUS_VALUES returns values of\n" );
  printf ( "  the Moebius function.\n" );
  printf ( "\n" );
  printf ( "     N         MU(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    moebius_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    } 
    printf ( "  %12d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test110 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST110 tests NEGATIVE_BINOMIAL_CDF_VALUES.

  Modified:

    13 June 2007

  Author:

    John Burkardt
*/
{
  double cdf;
  int f;
  int n_data;
  double p;
  int s;

  printf ( "\n" );
  printf ( "TEST110:\n" );
  printf ( "  NEGATIVE_BINOMIAL_CDF_VALUES stores values of\n" );
  printf ( "  the Negative Binomial Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     F     S         P         CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    negative_binomial_cdf_values ( &n_data, &f, &s, &p, &cdf );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12f  %24.16e\n", f, s, p, cdf );
  }
  return;
}
/******************************************************************************/

void test1105 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1105 demonstrates NINE_J_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double j1;
  double j2;
  double j3;
  double j4;
  double j5;
  double j6;
  double j7;
  double j8;
  double j9;
  int n_data;

  printf ( "\n" );
  printf ( "TEST1105:\n" );
  printf ( "  GSL_SF_COUPLING_9J returns values of\n" );
  printf ( "  the Wigner 9J coefficient.\n" );
  printf ( "\n" );
  printf ( "      J1      J2      J3      J4      J5      J6" );
  printf ( "      J7      J8      J9        NINE_J\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    nine_j_values ( &n_data, &j1, &j2, &j3, &j4, &j5, &j6, &j7, &j8, &j9, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f  %6f  %6f  %6f  %24.16e\n",
    j1, j2, j3, j4, j5, j6, j7, j8, j9, fx );
  }
  return;
}
/******************************************************************************/

void test111 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST111 tests NORMAL_CDF_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double mu;
  int n_data;
  double sigma;
  double x;

  printf ( "\n" );
  printf ( "TEST111:\n" );
  printf ( "  NORMAL_CDF_VALUES stores values of\n" );
  printf ( "  the Normal Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "            X                   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_cdf_values ( &n_data, &mu, &sigma, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", mu, sigma, x, fx );
  }
  return;
}
/******************************************************************************/

void test112 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST112 tests NORMAL_01_CDF_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST112:\n" );
  printf ( "  NORMAL_01_CDF_VALUES stores values of\n" );
  printf ( "  the Normal 01 Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "            X                   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test113 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST113 tests OMEGA_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST113:\n" );
  printf ( "  OMEGA_VALUES returns values of\n" );
  printf ( "  the Omega function.\n" );
  printf ( "\n" );
  printf ( "     N           OMEGA(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    omega_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test1135 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST1135 tests OWEN_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double a;
  double h;
  int n_data;
  double t;

  printf ( "\n" );
  printf ( "TEST1135\n" );
  printf ( "  OWEN_VALUES stores values of\n" );
  printf ( "  Owen's T function.\n" );
  printf ( "\n" );
  printf ( "          H            A            T\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( &n_data, &h, &a, &t );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f %12f  %24.16e\n", h, a, t );
  }
  return;
}
/******************************************************************************/

void test114 ( void )
 
/******************************************************************************/
/*
  Purpose:

    TEST114 tests PARTITION_COUNT_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST114:\n" );
  printf ( "  PARTITION_COUNT_VALUES returns values of \n" );
  printf ( "  the integer partition count function.\n" );
  printf ( "\n" );
  printf ( "     N         P(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    partition_count_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test115 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST115 tests PARTITION_DISTINCT_COUNT_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST115:\n" );
  printf ( "  PARTITION_DISTINCT_COUNT_VALUES returns values of \n" );
  printf ( "  the integer distinct partition count function.\n" );
  printf ( "\n" );
  printf ( "     N         Q(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    partition_distinct_count_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test116 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST116 tests PHI_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST116:\n" );
  printf ( "  PHI_VALUES returns values of\n" );
  printf ( "  the PHI function.\n" );
  printf ( "\n" );
  printf ( "     N         PHI(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    phi_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test117 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST117 tests PI_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST117:\n" );
  printf ( "  PI_VALUES returns values of\n" );
  printf ( "  the PI function.\n" );
  printf ( "\n" );
  printf ( "     N         PI(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    pi_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test118 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST118 tests POISSON_CDF_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double a;
  double fx;
  int n_data;
  int x;

  printf ( "\n" );
  printf ( "TEST118:\n" );
  printf ( "  POISSON_CDF_VALUES returns values of\n" );
  printf ( "  the Poisson Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "      A     X       CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %4d  %24.16e\n", a, x, fx );
  }
  return;
}
/******************************************************************************/

void test1185 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1185 tests POLYLOGARITHM_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n;
  int n_data;
  double z;

  printf ( "\n" );
  printf ( "TEST1185:\n" );
  printf ( "  POLYLOGARITHM_VALUES returns values of \n" );
  printf ( "  the polylogarithm function.\n" );
  printf ( "\n" );
  printf ( "     N      Z          Fx\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    polylogarithm_values ( &n_data, &n, &z, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24.16e  %24.16e\n", n, z, fx );
  }
  return;
}
/******************************************************************************/

void test119 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST119 tests PRANDTL_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double p;
  double pr;
  double tc;

  printf ( "\n" );
  printf ( "TEST119:\n" );
  printf ( "  PRANDTL_VALUES stores values of\n" );
  printf ( "  the Prandtl number of water\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            Pr(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    prandtl_values ( &n_data, &tc, &p, &pr );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, pr );
  }
  return;
}
/******************************************************************************/

void test120 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST120 tests PRIME_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int n;
  int n_data;
  int p;

  printf ( "\n" );
  printf ( "TEST120:\n" );
  printf ( "  PRIME_VALUES returns values of\n" );
  printf ( "  the prime function.\n" );
  printf ( "\n" );
  printf ( "           N          P[N]\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    prime_values ( &n_data, &n, &p );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %12d\n", n, p );
  }

  return;
}
/******************************************************************************/

void test121 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST121 tests PSAT_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double psat;
  double tc;

  printf ( "\n" );
  printf ( "TEST121:\n" );
  printf ( "  PSAT_VALUES stores values of\n" );
  printf ( "  the saturation pressure of water\n" );
  printf ( "  as a function of temperature.\n" );
  printf ( "\n" );
  printf ( "      T            PSAT(T)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psat_values ( &n_data, &tc, &psat );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", tc, psat );
  }
  return;
}
/******************************************************************************/

void test122 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST122 tests PSI_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST122\n" );
  printf ( "  PSI_VALUES stores values of\n" );
  printf ( "  the PSI function.\n" );
  printf ( "\n" );
  printf ( "      X            PSI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test123 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST123 tests R8_FACTORIAL_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST123:\n" );
  printf ( "  R8_FACTORIAL_VALUES stores values of\n" );
  printf ( "  the factorial function (using double arithmetic).\n" );
  printf ( "\n" );
  printf ( "      N       Factorial(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12d  %24.16e\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test124 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST124 tests R8_FACTORIAL_LOG_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST124:\n" );
  printf ( "  R8_FACTORIAL_LOG_VALUES stores values of\n" );
  printf ( "  the logarithm of the factorial function\n" );
  printf ( "  (using real arithmetic).\n" );
  printf ( "\n" );
  printf ( "      N       Log(Factorial(N))\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_log_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24.16e\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test125 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST125 tests SECVIR_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double tc;
  double vir;

  printf ( "\n" );
  printf ( "TEST125:\n" );
  printf ( "  SECVIR_VALUES stores values of\n" );
  printf ( "  the second virial coefficient of water\n" );
  printf ( "  as a function of temperature.\n" );
  printf ( "\n" );
  printf ( "      T            VIR(T)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   secvir_values ( &n_data, &tc, &vir );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", tc, vir );
  }
  return;
}
/******************************************************************************/

void test1255 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1255 tests SHI_VALUES.

  Modified:

    11 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1255:\n" );
  printf ( "  SHI_VALUES stores values of\n" );
  printf ( "  the hyperbolic sine integral function.\n" );
  printf ( "\n" );
  printf ( "      X            SHI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   shi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test126 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST126 tests SI_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST126:\n" );
  printf ( "  SI_VALUES stores values of\n" );
  printf ( "  the sine integral function.\n" );
  printf ( "\n" );
  printf ( "      X            SI(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   si_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test127 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST127 tests SIGMA_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST127:\n" );
  printf ( "  SIGMA_VALUES returns values of\n" );
  printf ( "  the SIGMA function.\n" );
  printf ( "\n" );
  printf ( "       N         SIGMA(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sigma_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test1275 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1275 tests SIN_VALUES.

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1275:\n" );
  printf ( "   SIN_VALUES stores values of the sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sin_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test128 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST128 tests SIN_POWER_INT_VALUES.

  Modified:

    09 February 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST128:\n" );
  printf ( "  SIN_POWER_INT_VALUES returns values of\n" );
  printf ( "  the integral of the N-th power of the sine function.\n" );
  printf ( "\n" );
  printf ( "         A         B       N        FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sin_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %6d  %24.16e\n", a, b, n, fx );
  }
  return;
}
/******************************************************************************/

void test1283 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1283 tests SINH_VALUES.

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1283:\n" );
  printf ( "   SINH_VALUES stores values of the hyperbolic sine function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sinh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test1285 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1285 tests SIX_J_VALUES.

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double j1;
  double j2;
  double j3;
  double j4;
  double j5;
  double j6;
  int n_data;

  printf ( "\n" );
  printf ( "TEST1285:\n" );
  printf ( "  SIX_J_VALUES returns values of \n" );
  printf ( "  the Wigner 6J coefficient.\n" );
  printf ( "\n" );
  printf ( "      J1      J2      J3      J4      J5      J6        SIX_J\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    six_j_values ( &n_data, &j1, &j2, &j3, &j4, &j5, &j6, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f  %24.16e\n", 
    j1, j2, j3, j4, j5, j6, fx );
  }

  return;
}
/******************************************************************************/

void test129 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST129 tests SOUND_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double c;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "TEST129:\n" );
  printf ( "  SOUND_VALUES stores values of\n" );
  printf ( "  the spead of sound in water\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            C(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sound_values ( &n_data, &tc, &p, &c );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, c );
  }
  return;
}
/******************************************************************************/

void test131 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST131 tests SPHERE_UNIT_AREA_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{  
  double fx;
  int n_data;
  int n;

  printf ( "\n" );
  printf ( "TEST131:\n" );
  printf ( "  SPHERE_UNIT_AREA_VALUES stores values of\n" );
  printf ( "  the area of the unit sphere in various dimensions.\n" );
  printf ( "\n" );
  printf ( "      N           AREA\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_area_values ( &n_data, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24.16e\n", n, fx );
  }
  return;
}
/******************************************************************************/

void test132 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST132 tests SPHERE_UNIT_VOLUME_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{  
  double fx;
  int n_data;
  int n;

  printf ( "\n" );
  printf ( "TEST132:\n" );
  printf ( "  SPHERE_UNIT_VOLUME_VALUES stores values of\n" );
  printf ( "  the volume of the unit sphere in various dimensions.\n" );
  printf ( "\n" );
  printf ( "      N           VOLUME\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    sphere_unit_volume_values ( &n_data, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24.16e\n", n, fx );
  }
  return;
}
/******************************************************************************/

void test1325 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1325 tests SPHERICAL_HARMONIC_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{

  int l;
  int m;
  int n_data;
  double phi;
  double theta;
  double yi;
  double yr;

  printf ( "\n" );
  printf ( "TEST1325:\n" );
  printf ( "  SPHERICAL_HARMONIC_VALUES stores values of\n" );
  printf ( "  the spherical harmonic function.\n" );
  printf ( "\n" );
  printf ( "   L   M    THETA       PHI           Yr                    Yi\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    spherical_harmonic_values ( &n_data, &l, &m, &theta, &phi, &yr, &yi );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %2d  %2d  %8f  %8f  %24.16e  %24.16e\n",
    l, m, theta, phi, yr, yi );
  }
  return;
}
/******************************************************************************/

void test130 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST130 tests SQRT_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST130:\n" );
  printf ( "  SQRT_VALUES returns some exact values.\n" );
  printf ( "\n" );
  printf ( "     X       Fx\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   sqrt_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test133 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST133 tests STIRLING1_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int n_data;
  int s1;

  printf ( "\n" );
  printf ( "TEST133:\n" );
  printf ( "  STIRLING1_VALUES returns values of\n" );
  printf ( "  the Stirling numbers of the first kind.\n" );
  printf ( "\n" );
  printf ( "     N     N        S1\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    stirling1_values ( &n_data, &n, &m, &s1 );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", n, m, s1 );
  }
  return;
}
/******************************************************************************/

void test134 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST134 tests STIRLING2_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int m;
  int n;
  int n_data;
  int s2;

  printf ( "\n" );
  printf ( "TEST134:\n" );
  printf ( "  STIRLING2_VALUES returns values of\n" );
  printf ( "  the Stirling numbers of the second kind.\n" );
  printf ( "\n" );
  printf ( "     N     N        S2\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    stirling1_values ( &n_data, &n, &m, &s2 );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %6d  %12d\n", n, m, s2 );
  }
  return;
}
/******************************************************************************/

void test135 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST135 tests STROMGEN_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST135:\n" );
  printf ( "  STROMGEN_VALUES stores values of \n" );
  printf ( "  the Stromgen function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    stromgen_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test136 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST136 tests STRUVE_H0_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{  
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST136:\n" );
  printf ( "  STRUVE_H0_VALUES stores values of\n" );
  printf ( "  the Struve H0 function.\n" );
  printf ( "\n" );
  printf ( "      X            H0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    struve_h0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test137 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST137 tests STRUVE_H1_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST137:\n" );
  printf ( "  STRUVE_H1_VALUES stores values of\n" );
  printf ( "  the Struve H1 function.\n" );
  printf ( "\n" );
  printf ( "      X            H1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   struve_h1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test138 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST138 tests STRUVE_L0_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{  
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST138:\n" );
  printf ( "  STRUVE_L0_VALUES stores values of\n" );
  printf ( "  the Struve L0 function.\n" );
  printf ( "\n" );
  printf ( "      X            L0(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    struve_l0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test139 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST139 tests STRUVE_L1_VALUES.

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST139:\n" );
  printf ( "  STRUVE_L1_VALUES stores values of\n" );
  printf ( "  the Struve L1 function.\n" );
  printf ( "\n" );
  printf ( "      X            L1(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   struve_l1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test140 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST140 tests STUDENT_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 November 2005

  Author:

    John Burkardt
*/
{
  double c;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST140:\n" );
  printf ( "  STUDENT_CDF_VALUES returns values of\n" );
  printf ( "  the Student T Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "      C     X       CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   student_cdf_values ( &n_data, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  16f  %16f  %24.16e\n", c, x, fx );
  }
  return;
}
/******************************************************************************/

void test141 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST141 tests STUDENT_NONCENTRAL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int df;
  double fx;
  double lambda;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST141:\n" );
  printf ( "  STUDENT_NONCENTRAL_CDF_VALUES returns values of\n" );
  printf ( "  the noncentral Student T Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "    DF     LAMBDA        X        CDF\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    student_noncentral_cdf_values ( &n_data, &df, &lambda, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %8f  %8f  %24.16e\n", df, lambda, x, fx );
  }
  return;
}
/******************************************************************************/

void test1415 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1415 tests SUBFACTORIAL_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 March 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( " TEST1415:\n" );
  printf ( "   SUBFACTORIAL_VALUES returns values of\n" );
  printf ( "   the subfactorial function.\n" );
  printf ( "\n" );
  printf ( "      N       Subfactorial[N]\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    subfactorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test142 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST142 tests SURTEN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double sigma;
  double tc;

  printf ( "\n" );
  printf ( "TEST142:\n" );
  printf ( "  SURTEN_VALUES stores values of\n" );
  printf ( "  the surface tension of water\n" );
  printf ( "  as a function of temperature.\n" );
  printf ( "\n" );
  printf ( "      T            SIGMA(T)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   surten_values ( &n_data, &tc, &sigma );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", tc, sigma );
  }
  return;
}
/******************************************************************************/

void test143 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST143 tests SYNCH1_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST143:\n" );
  printf ( "  SYNCH1_VALUES stores values of \n" );
  printf ( "  the Synchrotron function of order 1.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    synch1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test144 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST144 tests SYNCH2_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST144:\n" );
  printf ( "  SYNCH2_VALUES stores values of \n" );
  printf ( "  the Synchrotron function of order 2.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    synch2_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test1445 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1445 tests TAN_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1445:\n" );
  printf ( "   TAN_VALUES stores values of the tangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tan_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test1447 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1447 tests TANH_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    23 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1447:\n" );
  printf ( "   TANH_VALUES stores values of the hyperbolic tangent function.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tanh_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test145 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST145 tests TAU_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  int fn;
  int n;
  int n_data;

  printf ( "\n" );
  printf ( "TEST145:\n" );
  printf ( "  TAU_VALUES returns values of\n" );
  printf ( "  the TAU function.\n" );
  printf ( "\n" );
  printf ( "     N         TAU(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   tau_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %12d\n", n, fn );
  }
  return;
}
/******************************************************************************/

void test146 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST146 tests THERCON_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    08 February 2007

  Author:

    John Burkardt
*/
{
  double lambda;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "TEST146:\n" );
  printf ( "  THERCON_VALUES stores values of\n" );
  printf ( "  the thermal conductivity of water\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            LAMBDA(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   thercon_values ( &n_data, &tc, &p, &lambda );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, lambda );
  }
  return;
}
/******************************************************************************/

void test1465 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1465 tests THREE_J_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 February 2007

  Author:

    John Burkardt
*/
{
  double fx;
  double j1;
  double j2;
  double j3;
  double m1;
  double m2;
  double m3;
  int n_data;

  printf ( "\n" );
  printf ( "TEST1465:\n" );
  printf ( "  THREE_J_VALUES returns values of\n" );
  printf ( "  the Wigner 3J coefficient.\n" );
  printf ( "\n" );
  printf ( "      J1      J2      J3      M1      M2      M3        THREE_J\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    three_j_values ( &n_data, &j1, &j2, &j3, &m1, &m2, &m3, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6f  %6f  %6f  %6f  %6f  %6f  %24.16e\n",
    j1, j2, j3, m1, m2, m3, fx );
  }
  return;
}
/******************************************************************************/

void test147 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST147 tests TRAN02_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST147:\n" );
  printf ( "  TRAN02_VALUES stores values of \n" );
  printf ( "  the Transport function of order 2.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran02_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test148 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST148 tests TRAN03_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST148:\n" );
  printf ( "  TRAN03_VALUES stores values of \n" );
  printf ( "  the Transport function of order 3.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran03_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test149 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST149 tests TRAN04_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST149:\n" );
  printf ( "  TRAN04_VALUES stores values of \n" );
  printf ( "  the Transport function of order 4.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran04_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test150 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST150 tests TRAN05_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST150:\n" );
  printf ( "  TRAN05_VALUES stores values of \n" );
  printf ( "  the Transport function of order 5.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran05_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test151 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST151 tests TRAN06_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST151:\n" );
  printf ( "  TRAN06_VALUES stores values of \n" );
  printf ( "  the Transport function of order 6.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran06_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test152 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST152 tests TRAN07_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST152:\n" );
  printf ( "  TRAN07_VALUES stores values of \n" );
  printf ( "  the Transport function of order 7.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran07_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test153 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST153 tests TRAN08_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST153:\n" );
  printf ( "  TRAN08_VALUES stores values of \n" );
  printf ( "  the Transport function of order 8.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran08_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test154 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST154 tests TRAN09_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST154:\n" );
  printf ( "  TRAN09_VALUES stores values of \n" );
  printf ( "  the Transport function of order 9.\n" );
  printf ( "\n" );
  printf ( "                X                     FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    tran09_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}
/******************************************************************************/

void test1545 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST1545 tests TRIGAMMA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1545\n" );
  printf ( "  TRIGAMMA_VALUES stores values of\n" );
  printf ( "  the TriGamma function.\n" );
  printf ( "\n" );
  printf ( "      X            FX\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    trigamma_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %24.16e  %24.16e\n", x, fx );
  }
  return;
}

/******************************************************************************/

void test155 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST155 tests TSAT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "TEST155:\n" );
  printf ( "  TSAT_VALUES stores values of\n" );
  printf ( "  the saturation temperature\n" );
  printf ( "  as a function of pressure.\n" );
  printf ( "\n" );
  printf ( "      P           Tsat(P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   tsat_values ( &n_data, &p, &tc );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f\n", p, tc );
  }
  return;
}
/******************************************************************************/

void test156 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST156 tests VAN_DER_CORPUT_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int base;
  int n_data;
  int seed;
  double value;

  printf ( "\n" );
  printf ( "TEST156:\n" );
  printf ( "  VAN_DER_CORPUT_VALUES stores values of\n" );
  printf ( "  the van der Corput sequence in a given base.\n" );
  printf ( "\n" );
  printf ( "      BASE      SEED    VDC(BASE,SEED)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    van_der_corput_values ( &n_data, &base, &seed, &value );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8d  %8d  %14f\n", base, seed, value );
  }

  return;
}
/******************************************************************************/

void test157 ( void )

/******************************************************************************/
/*
  Purpose: 

    TEST157 tests VISCOSITY_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double eta;
  int n_data;
  double p;
  double tc;

  printf ( "\n" );
  printf ( "TEST157:\n" );
  printf ( "  VISCOSITY_VALUES stores values of\n" );
  printf ( "  the viscosity of water\n" );
  printf ( "  as a function of temperature and pressure.\n" );
  printf ( "\n" );
  printf ( "      T            P            ETA(T,P)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
   viscosity_values ( &n_data, &tc, &p, &eta );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f\n", tc, p, eta );
  }
  return;
}
/******************************************************************************/

void test1575 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST1575 tests VON_MISES_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST1575:\n" );
  printf ( "  VON_MISES_CDF_VALUES stores values of\n" );
  printf ( "  the von Mises CDF.\n" );
  printf ( "\n" );
  printf ( "      A            B            X            CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    von_mises_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %12f  %12f  %12f  %24.16e\n", a, b, x, fx );
  }
  return;
}
/******************************************************************************/

void test158 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST158 tests WEIBULL_CDF_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  double alpha;
  double beta;
  double fx;
  int n_data;
  double x;

  printf ( "\n" );
  printf ( "TEST158:\n" );
  printf ( "  WEIBULL_CDF_VALUES returns values of \n" );
  printf ( "  the Weibull Cumulative Density Function.\n" );
  printf ( "\n" );
  printf ( "     Alpha   Beta        X   CDF(X)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    weibull_cdf_values ( &n_data, &alpha, &beta, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %8f  %8f  %8f  %24.16e\n", alpha, beta, x, fx );
  }
  return;
}
/******************************************************************************/

void test159 ( void )

/******************************************************************************/
/*
  Purpose:

    TEST159 tests ZETA_VALUES.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    12 June 2007

  Author:

    John Burkardt
*/
{
  int n;
  int n_data;
  double zeta;

  printf ( "\n" );
  printf ( "TEST159:\n" );
  printf ( "  ZETA_VALUES returns values of \n" );
  printf ( "  the Riemann Zeta function.\n" );
  printf ( "\n" );
  printf ( "     N        ZETA(N)\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    zeta_values ( &n_data, &n, &zeta );

    if ( n_data == 0 )
    {
      break;
    }
    printf ( "  %6d  %24e\n", n, zeta );
  }
  return;
}
