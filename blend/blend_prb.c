# include <stdlib.h>
# include <stdio.h>

# include "blend.h"

int main ( );
void cubic_rs ( double r, double s, int i, double *xi );
void identity_r ( double r, int i, double *xi );
void identity_rs ( double r, double s, int i, double *xi );
void identity_rst ( double r, double s, double t, int i, double *xi );
void quad_rst ( double r, double s, double t, int i, double *xi );
void stretch_r ( double r, int i, double *xi );
void stretch_rs ( double r, double s, int i, double *xi );
void stretch_rst ( double r, double s, double t, int i, double *xi );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for BLEND_PRB.

  Discussion:

    BLEND_PRB tests the BLEND library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  timestamp ( );
  printf ( "\n" );
  printf ( "BLEND_PRB\n" );
  printf ( "  C version\n" );
  printf ( "  Test the BLEND library.\n" );

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "BLEND_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void test01 ( )

/******************************************************************************/
/*
  Purpose:

    TEST01 checks out BLEND_R_0DN on the identity map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double x[1];

  printf ( "\n" );
  printf ( "TEST01\n" );
  printf ( "  Identity test on BLEND_R_0DN.\n" );

  n = 1;

  r = 0.0;
  blend_r_0dn ( r, x, n, identity_r );
  printf ( "  %8g  %8g\n", r, x[0] );

  r = 1.0;
  blend_r_0dn ( r, x, n, identity_r );
  printf ( "  %8g  %8g\n", r, x[0] );

  r = 0.5;
  blend_r_0dn ( r, x, n, identity_r );
  printf ( "  %8g  %8g\n", r, x[0] );

  return;
}
/******************************************************************************/

void test02 ( )

/******************************************************************************/
/*
  Purpose:

    TEST02 checks out BLEND_RS_0DN on the identity map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double x[2];

  printf ( "\n" );
  printf ( "TEST02\n" );
  printf ( "  Identity test on BLEND_RS_0DN.\n" );

  n = 2;

  r = 0.0;
  s = 0.0;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 1.0;
  s = 0.0;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 0.0;
  s = 1.0;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 1.0;
  s = 1.0;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 0.5;
  s = 0.5;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  return;
}
/******************************************************************************/

void test03 ( )

/******************************************************************************/
/*
  Purpose:

    TEST03 checks out BLEND_RS_1DN on the identity map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double x[2];

  printf ( "\n" );
  printf ( "TEST03\n" );
  printf ( "  Identity test on BLEND_RS_1DN.\n" );

  n = 2;

  r = 0.0;
  s = 0.0;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 1.0;
  s = 0.0;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 0.0;
  s = 1.0;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 1.0;
  s = 1.0;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 0.5;
  s = 0.5;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  return;
}
/******************************************************************************/

void test04 ( )

/******************************************************************************/
/*
  Purpose:

    TEST04 checks out BLEND_RST_0DN on the identity map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  printf ( "\n" );
  printf ( "TEST04\n" );
  printf ( "  Identity test on BLEND_RST_0DN.\n" );

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  return;
}
/******************************************************************************/

void test05 ( )

/******************************************************************************/
/*
  Purpose:

    TEST05 checks out BLEND_RST_1DN on the identity map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  printf ( "\n" );
  printf ( "TEST05\n" );
  printf ( "  Identity test on BLEND_RST_1DN.\n" );

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  return;
}
/******************************************************************************/

void test06 ( )

/******************************************************************************/
/*
  Purpose:

    TEST06 checks out BLEND_RST_2DN on the identity map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  printf ( "\n" );
  printf ( "TEST06\n" );
  printf ( "  Identity test on BLEND_RST_2DN.\n" );

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  return;
}
/******************************************************************************/

void test07 ( )

/******************************************************************************/
/*
  Purpose:

    TEST07 checks out BLEND_R_0DN on the stretch map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double x[1];

  printf ( "\n" );
  printf ( "TEST07\n" );
  printf ( "  Stretch test on BLEND_R_0DN.\n" );

  n = 1;

  r = 0.0;
  blend_r_0dn ( r, x, n, stretch_r );
  printf ( "  %8g  %8g\n", r, x[0] );

  r = 1.0;
  blend_r_0dn ( r, x, n, stretch_r );
  printf ( "  %8g  %8g\n", r, x[0] );

  r = 0.5;
  blend_r_0dn ( r, x, n, stretch_r );
  printf ( "  %8g  %8g\n", r, x[0] );

  return;
}
/******************************************************************************/

void test08 ( )

/******************************************************************************/
/*
  Purpose:

    TEST08 checks out BLEND_RS_0DN on the stretch map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double x[2];

  printf ( "\n" );
  printf ( "TEST08\n" );
  printf ( "  Stretch test on BLEND_RS_0DN.\n" );

  n = 2;

  r = 0.0;
  s = 0.0;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 1.0;
  s = 0.0;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 0.0;
  s = 1.0;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 1.0;
  s = 1.0;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 0.5;
  s = 0.5;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  return;
}
/******************************************************************************/

void test09 ( )

/******************************************************************************/
/*
  Purpose:

    TEST09 checks out BLEND_RS_1DN on the stretch map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double x[2];

  printf ( "\n" );
  printf ( "TEST09\n" );
  printf ( "  Stretch test on BLEND_RS_1DN.\n" );

  n = 2;

  r = 0.0;
  s = 0.0;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 1.0;
  s = 0.0;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 0.0;
  s = 1.0;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 1.0;
  s = 1.0;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  r = 0.5;
  s = 0.5;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  printf ( "  %8g  %8g  %8g  %8g\n", r, s, x[0], x[1] );

  return;
}
/******************************************************************************/

void test10 ( )

/******************************************************************************/
/*
  Purpose:

    TEST10 checks out BLEND_RST_0DN on the stretch map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  printf ( "\n" );
  printf ( "TEST10\n" );
  printf ( "  Stretch test on BLEND_RST_0DN.\n" );

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  return;
}
/******************************************************************************/

void test11 ( )

/******************************************************************************/
/*
  Purpose:

    TEST11 checks out BLEND_RST_1DN on the stretch map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  printf ( "\n" );
  printf ( "TEST11\n" );
  printf ( "  Stretch test on BLEND_RST_1DN.\n" );

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  return;
}
/******************************************************************************/

void test12 ( )

/******************************************************************************/
/*
  Purpose:

    TEST12 checks out BLEND_RST_2DN on the stretch map.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  printf ( "\n" );
  printf ( "TEST12\n" );
  printf ( "  Stretch test on BLEND_RST_2DN.\n" );

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  printf ( "  %8g  %8g  %8g  %8g  %8g  %8g\n", r, s, t, x[0], x[1], x[2] );

  return;
}
/******************************************************************************/

void test13 ( )

/******************************************************************************/
/*
  Purpose:

    TEST13 checks out BLEND_I_0D1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int i;
  int m;
  double x[5];

  m = 5;
  x[0] = 100.0;
  x[m-1] = 100.0 + ( m - 1 ) * 5;

  printf ( "\n" );
  printf ( "TEST13\n" );
  printf ( "  BLEND_I_0D1 interpolates data in a vector.\n" );
  printf ( "\n" );
  printf ( "  X[0] = %g\n", x[0] );
  printf ( "  X(%d)= %g\n", m - 1, x[m-1] );
  printf ( "\n" );
  printf ( "  Interpolated values:\n" );
  printf ( "\n" );

  blend_i_0d1 ( x, m );

  for ( i = 0; i < m; i++ )
  {
    printf ( "  %6d  %8g\n", i, x[i] );
  }
  return;
}
/******************************************************************************/

void test14 ( )

/******************************************************************************/
/*
  Purpose:

    TEST14 checks out BLEND_IJ_0D1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int m1 = 5;
  int m2 = 4;
  double r;
  double s;
  double temp;
  double x[20];

  printf ( "\n" );
  printf ( "TEST14\n" );
  printf ( "  BLEND_IJ_0D1 interpolates data in a table,\n" );
  printf ( "  from corner data.\n" );
  printf ( "\n" );
  printf ( "  The table is %d rows by %d columns.\n", m1, m2 );
/*
  Load data in the corners only.
*/
  i = 0;
  j = 0;
  r = ( double ) i / ( double ) ( m1 - 1 );
  s = ( double ) j / ( double ) ( m2 - 1 );
  cubic_rs ( r, s, 1, &temp );
  x[i*m2+j] = temp;

  i = m1 - 1;
  j = 0;
  r = ( double ) i / ( double ) ( m1 - 1 );
  s = ( double ) j / ( double ) ( m2 - 1 );
  cubic_rs ( r, s, 1, &temp );
  x[i*m2+j] = temp;

  i = 0;
  j = m2 - 1;
  r = ( double ) i / ( double ) ( m1 - 1 );
  s = ( double ) j / ( double ) ( m2 - 1 );
  cubic_rs ( r, s, 1, &temp );
  x[i*m2+j] = temp;

  i = m1 - 1;
  j = m2 - 1;
  r = ( double ) i / ( double ) ( m1 - 1 );
  s = ( double ) j / ( double ) ( m2 - 1 );
  cubic_rs ( r, s, 1, &temp );
  x[i*m2+j] = temp;

  blend_ij_0d1 ( x, m1, m2 );

  printf ( "\n" );
  printf ( "  Values interpolated by BLEND_IJ_0D1:\n" );
  printf ( "\n" );

  for ( i = 0; i < m1; i++ )
  {
    printf ( "  %8g  %8g  %8g  %8g\n", x[i*m2], x[i*m2+1], x[i*m2+2], x[i*m2+3] );
  }

  return;
}
/******************************************************************************/

void test15 ( )

/******************************************************************************/
/*
  Purpose:

    TEST15 checks out BLEND_IJ_1D1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int m1 = 5;
  int m2 = 4;
  double r;
  double s;
  double temp;
  double x[20];

  printf ( "\n" );
  printf ( "TEST15\n" );
  printf ( "  BLEND_IJ_1D1 interpolates data in a table,\n" );
  printf ( "  from edge data.\n" );
  printf ( "\n" );
  printf ( "  The table is %d rows by %d columns.\n", m1, m2 );
/*
  Load data in the edges only.
*/
  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );

    j = 0;
    s = ( double ) j / ( double ) ( m2 - 1 );
    cubic_rs ( r, s, 1, &temp );
    x[i*m2+j] = temp;

    j = m2 - 1;
    s = ( double ) j / ( double ) ( m2 - 1 );
    cubic_rs ( r, s, 1, &temp );
    x[i*m2+j] = temp;
  }

  for ( j = 0; j < m2; j++ )
  {
    s = ( double ) j / ( double ) ( m2 - 1 );

    i = 0;
    r = ( double ) i / ( double ) ( m1 - 1 );
    cubic_rs ( r, s, 1, &temp );
    x[i*m2+j] = temp;

    i = m1 - 1;
    r = ( double ) i / ( double ) ( m1 - 1 );
    cubic_rs ( r, s, 1, &temp );
    x[i*m2+j] = temp;
  }

  blend_ij_1d1 ( x, m1, m2 );

  printf ( "\n" );
  printf ( "  Values interpolated by BLEND_IJ_1D1:\n" );
  printf ( "\n" );

  for ( i = 0; i < m1; i++ )
  {
    printf ( "  %8g  %8g  %8g  %8g\n", x[i*m2], x[i*m2+1], x[i*m2+2], x[i*m2+3] );
  }

  return;
}
/******************************************************************************/

void test16 ( )

/******************************************************************************/
/*
  Purpose:

    TEST16 checks out BLEND_IJK_0D1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;
  int m1 = 4;
  int m2 = 3;
  int m3 = 3;
  int num_extreme;
  double r;
  double s;
  double t;
  double temp;
  double x[36];

  printf ( "\n" );
  printf ( "TEST16\n" );
  printf ( "  BLEND_IJK_0D1 interpolates data in a 3D table,\n" );
  printf ( "  from corner data.\n" );
  printf ( "\n" );
  printf ( "  The table is %d rows by %d columns by %d layers.\n", m1, m2, m3 );
/*
  Load data in the faces only.
*/
  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        num_extreme = 0;
        if ( i == 0 || i == m1 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( j == 0 || j == m2 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( k == 0 || k == m3 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( num_extreme >= 3 )
        {
          quad_rst ( r, s, t, 1, &temp );
        }
        else
        {
          temp = 0.0;
        }
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  printf ( "\n" );
  printf ( "  Data given to BLEND_IJK_0D1:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }
  blend_ijk_0d1 ( x, m1, m2, m3 );

  printf ( "\n" );
  printf ( "  Values interpolated by BLEND_IJK_0D1:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }

  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        quad_rst ( r, s, t, 1, &temp );
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  printf ( "\n" );
  printf ( "  Exact data:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }
  return;
}
/******************************************************************************/

void test17 ( )

/******************************************************************************/
/*
  Purpose:

    TEST17 checks out BLEND_IJK_1D1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;
  int m1 = 4;
  int m2 = 3;
  int m3 = 3;
  int num_extreme;
  double r;
  double s;
  double t;
  double temp;
  double x[36];

  printf ( "\n" );
  printf ( "TEST17\n" );
  printf ( "  BLEND_IJK_1D1 interpolates data in a 3D table,\n" );
  printf ( "  from edge data.\n" );
  printf ( "\n" );
  printf ( "  The table is %d rows by %d columns by %d layers.\n", m1, m2, m3 );
/*
  Load data in the faces only.
*/
  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        num_extreme = 0;
        if ( i == 0 || i == m1 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( j == 0 || j == m2 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( k == 0 || k == m3 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( num_extreme >= 2 )
        {
          quad_rst ( r, s, t, 1, &temp );
        }
        else
        {
          temp = 0.0;
        }
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  printf ( "\n" );
  printf ( "  Data given to BLEND_IJK_1D1:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }
  blend_ijk_1d1 ( x, m1, m2, m3 );

  printf ( "\n" );
  printf ( "  Values interpolated by BLEND_IJK_1D1:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }

  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        quad_rst ( r, s, t, 1, &temp );
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  printf ( "\n" );
  printf ( "  Exact data:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }
  return;
}
/******************************************************************************/

void test18 ( )

/******************************************************************************/
/*
  Purpose:

    TEST18 checks out BLEND_IJK_2D1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  int i;
  int j;
  int k;
  int m1 = 4;
  int m2 = 3;
  int m3 = 3;
  double r;
  double s;
  double t;
  double temp;
  double x[36];

  printf ( "\n" );
  printf ( "TEST18\n" );
  printf ( "  BLEND_IJK_2D1 interpolates data in a 3D table,\n" );
  printf ( "  from face data.\n" );
  printf ( "\n" );
  printf ( "  The table is %d rows by %d columns by %d layers.\n", m1, m2, m3 );
/*
  Load data in the faces only.
*/
  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        if ( i == 0 || i == m1 - 1 ||
             j == 0 || j == m2 - 1 ||
             k == 0 || k == m3 - 1 )
        {
          quad_rst ( r, s, t, 1, &temp );
        }
        else
        {
          temp = 0.0;
        }
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  printf ( "\n" );
  printf ( "  Data given to BLEND_IJK_2D1:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }
  blend_ijk_2d1 ( x, m1, m2, m3 );

  printf ( "\n" );
  printf ( "  Values interpolated by BLEND_IJK_2D1:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }

  for ( i = 0; i < m1; i++ )
   {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        quad_rst ( r, s, t, 1, &temp );
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  printf ( "\n" );
  printf ( "  Exact data:\n" );
  printf ( "\n" );

  for ( k = 0; k < m3; k++ )
  {
    printf ( "\n" );
    printf ( "  Layer K = %d\n", k );
    printf ( "\n" );

    for ( i = 0; i < m1; i++ )
    {
      printf ( "  %8g  %8g  %8g\n", x[(i*m3+0)*m2+k], x[(i*m3+1)*m2+k], x[(i*m3+2)*m2+k] );
    }
  }
  return;
}
/******************************************************************************/

void cubic_rs ( double r, double s, int i, double *xi )

/******************************************************************************/
/*
  Purpose:

    CUBIC_RS evaluates a function of R and S used for some tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  *xi = 20.0 * ( r * r * s * s * s );

  return;
}
/******************************************************************************/

void quad_rst ( double r, double s, double t, int i, double *xi )

/******************************************************************************/
/*
  Purpose:

    QUAD_RST evaluates a function of R, S and T used for some tests.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt
*/
{
  *xi = 18.0 * ( r * r + s + t );

  return;
}
/******************************************************************************/

void identity_r ( double r, int i, double *xi )

/******************************************************************************/
/*
  Purpose:

    IDENTITY_R returns a data component given (R).

  Discussion:

    This is a dummy routine, which simply returns (R).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, the coordinate of a point that lies on the
    boundary of the cube.

    Input, int I, the component of the data which is to be returned.

    Output, double *XI, the I-th component of the data vector X, evaluated
    at the point (R), which is on an endpoint of the unit line segment.
*/
{
  if ( i == 0 )
  {
    *xi = r;
  }
  else
  {
    printf ( "\n" );
    printf ( "IDENTITY_R - Fatal error!\n" );
    printf ( "  Illegal component index I = %d\n", i );
    *xi = 0.0;
  }

  return;
}
/******************************************************************************/

void identity_rs ( double r, double s, int i, double *xi )

/******************************************************************************/
/*
  Purpose:

    IDENTITY_RS returns a data component given (R,S).

  Discussion:

    This is a dummy routine, which simply returns (R,S).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the coordinates of a point that lies on the
    boundary of the square.

    Input, int I, the component of the data which is to be returned.

    Output, double *XI, the I-th component of the data vector X, evaluated
    at the point (R,S), which is on a corner, or edge, of the unit square.
*/
{
  if ( i == 0 )
  {
    *xi = r;
  }
  else if ( i == 1 )
  {
    *xi = s;
  }
  else
  {
    printf ( "\n" );
    printf ( "IDENTITY_RS - Fatal error!\n" );
    printf ( "  Illegal component index I = %d\n", i );
    *xi = 0.0;
  }

  return;
}
/******************************************************************************/

void identity_rst ( double r, double s, double t, int i, double *xi )

/******************************************************************************/
/*
  Purpose:

    IDENTITY_RST returns a data component given (R,S,T).

  Discussion:

    This is a dummy routine, which simply returns (R,S,T).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, T, the coordinates of a point that lies on the
    boundary of the cube.

    Input, int I, the component of the data which is to be returned.

    Output, double *XI, the I-th component of the data vector X, evaluated
    at the point (R,S), which is on a corner, edge or face of the unit cube.
*/
{
  if ( i == 0 )
  {
    *xi = r;
  }
  else if ( i == 1 )
  {
    *xi = s;
  }
  else if ( i == 2 )
  {
    *xi = t;
  }
  else
  {
    printf ( "\n" );
    printf ( "IDENTITY_RST - Fatal error!\n" );
    printf ( "  Illegal component index I = %d\n", i );
    *xi = 0.0;
  }

  return;
}
/******************************************************************************/

void stretch_r ( double r, int i, double *xi )

/******************************************************************************/
/*
  Purpose:

    STRETCH_R returns a data component given (R).

  Discussion:

    This routine shifts by 1 and stretches by 2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, the coordinate of a point that lies on the
    boundary of the cube.

    Input, int I, the component of the data which is to be returned.

    Output, double *XI, the I-th component of the data vector X, evaluated
    at the point (R), which is on an endpoint of the unit line segment.
*/
{
  if ( i == 0 )
  {
    *xi = 2.0 * r + 1.0;
  }
  else
  {
    printf ( "\n" );
    printf ( "STRETCH_R - Fatal error\n" );
    printf ( "  Illegal component index I = %d\n", i );
    *xi = 0.0;
  }

  return;
}
/******************************************************************************/

void stretch_rs ( double r, double s, int i, double *xi )

/******************************************************************************/
/*
  Purpose:

    STRETCH_RS returns a data component given (R,S).

  Discussion:

    This routine shifts by (1,2) and stretches by (3,4).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, the coordinates of a point that lies on the
    boundary of the square.

    Input, int I, the component of the data which is to be returned.

    Output, double *XI, the I-th component of the data vector X, evaluated
    at the point (R,S), which is on a corner, or edge, of the unit square.
*/
{
  if ( i == 0 )
  {
    *xi = 3.0 * r + 1.0;
  }
  else if ( i == 1 )
  {
    *xi = 4.0 * s + 2.0;
  }
  else
  {
    printf ( "\n" );
    printf ( "STRETCH_RS - Fatal error!\n" );
    printf ( "  Illegal component index I = %d\n", i );
    *xi = 0.0;
  }

  return;
}
/******************************************************************************/

void stretch_rst ( double r, double s, double t, int i, double *xi )

/******************************************************************************/
/*
  Purpose:

    STRETCH_RST returns a data component given (R,S,T).

  Discussion:

    This routine shifts by (1,2,3) and stretches by (4,5,6)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 December 2013

  Author:

    John Burkardt

  Parameters:

    Input, double R, S, T, the coordinates of a point that lies on the
    boundary of the cube.

    Input, int I, the component of the data which is to be returned.

    Output, double *XI, the I-th component of the data vector X, evaluated
    at the point (R,S), which is on a corner, edge or face of the unit cube.
*/
{
  if ( i == 0 )
  {
    *xi = 4.0 * r + 1.0;
  }
  else if ( i == 1 )
  {
    *xi = 5.0 * s + 2.0;
  }
  else if ( i == 2 )
  {
    *xi = 6.0 * t + 3.0;
  }
  else
  {
    printf ( "\n" );
    printf ( "STRETCH_RST - Fatal error\n" );
    printf ( "  Illegal component index I = %d\n", i );
    *xi = 0.0;
  }

  return;
}
