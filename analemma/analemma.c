# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>
# include <unistd.h>

int main ( int argc, char **argv );
void timestamp ( void );

/******************************************************************************/

int main ( int argc, char **argv )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ANALEMMA.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 January 2013

  Author:

    Original C version by Brian Tung.
    This C version by John Burkardt.

  Local parameters:

    Local, double ECC, the orbital eccentricity.

    Local, double LON, the longitude of the perihelion in radians.

    Local, double OBLIQ, the obliquity in radians.
*/
{
  char c;
  char command_filename[] = "analemma_commands.txt";
  FILE *command_unit;
  char data_filename[] = "analemma_data.txt";
  FILE *data_unit;
  double days = 365.242;
  double dec;
  double degrees = ( 3.141592653589793 / 180.0 );
  double ecc = 0.01671;
  double eot;
  double f;
  double lon = 1.347;
  double obliq = 0.4091;
  extern char *optarg;
  double pi = 3.141592653589793;
  double t;
  double tau;
  double theta;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;
  double z1;
  double z2;
  double z3;

  timestamp ( );
  printf ( "\n" );
  printf ( "ANALEMMA\n" );
  printf ( "  C version\n" );
  printf ( "  Compute and plot the analemma, equation of time, and declination.\n" );
  printf ( "  This program is based on a C program by Brian Tung.\n" );
/* 
  Parse the arguments 
*/
  while ( ( c = getopt ( argc, argv, "e:l:o:h" ) ) >= 0)
  {
    switch ( c ) 
    {
      case 'e':
        ecc = atof ( optarg );
        break;
      case 'l':
        lon = atof ( optarg ) * degrees;
        break;
      case 'o':
        obliq = atof ( optarg ) * degrees;
        break;
      default:
        fprintf ( stderr, "Usage: analemma [options]\n" );
        fprintf ( stderr, "    -e <ecc>    eccentricity " );
        fprintf ( stderr, "(default value: %.5f)\n", ecc );
        fprintf ( stderr, "    -l <lon>    longitude of perihelion in deg " );
        fprintf ( stderr, "(default value: %.2f)\n", lon / degrees );
        fprintf ( stderr, "    -o <obliq>  axial obliquity in deg " );
        fprintf ( stderr, "(default value: %.2f)\n", obliq / degrees );
        fprintf ( stderr, "    -h          print this page\n" );
        exit ( 0 );
    }
  }
/*
  Compute the data.
*/
  data_unit = fopen ( data_filename, "wt" );

  if ( !data_unit )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "ANALEMMA - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the data file \"%s\"\n", data_filename );
    exit ( 1 );
  }

  for ( f = 0.0; f <= 1.0; f = f + 0.0001 ) 
  {
    tau = 2.0 * pi * f;
/* 
  Set theta to the current longitude. 
*/
    theta = atan2 ( sqrt ( 1.0 - ecc * ecc ) * sin ( tau ), cos ( tau ) - ecc );
/* 
  Rotate clockwise in XY plane by theta, corrected by lon.
*/
    x1 = cos ( theta - ( lon - pi / 2.0 ) );
    y1 = sin ( theta - ( lon - pi / 2.0 ) );
    z1 = 0.0;
/* 
  Rotate counter-clockwise in XZ plane by obliq.
*/
    x2 = cos ( obliq ) * x1 + sin ( obliq ) * z1;
    y2 = y1;
    z2 = - sin ( obliq ) * x1 + cos ( obliq ) * z1;
/* 
  Set t equal to real time from tau and
  rotate counter-clockwise by t, corrected by lon 
*/
    t = tau - ecc * sin ( tau );
    x3 =   cos ( t - ( lon - pi / 2.0 ) ) * x2 + sin ( t - ( lon - pi / 2.0 ) ) * y2;
    y3 = - sin ( t - ( lon - pi / 2.0 ) ) * x2 + cos ( t - ( lon - pi / 2.0 ) ) * y2;
    z3 = z2;

    eot = - atan2 ( y3, x3 ) * 4.0 / degrees * days / ( days + 1.0 );
    dec = asin ( z3 ) / degrees;
/* 
  Print results in minutes early/late and degrees north/south 
*/
    fprintf( data_unit, "%.9f  %.9f  %.9f\n", t / ( 2.0 * pi ), eot, dec );
  }

  fclose ( data_unit );

  printf ( "\n" );
  printf ( "  Created data file \"%s\".\n", data_filename );
/*
  Create the command file.
*/
  command_unit = fopen ( command_filename, "wt" );
  fprintf ( command_unit, "set term png\n" );
  fprintf ( command_unit, "set grid\n" );
  fprintf ( command_unit, "set style data lines\n" );
  fprintf ( command_unit, "unset key\n" );
  fprintf ( command_unit, "set output \"eot.png\"\n" );
  fprintf ( command_unit, "set xlabel '<---Normalized Date--->'\n" );
  fprintf ( command_unit, "set ylabel '<---Minutes Early/Late--->'\n" );
  fprintf ( command_unit, "set title 'The equation of time'\n" );
  fprintf ( command_unit, "plot '%s' using 1:2 with lines\n", data_filename );
  fprintf ( command_unit, "set output \"declination.png\"\n" );
  fprintf ( command_unit, "set xlabel '<---Normalized Date--->'\n" );
  fprintf ( command_unit, "set ylabel '<---Degrees North/South--->'\n" );
  fprintf ( command_unit, "set title 'Declination'\n" );
  fprintf ( command_unit, "plot '%s' using 1:3 with lines\n", data_filename );
  fprintf ( command_unit, "set output \"analemma.png\"\n" );
  fprintf ( command_unit, "set xlabel '<---Minutes Early/Late--->'\n" );
  fprintf ( command_unit, "set ylabel '<---Degrees North/South--->'\n" );
  fprintf ( command_unit, "set title 'The analemma'\n" );
  fprintf ( command_unit, "plot '%s' using 2:3 with lines\n", data_filename );
  fprintf ( command_unit, "quit\n" );
  fprintf ( command_unit, "\n" );

  fclose ( command_unit );

  printf ( "  Created command file \"%s\".\n", command_filename );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "ANALEMMA\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
