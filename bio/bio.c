/*
 *	Binary I/O utility functions
 *
 *	Copyright (c) 1996,1999  Pittsburgh Supercomputing Center
 *                                                          *
 *  This program is distributed in the hope that it will    *
 *  be useful, but WITHOUT ANY WARRANTY; without even the   *
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
 *  nor any of the authors assume any liability for         *
 *  damages, incidental or otherwise, caused by the         *
 *  installation or use of this software.                   *
 *                                                          *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
 *  FDA FOR ANY CLINICAL USE.                               *
 *                                                          *
 *
 *	HISTORY
 *		1/96 Written by Greg Hood (ghood@psc.edu)
 *		3/96 Minor modifications to use with libmri (Hood)
 *		9/99 Made array sizes into long longs to
 *			accommodate large datasets (Hood)
 *		4/00 Added InitBIO routine to set endianness variables
 *			according to machine type (Hood)
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "bio.h"

/* 
  Unless otherwise specified, we do all I/O in native-endian format 
*/

static char rcsid[] = "$Id: libbio.c,v 1.3 1999/07/07 20:14:50 welling Exp $";

#if defined(LITTLE_ENDIAN)
int bio_big_endian_machine = 0;
int bio_big_endian_input = 0;
int bio_big_endian_output = 0;
#else
int bio_big_endian_machine = 1;
int bio_big_endian_input = 1;
int bio_big_endian_output = 1;
#endif

/* 
  This flag will be set to 1 upon detection of any I/O error 
*/

int bio_error = 0;

/******************************************************************************/

void InitBIO ( )

/******************************************************************************/
{
  union {
    long l;
    unsigned char c[4];
  } test4;
  union {
    long l;
    unsigned char c[8];
  } test8;
  int big = 1;
  
  if (sizeof(long) == 4)
    {
      test4.l = 1;
      if (test4.c[0] == 1)
	big = 0;
      else if (test4.c[3] == 1)
	big = 1;
      else
	bio_error = 1;
    }
  else if (sizeof(long) == 8)
    {
      test8.l = 1;
      if (test8.c[0] == 1)
	big = 0;
      else if (test8.c[7] == 1)
	big = 1;
      else
	bio_error = 1;
    }
  else
    bio_error = 1;
  bio_big_endian_machine = big;
  bio_big_endian_input = big;
  bio_big_endian_output = big;

  return;
}
/******************************************************************************/

int BRdUInt8 ( unsigned char *addr )

/******************************************************************************/
{
  return(*((unsigned char *) addr));
}

/******************************************************************************/

void BWrUInt8 ( unsigned char *addr, int v )

/******************************************************************************/
{
  *((unsigned char *) addr) = v;
}
/******************************************************************************/

int BRdInt8 ( unsigned char *addr )

/******************************************************************************/
{
  return(*((signed char *) addr));
}
/******************************************************************************/

void BWrInt8 ( unsigned char *addr, int v )

/******************************************************************************/
{
  *((signed char *) addr) = v;
}
/******************************************************************************/

int BRdInt16 ( unsigned char *addr )

/******************************************************************************/
{
  if (bio_big_endian_input)
    return((((signed char) addr[0]) << 8) |
	   addr[1]);
  else
    return((((signed char) addr[1]) << 8) |
	   addr[0]);
}
/******************************************************************************/

void BWrInt16 ( unsigned char *addr, int v )

/******************************************************************************/
{
  if (bio_big_endian_output)
    {
      addr[0] = (v >> 8) & 0xff;
      addr[1] = v & 0xff;
    }
  else
    {
      addr[0] = v & 0xff;
      addr[1] = (v >> 8) & 0xff;
    }
  return;
}
/******************************************************************************/

int BRdInt32 ( unsigned char *addr )

/******************************************************************************/
{
  if (bio_big_endian_input)
    return((((signed char) addr[0]) << 24) |
	   (addr[1] << 16) |
	   (addr[2] << 8) |
	   addr[3]);
  else
    return((((signed char) addr[3]) << 24) |
	   (addr[2] << 16) |
	   (addr[1] << 8) |
	   addr[0]);
}
/******************************************************************************/

void BWrInt32 ( unsigned char *addr, int v )

/******************************************************************************/
{
  if (bio_big_endian_output)
    {
      addr[0] = (v >> 24) & 0xff;
      addr[1] = (v >> 16) & 0xff;
      addr[2] = (v >> 8) & 0xff;
      addr[3] = v & 0xff;
    }
  else
    {
      addr[0] = v & 0xff;
      addr[1] = (v >> 8) & 0xff;
      addr[2] = (v >> 16) & 0xff;
      addr[3] = (v >> 24) & 0xff;
    }
}
/******************************************************************************/

long long BRdInt64 ( unsigned char *addr )

/******************************************************************************/
{
  if (bio_big_endian_input)
    return((((long long) ((signed char) addr[0])) << 56) |
	   (((long long) addr[1]) << 48) |
	   (((long long) addr[2]) << 40) |
	   (((long long) addr[3]) << 32) |
	   (((long long) addr[4]) << 24) |
	   (((long long) addr[5]) << 16) |
	   (((long long) addr[6]) << 8) |
	   ((long long) addr[7]));
  else
    return((((long long) ((signed char) addr[7])) << 56) |
	   (((long long) addr[6]) << 48) |
	   (((long long) addr[5]) << 40) |
	   (((long long) addr[4]) << 32) |
	   (((long long) addr[3]) << 24) |
	   (((long long) addr[2]) << 16) |
	   (((long long) addr[1]) << 8) |
	   ((long long) addr[0]));
}
/******************************************************************************/

void BWrInt64 ( unsigned char *addr, long long v )

/******************************************************************************/
{
  if (bio_big_endian_output)
    {
      addr[0] = (v >> 56) & 0xff;
      addr[1] = (v >> 48) & 0xff;
      addr[2] = (v >> 40) & 0xff;
      addr[3] = (v >> 32) & 0xff;
      addr[4] = (v >> 24) & 0xff;
      addr[5] = (v >> 16) & 0xff;
      addr[6] = (v >> 8) & 0xff;
      addr[7] = v & 0xff;
    }
  else
    {
      addr[0] = v & 0xff;
      addr[1] = (v >> 8) & 0xff;
      addr[2] = (v >> 16) & 0xff;
      addr[3] = (v >> 24) & 0xff;
      addr[4] = (v >> 32) & 0xff;
      addr[5] = (v >> 40) & 0xff;
      addr[6] = (v >> 48) & 0xff;
      addr[7] = (v >> 56) & 0xff;
    }
}
/******************************************************************************/

float BRdFloat32 ( unsigned char *addr )

/******************************************************************************/
{
  union {
    unsigned char temp[4];
    float val;
  } u;

  if (bio_big_endian_machine ^ bio_big_endian_input)
    {
      u.temp[0] = addr[3];
      u.temp[1] = addr[2];
      u.temp[2] = addr[1];
      u.temp[3] = addr[0];
    }
  else
    {
      u.temp[0] = addr[0];
      u.temp[1] = addr[1];
      u.temp[2] = addr[2];
      u.temp[3] = addr[3];
    }
  return (u.val);
}
/******************************************************************************/

void BWrFloat32 ( unsigned char *addr, float v )

/******************************************************************************/
{
  unsigned char *temp;

  temp = (unsigned char *) &v;
  if (bio_big_endian_machine ^ bio_big_endian_output)
    {
      addr[0] = temp[3];
      addr[1] = temp[2];
      addr[2] = temp[1];
      addr[3] = temp[0];
    }
  else
    {
      addr[0] = temp[0];
      addr[1] = temp[1];
      addr[2] = temp[2];
      addr[3] = temp[3];
    }
}
/******************************************************************************/

double BRdFloat64 ( unsigned char *addr )

/******************************************************************************/
{
  union {
    unsigned char temp[8];
    double val;
  } u;

  if (bio_big_endian_machine ^ bio_big_endian_input)
    {
      u.temp[0] = addr[7];
      u.temp[1] = addr[6];
      u.temp[2] = addr[5];
      u.temp[3] = addr[4];
      u.temp[4] = addr[3];
      u.temp[5] = addr[2];
      u.temp[6] = addr[1];
      u.temp[7] = addr[0];
    }
  else
  {
    memcpy(u.temp, addr, 8);
  }
  return (u.val);
}
/******************************************************************************/

void BWrFloat64 ( unsigned char *addr, double v )

/******************************************************************************/
{
  unsigned char *temp;

  temp = (unsigned char *) &v;
  if (bio_big_endian_machine ^ bio_big_endian_output)
    {
      addr[0] = temp[7];
      addr[1] = temp[6];
      addr[2] = temp[5];
      addr[3] = temp[4];
      addr[4] = temp[3];
      addr[5] = temp[2];
      addr[6] = temp[1];
      addr[7] = temp[0];
    }
  else
    memcpy(addr, temp, 8);
}
/******************************************************************************/

void BRdUInt8Array ( unsigned char *buf, unsigned char *a, long long n )

/******************************************************************************/
{
  memcpy(a, buf, n);
}
/******************************************************************************/

void BWrUInt8Array ( unsigned char *buf, unsigned char *a, long long n )

/******************************************************************************/
{
  memcpy(buf, a, n);
}
/******************************************************************************/

void BRdInt8Array ( unsigned char *buf, signed char *a, long long n )

/******************************************************************************/
{
  memcpy(a, buf, n);
}
/******************************************************************************/

void BWrInt8Array ( unsigned char *buf, signed char *a, long long n)

/******************************************************************************/
{
  memcpy(buf, a, n);
}
/******************************************************************************/

void BRdInt16Array ( unsigned char *buf, short *a, long long n)

/******************************************************************************/
{
  long long i, len;
  unsigned char *dest, temp;

  if (sizeof(short) == 2)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 2*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* swap bytes */
	  dest = (unsigned char *) a;
	  len = 2*n;
	  for (i = 0; i < len; i+=2)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+1];
	      dest[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdInt16(&buf[2*i]);
}
/******************************************************************************/

void BWrInt16Array ( unsigned char *buf, short *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char temp;

  if (sizeof(short) == 2)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 2*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* swap bytes */
	  len = 2*n;
	  for (i = 0; i < len; i += 2)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+1];
	      buf[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrInt16(&buf[2*i], a[i]);
}
/******************************************************************************/

void BRdInt32Array ( unsigned char *buf, int *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char *dest, temp;

  if (sizeof(int) == 4)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 4*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  dest = (unsigned char *) a;
	  len = 4*n;
	  for (i = 0; i < len; i+=4)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+3];
	      dest[i+3] = temp;
	      temp = dest[i+2];
	      dest[i+2] = dest[i+1];
	      dest[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdInt32(&buf[4*i]);
}
/******************************************************************************/

void BWrInt32Array ( unsigned char *buf, int *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char temp;

  if (sizeof(int) == 4)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 4*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* reverse byte order */
	  len = 4*n;
	  for (i = 0; i < len; i += 4)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+3];
	      buf[i+3] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+1];
	      buf[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrInt32(&buf[4*i], a[i]);
}
/******************************************************************************/

void BRdInt64Array ( unsigned char *buf, long long *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char *dest, temp;

  if (sizeof(long long) == 8)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 8*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  dest = (unsigned char *) a;
	  len = 8*n;
	  for (i = 0; i < len; i+=8)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+7];
	      dest[i+7] = temp;
	      temp = dest[i+1];
	      dest[i+1] = dest[i+6];
	      dest[i+6] = temp;
	      temp = dest[i+2];
	      dest[i+2] = dest[i+5];
	      dest[i+5] = temp;
	      temp = dest[i+3];
	      dest[i+3] = dest[i+4];
	      dest[i+4] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdInt64(&buf[8*i]);
}
/******************************************************************************/

void BWrInt64Array ( unsigned char *buf, long long *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char temp;

  if (sizeof(long long) == 8)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 8*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* reverse byte order */
	  len = 8*n;
	  for (i = 0; i < len; i += 8)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+7];
	      buf[i+7] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+6];
	      buf[i+6] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+5];
	      buf[i+5] = temp;
	      temp = buf[i+3];
	      buf[i+3] = buf[i+4];
	      buf[i+4] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrInt64(&buf[8*i], a[i]);
}
/******************************************************************************/

void BRdFloat32Array ( unsigned char *buf, float *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char *dest, temp;

  if (sizeof(float) == 4)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 4*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  dest = (unsigned char *) a;
	  len = 4*n;
	  for (i = 0; i < len; i+=4)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+3];
	      dest[i+3] = temp;
	      temp = dest[i+2];
	      dest[i+2] = dest[i+1];
	      dest[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdFloat32(&buf[4*i]);
}
/******************************************************************************/

void BWrFloat32Array ( unsigned char *buf, float *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char temp;

  if (sizeof(int) == 4)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 4*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* reverse byte order */
	  len = 4*n;
	  for (i = 0; i < len; i += 4)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+3];
	      buf[i+3] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+1];
	      buf[i+1] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrFloat32(&buf[4*i], a[i]);
}
/******************************************************************************/

void BRdFloat64Array ( unsigned char *buf, double *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char *dest, temp;

  if (sizeof(double) == 8)
    {
      /* copy directly into the destination */
      memcpy(a, buf, 8*n);
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  dest = (unsigned char *) a;
	  len = 8*n;
	  for (i = 0; i < len; i+=8)
	    {
	      temp = dest[i];
	      dest[i] = dest[i+7];
	      dest[i+7] = temp;
	      temp = dest[i+1];
	      dest[i+1] = dest[i+6];
	      dest[i+6] = temp;
	      temp = dest[i+2];
	      dest[i+2] = dest[i+5];
	      dest[i+5] = temp;
	      temp = dest[i+3];
	      dest[i+3] = dest[i+4];
	      dest[i+4] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      a[i] = BRdFloat64(&buf[8*i]);
}
/******************************************************************************/

void BWrFloat64Array ( unsigned char *buf, double *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char temp;

  if (sizeof(double) == 8)
    {
      /* copy directly into the destination */
      memcpy(buf, a, 8*n);
      if (bio_big_endian_machine ^ bio_big_endian_output)
	{
	  /* reverse byte order */
	  len = 8*n;
	  for (i = 0; i < len; i += 8)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+7];
	      buf[i+7] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+6];
	      buf[i+6] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+5];
	      buf[i+5] = temp;
	      temp = buf[i+3];
	      buf[i+3] = buf[i+4];
	      buf[i+4] = temp;
	    }
	}
    }
  else
    for (i = 0; i < n; ++i)
      BWrFloat64(&buf[8*i], a[i]);
}
/******************************************************************************/

int FRdUInt8 ( FILE *stream )

/******************************************************************************/
{
  int c;

  if ((c = fgetc(stream)) == EOF)
    {
      bio_error = 1;
      return(0);
    }
  return(c & 0xff);
}
/******************************************************************************/

void FWrUInt8 ( FILE *stream, int v )

/******************************************************************************/
{
  if (fputc(v & 0xff, stream) == EOF)
    bio_error = 1;
}
/******************************************************************************/

int FRdInt8 ( FILE *stream )

/******************************************************************************/
{
  int c;

  if ((c = fgetc(stream)) == EOF)
    {
      bio_error = 1;
      return(0);
    }
  return((int) ((signed char) (c & 0xff)));
}
/******************************************************************************/

void FWrInt8 ( FILE *stream, int v)

/******************************************************************************/
{
  if (fputc(v & 0xff, stream) == EOF)
    bio_error = 1;
}
/******************************************************************************/

int FRdInt16 ( FILE *stream)

/******************************************************************************/
{
  unsigned char buf[2];

  if (fread(buf, 2, 1, stream) != 1)
    {
      bio_error = 1;
      return(0);
    }
  return(BRdInt16(buf));
}
/******************************************************************************/

void FWrInt16 ( FILE *stream, int v)

/******************************************************************************/
{
  unsigned char buf[2];

  BWrInt16(buf, v);
  if (fwrite(buf, 2, 1, stream) != 1)
    bio_error = 1;
}
/******************************************************************************/

int FRdInt32 ( FILE *stream)

/******************************************************************************/
{
  unsigned char buf[4];

  if (fread(buf, 4, 1, stream) != 1)
    {
      bio_error = 1;
      return(0);
    }
  return(BRdInt32(buf));
}
/******************************************************************************/

void FWrInt32 ( FILE *stream, int v )

/******************************************************************************/
{
  unsigned char buf[4];

  BWrInt32(buf, v);
  if (fwrite(buf, 4, 1, stream) != 1)
    bio_error = 1;
}
/******************************************************************************/

long long FRdInt64 ( FILE *stream )

/******************************************************************************/
{
  unsigned char buf[8];

  if (fread(buf, 8, 1, stream) != 1)
    {
      bio_error = 1;
      return(0);
    }
  return(BRdInt64(buf));
}
/******************************************************************************/

void FWrInt64 ( FILE *stream, long long v )

/******************************************************************************/
{
  unsigned char buf[8];

  BWrInt64(buf, v);
  if (fwrite(buf, 8, 1, stream) != 1)
    bio_error = 1;
}
/******************************************************************************/

float FRdFloat32 ( FILE *stream )

/******************************************************************************/
{
  unsigned char buf[4];

  if (fread(buf, 4, 1, stream) != 1)
    {
      bio_error = 1;
      return(0.0);
    }
  return(BRdFloat32(buf));
}
/******************************************************************************/

void FWrFloat32 ( FILE *stream, float v )

/******************************************************************************/
{
  unsigned char buf[4];

  BWrFloat32(buf, v);
  if (fwrite(buf, 4, 1, stream) != 1)
    bio_error = 1;
}
/******************************************************************************/

double FRdFloat64 ( FILE *stream )

/******************************************************************************/
{
  unsigned char buf[8];

  if (fread(buf, 1, 8, stream) != 8)
    {
      bio_error = 1;
      return(0);
    }
  return(BRdFloat64(buf));
}
/******************************************************************************/

void FWrFloat64 ( FILE *stream, double v )

/******************************************************************************/
{
  unsigned char buf[8];

  BWrFloat64(buf, v);
  if (fwrite(buf, 8, 1, stream) != 1)
    bio_error = 1;
}
/******************************************************************************/

void FRdUInt8Array ( FILE *stream, unsigned char *a, long long n )

/******************************************************************************/
{
  if (fread(a, 1, n, stream) != n)
    bio_error = 1;
}
/******************************************************************************/

void FWrUInt8Array ( FILE *stream, unsigned char *a, long long n )

/******************************************************************************/
{
  if (fwrite(a, 1, n, stream) != n)
    bio_error = 1;
}
/******************************************************************************/

void FRdInt8Array ( FILE *stream, signed char *a, long long n )

/******************************************************************************/
{
  if (fread(a, 1, n, stream) != n)
    bio_error = 1;
}
/******************************************************************************/

void FWrInt8Array ( FILE *stream, signed char *a, long long n )

/******************************************************************************/
{
  if (fwrite(a, 1, n, stream) != n)
    bio_error = 1;
}
/******************************************************************************/

void FRdInt16Array ( FILE *stream, short *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char *buf, temp;

  if (sizeof(short) == 2)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 2, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* swap bytes */
	  len = 2*n;
	  for (i = 0; i < len; i+=2)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+1];
	      buf[i+1] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(2*n);
      if (fread(buf, 2, n, stream) != n)
	bio_error = 1;
      BRdInt16Array(buf, a, n);
      free(buf);
    }
}
/******************************************************************************/

void FWrInt16Array ( FILE *stream, short *a, long long n )

/******************************************************************************/
{
  unsigned char *buf;

  if (sizeof(short) == 2 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 2, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(2*n);
      BWrInt16Array(buf, a, n);
      if (fwrite(buf, 2, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}
/******************************************************************************/

void FRdInt32Array ( FILE *stream, int *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char *buf, temp;

  if (sizeof(int) == 4)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 4, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  len = 4*n;
	  for (i = 0; i < len; i+=4)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+3];
	      buf[i+3] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+2];
	      buf[i+2] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(4*n);
      if (fread(buf, 4, n, stream) != n)
	bio_error = 1;
      BRdInt32Array(buf, a, n);
      free(buf);
    }
}
/******************************************************************************/

void FWrInt32Array ( FILE *stream, int *a, long long n )

/******************************************************************************/
{
  unsigned char *buf;

  if (sizeof(int) == 4 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 4, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(4*n);
      BWrInt32Array(buf, a, n);
      if (fwrite(buf, 4, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}
/******************************************************************************/

void FRdInt64Array ( FILE *stream, long long *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char *buf, temp;

  if (sizeof(long long) == 8)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 8, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  len = 8*n;
	  for (i = 0; i < len; i+=8)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+7];
	      buf[i+7] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+6];
	      buf[i+6] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+5];
	      buf[i+5] = temp;
	      temp = buf[i+3];
	      buf[i+3] = buf[i+4];
	      buf[i+4] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(8*n);
      if (fread(buf, 8, n, stream) != n)
	bio_error = 1;
      BRdInt64Array(buf, a, n);
      free(buf);
    }
}
/******************************************************************************/

void FWrInt64Array ( FILE *stream, long long *a, long long n )

/******************************************************************************/
{
  unsigned char *buf;

  if (sizeof(long long) == 8 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 8, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(8*n);
      BWrInt64Array(buf, a, n);
      if (fwrite(buf, 8, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}
/******************************************************************************/

void FRdFloat32Array ( FILE *stream, float *a, long long n ) 

/******************************************************************************/
{
  long long i, len;
  unsigned char *buf, temp;

  if (sizeof(float) == 4)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 4, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  len = 4*n;
	  for (i = 0; i < len; i+=4)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+3];
	      buf[i+3] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+2];
	      buf[i+2] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(4*n);
      if (fread(buf, 4, n, stream) != n)
	bio_error = 1;
      BRdFloat32Array(buf, a, n);
      free(buf);
    }
}
/******************************************************************************/

void FWrFloat32Array ( FILE *stream, float *a, long long n )

/******************************************************************************/
{
  unsigned char *buf;

  if (sizeof(float) == 4 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 4, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(4*n);
      BWrFloat32Array(buf, a, n);
      if (fwrite(buf, 4, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}
/******************************************************************************/

void FRdFloat64Array ( FILE *stream, double *a, long long n )

/******************************************************************************/
{
  long long i, len;
  unsigned char *buf, temp;

  if (sizeof(double) == 8)
    {
      buf = (unsigned char *) a;
      if (fread(buf, 8, n, stream) != n)
	bio_error = 1;
      if (bio_big_endian_machine ^ bio_big_endian_input)
	{
	  /* reverse byte order */
	  len = 8*n;
	  for (i = 0; i < len; i+=8)
	    {
	      temp = buf[i];
	      buf[i] = buf[i+7];
	      buf[i+7] = temp;
	      temp = buf[i+1];
	      buf[i+1] = buf[i+6];
	      buf[i+6] = temp;
	      temp = buf[i+2];
	      buf[i+2] = buf[i+5];
	      buf[i+5] = temp;
	      temp = buf[i+3];
	      buf[i+3] = buf[i+4];
	      buf[i+4] = temp;
	    }
	}
    }
  else
    {
      buf = (unsigned char *) malloc(8*n);
      if (fread(buf, 8, n, stream) != n)
	bio_error = 1;
      BRdFloat64Array(buf, a, n);
      free(buf);
    }
}
/******************************************************************************/

void FWrFloat64Array ( FILE *stream, double *a, long long n )

/******************************************************************************/
{
  unsigned char *buf;

  if (sizeof(double) == 8 &&
      !(bio_big_endian_machine ^ bio_big_endian_output))
    {
      /* directly write out to the file */
      if (fwrite(a, 8, n, stream) != n)
	bio_error = 1;
    }
  else
    {
      buf = (unsigned char *) malloc(8*n);
      BWrFloat64Array(buf, a, n);
      if (fwrite(buf, 8, n, stream) != n)
	bio_error = 1;
      free(buf);
    }
}
