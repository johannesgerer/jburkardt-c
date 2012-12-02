# include <stdlib.h>
#include <math.h>

#include "vec.h"

#define max(a,b) ((a>b)?a:b)

Vec zero_vec(void)
{
   Vec vec;

   vec = (Vec)malloc(sizeof(_Vec));
   vec->x = vec->y = vec->z = 0;
   return vec;
}

Vec new_vec(double X, double Y, double Z)
{
   Vec vec;

   vec = zero_vec();
   vec->x = X; vec->y = Y; vec->z = Z;
   return vec;
}

Vec vec_copy(Vec v1)
{
   Vec vec;

   vec = zero_vec();
   vec->x = v1->x; vec->y = v1->y; vec->z = v1->z;
   return vec;
}

void vec_clear(Vec vec)
{
   vec->x = vec->y = vec->z = 0;
}

void vec_diff(Vec vec1, Vec vec2, Vec vec3)
{
   vec3->x = vec1->x - vec2->x;
   vec3->y = vec1->y - vec2->y;
   vec3->z = vec1->z - vec2->z;
}

void vec_add(Vec vec1, Vec vec2)
{
   vec1->x += vec2->x;
   vec1->y += vec2->y;
   vec1->z += vec2->z;
}

void vec_smul(Vec vec, double scalar)
{
   vec->x *= scalar;
   vec->y *= scalar;
   vec->z *= scalar;
}

void vec_sdiv(Vec vec, double scalar)
{
   vec->x /= scalar;
   vec->y /= scalar;
   vec->z /= scalar;
}

void vec_rshift(Vec vec, int n)
{
   vec->x = (int)vec->x >>n;
   vec->y = (int)vec->y >>n;
   vec->z = (int)vec->z >>n;
}

void vec_lshift(Vec vec, int n)
{
   vec->x = (int)vec->x <<n;
   vec->y = (int)vec->y <<n;
   vec->z = (int)vec->z <<n;
}

/*
 * Limit the length of the longest
 * component to lim, while keeping others
 * in proportion
 */
void vec_limit(Vec vec, double lim)
{
   double m,f;

   m = max(fabs(vec->x), fabs(vec->y));
   m = max(m, fabs(vec->z));

   if(m <= lim) return;

   f = lim/m;
   vec_smul(vec, f);
}

/*
 * Set the magnitude of the vector to a
 * particular value
 */
void vec_setmag(Vec vec, double mag)
{
   double m,f;

   m = max(fabs(vec->x), fabs(vec->y));
   m = max(m, fabs(vec->z));

   f = mag/m;
   vec_smul(vec, f);
}

/*
 * Rectangular (ie. min component) distance
 * between this and vec2
 */
double vec_rdist(Vec vec1, Vec vec2)
{
   double dx,dy,dz, dm;

   dx = vec1->x - vec2->x;
   dy = vec1->y - vec2->y;
   dz = vec1->z - vec2->z;

   dm = max(fabs(dx), fabs(dy));
   dm = max(dm, fabs(dz));

   return dm;
}
