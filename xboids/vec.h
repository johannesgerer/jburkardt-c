#ifndef __VEC_H__
#define __VEC_H__

typedef struct {
   double x, y, z;
} _Vec, *Vec;

Vec zero_vec(void);
Vec new_vec(double X, double Y, double Z);
Vec vec_copy(Vec v1);
void vec_clear(Vec vec);
void vec_diff(Vec vec1, Vec vec2, Vec vec3);
void vec_add(Vec vec1, Vec vec2);
void vec_smul(Vec vec, double scalar);
void vec_sdiv(Vec vec, double scalar);
void vec_rshift(Vec vec, int n);
void vec_lshift(Vec vec, int n);

/*
 * Limit the length of the longest
 * component to lim, while keeping others
 * in proportion
 */
void vec_limit(Vec vec, double lim);

/*
 * Set the magnitude of the vector to a
 * particular value
 */
void vec_setmag(Vec vec, double mag);

/*
 * Rectangular (ie. min component) distance
 * between this and vec2
 */
double vec_rdist(Vec vec1, Vec vec2);


#endif  /*  __VEC_H__  */
