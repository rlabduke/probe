/* name: geom3d.h                                 */
/* author: J. Michael Word   date written: 2/9/96 */
/* purpose: collection of geometric primatives    */
/*          from Graphics Gems book.              */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef GEOM3D_H
#define GEOM3D_H 1

typedef struct {
	double x, y, z;
} point3d, vector3d;

typedef struct {
	double element[4][4];
} matrix4;

/* linear interpolation from l (when a=0) to h (when a=1) */
#define LERP(a,l,h) ((l)+(((h)-(l))*(a)))

#define DEG2RAD 0.017453293 /* convert degrees to radians */
#define RAD2DEG 57.29578    /* convert radians to degrees */

double    v3squaredLength(vector3d* a);
double    v3length(vector3d* a);
vector3d* v3negate(vector3d* a);
vector3d* v3normalize(vector3d* a);
vector3d* v3scale(vector3d* a, double newlen);
vector3d* v3add(vector3d* a, vector3d* b, vector3d* c);
vector3d* v3sub(vector3d* a, vector3d* b, vector3d* c);
double    v3dot(vector3d* a, vector3d* b);
vector3d* v3lerp(vector3d* lo, vector3d* hi, double alpha, vector3d* rslt);
vector3d* v3cross(vector3d* a, vector3d* b, vector3d* c);
double    v3distanceSq(point3d* a, point3d* b);
double    v3distance(point3d* a, point3d* b);
vector3d* v3makeVec(point3d* a, point3d* b, vector3d* v);
matrix4*  v3identityMat(matrix4* m);
point3d*  v3mulPointMat(point3d* p, matrix4* m);
matrix4*  v3matMul(matrix4* a, matrix4* b, matrix4* c);
matrix4*  v3rotationMat(vector3d* a, double theta, matrix4* m);
matrix4*  v3translationMat(vector3d* a, matrix4* m);
point3d*  v3rotate(point3d* a,double theta,point3d* b,point3d* c);
double    v3angle(point3d* p1, point3d* p2, point3d* p3);
double    v3dihedral(point3d* p1, point3d* p2, point3d* p3, point3d* p4);

#endif

