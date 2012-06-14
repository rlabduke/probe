/* name: geom3d.c                                 */
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

#include <math.h>
#include "geom3d.h"

double    v3squaredLength(vector3d* a) {
	return ((a->x * a->x)+(a->y * a->y)+(a->z * a->z));
}

double    v3length(vector3d* a) {
	return(sqrt(v3squaredLength(a)));
}

vector3d* v3negate(vector3d* v) {
	v->x = -v->x; v->y = -v->y; v->z = -v->z;
	return(v);
}

vector3d* v3normalize(vector3d* v) {
	double len = v3length(v);
	if (len != 0.0) { v->x /= len; v->y /= len; v->z /= len; }
	return(v);
}

vector3d* v3scale(vector3d* v, double newlen) {
	double len = v3length(v), scale;
	if (len != 0.0) {
		scale = newlen/len;
		v->x *= scale; v->y *= scale; v->z *= scale;
	}
	return(v);
}

vector3d* v3add(vector3d* a, vector3d* b, vector3d* c) {
	c->x = a->x+b->x; c->y = a->y+b->y; c->z = a->z+b->z;
	return(c);
}

vector3d* v3sub(vector3d* a, vector3d* b, vector3d* c) {
	c->x = a->x-b->x; c->y = a->y-b->y; c->z = a->z-b->z;
	return(c);
}

double    v3dot(vector3d* a, vector3d* b) {
	return((a->x*b->x) + (a->y*b->y) + (a->z*b->z));
}

vector3d* v3lerp(vector3d* lo, vector3d* hi, double alpha, vector3d* rslt) {

	rslt->x = LERP(alpha, lo->x, hi->x);
	rslt->y = LERP(alpha, lo->y, hi->y);
	rslt->z = LERP(alpha, lo->z, hi->z);
	return(rslt);
}

vector3d* v3cross(vector3d* a, vector3d* b, vector3d* c) {
	c->x = (a->y*b->z) - (a->z*b->y);
	c->y = (a->z*b->x) - (a->x*b->z);
	c->z = (a->x*b->y) - (a->y*b->x);
	return(c);
}

double    v3distanceSq(point3d* a, point3d* b) {
	double dx = a->x-b->x;
	double dy = a->y-b->y;
	double dz = a->z-b->z;
	return((dx*dx)+(dy*dy)+(dz*dz));
}

double    v3distance(point3d* a, point3d* b) {
	return(sqrt(v3distanceSq(a, b)));
}

vector3d*    v3makeVec(point3d* a, point3d* b, vector3d* v) {
	v->x = a->x - b->x;
	v->y = a->y - b->y;
	v->z = a->z - b->z;
	return(v3normalize(v));
}

matrix4*    v3identityMat(matrix4* m) {
	m->element[0][0] = 1.0;
	m->element[0][1] = 0.0;
	m->element[0][2] = 0.0;
	m->element[0][3] = 0.0;

	m->element[1][0] = 0.0;
	m->element[1][1] = 1.0;
	m->element[1][2] = 0.0;
	m->element[1][3] = 0.0;

	m->element[2][0] = 0.0;
	m->element[2][1] = 0.0;
	m->element[2][2] = 1.0;
	m->element[2][3] = 0.0;

	m->element[3][0] = 0.0;
	m->element[3][1] = 0.0;
	m->element[3][2] = 0.0;
	m->element[3][3] = 1.0;

	return(m);
}

point3d*    v3mulPointMat(point3d* p, matrix4* m) {
	double x, y, z, w;
	x = (p->x * m->element[0][0]) +
	    (p->y * m->element[1][0]) +
	    (p->z * m->element[2][0]) +
	            m->element[3][0]; 
	y = (p->x * m->element[0][1]) +
	    (p->y * m->element[1][1]) +
	    (p->z * m->element[2][1]) +
	            m->element[3][1]; 
	z = (p->x * m->element[0][2]) +
	    (p->y * m->element[1][2]) +
	    (p->z * m->element[2][2]) +
	            m->element[3][2]; 
	w = (p->x * m->element[0][3]) +
	    (p->y * m->element[1][3]) +
	    (p->z * m->element[2][3]) +
	            m->element[3][3]; 
	p->x = x;
	p->y = y;
	p->z = z;

	if (w != 0.0) { p->x /= w; p->y /= w; p->z /= w; }
	return(p);
}

matrix4*    v3matMul(matrix4* a, matrix4* b, matrix4* c) {
	int i, j, k;
	for (i=0; i < 4; i++) {
		for (j=0; j < 4; j++) {
			c->element[i][j] = 0.0;
			for (k=0; k < 4; k++) {
				c->element[i][j] +=
					a->element[i][k] * b->element[k][j];
			}
		}
	}
	return(c);
}

/* v3rotationMat() - build a rotation matrix */

matrix4*  v3rotationMat(vector3d* axis, double theta, matrix4* m) {
	double c, s, t, x, y, z;

	x = axis->x;
	y = axis->y;
	z = axis->z;

	c = cos(theta*DEG2RAD);
	s = sin(theta*DEG2RAD);
	t = 1.0 - c;

	v3identityMat(m);

	m->element[0][0] = t * x * x    +        c;
	m->element[0][1] = t * x * y    +    z * s;
	m->element[0][2] = t * x * z    -    y * s;

	m->element[1][0] = t * y * x    -    z * s;
	m->element[1][1] = t * y * y    +        c;
	m->element[1][2] = t * y * z    +    x * s;

	m->element[2][0] = t * z * x    +    y * s;
	m->element[2][1] = t * z * y    -    x * s;
	m->element[2][2] = t * z * z    +        c;

	return(m);
}

/* v3translationMat() - build a translation matrix */

matrix4*  v3translationMat(vector3d* axis, matrix4* m) {

	v3identityMat(m);

	m->element[3][0] = axis->x;
	m->element[3][1] = axis->y;
	m->element[3][2] = axis->z;

	return(m);
}

/* v3rotate() - rotate p around the a->b axis */

point3d*  v3rotate(point3d* p, double theta, point3d* a, point3d* b) {
	vector3d axis;
	point3d  q;
	matrix4  rotmat;

	v3makeVec(b, a, &axis);
	v3rotationMat(&axis, theta, &rotmat);

	v3sub(p, b, &q);
	v3mulPointMat(&q, &rotmat);
	v3add(&q, b, p);

	return(p);
}

/* v3angle() - calculate the angle (radians) between 3 points */

double  v3angle(point3d* p1, point3d* p2, point3d* p3) {
	vector3d a, b;
	double amag, bmag, dotval, theta;

	dotval = v3dot(v3sub(p1, p2, &a), v3sub(p3, p2, &b));
	amag = v3length(&a);
	bmag = v3length(&b);

	if (amag*bmag < 0.0001) { theta = 0.0; }
	else { theta = acos(dotval/(amag*bmag)); }

	return(theta);
}

/* v3dihedral() - calculate the dihedral angle (radians) given 4 points */

double  v3dihedral(point3d* p1, point3d* p2, point3d* p3, point3d* p4) {
	vector3d a, b, c, d, e, f;
	double dmag, emag, fmag, theta, phi;

	v3cross(v3sub(p1, p2, &a), v3sub(p3, p2, &b), &d);
	v3cross(v3sub(p2, p3, &b), v3sub(p4, p3, &c), &e);
	dmag = v3length(&d);
	emag = v3length(&e);

	if (dmag*emag < 0.0001) { theta = 0.0; }
	else { theta = acos(v3dot(&d, &e)/(dmag*emag)); }

	v3cross(&d, &b, &f); /* this part sets the correct handedness */
	fmag = v3length(&f);
	if (fmag*emag < 0.0001) { phi = 0.0; }
	else { phi = acos(v3dot(&f, &e)/(fmag*emag)); }

	if (phi*RAD2DEG > 90.0) { theta = - theta; };

	return(theta);
}
