/* name: dots.h                                     */
/* purpose: routines to generate points on a sphere */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef DOTS_H
#define DOTS_H 1

#include "geom3d.h"

typedef struct {
	int n;      /* number of three-dimensional points */
	point3d *p; /* array of points */
} pointSet;

void dotSphere(pointSet *set, float radius, float density);
void freeDotSphere(pointSet *set);
int estNumDots(float radius, float density);
int makeDots(float radius, point3d points[], int maxpnts);

#endif
