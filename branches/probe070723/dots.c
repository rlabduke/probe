/* name: dots.c                                                 */
/* author: J. Michael Word (port from dcr and mez fortran code) */
/* date written: 2/20/96                                        */
/* purpose: generate points on a sphere                         */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#include <stdlib.h>
#include <math.h>
#include "dots.h"

#define PI 3.14159265359

void dotSphere(pointSet *set, float radius, float density) {
   int m;

   if (set) {
      m = estNumDots(radius, density);
      set->n = 0;
      set->p = (point3d *)malloc(m*sizeof(point3d));
      if (set->p) {
	 set->n = makeDots(radius, set->p, m);
      }
   }
}

void freeDotSphere(pointSet *set){
   if (set) {
      set->n = 0;
      if (set->p) {
	 free(set->p);
	 set->p = NULL;
      }
   }
}

int estNumDots(float radius, float density) {
   float sizefact = 1.0;

   /* overestimate of the number of dots */
   return (int)floor(4.0 * PI * density * sizefact * (radius * radius));
}

int makeDots(float radius, point3d points[], int maxpnts) {
   float offset = 0.2;
   double ang, cosang, sinang, phi, theta, xy0, x0, y0, z0;
   int i, j, k, odd, nequator, nvert, nhoriz;

   nequator = (int)floor(sqrt(maxpnts * PI));

   odd = 1;
   ang = 5.0 * PI / 360.0;
   cosang = cos(ang);
   sinang = sin(ang);

   i = 0;
   nvert = nequator / 2;
   for (j = 0; j <= nvert; j++) {
      phi = (PI * j) / nvert;
      z0 = cos(phi) * radius;
      xy0= sin(phi) * radius;

      nhoriz = (int)floor(nequator * sin(phi));
      if (nhoriz < 1) nhoriz = 1;
      for (k = 0; k < nhoriz; k++) {
	 if(odd) {theta = (2.0 * PI * k + offset)/nhoriz; }
	 else    {theta = (2.0 * PI * k         )/nhoriz; }
	 x0 = cos(theta) * xy0;
	 y0 = sin(theta) * xy0;

	 if (i >= maxpnts) return i;
	 points[i].x = x0;
	 points[i].y = y0*cosang - z0*sinang;
	 points[i].z = y0*sinang + z0*cosang;
	 i++;
      }
      odd = !odd;
   }
   return i;
}
