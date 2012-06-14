/* name: abin.c                                                 */
/* author: J. Michael Word (port from dcr and mez fortran code) */
/* date written: 2/20/96                                        */
/* purpose: organize atoms into bins to aid finding neighbors   */

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
#include <stdio.h>
#include "utility.h"
#include "abin.h"

atomBins* initBins(char serialNum, region *boundingBox, float delta) {
	atomBins *bp;
	int i, j, k, nx, ny, nz;
	int n1, n2, n3;
	atomPtr ***b1, **b2, *b3;

	if (delta < 0.1) delta = 5.0;

	nx = OFFSET(boundingBox->max.x, boundingBox->min.x, delta) + 2;
	ny = OFFSET(boundingBox->max.y, boundingBox->min.y, delta) + 2;
	nz = OFFSET(boundingBox->max.z, boundingBox->min.z, delta) + 2;

	bp = (atomBins *)malloc(sizeof(atomBins));

	bp->delta = delta;

	bp->min = boundingBox->min;
	bp->max = boundingBox->max;
	v3sub(&(bp->max), &(bp->min), &(bp->sz));
	bp->nx = nx;
	bp->ny = ny;
	bp->nz = nz;
	bp->binSerialNum = serialNum;

	n1 = nx;
	n2 = nx*ny;
	n3 = nx*ny*nz;

	b1 = (atomPtr ***)malloc(sizeof(atomPtr **)*n1);
	b2 = (atomPtr  **)malloc(sizeof(atomPtr  *)*n2);
	b3 = (atomPtr   *)malloc(sizeof(atomPtr   )*n3);
	if (!bp || !b1 || !b2 || !b3) {
		warn("could not create atom bins");
		if (bp) free(bp);
		if (b1) free(b1);
		if (b2) free(b2);
		if (b3) free(b3);
		return NULL;
	}

	bp->list = b1;

	for(i = 0; i < nx; i++) {
		bp->list[i] = b2; b2 += ny;

		for(j = 0; j < ny; j++) {
			bp->list[i][j] = b3; b3 += nz;

			for(k = 0; k < nz; k++) {
				bp->list[i][j][k] = 0;
			}
		}
	}

	return bp;
}

void disposeBins(atomBins* bins) {
	if (bins) {
		free(bins->list[0][0]);
		free(bins->list[0]);
		free(bins->list);
		free(bins);
	}
}

void updateBoundingBox(point3d *loc, region *bb) {
	if (loc->x < bb->min.x) {bb->min.x = loc->x;}
	if (loc->y < bb->min.y) {bb->min.y = loc->y;}
	if (loc->z < bb->min.z) {bb->min.z = loc->z;}

	if (loc->x > bb->max.x) {bb->max.x = loc->x;}
	if (loc->y > bb->max.y) {bb->max.y = loc->y;}
	if (loc->z > bb->max.z) {bb->max.z = loc->z;}
}

void growBoundingBox(float margin, region *bb) {
  bb->min.x -= margin; bb->min.y -= margin; bb->min.z -= margin;
  bb->max.x += margin; bb->max.y += margin; bb->max.z += margin;
}

void addBBox2BBox(region *abb, region *bb) {
	if (abb->min.x < bb->min.x) {bb->min.x = abb->min.x;}
	if (abb->min.y < bb->min.y) {bb->min.y = abb->min.y;}
	if (abb->min.z < bb->min.z) {bb->min.z = abb->min.z;}

	if (abb->max.x > bb->max.x) {bb->max.x = abb->max.x;}
	if (abb->max.y > bb->max.y) {bb->max.y = abb->max.y;}
	if (abb->max.z > bb->max.z) {bb->max.z = abb->max.z;}
}

void addNeighbor(atom *a, atomBins *bins) {

	binLoc(&(a->loc), &(a->ix), &(a->iy), &(a->iz), bins);

	/* record bin identifier */
	a->binSerialNum = bins->binSerialNum;

	/* add atom in as first link in neighbor list */
	a->nextInBin = bins->list[a->ix][a->iy][a->iz];

	bins->list[a->ix][a->iy][a->iz] = a;
}

void binLoc(point3d *loc, int *ix, int *iy, int *iz, atomBins *bins) {

	*ix = OFFSET(loc->x, bins->min.x, bins->delta);	
	*iy = OFFSET(loc->y, bins->min.y, bins->delta);	
	*iz = OFFSET(loc->z, bins->min.z, bins->delta);	
}

