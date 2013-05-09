/* name: abin.h                                                 */
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

#ifndef ABIN_H
#define ABIN_H 1

#include <math.h>
#include "geom3d.h"

typedef struct {  /* 3D bounding region */
	point3d min, max;
} region;

#define FPRINT_ATOMNAME(outf, aptr) \
              fprintf((outf), "%4.4s%c%3.3s%2s%4d%c",            \
                           (aptr)->atomname, (aptr)->altConf,    \
                           (aptr)->r->resname, (aptr)->r->chain,       \
                           (aptr)->r->resid, (aptr)->r->resInsCode)

struct residue_t;

typedef struct atom_t {
	struct atom_t *next;      /* primary list of atoms         */
	struct atom_t *nextInBin; /* list of atoms in same bin     */
	struct atom_t *scratch;   /* used for arbitrary atom lists */

	char* bondedto; /* for H atoms where the parent atom is specified  */
	                /* for heavy atoms, the list of bonded heavy atoms */
			/* are listed (for standard residues)              */

	int mark; /* general purpose (e.g. used for marking bonded atoms) */

	int  flags;      /* selection mechanism */
	int  props;      /* property flags      */

	int  elem;       /* elememt number  */
	int  atomclass;  /* either elem or (for H) the parent elem or -1 (for unknown) */

	point3d loc;     /* xyz position    */
	int ix, iy, iz;  /* address in bins */

	float radius;    /* VDW radius      */
	float covRad;    /* covalent radius */

	float occ;       /* fractional occupancy */
	float bval;      /* temperature factor   */

	struct residue_t *r;
	char   atomname[5];

	char altConf;    /* atom alternate conformation code */
	char binSerialNum; /* identifies bin where (ix,iy,iz) applies */
} atom, *atomPtr;

typedef struct ringInfo_t {
	struct ringInfo_t *nextRing; /* ring info list */

	point3d ringNormal; /* direction for aromatic ring */
	point3d ringCenter; /* center of aromatic ring */
} ringInfo, *ringInfoPtr;

typedef struct residue_t {
	struct residue_t *nextRes; /* list of residues */
	atom *a;         /* pointer to first atom */

	ringInfo  *ring; /* pointer to ring info list */

	int  file;       /* which file did atom come from? */
	int  model;      /* which model is the atom for? */

	int  resid;        /* residue number */
        char Hy36resno[5]; /* Hybrid 36 residue number */
	int  rescnt;       /* residue order 1, 2, 3, ... */

	char  segid[5];  /* segment identifier   */
	char resname[5]; /* residue name */
	char resInsCode; /* insertion code */

	char chain[5];      /* peptide chain code */
} residue, *residuePtr;

typedef struct {
	point3d min, max, sz; /* bounding box and span */
	int nx, ny, nz;       /* bin size in each dimension */
	float  delta;         /* size bin edges in angstroms */
	char   binSerialNum;  /* identifier */
	atomPtr ***list;      /* 3d matrix of atom lists */
} atomBins;

#define OFFSET(val, min, delta) floor(1.0 + ((val)-(min))/(delta))

#define isDummyAtom(a) ((a).loc.x >= 9998.0 \
                     && (a).loc.y >= 9998.0 \
		     && (a).loc.z >= 9998.0)

atomBins* initBins(char serialNum, region *boundingBox, float delta);
void disposeBins(atomBins* bins);
void updateBoundingBox(point3d *loc, region *boundingBox);
void growBoundingBox(float margin, region *bb);
void addBBox2BBox(region *abb, region *bb);
void addNeighbor(atom *a, atomBins *bins);
void binLoc(point3d *loc, int *ix, int *iy, int *iz, atomBins *bins);

#endif

