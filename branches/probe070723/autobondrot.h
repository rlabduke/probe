/* name: autobondrot.h                                          */
/* author: J. Michael Word                                      */
/* date written: 7/ 9/99                                        */
/* purpose: manipulations of moving atoms                       */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef AUTOBONDROT_H
#define AUTOBONDROT_H 1

#include <stdio.h>
#include <math.h>
#include "abin.h"
#include "geom3d.h"

#define XFORMREC_DELIM ":"

#define MAX_XFORM_LEVELS       8
#define MAX_XFORM_STACK_DEPTH  8

#define NULLxform  0
#define ROTxform   1
#define TRANSxform 2
#define SAVExform  3
#define RESTORExform  4

#define CONSTbiasfunc  0
#define COSbiasfunc    1
#define POLYbiasfunc   2

typedef atomPtr (abrMkAtomProc)(char*, void*);    /* make an atom from a record */
typedef void (abrAtomListProc)(atom*, void*);     /* process a list of atoms */
typedef void (abrProbeProc)(char*, double, atom*, void*); /* process a list of atoms */

typedef struct xformAtomRecords_t {
   struct xformAtomRecords_t *next; /* next record in list */
   atom    *a;                      /* atom being worked on */
   point3d loc;                     /* original location (at "phase") */
} xformAtomRecords;

typedef struct goToRec_t {
   struct goToRec_t *next; /* next record in list */
   double level[MAX_XFORM_LEVELS]; /* value of each */
} goToRec;

typedef struct biasFunction_t {
   struct biasFunction_t *next; /* next record in list */
   int type;       /* bias function type */
   double v;       /* scale */
   double ph;      /* phase angle or zeropoint */
   double freq;    /* cyclical frequency */
   double shift;   /* [shift - cos()] usually 1 */
} biasFunction;

typedef struct transformdata_t {
   struct transformdata_t *next; /* list of transforms */

   xformAtomRecords* recs; /* for which the current transformations apply */
   xformAtomRecords* lastrec; /* last atom record in the list */

   char *name;      /* used to show the name */
   int num;         /* level number (helps when debugging) */

   int type;        /* ROTxform, etc. */

   double phase;    /* e.g. original angle in the starting structure */
   double startVal; /* starting value  */
   double endVal;   /* ending value    */
   double stepVal;  /* change in value */

   double currVal;  /* working value (changes during processing) */

   point3d a1, a2;  /* axis vector end points */

   matrix4 tmat;    /* working transformation matrix (changes) */

   biasFunction* funcs; /* list of bias functions */
} transformData;

typedef struct {
   transformData *xforms; /* list of transformations */
   transformData *last;   /* last transformation in list */
   
   int nlevs;       /* number of adjustable levels   */
   transformData* level[MAX_XFORM_LEVELS]; /* index of each */

   int nstack; /* how many CTMs on the stack */
   matrix4* ctmstack[MAX_XFORM_STACK_DEPTH];

   goToRec* golist; /* if non null, we do these instead of loop*/
   goToRec* glast;  /* last goTo record in list */

   xformAtomRecords** sortedByRes; /* array of atom recs sorted by residue */
} xformDatabase;

xformDatabase* readTransformationDatabase(FILE *inf, FILE *outf,
	    abrMkAtomProc mkAtom, void *atomstuff,
	    abrAtomListProc inputListProc, void *liststuff,
	    char* cch);
void discardXformDB(xformDatabase* db,
	       abrAtomListProc delAtomProc, void *deletestuff);
void autobondrot(FILE *outf, xformDatabase* xdb,
	       abrProbeProc    probeProc,   void *probestuff,
	       abrAtomListProc delAtomProc, void *deletestuff,
	       int dumpGoToAtoms);
#endif
