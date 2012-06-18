/* name: select.h                                */
/* author: J. Michael Word                       */
/* date written: 2/20/96                         */
/* purpose: parse atom selections                */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef SELECT_H
#define SELECT_H 1

#include "abin.h"
#include "parse.h"

typedef struct {char * rlist; char * alist; int bits;} ResidueAndAtomPair;

pattern* getPat(char *line, char *which, int verbose);
void setProperties(atom *a, int hetflag, int hb2aromaticFace, int chohb);
int matchPat(atom *a, pattern *pat);
void setAromaticProp(atom *a);
void setMethylProp(atom *a);
void setMethyleneProp(atom *a); /*20111211dcr*/
int isCarbonylAtom(atom *a);
int atomWithinDistance(atom *a, float *fvec);
char * setHydrogenParentName(char *rname, char *aname);
char * setMainchainBonding(char *rname, char *aname);

int naBaseCategory(atom *a);
#endif
