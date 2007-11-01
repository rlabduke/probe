/* name: readPDBrecs.h                               */
/* author: J. Michael Word     date written: 2/16/96 */
/* purpose: read records in a pdb file               */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef READPDBRECS_H
#define READPDBRECS_H 1

#include <stdio.h>
#include "utility.h"
#include "geom3d.h"
#include "hybrid_36_c.h"

/* while parsing, the last rec is kept in globPDBrec */
#define PDBRECSIZE 80
extern char globPDBrec[];

char * getPDBrecord(FILE *inf);
int readRecord(FILE *inf, char buffer[], int maxChars);
int isAtom(char *line);
int isHet(char *line);
int isPseudoAtom(char *line);
int isTer(char *line);
int isModel(char *line);
int parseModel(char *line);
/*int parseAtomNumber(char *line);
int parseResidueNumber(char *line);*/
int parseAtomNumber(char *line);
int parseResidueNumber(char *line); 
void parseResidueHy36Num(char *line, char Hy36resno[]); 
void parseChain(char *line, char chain[]); 
/*char parseChain(char *line);*/
float parseOccupancy(char *line);
float parseTempFactor(char *line);
char parseResidueInsertionCode(char *line);
char parseAltLocCode(char *line);
void parseResidueName(char *line, char name[]);
void parseAtomName(char *line, char name[]);
void parseXYZ(char *line, point3d *loc);
void parseSegID(char *line, char id[]);
#endif

