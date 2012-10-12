/* name: stdconntable.h                            */
/* author: J. Michael Word Date Written: 11/15/98  */
/* purpose: describes the bonding connectivies for */
/*          each std amino acid residue or         */
/*          nucleic acid base                      */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef STDCONNTABLE_H
#define STDCONNTABLE_H 1

void initStdConnTable();
char * searchForStdBondingPartner(char *resname, char *atomname, int isAhydrogen);
void dumpStdConnTable(FILE * outf);

#endif
