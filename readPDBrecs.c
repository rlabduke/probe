/* name: readPDBrecs.c                               */
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

#include <string.h>
#include <sys/types.h>
#include <ctype.h>
#include "hybrid_36_c.h"
#include "readPDBrecs.h"

char globPDBrec[PDBRECSIZE + 1];
char globPDBrecLen = 0;

char * getPDBrecord(FILE *inf) {
   int rlen;

   rlen = readRecord(inf, globPDBrec, PDBRECSIZE );
   if (rlen < 0) {
      globPDBrecLen = 0;
      return NULL;
   }
   else {
      globPDBrecLen = rlen;
      return globPDBrec;
   }
}

/* readRecord() - read characters until the end of the line  */
/*                putting the first maxChars characters into */
/*                buffer and adding a trailing end-of-string.*/
/*                Returns the length of the buffer string.   */
int readRecord(FILE *inf, char buffer[], int maxChars) {
   register int ch;
   register int count = 0;
   register int done  = 0;
   int makeUpperCase  = 1;
   int setCase = 0;
   char linecard[5];

   if (!inf || feof(inf)) {
      count = -1; /* signal end-of-file */
      buffer[0] = '\0';
   }
   else {
      while(!done) {
	 ch = getc(inf);
	 // This block added by cjw, 01-21-2020
	 // Forces uppercase on record type, resname, and atomname fields
	 // Preserves existing case for all other fields, assumes internal file
	 //   consistency.
	 if (count < 6){
	   //set 6-char line identifier (REMARK, ATOM) to uppercase
	   ch = toupper(ch);
	   if (count < 4){
	     //build short record type from first 4 chars
	     linecard[count] = ch;
	     }
	   else if (count ==4){
	     linecard[count] = '\0'; //finish linecard string
	     //setCase if line is applicable record
	     if ((strncmp(linecard, "ATOM", 4) == 0) || (strncmp(linecard, "HETA", 4) == 0) || (strncmp(linecard, "TER", 3) == 0)){
	       setCase = 1;
	       }
	     }
	   }
//0123456789*123456789*123456789*123456789*123456789*123456789*123456789*1234567
//ATOM   1568  CE2 TYR B  90       0.565  40.044  57.767  1.00 16.94           C
	 else if (setCase && count >= 12 && count < 20 && count != 16){
	   //fix case for atom (12-15) and resname (17-19)
	   ch = toupper(ch);
	   }

	 // This block removed by cjw, 12-09-2019
	 // It appears to convert PDB lines to all-uppercase on read-in
	 // Some large files use lowercase chain IDs, which *must* be preserved
	 // Be warned that its removal allows past unsupported lower case chars in
	 //   atomnames, etc
	 //if ((ch == '@') || (ch == '#')) {
	 //   /* special codes turn off forced uppercase */
	 //   /* "at sign" used for commands */
	 //   /* "hash" used for comments */
	 //   makeUpperCase  = 0;
	 //}
	 //if (makeUpperCase) {
	 //   ch = toupper(ch);
	 //}
	 if ((ch == EOF) || (ch == '\n')) {
	    done = 1;
	 }
	 else if (ch == '\r') {
	    ch = getc(inf);
	    if (ch != '\n') { ungetc(ch, inf); }
	    done = 1;
	 }
	 else if (count < maxChars) {
	    buffer[count++] = ch;
	 }
	 /*   otherwise we drop ch in the bit bucket */
      }
      buffer[count] = '\0';
   }

   return count;
}

/* isAtom() - basic validation for an atom record         */
/*                  must be long enough for 1 digit temp factor */
int isAtom(char *line) {
   if (line) {
      return (strncmp(line, "ATOM", 4) == 0) && (strlen(line) >= 47);
   }
   else return 0;
}

/* isHet() - basic validation for an het record         */
int isHet(char *line) {
   if (line) {
      return (strncmp(line, "HETA", 4) == 0) && (strlen(line) >= 47);
   }
   else return 0;
}

int isPseudoAtom(char *line) {
   if (line) {
      return (line[13] == 'Q');
   }
   else return 0;
}

/* isTer() - is this a ter record? */
int isTer(char *line) {
   if (line) {
      return (strncmp(line, "TER", 3) == 0) && (strlen(line) >= 27);
   }
   else return 0;
}

/* isModel() - is this a Model record? */
int isModel(char *line) {
   if (line) {
      return (strncmp(line, "MODEL", 5) == 0) && (strlen(line) >= 14);
   }
   else return 0;
}

int parseModel(char *line) {

   return parseInteger(line, 6, 8);
}

/*
int parseAtomNumber(char *line) {

   return parseInteger(line, 6, 5);
}

int parseResidueNumber(char *line) {

   return parseInteger(line, 22, 4);
}
*/

int parseAtomNumber(char *line) {
   int atomno;
   const char* errmsg = hy36decode(5, &line[6], 5, &atomno);
/*   if (errmsg) throw std::runtime_error(errmsg);
     fprintf(stderr, "ATOM NUMBER   %d\n", atomno); */

   return atomno;
}

int parseResidueNumber(char *line) {
   int resid;
   const char* errmsg = hy36decode(4, &line[22], 4, &resid);
/*   if (errmsg) throw std::runtime_error(errmsg);
   fprintf(stderr, "RESIDUE NUMBER   %d\n", resid); */

   return resid;
}

void parseResidueHy36Num(char *line, char Hy36resno[]) {
   Hy36resno[0] = line[22];
   Hy36resno[1] = line[23];
   Hy36resno[2] = line[24];
   Hy36resno[3] = line[25];
   Hy36resno[4] = '\0';
}

void parseChain(char *line, char chain[]) {
   chain[0] = line[20];
   chain[1] = line[21];
   chain[2] = '\0';
}

/*  Now using two character chain ids
char parseChain(char *line) {

   return line[21];
}
*/

float parseOccupancy(char *line) {

   return nonblankrange(line, 54, 6) ? parseReal(line, 54, 6) : 1.0;
}

float parseTempFactor(char *line) {

   return nonblankrange(line, 60, 6) ? parseReal(line, 60, 6) : 0.0;
}

char parseResidueInsertionCode(char *line) {

   return line[26];
}

char parseAltLocCode(char *line) {

   return line[16];
}

void parseResidueName(char *line, char name[]) {
   name[0] = line[17];
   name[1] = line[18];
   name[2] = line[19];
   name[3] = '\0';
}

void parseAtomName(char *line, char name[]) {
   name[0] = line[12];
   name[1] = line[13];
   name[2] = line[14];
   name[3] = line[15];
   name[4] = '\0';
}

void parseXYZ(char *line, point3d *loc) {
   loc->x = parseReal(line, 30, 8);
   loc->y = parseReal(line, 38, 8);
   loc->z = parseReal(line, 46, 8);
}


void parseSegID(char *line, char id[]) {
   int reclen = strlen(line);
   id[0] = (reclen >= 73) ? line[72] : ' ';
   id[1] = (reclen >= 74) ? line[73] : ' ';
   id[2] = (reclen >= 75) ? line[74] : ' ';
   id[3] = (reclen >= 76) ? line[75] : ' ';
   id[4] = '\0';
}
