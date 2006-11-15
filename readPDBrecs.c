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

   if (!inf || feof(inf)) {
      count = -1; /* signal end-of-file */
      buffer[0] = '\0';
   }
   else {
      while(!done) {
	 ch = getc(inf);
	 if ((ch == '@') || (ch == '#')) {
	    /* special codes turn off forced uppercase */
	    /* "at sign" used for commands */
	    /* "hash" used for comments */
	    makeUpperCase  = 0;
	 }
	 if (makeUpperCase) {
	    ch = toupper(ch);
	 }
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

int parseAtomNumber(char *line) {
	
   return parseInteger(line, 6, 5);
}

int parseResidueNumber(char *line) {
	
   return parseInteger(line, 22, 4);
}

char parseChain(char *line) {
	
   return line[21];
}

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
