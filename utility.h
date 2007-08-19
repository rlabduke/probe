/* name: utility.h              */
/* author: J. Michael Word      */
/* date written: 2/26/96        */
/* purpose: utility functions   */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef UTILITY_H
#define UTILITY_H 1

#define TRUE  1
#define FALSE 0

#define FABS(x)  (((x)<0.0)?-(x):(x))
#define MIN(x,y) (((x)<(y))? (x):(y))
#define MAX(x,y) (((x)>(y))? (x):(y))

#if defined(__DECCXX_VER) || defined(_MSC_VER)
#define NEEDSTRCASECMP
#endif
#ifdef NEEDSTRCASECMP
int strncasecmp(const char *buf, const char *pat, int sz);
int strcasecmp(const char *buf, const char *pat);
#endif

void note(char *message);
void warn(char *message);
void errmsg(char *message);
void halt(char *message);

int parseInteger(char *str, int start, int len);
float parseReal(char *str, int start, int len);
void copyChars(char *to, char *from, int n);
int compArgStr(char *str, char *arg, int n);
int nonblankstr(char *str);
int nonblankrange(char *str, int start, int len);
#endif
