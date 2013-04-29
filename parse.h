/* name: parse.h                                         */
/* author: J. Michael Word     date written: 2/29/96     */
/* purpose: build a parse tree from an input string.     */
/*          based on "Compilers", Aho, Sethi, and Ullman */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef PARSE_H
#define PARSE_H 1

#include <stdio.h>

/* types of pattern entries */
#define         OR_NODE 300
#define        AND_NODE 301
#define        NOT_NODE 302
#define      RANGE_NODE 303
#define       FILE_NODE 304
#define      MODEL_NODE 305
#define      CHAIN_NODE 306
#define        ALT_NODE 307
#define        RES_NODE 308
#define      RTYPE_NODE 309
#define       TRUE_NODE 310
#define      FALSE_NODE 311
#define       PROP_NODE 312
#define      ANAME_NODE 313
#define     OCC_LT_NODE 314
#define     OCC_GT_NODE 315
#define       B_LT_NODE 316
#define       B_GT_NODE 317
#define       DIST_NODE 318
#define      SEGID_NODE 319
#define        INS_NODE 320
#define  INS_RANGE_NODE 321

/* lexical analyzer tokens */
#define      NOT_TOK 400
#define      NUM_TOK 401
#define  REALNUM_TOK 402
#define       ID_TOK 403
#define   WITHIN_TOK 404
#define     DONE_TOK 405

/* flags for atom properties */
#define         MC_PROP (1 <<  0)
#define         SC_PROP (1 <<  1)
#define      ALPHA_PROP (1 <<  2)
#define       BETA_PROP (1 <<  3)
#define          N_PROP (1 <<  4)
#define          C_PROP (1 <<  5)
#define          O_PROP (1 <<  6)
#define          S_PROP (1 <<  7)
#define          P_PROP (1 <<  8)
#define          H_PROP (1 <<  9)
#define      WATER_PROP (1 << 10)
#define        HET_PROP (1 << 11)
#define       PROT_PROP (1 << 12)
#define        DNA_PROP (1 << 13)
#define      METAL_PROP (1 << 14) /*re atomprops.h METALIC_ATOM_FLAG */
#define     METHYL_PROP (1 << 15)
#define      DONOR_PROP (1 << 16)
#define   ACCEPTOR_PROP (1 << 17)
#define   AROMATIC_PROP (1 << 18)
#define   CH_DONOR_PROP (1 << 19)
#define TEST_ACCEPT_ANGLE_PROP (1 << 20)

#define   NEGATIVE_PROP (1 << 21)
#define   POSITIVE_PROP (1 << 22)

#define   RHPHOBIC_PROP (1 << 23)
#define   RHPHILIC_PROP (1 << 24)
#define   RCHARGED_PROP (1 << 25)

#define   MAYBECHG_PROP (1 << 26)
#define CHECK_ENDS_PROP (1 << 27)
#define AMBIGWATER_PROP (1 << 28)
#define        ION_PROP (1 << 29) /*re atomprops.h IONIC_ATOM_FLAG dcr041007*/
#define  METHYLENE_PROP (1 << 30) /*re atomprops.h dcr20111210 */
typedef struct patternType {
   struct patternType *lhs; /* link to sub-pattern */
   struct patternType *rhs; /*       ditto         */
   int type; /* pattern node type    */
   int val;  /* value for leaf nodes */
   float fvec[4];  /* real vector leaf data */
} pattern;

/* lexical analyzer symbol table entry */
typedef struct {
   char *lexptr; /* string */
   int token;    /* type   */
} symbolEntry;

#define LXBUFSZ       100
#define LEXSTRINGMAX 1000
#define LEXSYMMAX     500

pattern* parseArg(char *s);
pattern* freePattern(pattern *p);
pattern* exprItem();
pattern* setItem();
pattern* limitItem();
pattern* featureItem();
pattern* numericItem(int signFactor);
pattern* idItem();
pattern* distItem();
pattern* makeNode(pattern *lhs, pattern *rhs, int type);
pattern* makeTerminal(int type, int val);
void matchTok(int t);

char* lexString(int p);
int lexan();
int lookup(char *s);
int insert(char *s, int tok);

void printPattern(FILE *outf, pattern *pat);
void recPrint(FILE *outf, pattern *pat);
char* describeLookahead(char formatstring[], char *s);

int  modelInPattern(pattern *pat);  /*041114*/

#endif

