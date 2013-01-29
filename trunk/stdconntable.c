/* name: stdconntable.c                            */
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "stdconntable.h"

/* the hash modulus M must be prime */
#define M 2617

typedef struct _StdResConnTableEntry {
   char *key;
   char *value;
   struct _StdResConnTableEntry *next;
} StdResConnTableEntry_t;

int HashInStdResTbl(char *s);
int InsertInStdResConnTable(StdResConnTableEntry_t *e);
char *SearchStdResConnTable(char *key);

static StdResConnTableEntry_t *StdResTblBucket[M] = {NULL};

/* For each of the standard amino acids and bases, we  */
/* list the parent atom for each hydrogen in the table */
/* below. We want this table to be as general as       */
/* possible so where alternative names are commonly    */
/* employed for an H atom, we list each of them. Where */
/* the parent atom name has several forms, we list     */
/* each of them separated by a colon                   */
/* Farther down the table is the connectivity for the  */
/* non hydrogens in each standard residue or base.     */

static StdResConnTableEntry_t StandardResAtomConnRec[] = {

/* fix for problem with output of truncated H names from O */
#define ALLOW_TRUNCATED_H_NAMES 1

/*---------------------------------------*/
{"ALA: H  ", " N  : NT ", NULL},
{"ALA: HN ", " N  : NT ", NULL},
{"ALA:1H  ", " N  : NT ", NULL},
{"ALA:2H  ", " N  : NT ", NULL},
{"ALA:3H  ", " N  : NT ", NULL},
{"ALA: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"ALA: H2 ", " N  : NT ", NULL}, /*  */
{"ALA: H3 ", " N  : NT ", NULL}, /*  */
{"ALA: HT1", " N  : NT ", NULL},
{"ALA: HT2", " N  : NT ", NULL},
{"ALA: HT3", " N  : NT ", NULL},
{"ALA: HA ", " CA ", NULL},
{"ALA:1HB ", " CB ", NULL},
{"ALA:2HB ", " CB ", NULL},
{"ALA:3HB ", " CB ", NULL},
{"ALA: HB1", " CB ", NULL},
{"ALA: HB2", " CB ", NULL},
{"ALA: HB3", " CB ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"ALA: HT ", " N  : NT ", NULL},
{"ALA: HB ", " CB ", NULL},
#endif

{"CYS: H  ", " N  : NT ", NULL},
{"CYS: HN ", " N  : NT ", NULL},
{"CYS:1H  ", " N  : NT ", NULL},
{"CYS:2H  ", " N  : NT ", NULL},
{"CYS:3H  ", " N  : NT ", NULL},
{"CYS: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"CYS: H2 ", " N  : NT ", NULL}, /*  */
{"CYS: H3 ", " N  : NT ", NULL}, /*  */
{"CYS: HT1", " N  : NT ", NULL},
{"CYS: HT2", " N  : NT ", NULL},
{"CYS: HT3", " N  : NT ", NULL},
{"CYS: HA ", " CA ", NULL},
{"CYS:1HB ", " CB ", NULL},
{"CYS:2HB ", " CB ", NULL},
{"CYS: HB1", " CB ", NULL},
{"CYS: HB2", " CB ", NULL},
/* {"CYS: HB2", " CB ", NULL},*/ /* remediated names RMI 070718 */
{"CYS: HB3", " CB ", NULL}, /*  */
{"CYS: HG ", " SG ", NULL}, /* check for SH2 ? reduce has HG2 and HG3...  */
#ifdef ALLOW_TRUNCATED_H_NAMES
{"CYS: HT ", " N  : NT ", NULL},
{"CYS: HB ", " CB ", NULL},
#endif

{"ASP: H  ", " N  : NT ", NULL},
{"ASP: HN ", " N  : NT ", NULL},
{"ASP:1H  ", " N  : NT ", NULL},
{"ASP:2H  ", " N  : NT ", NULL},
{"ASP:3H  ", " N  : NT ", NULL},
{"ASP: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"ASP: H2 ", " N  : NT ", NULL}, /*  */
{"ASP: H3 ", " N  : NT ", NULL}, /*  */
{"ASP: HT1", " N  : NT ", NULL},
{"ASP: HT2", " N  : NT ", NULL},
{"ASP: HT3", " N  : NT ", NULL},
{"ASP: HA ", " CA ", NULL},
{"ASP:1HB ", " CB ", NULL},
{"ASP:2HB ", " CB ", NULL},
{"ASP: HB1", " CB ", NULL},
{"ASP: HB2", " CB ", NULL},
/* {"ASP: HB2", " CB ", NULL},*/ /* remediated names RMI 070718 */
{"ASP: HB3", " CB ", NULL}, /*  */
#ifdef ALLOW_TRUNCATED_H_NAMES
{"ASP: HT ", " N  : NT ", NULL},
{"ASP: HB ", " CB ", NULL},
#endif

{"GLU: H  ", " N  : NT ", NULL},
{"GLU: HN ", " N  : NT ", NULL},
{"GLU:1H  ", " N  : NT ", NULL},
{"GLU:2H  ", " N  : NT ", NULL},
{"GLU:3H  ", " N  : NT ", NULL},
{"GLU: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"GLU: H2 ", " N  : NT ", NULL}, /*  */
{"GLU: H3 ", " N  : NT ", NULL}, /*  */
{"GLU: HT1", " N  : NT ", NULL},
{"GLU: HT2", " N  : NT ", NULL},
{"GLU: HT3", " N  : NT ", NULL},
{"GLU: HA ", " CA ", NULL},
{"GLU:1HB ", " CB ", NULL},
{"GLU:2HB ", " CB ", NULL},
{"GLU: HB1", " CB ", NULL},
{"GLU: HB2", " CB ", NULL},
/* {"GLU: HB2", " CB ", NULL}, */ /* remediated names RMI 070718 */
{"GLU: HB3", " CB ", NULL}, /*  */
{"GLU:1HG ", " CG ", NULL},
{"GLU:2HG ", " CG ", NULL},
{"GLU: HG1", " CG ", NULL},
{"GLU: HG2", " CG ", NULL},
/* {"GLU: HG2", " CG ", NULL},  */ /* remediated names RMI 070718 */
{"GLU: HG3", " CG ", NULL}, /*  */
#ifdef ALLOW_TRUNCATED_H_NAMES
{"GLU: HT ", " N  : NT ", NULL},
{"GLU: HB ", " CB ", NULL},
{"GLU: HG ", " CG ", NULL},
#endif

{"PHE: H  ", " N  : NT ", NULL},
{"PHE: HN ", " N  : NT ", NULL},
{"PHE:1H  ", " N  : NT ", NULL},
{"PHE:2H  ", " N  : NT ", NULL},
{"PHE:3H  ", " N  : NT ", NULL},
{"PHE: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"PHE: H2 ", " N  : NT ", NULL}, /*  */
{"PHE: H3 ", " N  : NT ", NULL}, /*  */
{"PHE: HT1", " N  : NT ", NULL},
{"PHE: HT2", " N  : NT ", NULL},
{"PHE: HT3", " N  : NT ", NULL},
{"PHE: HA ", " CA ", NULL},
{"PHE:1HB ", " CB ", NULL},
{"PHE:2HB ", " CB ", NULL},
{"PHE: HB1", " CB ", NULL},
{"PHE: HB2", " CB ", NULL},
/* {"PHE: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"PHE: HB3", " CB ", NULL}, /*  */
{"PHE: HD1", " CD1", NULL},
{"PHE: HD2", " CD2", NULL},
{"PHE: HE1", " CE1", NULL},
{"PHE: HE2", " CE2", NULL},
{"PHE: HZ ", " CZ ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"PHE: HT ", " N  : NT ", NULL},
{"PHE: HB ", " CB ", NULL},
{"PHE: HD ", " CD1: CD2", NULL},
{"PHE: HE ", " CE1: CE2", NULL},
#endif

{"GLY: H  ", " N  : NT ", NULL},
{"GLY: HN ", " N  : NT ", NULL},
{"GLY:1H  ", " N  : NT ", NULL},
{"GLY:2H  ", " N  : NT ", NULL},
{"GLY:3H  ", " N  : NT ", NULL},
{"GLY: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"GLY: H2 ", " N  : NT ", NULL}, /*  */
{"GLY: H3 ", " N  : NT ", NULL}, /*  */
{"GLY: HT1", " N  : NT ", NULL},
{"GLY: HT2", " N  : NT ", NULL},
{"GLY: HT3", " N  : NT ", NULL},
{"GLY:1HA ", " CA ", NULL},
{"GLY:2HA ", " CA ", NULL},
{"GLY: HA1", " CA ", NULL},
{"GLY: HA2", " CA ", NULL},
/* {"GLY: HA2", " CA ", NULL},  */ /* remediated names RMI 070718 */
{"GLY: HA3", " CA ", NULL}, /*  */
#ifdef ALLOW_TRUNCATED_H_NAMES
{"GLY: HT ", " N  : NT ", NULL},
{"GLY: HA ", " CA ", NULL},
#endif

{"HIS: H  ", " N  : NT ", NULL},
{"HIS: HN ", " N  : NT ", NULL},
{"HIS:1H  ", " N  : NT ", NULL},
{"HIS:2H  ", " N  : NT ", NULL},
{"HIS:3H  ", " N  : NT ", NULL},
{"HIS: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"HIS: H2 ", " N  : NT ", NULL}, /*  */
{"HIS: H3 ", " N  : NT ", NULL}, /*  */
{"HIS: HT1", " N  : NT ", NULL},
{"HIS: HT2", " N  : NT ", NULL},
{"HIS: HT3", " N  : NT ", NULL},
{"HIS: HA ", " CA ", NULL},
{"HIS:1HB ", " CB ", NULL},
{"HIS:2HB ", " CB ", NULL},
{"HIS: HB1", " CB ", NULL},
{"HIS: HB2", " CB ", NULL},
/* {"HIS: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"HIS: HB3", " CB ", NULL}, /*  */
{"HIS: HD1", " ND1", NULL},
{"HIS: HD2", " CD2", NULL},
{"HIS: HE1", " CE1", NULL},
{"HIS: HE2", " NE2", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"HIS: HT ", " N  : NT ", NULL},
{"HIS: HB ", " CB ", NULL},
{"HIS: HD ", " ND1: CD2", NULL},
{"HIS: HE ", " CE1: NE2", NULL},
#endif

{"ILE: H  ", " N  : NT ", NULL},
{"ILE: HN ", " N  : NT ", NULL},
{"ILE:1H  ", " N  : NT ", NULL},
{"ILE:2H  ", " N  : NT ", NULL},
{"ILE:3H  ", " N  : NT ", NULL},
{"ILE: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"ILE: H2 ", " N  : NT ", NULL}, /*  */
{"ILE: H3 ", " N  : NT ", NULL}, /*  */
{"ILE: HT1", " N  : NT ", NULL},
{"ILE: HT2", " N  : NT ", NULL},
{"ILE: HT3", " N  : NT ", NULL},
{"ILE: HA ", " CA ", NULL},
{"ILE: HB ", " CB ", NULL},
{"ILE:1HG1", " CG1", NULL},
{"ILE:2HG1", " CG1", NULL},
{"ILE:HG11", " CG1", NULL},
{"ILE:HG12", " CG1", NULL},
{"ILE:1HG2", " CG2", NULL},
{"ILE:2HG2", " CG2", NULL},
{"ILE:3HG2", " CG2", NULL},
{"ILE:HG21", " CG2", NULL},
{"ILE:HG22", " CG2", NULL},
{"ILE:HG23", " CG2", NULL},
{"ILE:1HD1", " CD1", NULL},
{"ILE:2HD1", " CD1", NULL},
{"ILE:3HD1", " CD1", NULL},
{"ILE:HD11", " CD1", NULL},
{"ILE:HD12", " CD1", NULL},
{"ILE:HD13", " CD1", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"ILE: HT ", " N  : NT ", NULL},
{"ILE:1HG ", " CG1: CG2", NULL},
{"ILE:2HG ", " CG1: CG2", NULL},
{"ILE:3HG ", " CG2", NULL},
{"ILE:HG1 ", " CG1", NULL},
{"ILE:HG2 ", " CG2", NULL},
{"ILE:1HD ", " CD1", NULL},
{"ILE:2HD ", " CD1", NULL},
{"ILE:3HD ", " CD1", NULL},
{"ILE:HD1 ", " CD1", NULL},
#endif

{"LYS: H  ", " N  : NT ", NULL},
{"LYS: HN ", " N  : NT ", NULL},
{"LYS:1H  ", " N  : NT ", NULL},
{"LYS:2H  ", " N  : NT ", NULL},
{"LYS:3H  ", " N  : NT ", NULL},
{"LYS: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"LYS: H2 ", " N  : NT ", NULL}, /*  */
{"LYS: H3 ", " N  : NT ", NULL}, /*  */
{"LYS: HT1", " N  : NT ", NULL},
{"LYS: HT2", " N  : NT ", NULL},
{"LYS: HT3", " N  : NT ", NULL},
{"LYS: HA ", " CA ", NULL},
{"LYS:1HB ", " CB ", NULL},
{"LYS:2HB ", " CB ", NULL},
{"LYS: HB1", " CB ", NULL},
{"LYS: HB2", " CB ", NULL},
/* {"LYS: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"LYS: HB3", " CB ", NULL}, /*  */
{"LYS:1HG ", " CG ", NULL},
{"LYS:2HG ", " CG ", NULL},
{"LYS: HG1", " CG ", NULL},
{"LYS: HG2", " CG ", NULL},
/* {"LYS: HG2", " CG ", NULL},  */ /* remediated names RMI 070718 */
{"LYS: HG3", " CG ", NULL}, /*  */
{"LYS:1HD ", " CD ", NULL},
{"LYS:2HD ", " CD ", NULL},
{"LYS: HD1", " CD ", NULL},
{"LYS: HD2", " CD ", NULL},
/* {"LYS: HD2", " CD ", NULL},  */ /* remediated names RMI 070718 */
{"LYS: HD3", " CD ", NULL}, /*  */
{"LYS:1HE ", " CE ", NULL},
{"LYS:2HE ", " CE ", NULL},
{"LYS: HE1", " CE ", NULL},
{"LYS: HE2", " CE ", NULL},
/* {"LYS: HE2", " CE ", NULL},  */ /* remediated names RMI 070718 */
{"LYS: HE3", " CE ", NULL}, /*  */
{"LYS:1HZ ", " NZ ", NULL},
{"LYS:2HZ ", " NZ ", NULL},
{"LYS:3HZ ", " NZ ", NULL},
{"LYS: HZ1", " NZ ", NULL},
{"LYS: HZ2", " NZ ", NULL},
{"LYS: HZ3", " NZ ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"LYS: HT ", " N  : NT ", NULL},
{"LYS: HB ", " CB ", NULL},
{"LYS: HG ", " CG ", NULL},
{"LYS: HD ", " CD ", NULL},
{"LYS: HE ", " CE ", NULL},
{"LYS: HZ ", " NZ ", NULL},
#endif

{"LEU: H  ", " N  : NT ", NULL},
{"LEU: HN ", " N  : NT ", NULL},
{"LEU:1H  ", " N  : NT ", NULL},
{"LEU:2H  ", " N  : NT ", NULL},
{"LEU:3H  ", " N  : NT ", NULL},
{"LEU: HT1", " N  : NT ", NULL},
{"LEU: HT2", " N  : NT ", NULL},
{"LEU: HT3", " N  : NT ", NULL},
{"LEU: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"LEU: H2 ", " N  : NT ", NULL}, /*  */
{"LEU: H3 ", " N  : NT ", NULL}, /*  */
{"LEU: HA ", " CA ", NULL},
{"LEU:1HB ", " CB ", NULL},
{"LEU:2HB ", " CB ", NULL},
{"LEU: HB1", " CB ", NULL},
{"LEU: HB2", " CB ", NULL},
/* {"LEU: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"LEU: HB3", " CB ", NULL}, /*  */
{"LEU: HG ", " CG ", NULL},
{"LEU:1HD1", " CD1", NULL},
{"LEU:2HD1", " CD1", NULL},
{"LEU:3HD1", " CD1", NULL},
{"LEU:HD11", " CD1", NULL},
{"LEU:HD12", " CD1", NULL},
{"LEU:HD13", " CD1", NULL},
{"LEU:1HD2", " CD2", NULL},
{"LEU:2HD2", " CD2", NULL},
{"LEU:3HD2", " CD2", NULL},
{"LEU:HD21", " CD2", NULL},
{"LEU:HD22", " CD2", NULL},
{"LEU:HD23", " CD2", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"LEU: HT ", " N  : NT ", NULL},
{"LEU: HB ", " CB ", NULL},
{"LEU:1HD ", " CD1: CD2", NULL},
{"LEU:2HD ", " CD1: CD2", NULL},
{"LEU:3HD ", " CD1: CD2", NULL},
{"LEU:HD1 ", " CD1", NULL},
{"LEU:HD2 ", " CD2", NULL},
#endif

{"MET: H  ", " N  : NT ", NULL},
{"MET: HN ", " N  : NT ", NULL},
{"MET:1H  ", " N  : NT ", NULL},
{"MET:2H  ", " N  : NT ", NULL},
{"MET:3H  ", " N  : NT ", NULL},
{"MET: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"MET: H2 ", " N  : NT ", NULL}, /*  */
{"MET: H3 ", " N  : NT ", NULL}, /*  */
{"MET: HT1", " N  : NT ", NULL},
{"MET: HT2", " N  : NT ", NULL},
{"MET: HT3", " N  : NT ", NULL},
{"MET: HA ", " CA ", NULL},
{"MET:1HB ", " CB ", NULL},
{"MET:2HB ", " CB ", NULL},
{"MET: HB1", " CB ", NULL},
{"MET: HB2", " CB ", NULL},
/* {"MET: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"MET: HB3", " CB ", NULL}, /*  */
{"MET:1HG ", " CG ", NULL},
{"MET:2HG ", " CG ", NULL},
{"MET: HG1", " CG ", NULL},
{"MET: HG2", " CG ", NULL},
/* {"MET: HG2", " CG ", NULL},  */ /* remediated names RMI 070718 */
{"MET: HG3", " CG ", NULL}, /*  */
{"MET:1HE ", " CE ", NULL},
{"MET:2HE ", " CE ", NULL},
{"MET:3HE ", " CE ", NULL},
{"MET: HE1", " CE ", NULL},
{"MET: HE2", " CE ", NULL},
{"MET: HE3", " CE ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"MET: HT ", " N  : NT ", NULL},
{"MET: HB ", " CB ", NULL},
{"MET: HG ", " CG ", NULL},
{"MET: HE ", " CE ", NULL},
#endif

{"ASN: H  ", " N  : NT ", NULL},
{"ASN: HN ", " N  : NT ", NULL},
{"ASN:1H  ", " N  : NT ", NULL},
{"ASN:2H  ", " N  : NT ", NULL},
{"ASN:3H  ", " N  : NT ", NULL},
{"ASN: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"ASN: H2 ", " N  : NT ", NULL}, /*  */
{"ASN: H3 ", " N  : NT ", NULL}, /*  */
{"ASN: HT1", " N  : NT ", NULL},
{"ASN: HT2", " N  : NT ", NULL},
{"ASN: HT3", " N  : NT ", NULL},
{"ASN: HA ", " CA ", NULL},
{"ASN:1HB ", " CB ", NULL},
{"ASN:2HB ", " CB ", NULL},
{"ASN: HB1", " CB ", NULL},
{"ASN: HB2", " CB ", NULL},
/* {"ASN: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"ASN: HB3", " CB ", NULL}, /*  */
{"ASN:1HD2", " ND2", NULL},
{"ASN:2HD2", " ND2", NULL},
{"ASN:HD21", " ND2", NULL},
{"ASN:HD22", " ND2", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"ASN: HT ", " N  : NT ", NULL},
{"ASN: HB ", " CB ", NULL},
{"ASN:1HD ", " ND2", NULL},
{"ASN:2HD ", " ND2", NULL},
{"ASN:HD2 ", " ND2", NULL},
#endif

{"PRO:1H  ", " N  : NT ", NULL},
{"PRO:2H  ", " N  : NT ", NULL},
{"PRO: HT1", " N  : NT ", NULL},
{"PRO: HT2", " N  : NT ", NULL},
{"PRO: HA ", " CA ", NULL},
{"PRO:1HB ", " CB ", NULL},
{"PRO:2HB ", " CB ", NULL},
{"PRO: HB1", " CB ", NULL},
{"PRO: HB2", " CB ", NULL},
/* {"PRO: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"PRO: HB3", " CB ", NULL}, /*  */
{"PRO:1HG ", " CG ", NULL},
{"PRO:2HG ", " CG ", NULL},
{"PRO: HG1", " CG ", NULL},
{"PRO: HG2", " CG ", NULL},
/* {"PRO: HG2", " CG ", NULL},  */ /* remediated names RMI 070718 */
{"PRO: HG3", " CG ", NULL}, /*  */
{"PRO:1HD ", " CD ", NULL},
{"PRO:2HD ", " CD ", NULL},
{"PRO: HD1", " CD ", NULL},
{"PRO: HD2", " CD ", NULL},
/* {"PRO: HD2", " CD ", NULL},  */ /* remediated names RMI 070718 */
{"PRO: HD3", " CD ", NULL}, /*  */
#ifdef ALLOW_TRUNCATED_H_NAMES
{"PRO: HT ", " N  : NT ", NULL},
{"PRO: HB ", " CB ", NULL},
{"PRO: HG ", " CG ", NULL},
{"PRO: HD ", " CD ", NULL},
#endif

{"GLN: H  ", " N  : NT ", NULL},
{"GLN: HN ", " N  : NT ", NULL},
{"GLN:1H  ", " N  : NT ", NULL},
{"GLN:2H  ", " N  : NT ", NULL},
{"GLN:3H  ", " N  : NT ", NULL},
{"GLN: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"GLN: H2 ", " N  : NT ", NULL}, /*  */
{"GLN: H3 ", " N  : NT ", NULL}, /*  */
{"GLN: HT1", " N  : NT ", NULL},
{"GLN: HT2", " N  : NT ", NULL},
{"GLN: HT3", " N  : NT ", NULL},
{"GLN: HA ", " CA ", NULL},
{"GLN:1HB ", " CB ", NULL},
{"GLN:2HB ", " CB ", NULL},
{"GLN: HB1", " CB ", NULL},
{"GLN: HB2", " CB ", NULL},
/* {"GLN: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"GLN: HB3", " CB ", NULL}, /*  */
{"GLN:1HG ", " CG ", NULL},
{"GLN:2HG ", " CG ", NULL},
{"GLN: HG1", " CG ", NULL},
{"GLN: HG2", " CG ", NULL},
/* {"GLN: HG2", " CG ", NULL},  */ /* remediated names RMI 070718 */
{"GLN: HG3", " CG ", NULL}, /*  */
{"GLN:1HE2", " NE2", NULL},
{"GLN:2HE2", " NE2", NULL},
{"GLN:HE21", " NE2", NULL},
{"GLN:HE22", " NE2", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"GLN: HT ", " N  : NT ", NULL},
{"GLN: HB ", " CB ", NULL},
{"GLN: HG ", " CG ", NULL},
{"GLN:1HE ", " NE2", NULL},
{"GLN:2HE ", " NE2", NULL},
{"GLN:HE2 ", " NE2", NULL},
#endif

{"ARG: H  ", " N  : NT ", NULL},
{"ARG: HN ", " N  : NT ", NULL},
{"ARG:1H  ", " N  : NT ", NULL},
{"ARG:2H  ", " N  : NT ", NULL},
{"ARG:3H  ", " N  : NT ", NULL},
{"ARG: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"ARG: H2 ", " N  : NT ", NULL}, /*  */
{"ARG: H3 ", " N  : NT ", NULL}, /*  */
{"ARG: HT1", " N  : NT ", NULL},
{"ARG: HT2", " N  : NT ", NULL},
{"ARG: HT3", " N  : NT ", NULL},
{"ARG: HA ", " CA ", NULL},
{"ARG:1HB ", " CB ", NULL},
{"ARG:2HB ", " CB ", NULL},
{"ARG: HB1", " CB ", NULL},
{"ARG: HB2", " CB ", NULL},
/* {"ARG: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"ARG: HB3", " CB ", NULL}, /*  */
{"ARG:1HG ", " CG ", NULL},
{"ARG:2HG ", " CG ", NULL},
{"ARG: HG1", " CG ", NULL},
{"ARG: HG2", " CG ", NULL},
/* {"ARG: HG2", " CG ", NULL},  */ /* remediated names RMI 070718 */
{"ARG: HG3", " CG ", NULL}, /*  */
{"ARG:1HD ", " CD ", NULL},
{"ARG:2HD ", " CD ", NULL},
{"ARG: HD1", " CD ", NULL},
{"ARG: HD2", " CD ", NULL},
/* {"ARG: HD2", " CD ", NULL},  */ /* remediated names RMI 070718 */
{"ARG: HD3", " CD ", NULL}, /*  */
{"ARG: HE ", " NE ", NULL},
{"ARG:1HH1", " NH1", NULL},
{"ARG:2HH1", " NH1", NULL},
{"ARG:1HH2", " NH2", NULL},
{"ARG:2HH2", " NH2", NULL},
{"ARG:HH11", " NH1", NULL},
{"ARG:HH12", " NH1", NULL},
{"ARG:HH21", " NH2", NULL},
{"ARG:HH22", " NH2", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"ARG: HT ", " N  : NT ", NULL},
{"ARG: HB ", " CB ", NULL},
{"ARG: HG ", " CG ", NULL},
{"ARG: HD ", " CD ", NULL},
{"ARG:1HH ", " NH1: NH2", NULL},
{"ARG:2HH ", " NH1: NH2", NULL},
{"ARG:HH1 ", " NH1", NULL},
{"ARG:HH2 ", " NH2", NULL},
#endif

{"SER: H  ", " N  : NT ", NULL},
{"SER: HN ", " N  : NT ", NULL},
{"SER:1H  ", " N  : NT ", NULL},
{"SER:2H  ", " N  : NT ", NULL},
{"SER:3H  ", " N  : NT ", NULL},
{"SER: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"SER: H2 ", " N  : NT ", NULL}, /*  */
{"SER: H3 ", " N  : NT ", NULL}, /*  */
{"SER: HT1", " N  : NT ", NULL},
{"SER: HT2", " N  : NT ", NULL},
{"SER: HT3", " N  : NT ", NULL},
{"SER: HA ", " CA ", NULL},
{"SER:1HB ", " CB ", NULL},
{"SER:2HB ", " CB ", NULL},
{"SER: HB1", " CB ", NULL},
{"SER: HB2", " CB ", NULL},
/* {"SER: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"SER: HB3", " CB ", NULL}, /*  */
{"SER: HG ", " OG ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"SER: HT ", " N  : NT ", NULL},
{"SER: HB ", " CB ", NULL},
#endif

{"THR: H  ", " N  : NT ", NULL},
{"THR: HN ", " N  : NT ", NULL},
{"THR:1H  ", " N  : NT ", NULL},
{"THR:2H  ", " N  : NT ", NULL},
{"THR:3H  ", " N  : NT ", NULL},
{"THR: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"THR: H2 ", " N  : NT ", NULL}, /*  */
{"THR: H3 ", " N  : NT ", NULL}, /*  */
{"THR: HT1", " N  : NT ", NULL},
{"THR: HT2", " N  : NT ", NULL},
{"THR: HT3", " N  : NT ", NULL},
{"THR: HA ", " CA ", NULL},
{"THR: HB ", " CB ", NULL},
{"THR: HG1", " OG1", NULL},
{"THR: HG ", " OG1", NULL},
{"THR:1HG2", " CG2", NULL},
{"THR:2HG2", " CG2", NULL},
{"THR:3HG2", " CG2", NULL},
{"THR:HG21", " CG2", NULL},
{"THR:HG22", " CG2", NULL},
{"THR:HG23", " CG2", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"THR: HT ", " N  : NT ", NULL},
{"THR:1HG ", " CG2", NULL},
{"THR:2HG ", " CG2", NULL},
{"THR:3HG ", " CG2", NULL},
{"THR:HG2 ", " CG2", NULL},
#endif

{"VAL: H  ", " N  : NT ", NULL},
{"VAL: HN ", " N  : NT ", NULL},
{"VAL:1H  ", " N  : NT ", NULL},
{"VAL:2H  ", " N  : NT ", NULL},
{"VAL:3H  ", " N  : NT ", NULL},
{"VAL: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"VAL: H2 ", " N  : NT ", NULL}, /*  */
{"VAL: H3 ", " N  : NT ", NULL}, /*  */
{"VAL: HT1", " N  : NT ", NULL},
{"VAL: HT2", " N  : NT ", NULL},
{"VAL: HT3", " N  : NT ", NULL},
{"VAL: HA ", " CA ", NULL},
{"VAL: HB ", " CB ", NULL},
{"VAL:1HG1", " CG1", NULL},
{"VAL:2HG1", " CG1", NULL},
{"VAL:3HG1", " CG1", NULL},
{"VAL:HG11", " CG1", NULL},
{"VAL:HG12", " CG1", NULL},
{"VAL:HG13", " CG1", NULL},
{"VAL:1HG2", " CG2", NULL},
{"VAL:2HG2", " CG2", NULL},
{"VAL:3HG2", " CG2", NULL},
{"VAL:HG21", " CG2", NULL},
{"VAL:HG22", " CG2", NULL},
{"VAL:HG23", " CG2", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"VAL: HT ", " N  : NT ", NULL},
{"VAL:1HG ", " CG1: CG2", NULL},
{"VAL:2HG ", " CG1: CG2", NULL},
{"VAL:3HG ", " CG1: CG2", NULL},
{"VAL:HG1 ", " CG1", NULL},
{"VAL:HG2 ", " CG2", NULL},
#endif

{"TRP: H  ", " N  : NT ", NULL},
{"TRP: HN ", " N  : NT ", NULL},
{"TRP:1H  ", " N  : NT ", NULL},
{"TRP:2H  ", " N  : NT ", NULL},
{"TRP:3H  ", " N  : NT ", NULL},
{"TRP: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"TRP: H2 ", " N  : NT ", NULL}, /*  */
{"TRP: H3 ", " N  : NT ", NULL}, /*  */
{"TRP: HT1", " N  : NT ", NULL},
{"TRP: HT2", " N  : NT ", NULL},
{"TRP: HT3", " N  : NT ", NULL},
{"TRP: HA ", " CA ", NULL},
{"TRP:1HB ", " CB ", NULL},
{"TRP:2HB ", " CB ", NULL},
{"TRP: HB1", " CB ", NULL},
{"TRP: HB2", " CB ", NULL},
/* {"TRP: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"TRP: HB3", " CB ", NULL}, /*  */
{"TRP: HD1", " CD1", NULL},
{"TRP: HE1", " NE1", NULL},
{"TRP: HE3", " CE3", NULL},
{"TRP: HZ2", " CZ2", NULL},
{"TRP: HZ3", " CZ3", NULL},
{"TRP: HH2", " CH2", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"TRP: HT ", " N  : NT ", NULL},
{"TRP: HB ", " CB ", NULL},
{"TRP: HD ", " CD1", NULL},
{"TRP: HE ", " NE1: CE3", NULL},
{"TRP: HZ ", " CZ2: CZ3", NULL},
{"TRP: HH ", " CH2", NULL},
#endif

{"TYR: H  ", " N  : NT ", NULL},
{"TYR: HN ", " N  : NT ", NULL},
{"TYR:1H  ", " N  : NT ", NULL},
{"TYR:2H  ", " N  : NT ", NULL},
{"TYR:3H  ", " N  : NT ", NULL},
{"TYR: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"TYR: H2 ", " N  : NT ", NULL}, /*  */
{"TYR: H3 ", " N  : NT ", NULL}, /*  */
{"TYR: HT1", " N  : NT ", NULL},
{"TYR: HT2", " N  : NT ", NULL},
{"TYR: HT3", " N  : NT ", NULL},
{"TYR: HA ", " CA ", NULL},
{"TYR:1HB ", " CB ", NULL},
{"TYR:2HB ", " CB ", NULL},
{"TYR: HB1", " CB ", NULL},
{"TYR: HB2", " CB ", NULL},
/* {"TYR: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"TYR: HB3", " CB ", NULL}, /*  */
{"TYR: HD1", " CD1", NULL},
{"TYR: HD2", " CD2", NULL},
{"TYR: HE1", " CE1", NULL},
{"TYR: HE2", " CE2", NULL},
{"TYR: HH ", " OH ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"TYR: HT ", " N  : NT ", NULL},
{"TYR: HB ", " CB ", NULL},
{"TYR: HD ", " CD1: CD2", NULL},
{"TYR: HE ", " CE1: CE2", NULL},
#endif

/*---------------------------------------*/

{"ASX: H  ", " N  : NT ", NULL},
{"ASX: HN ", " N  : NT ", NULL},
{"ASX:1H  ", " N  : NT ", NULL},
{"ASX:2H  ", " N  : NT ", NULL},
{"ASX:3H  ", " N  : NT ", NULL},
{"ASX: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"ASX: H2 ", " N  : NT ", NULL}, /*  */
{"ASX: H3 ", " N  : NT ", NULL}, /*  */
{"ASX: HT1", " N  : NT ", NULL},
{"ASX: HT2", " N  : NT ", NULL},
{"ASX: HT3", " N  : NT ", NULL},
{"ASX: HA ", " CA ", NULL},
{"ASX:1HB ", " CB ", NULL},
{"ASX:2HB ", " CB ", NULL},
{"ASX: HB1", " CB ", NULL},
{"ASX: HB2", " CB ", NULL},
/* {"ASX: HB2", " CB ", NULL},  */ /* remediated names RMI 070718  */
{"ASX: HB3", " CB ", NULL}, /*  */
#ifdef ALLOW_TRUNCATED_H_NAMES
{"ASX: HT ", " N  : NT ", NULL},
{"ASX: HB ", " CB ", NULL},
#endif

{"GLX: H  ", " N  : NT ", NULL},
{"GLX: HN ", " N  : NT ", NULL},
{"GLX:1H  ", " N  : NT ", NULL},
{"GLX:2H  ", " N  : NT ", NULL},
{"GLX:3H  ", " N  : NT ", NULL},
{"GLX: H1 ", " N  : NT ", NULL}, /* remediated names RMI 070718 */
{"GLX: H2 ", " N  : NT ", NULL}, /*  */
{"GLX: H3 ", " N  : NT ", NULL}, /*  */
{"GLX: HT1", " N  : NT ", NULL},
{"GLX: HT2", " N  : NT ", NULL},
{"GLX: HT3", " N  : NT ", NULL},
{"GLX: HA ", " CA ", NULL},
{"GLX:1HB ", " CB ", NULL},
{"GLX:2HB ", " CB ", NULL},
{"GLX: HB1", " CB ", NULL},
{"GLX: HB2", " CB ", NULL},
/* {"GLX: HB2", " CB ", NULL},  */ /* remediated names RMI 070718 */
{"GLX: HB3", " CB ", NULL}, /*  */
{"GLX:1HG ", " CG ", NULL},
{"GLX:2HG ", " CG ", NULL},
{"GLX: HG1", " CG ", NULL},
{"GLX: HG2", " CG ", NULL},
/* {"GLX: HG2", " CG ", NULL},  */ /* remediated names RMI 070718 */
{"GLX: HG3", " CG ", NULL}, /*  */
#ifdef ALLOW_TRUNCATED_H_NAMES
{"GLX: HT ", " N  : NT ", NULL},
{"GLX: HB ", " CB ", NULL},
{"GLX: HG ", " CG ", NULL},
#endif

/*---------------------------------------*/

{"  U: H1*", " C1*: C1'", NULL},
{"  U: H1'", " C1*: C1'", NULL},
{"  U: H3*", " C3*: C3'", NULL},
{"  U: H3'", " C3*: C3'", NULL},
{"  U: H4*", " C4*: C4'", NULL},
{"  U: H4'", " C4*: C4'", NULL},
{"  U:1H5*", " C5*: C5'", NULL},
{"  U:2H5*", " C5*: C5'", NULL},
{"  U:*H51", " C5*: C5'", NULL},
{"  U:*H52", " C5*: C5'", NULL},
{"  U:H5''", " C5*: C5'", NULL},
{"  U: H5'", " C5*: C5'", NULL},
{"  U:1H2*", " C2*: C2'", NULL},
{"  U:2H2*", " C2*: C2'", NULL},
{"  U: H2*", " C2*: C2'", NULL},
{"  U:H2''", " C2*: C2'", NULL},
{"  U: H2'", " C2*: C2'", NULL},
{"  U:2HO*", " O2*: O2'", NULL},
{"  U:*HO2", " O2*: O2'", NULL},
{"  U:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{"  U: H3T", " O3*: O3'", NULL},
{"  U: H5T", " O5*: O5'", NULL},
{"  U:3HO*", " O3*: O3'", NULL},
{"  U:*HO3", " O3*: O3'", NULL},
{"  U:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{"  U:5HO*", " O5*: O5'", NULL},
{"  U:*HO5", " O5*: O5'", NULL},
{"  U:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{"  U: HN3", " N3 ", NULL},
/* {"  U: H3 ", " N3 ", NULL},  */ /* remediated names RMI 070718 */
{"  U: H6 ", " C6 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"  U: H1 ", " C1*: C1'", NULL},
{"  U: H4 ", " C4*: C4'", NULL},
{"  U:1H5 ", " C5*: C5'", NULL},
{"  U:2H5 ", " C5*: C5'", NULL},
{"  U:*H5 ", " C5*", NULL},
{"  U:H5' ", " C5'", NULL},
{"  U:1H2 ", " C2*: C2'", NULL},
{"  U:2H2 ", " C2*: C2'", NULL},
{"  U: H2 ", " C2*: C2': O2*: O2'", NULL},
{"  U:H2' ", " C2'", NULL},
{"  U:2HO ", " O2*: O2'", NULL},
{"  U:*HO ", " O2*: O3*: O5*", NULL},
{"  U:3HO ", " O3*: O3'", NULL},
{"  U:5HO ", " O5*: O5'", NULL},

{"  U: HN ", " N3 ", NULL},

{"  U: H3 ", " N3 : C3*: C3': O3*: O3'", NULL},
{"  U: H5 ", " C5 : C5': C5*: O5*: O5'", NULL},
#else
{"  U: H3 ", " N3 ", NULL},
{"  U: H5 ", " C5 ", NULL},
#endif

{"  T: H1*", " C1*: C1'", NULL},
{"  T: H1'", " C1*: C1'", NULL},
{"  T: H3*", " C3*: C3'", NULL},
{"  T: H3'", " C3*: C3'", NULL},
{"  T: H4*", " C4*: C4'", NULL},
{"  T: H4'", " C4*: C4'", NULL},
{"  T:1H5*", " C5*: C5'", NULL},
{"  T:2H5*", " C5*: C5'", NULL},
{"  T:*H51", " C5*: C5'", NULL},
{"  T:*H52", " C5*: C5'", NULL},
{"  T:H5''", " C5*: C5'", NULL},
{"  T: H5'", " C5*: C5'", NULL},
{"  T:1H2*", " C2*: C2'", NULL},
{"  T:2H2*", " C2*: C2'", NULL},
{"  T: H2*", " C2*: C2'", NULL},
{"  T:H2''", " C2*: C2'", NULL},
{"  T: H2'", " C2*: C2'", NULL},
{"  T:2HO*", " O2*: O2'", NULL},
{"  T:*HO2", " O2*: O2'", NULL},
{"  T:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{"  T: H3T", " O3*: O3'", NULL},
{"  T: H5T", " O5*: O5'", NULL},
{"  T:3HO*", " O3*: O3'", NULL},
{"  T:*HO3", " O3*: O3'", NULL},
{"  T:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{"  T:5HO*", " O5*: O5'", NULL},
{"  T:*HO5", " O5*: O5'", NULL},
{"  T:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{"  T: HN3", " N3 ", NULL},
{"  T:1H5M", " C7 : C5A: CA5: C5M", NULL},
{"  T:2H5M", " C7 : C5A: CA5: C5M", NULL},
{"  T:3H5M", " C7 : C5A: CA5: C5M", NULL},
{"  T:1HM5", " C7 : C5A: CA5: C5M", NULL},
{"  T:2HM5", " C7 : C5A: CA5: C5M", NULL},
{"  T:3HM5", " C7 : C5A: CA5: C5M", NULL},
{"  T:H5M1", " C7 : C5A: CA5: C5M", NULL},
{"  T:H5M2", " C7 : C5A: CA5: C5M", NULL},
{"  T:H5M3", " C7 : C5A: CA5: C5M", NULL},
{"  T:1H5A", " C7 : C5A: CA5: C5M", NULL},
{"  T:2H5A", " C7 : C5A: CA5: C5M", NULL},
{"  T:3H5A", " C7 : C5A: CA5: C5M", NULL},
{"  T:1HA5", " C7 : C5A: CA5: C5M", NULL},
{"  T:2HA5", " C7 : C5A: CA5: C5M", NULL},
{"  T:3HA5", " C7 : C5A: CA5: C5M", NULL},
{"  T:H5A1", " C7 : C5A: CA5: C5M", NULL},
{"  T:H5A2", " C7 : C5A: CA5: C5M", NULL},
{"  T:H5A3", " C7 : C5A: CA5: C5M", NULL},
{"  T: H71", " C7 : C5A: CA5: C5M", NULL}, /* remediated names RMI 070718 */
{"  T: H72", " C7 : C5A: CA5: C5M", NULL}, /*  */
{"  T: H73", " C7 : C5A: CA5: C5M", NULL}, /*  */
{"  T: H6 ", " C6 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"  T: H1 ", " C1*: C1'", NULL},
{"  T: H4 ", " C4*: C4'", NULL},
{"  T:*H5 ", " C5*", NULL},
{"  T:H5' ", " C5'", NULL},
{"  T: H5 ", " C5': C5*: O5*: O5'", NULL},
{"  T:1H2 ", " C2*: C2'", NULL},
{"  T:2H2 ", " C2*: C2'", NULL},
{"  T: H2 ", " C2*: C2': O2*: O2'", NULL},
{"  T:H2' ", " C2'", NULL},
{"  T:2HO ", " O2*: O2'", NULL},
{"  T:*HO ", " O2*: O3*: O5*", NULL},
{"  T:3HO ", " O3*: O3'", NULL},
{"  T:5HO ", " O5*: O5'", NULL},

{"  T: HN ", " N3 ", NULL},
{"  T:1H5 ", " C5M: C5A: C5*: C5'", NULL},
{"  T:2H5 ", " C5M: C5A: C5*: C5'", NULL},
{"  T:3H5 ", " C5M: C5A", NULL},
{"  T:1HM ", " C5M: CM5", NULL},
{"  T:2HM ", " C5M: CM5", NULL},
{"  T:3HM ", " C5M: CM5", NULL},
{"  T:1HA ", " C5A: CA5", NULL},
{"  T:2HA ", " C5A: CA5", NULL},
{"  T:3HA ", " C5A: CA5", NULL},
{"  T:H5M ", " C5M", NULL},
{"  T:H5A ", " C5A", NULL},

{"  T: H3 ", " N3 : C3*: C3': O3*: O3'", NULL},
#else
{"  T: H3 ", " N3 ", NULL},
#endif

{"  A: H1*", " C1*: C1'", NULL},
{"  A: H1'", " C1*: C1'", NULL},
{"  A: H3*", " C3*: C3'", NULL},
{"  A: H3'", " C3*: C3'", NULL},
{"  A: H4*", " C4*: C4'", NULL},
{"  A: H4'", " C4*: C4'", NULL},
{"  A:1H5*", " C5*: C5'", NULL},
{"  A:2H5*", " C5*: C5'", NULL},
{"  A:*H51", " C5*: C5'", NULL},
{"  A:*H52", " C5*: C5'", NULL},
{"  A:H5''", " C5*: C5'", NULL},
{"  A: H5'", " C5*: C5'", NULL},
{"  A:1H2*", " C2*: C2'", NULL},
{"  A:2H2*", " C2*: C2'", NULL},
{"  A: H2*", " C2*: C2'", NULL},
{"  A:H2''", " C2*: C2'", NULL},
{"  A: H2'", " C2*: C2'", NULL},
{"  A:2HO*", " O2*: O2'", NULL},
{"  A:*HO2", " O2*: O2'", NULL},
{"  A:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{"  A: H3T", " O3*: O3'", NULL},
{"  A: H5T", " O5*: O5'", NULL},
{"  A:3HO*", " O3*: O3'", NULL},
{"  A:*HO3", " O3*: O3'", NULL},
{"  A:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{"  A:5HO*", " O5*: O5'", NULL},
{"  A:*HO5", " O5*: O5'", NULL},
{"  A:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{"  A:1H6 ", " N6 ", NULL},
{"  A:2H6 ", " N6 ", NULL},
{"  A:1HN6", " N6 ", NULL},
{"  A:2HN6", " N6 ", NULL},
{"  A: H61", " N6 ", NULL},
{"  A: H62", " N6 ", NULL},
{"  A: H8 ", " C8 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"  A: H1 ", " C1*: C1'", NULL},
{"  A: H3 ", " C3*: C3': O3*: O3'", NULL},
{"  A: H4 ", " C4*: C4'", NULL},
{"  A:1H5 ", " C5*: C5'", NULL},
{"  A:2H5 ", " C5*: C5'", NULL},
{"  A:*H5 ", " C5*: C5'", NULL},
{"  A:H5' ", " C5*: C5*", NULL},
{"  A: H5 ", " C5': C5*: O5*: O5'", NULL},
{"  A:1H2 ", " C2*: C2'", NULL},
{"  A:2H2 ", " C2*: C2'", NULL},
{"  A:H2' ", " C2*: C2'", NULL},
{"  A:2HO ", " O2*: O2'", NULL},
{"  A:*HO ", " O2*: O3*: O5*", NULL},
{"  A:3HO ", " O3*: O3'", NULL},
{"  A:5HO ", " O5*: O5'", NULL},

{"  A:1HN ", " N6 ", NULL},
{"  A:2HN ", " N6 ", NULL},
{"  A: H6 ", " N6 ", NULL},

{"  A: H2 ", " C2 : C2*: C2': O2*: O2'", NULL},
#else
{"  A: H2 ", " C2 ", NULL},
#endif

{"  C: H1*", " C1*: C1'", NULL},
{"  C: H1'", " C1*: C1'", NULL},
{"  C: H3*", " C3*: C3'", NULL},
{"  C: H3'", " C3*: C3'", NULL},
{"  C: H4*", " C4*: C4'", NULL},
{"  C: H4'", " C4*: C4'", NULL},
{"  C:1H5*", " C5*: C5'", NULL},
{"  C:2H5*", " C5*: C5'", NULL},
{"  C:*H51", " C5*: C5'", NULL},
{"  C:*H52", " C5*: C5'", NULL},
{"  C:H5''", " C5*: C5'", NULL},
{"  C: H5'", " C5*: C5'", NULL},
{"  C:1H2*", " C2*: C2'", NULL},
{"  C:2H2*", " C2*: C2'", NULL},
{"  C: H2*", " C2*: C2'", NULL},
{"  C:H2''", " C2*: C2'", NULL},
{"  C: H2'", " C2*: C2'", NULL},
{"  C:2HO*", " O2*: O2'", NULL},
{"  C:*HO2", " O2*: O2'", NULL},
{"  C:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{"  C: H3T", " O3*: O3'", NULL},
{"  C: H5T", " O5*: O5'", NULL},
{"  C:3HO*", " O3*: O3'", NULL},
{"  C:*HO3", " O3*: O3'", NULL},
{"  C:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{"  C:5HO*", " O5*: O5'", NULL},
{"  C:*HO5", " O5*: O5'", NULL},
{"  C:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{"  C:1H4 ", " N4 ", NULL},
{"  C:2H4 ", " N4 ", NULL},
{"  C:1HN4", " N4 ", NULL},
{"  C:2HN4", " N4 ", NULL},
{"  C: H41", " N4 ", NULL},
{"  C: H42", " N4 ", NULL},
{"  C: H6 ", " C6 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"  C: H1 ", " C1*: C1'", NULL},
{"  C: H3 ", " C3*: C3': O3*: O3'", NULL},
{"  C:1H5 ", " C5*: C5'", NULL},
{"  C:2H5 ", " C5*: C5'", NULL},
{"  C:*H5 ", " C5*", NULL},
{"  C:H5' ", " C5'", NULL},
{"  C:1H2 ", " C2*: C2'", NULL},
{"  C:2H2 ", " C2*: C2'", NULL},
{"  C: H2 ", " C2*: C2': O2*: O2'", NULL},
{"  C:H2' ", " C2'", NULL},
{"  C:2HO ", " O2*: O2'", NULL},
{"  C:*HO ", " O2*: O3*: O5*", NULL},
{"  C:3HO ", " O3*: O3'", NULL},
{"  C:5HO ", " O5*: O5'", NULL},

{"  C:1HN ", " N4 ", NULL},
{"  C:2HN ", " N4 ", NULL},
{"  C: H4 ", " N4 : C4*: C4'", NULL},

{"  C: H5 ", " C5 : C5': C5*: O5*: O5'", NULL},
#else
{"  C: H5 ", " C5 ", NULL},
#endif

{"  G: H1*", " C1*: C1'", NULL},
{"  G: H1'", " C1*: C1'", NULL},
{"  G: H3*", " C3*: C3'", NULL},
{"  G: H3'", " C3*: C3'", NULL},
{"  G: H4*", " C4*: C4'", NULL},
{"  G: H4'", " C4*: C4'", NULL},
{"  G:1H5*", " C5*: C5'", NULL},
{"  G:2H5*", " C5*: C5'", NULL},
{"  G:*H51", " C5*: C5'", NULL},
{"  G:*H52", " C5*: C5'", NULL},
{"  G:H5''", " C5*: C5'", NULL},
{"  G: H5'", " C5*: C5'", NULL},
{"  G:1H2*", " C2*: C2'", NULL},
{"  G:2H2*", " C2*: C2'", NULL},
{"  G: H2*", " C2*: C2'", NULL},
{"  G:H2''", " C2*: C2'", NULL},
{"  G: H2'", " C2*: C2'", NULL},
{"  G:2HO*", " O2*: O2'", NULL},
{"  G:*HO2", " O2*: O2'", NULL},
{"  G:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{"  G: H3T", " O3*: O3'", NULL},
{"  G: H5T", " O5*: O5'", NULL},
{"  G:3HO*", " O3*: O3'", NULL},
{"  G:*HO3", " O3*: O3'", NULL},
{"  G:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{"  G:5HO*", " O5*: O5'", NULL},
{"  G:*HO5", " O5*: O5'", NULL},
{"  G:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{"  G: HN1", " N1 ", NULL},
{"  G:1HN2", " N2 ", NULL},
{"  G:2HN2", " N2 ", NULL},
{"  G: H21", " N2 ", NULL},
{"  G: H22", " N2 ", NULL},
{"  G: H8 ", " C8 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{"  G: H3 ", " C3*: C3': O3*: O3'", NULL},
{"  G: H4 ", " C4*: C4'", NULL},
{"  G:1H5 ", " C5*: C5'", NULL},
{"  G:2H5 ", " C5*: C5'", NULL},
{"  G:*H5 ", " C5*", NULL},
{"  G:H5' ", " C5'", NULL},
{"  G: H5 ", " C5': C5*: O5*: O5'", NULL},
{"  G:H2' ", " C2'", NULL},
{"  G:2HO ", " O2*: O2'", NULL},
{"  G:*HO ", " O2*: O3*: O5*", NULL},
{"  G:3HO ", " O3*: O3'", NULL},
{"  G:5HO ", " O5*: O5'", NULL},

{"  G: HN ", " N1 ", NULL},
{"  G:1HN ", " N2 ", NULL},
{"  G:2HN ", " N2 ", NULL},
{"  G: H2 ", " N2 : C2*: C2': O2*: O2'", NULL},

{"  G: H1 ", " N1 : C1*: C1'", NULL},
{"  G:1H2 ", " N2 : C2*: C2'", NULL},
{"  G:2H2 ", " N2 : C2*: C2'", NULL},
#else
{"  G: H1 ", " N1 ", NULL},
{"  G:1H2 ", " N2 ", NULL},
{"  G:2H2 ", " N2 ", NULL},
#endif

/* RNA for Coot and CCP4 Ar, Ur, Tr, Cr, Gr */

{" UR: H1*", " C1*: C1'", NULL},
{" UR: H1'", " C1*: C1'", NULL},
{" UR: H3*", " C3*: C3'", NULL},
{" UR: H3'", " C3*: C3'", NULL},
{" UR: H4*", " C4*: C4'", NULL},
{" UR: H4'", " C4*: C4'", NULL},
{" UR:1H5*", " C5*: C5'", NULL},
{" UR:2H5*", " C5*: C5'", NULL},
{" UR:*H51", " C5*: C5'", NULL},
{" UR:*H52", " C5*: C5'", NULL},
{" UR:H5''", " C5*: C5'", NULL},
{" UR: H5'", " C5*: C5'", NULL},
{" UR:1H2*", " C2*: C2'", NULL},
{" UR:2H2*", " C2*: C2'", NULL},
{" UR: H2*", " C2*: C2'", NULL},
{" UR:H2''", " C2*: C2'", NULL},
{" UR: H2'", " C2*: C2'", NULL},
{" UR:2HO*", " O2*: O2'", NULL},
{" UR:*HO2", " O2*: O2'", NULL},
{" UR:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{" UR: H3T", " O3*: O3'", NULL},
{" UR: H5T", " O5*: O5'", NULL},
{" UR:3HO*", " O3*: O3'", NULL},
{" UR:*HO3", " O3*: O3'", NULL},
{" UR:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{" UR:5HO*", " O5*: O5'", NULL},
{" UR:*HO5", " O5*: O5'", NULL},
{" UR:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{" UR: HN3", " N3 ", NULL},
/* {"  U: H3 ", " N3 ", NULL},  */ /* remediated names RMI 070718 */
{" UR: H6 ", " C6 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" UR: H1 ", " C1*: C1'", NULL},
{" UR: H4 ", " C4*: C4'", NULL},
{" UR:1H5 ", " C5*: C5'", NULL},
{" UR:2H5 ", " C5*: C5'", NULL},
{" UR:*H5 ", " C5*", NULL},
{" UR:H5' ", " C5'", NULL},
{" UR:1H2 ", " C2*: C2'", NULL},
{" UR:2H2 ", " C2*: C2'", NULL},
{" UR: H2 ", " C2*: C2': O2*: O2'", NULL},
{" UR:H2' ", " C2'", NULL},
{" UR:2HO ", " O2*: O2'", NULL},
{" UR:*HO ", " O2*: O3*: O5*", NULL},
{" UR:3HO ", " O3*: O3'", NULL},
{" UR:5HO ", " O5*: O5'", NULL},

{" UR: HN ", " N3 ", NULL},

{" UR: H3 ", " N3 : C3*: C3': O3*: O3'", NULL},
{" UR: H5 ", " C5 : C5': C5*: O5*: O5'", NULL},
#else
{" UR: H3 ", " N3 ", NULL},
{" UR: H5 ", " C5 ", NULL},
#endif

{" TR: H1*", " C1*: C1'", NULL},
{" TR: H1'", " C1*: C1'", NULL},
{" TR: H3*", " C3*: C3'", NULL},
{" TR: H3'", " C3*: C3'", NULL},
{" TR: H4*", " C4*: C4'", NULL},
{" TR: H4'", " C4*: C4'", NULL},
{" TR:1H5*", " C5*: C5'", NULL},
{" TR:2H5*", " C5*: C5'", NULL},
{" TR:*H51", " C5*: C5'", NULL},
{" TR:*H52", " C5*: C5'", NULL},
{" TR:H5''", " C5*: C5'", NULL},
{" TR: H5'", " C5*: C5'", NULL},
{" TR:1H2*", " C2*: C2'", NULL},
{" TR:2H2*", " C2*: C2'", NULL},
{" TR: H2*", " C2*: C2'", NULL},
{" TR:H2''", " C2*: C2'", NULL},
{" TR: H2'", " C2*: C2'", NULL},
{" TR:2HO*", " O2*: O2'", NULL},
{" TR:*HO2", " O2*: O2'", NULL},
{" TR:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{" TR: H3T", " O3*: O3'", NULL},
{" TR: H5T", " O5*: O5'", NULL},
{" TR:3HO*", " O3*: O3'", NULL},
{" TR:*HO3", " O3*: O3'", NULL},
{" TR:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{" TR:5HO*", " O5*: O5'", NULL},
{" TR:*HO5", " O5*: O5'", NULL},
{" TR:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{" TR: HN3", " N3 ", NULL},
{" TR:1H5M", " C7 : C5A: CA5: C5M", NULL},
{" TR:2H5M", " C7 : C5A: CA5: C5M", NULL},
{" TR:3H5M", " C7 : C5A: CA5: C5M", NULL},
{" TR:1HM5", " C7 : C5A: CA5: C5M", NULL},
{" TR:2HM5", " C7 : C5A: CA5: C5M", NULL},
{" TR:3HM5", " C7 : C5A: CA5: C5M", NULL},
{" TR:H5M1", " C7 : C5A: CA5: C5M", NULL},
{" TR:H5M2", " C7 : C5A: CA5: C5M", NULL},
{" TR:H5M3", " C7 : C5A: CA5: C5M", NULL},
{" TR:1H5A", " C7 : C5A: CA5: C5M", NULL},
{" TR:2H5A", " C7 : C5A: CA5: C5M", NULL},
{" TR:3H5A", " C7 : C5A: CA5: C5M", NULL},
{" TR:1HA5", " C7 : C5A: CA5: C5M", NULL},
{" TR:2HA5", " C7 : C5A: CA5: C5M", NULL},
{" TR:3HA5", " C7 : C5A: CA5: C5M", NULL},
{" TR:H5A1", " C7 : C5A: CA5: C5M", NULL},
{" TR:H5A2", " C7 : C5A: CA5: C5M", NULL},
{" TR:H5A3", " C7 : C5A: CA5: C5M", NULL},
{" TR: H71", " C7 : C5A: CA5: C5M", NULL}, /* remediated names RMI 070718 */
{" TR: H72", " C7 : C5A: CA5: C5M", NULL}, /*  */
{" TR: H73", " C7 : C5A: CA5: C5M", NULL}, /*  */
{" TR: H6 ", " C6 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" TR: H1 ", " C1*: C1'", NULL},
{" TR: H4 ", " C4*: C4'", NULL},
{" TR:*H5 ", " C5*", NULL},
{" TR:H5' ", " C5'", NULL},
{" TR: H5 ", " C5': C5*: O5*: O5'", NULL},
{" TR:1H2 ", " C2*: C2'", NULL},
{" TR:2H2 ", " C2*: C2'", NULL},
{" TR: H2 ", " C2*: C2': O2*: O2'", NULL},
{" TR:H2' ", " C2'", NULL},
{" TR:2HO ", " O2*: O2'", NULL},
{" TR:*HO ", " O2*: O3*: O5*", NULL},
{" TR:3HO ", " O3*: O3'", NULL},
{" TR:5HO ", " O5*: O5'", NULL},

{" TR: HN ", " N3 ", NULL},
{" TR:1H5 ", " C5M: C5A: C5*: C5'", NULL},
{" TR:2H5 ", " C5M: C5A: C5*: C5'", NULL},
{" TR:3H5 ", " C5M: C5A", NULL},
{" TR:1HM ", " C5M: CM5", NULL},
{" TR:2HM ", " C5M: CM5", NULL},
{" TR:3HM ", " C5M: CM5", NULL},
{" TR:1HA ", " C5A: CA5", NULL},
{" TR:2HA ", " C5A: CA5", NULL},
{" TR:3HA ", " C5A: CA5", NULL},
{" TR:H5M ", " C5M", NULL},
{" TR:H5A ", " C5A", NULL},

{" TR: H3 ", " N3 : C3*: C3': O3*: O3'", NULL},
#else
{" TR: H3 ", " N3 ", NULL},
#endif

{" AR: H1*", " C1*: C1'", NULL},
{" AR: H1'", " C1*: C1'", NULL},
{" AR: H3*", " C3*: C3'", NULL},
{" AR: H3'", " C3*: C3'", NULL},
{" AR: H4*", " C4*: C4'", NULL},
{" AR: H4'", " C4*: C4'", NULL},
{" AR:1H5*", " C5*: C5'", NULL},
{" AR:2H5*", " C5*: C5'", NULL},
{" AR:*H51", " C5*: C5'", NULL},
{" AR:*H52", " C5*: C5'", NULL},
{" AR:H5''", " C5*: C5'", NULL},
{" AR: H5'", " C5*: C5'", NULL},
{" AR:1H2*", " C2*: C2'", NULL},
{" AR:2H2*", " C2*: C2'", NULL},
{" AR: H2*", " C2*: C2'", NULL},
{" AR:H2''", " C2*: C2'", NULL},
{" AR: H2'", " C2*: C2'", NULL},
{" AR:2HO*", " O2*: O2'", NULL},
{" AR:*HO2", " O2*: O2'", NULL},
{" AR:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{" AR: H3T", " O3*: O3'", NULL},
{" AR: H5T", " O5*: O5'", NULL},
{" AR:3HO*", " O3*: O3'", NULL},
{" AR:*HO3", " O3*: O3'", NULL},
{" AR:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{" AR:5HO*", " O5*: O5'", NULL},
{" AR:*HO5", " O5*: O5'", NULL},
{" AR:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{" AR:1H6 ", " N6 ", NULL},
{" AR:2H6 ", " N6 ", NULL},
{" AR:1HN6", " N6 ", NULL},
{" AR:2HN6", " N6 ", NULL},
{" AR: H61", " N6 ", NULL},
{" AR: H62", " N6 ", NULL},
{" AR: H8 ", " C8 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" AR: H1 ", " C1*: C1'", NULL},
{" AR: H3 ", " C3*: C3': O3*: O3'", NULL},
{" AR: H4 ", " C4*: C4'", NULL},
{" AR:1H5 ", " C5*: C5'", NULL},
{" AR:2H5 ", " C5*: C5'", NULL},
{" AR:*H5 ", " C5*: C5'", NULL},
{" AR:H5' ", " C5*: C5*", NULL},
{" AR: H5 ", " C5': C5*: O5*: O5'", NULL},
{" AR:1H2 ", " C2*: C2'", NULL},
{" AR:2H2 ", " C2*: C2'", NULL},
{" AR:H2' ", " C2*: C2'", NULL},
{" AR:2HO ", " O2*: O2'", NULL},
{" AR:*HO ", " O2*: O3*: O5*", NULL},
{" AR:3HO ", " O3*: O3'", NULL},
{" AR:5HO ", " O5*: O5'", NULL},

{" AR:1HN ", " N6 ", NULL},
{" AR:2HN ", " N6 ", NULL},
{" AR: H6 ", " N6 ", NULL},

{" AR: H2 ", " C2 : C2*: C2': O2*: O2'", NULL},
#else
{" AR: H2 ", " C2 ", NULL},
#endif

{" CR: H1*", " C1*: C1'", NULL},
{" CR: H1'", " C1*: C1'", NULL},
{" CR: H3*", " C3*: C3'", NULL},
{" CR: H3'", " C3*: C3'", NULL},
{" CR: H4*", " C4*: C4'", NULL},
{" CR: H4'", " C4*: C4'", NULL},
{" CR:1H5*", " C5*: C5'", NULL},
{" CR:2H5*", " C5*: C5'", NULL},
{" CR:*H51", " C5*: C5'", NULL},
{" CR:*H52", " C5*: C5'", NULL},
{" CR:H5''", " C5*: C5'", NULL},
{" CR: H5'", " C5*: C5'", NULL},
{" CR:1H2*", " C2*: C2'", NULL},
{" CR:2H2*", " C2*: C2'", NULL},
{" CR: H2*", " C2*: C2'", NULL},
{" CR:H2''", " C2*: C2'", NULL},
{" CR: H2'", " C2*: C2'", NULL},
{" CR:2HO*", " O2*: O2'", NULL},
{" CR:*HO2", " O2*: O2'", NULL},
{" CR:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{" CR: H3T", " O3*: O3'", NULL},
{" CR: H5T", " O5*: O5'", NULL},
{" CR:3HO*", " O3*: O3'", NULL},
{" CR:*HO3", " O3*: O3'", NULL},
{" CR:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{" CR:5HO*", " O5*: O5'", NULL},
{" CR:*HO5", " O5*: O5'", NULL},
{" CR:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{" CR:1H4 ", " N4 ", NULL},
{" CR:2H4 ", " N4 ", NULL},
{" CR:1HN4", " N4 ", NULL},
{" CR:2HN4", " N4 ", NULL},
{" CR: H41", " N4 ", NULL},
{" CR: H42", " N4 ", NULL},
{" CR: H6 ", " C6 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" CR: H1 ", " C1*: C1'", NULL},
{" CR: H3 ", " C3*: C3': O3*: O3'", NULL},
{" CR:1H5 ", " C5*: C5'", NULL},
{" CR:2H5 ", " C5*: C5'", NULL},
{" CR:*H5 ", " C5*", NULL},
{" CR:H5' ", " C5'", NULL},
{" CR:1H2 ", " C2*: C2'", NULL},
{" CR:2H2 ", " C2*: C2'", NULL},
{" CR: H2 ", " C2*: C2': O2*: O2'", NULL},
{" CR:H2' ", " C2'", NULL},
{" CR:2HO ", " O2*: O2'", NULL},
{" CR:*HO ", " O2*: O3*: O5*", NULL},
{" CR:3HO ", " O3*: O3'", NULL},
{" CR:5HO ", " O5*: O5'", NULL},

{" CR:1HN ", " N4 ", NULL},
{" CR:2HN ", " N4 ", NULL},
{" CR: H4 ", " N4 : C4*: C4'", NULL},

{" CR: H5 ", " C5 : C5': C5*: O5*: O5'", NULL},
#else
{" CR: H5 ", " C5 ", NULL},
#endif

{" GR: H1*", " C1*: C1'", NULL},
{" GR: H1'", " C1*: C1'", NULL},
{" GR: H3*", " C3*: C3'", NULL},
{" GR: H3'", " C3*: C3'", NULL},
{" GR: H4*", " C4*: C4'", NULL},
{" GR: H4'", " C4*: C4'", NULL},
{" GR:1H5*", " C5*: C5'", NULL},
{" GR:2H5*", " C5*: C5'", NULL},
{" GR:*H51", " C5*: C5'", NULL},
{" GR:*H52", " C5*: C5'", NULL},
{" GR:H5''", " C5*: C5'", NULL},
{" GR: H5'", " C5*: C5'", NULL},
{" GR:1H2*", " C2*: C2'", NULL},
{" GR:2H2*", " C2*: C2'", NULL},
{" GR: H2*", " C2*: C2'", NULL},
{" GR:H2''", " C2*: C2'", NULL},
{" GR: H2'", " C2*: C2'", NULL},
{" GR:2HO*", " O2*: O2'", NULL},
{" GR:*HO2", " O2*: O2'", NULL},
{" GR:HO2'", " O2*: O2'", NULL}, /* remediated names RMI 070718 */
{" GR: H3T", " O3*: O3'", NULL},
{" GR: H5T", " O5*: O5'", NULL},
{" GR:3HO*", " O3*: O3'", NULL},
{" GR:*HO3", " O3*: O3'", NULL},
{" GR:HO3'", " O3*: O3'", NULL}, /* remediated names RMI 070718 */
{" GR:5HO*", " O5*: O5'", NULL},
{" GR:*HO5", " O5*: O5'", NULL},
{" GR:HO5'", " O5*: O5'", NULL}, /* remediated names RMI 070718 */
{" GR: HN1", " N1 ", NULL},
{" GR:1HN2", " N2 ", NULL},
{" GR:2HN2", " N2 ", NULL},
{" GR: H21", " N2 ", NULL},
{" GR: H22", " N2 ", NULL},
{" GR: H8 ", " C8 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" GR: H3 ", " C3*: C3': O3*: O3'", NULL},
{" GR: H4 ", " C4*: C4'", NULL},
{" GR:1H5 ", " C5*: C5'", NULL},
{" GR:2H5 ", " C5*: C5'", NULL},
{" GR:*H5 ", " C5*", NULL},
{" GR:H5' ", " C5'", NULL},
{" GR: H5 ", " C5': C5*: O5*: O5'", NULL},
{" GR:H2' ", " C2'", NULL},
{" GR:2HO ", " O2*: O2'", NULL},
{" GR:*HO ", " O2*: O3*: O5*", NULL},
{" GR:3HO ", " O3*: O3'", NULL},
{" GR:5HO ", " O5*: O5'", NULL},

{" GR: HN ", " N1 ", NULL},
{" GR:1HN ", " N2 ", NULL},
{" GR:2HN ", " N2 ", NULL},
{" GR: H2 ", " N2 : C2*: C2': O2*: O2'", NULL},

{" GR: H1 ", " N1 : C1*: C1'", NULL},
{" GR:1H2 ", " N2 : C2*: C2'", NULL},
{" GR:2H2 ", " N2 : C2*: C2'", NULL},
#else
{" GR: H1 ", " N1 ", NULL},
{" GR:1H2 ", " N2 ", NULL},
{" GR:2H2 ", " N2 ", NULL},
#endif

/* DNA in the remediated system is DA, DT, DG, DC */

{" DT: H1*", " C1*", NULL},
{" DT: H1'", " C1'", NULL},
{" DT: H3*", " C3*", NULL},
{" DT: H3'", " C3'", NULL},
{" DT: H4*", " C4*", NULL},
{" DT: H4'", " C4'", NULL},
{" DT:1H5*", " C5*", NULL},
{" DT:2H5*", " C5*", NULL},
{" DT:*H51", " C5*", NULL},
{" DT:*H52", " C5*", NULL},
{" DT:H5''", " C5'", NULL},
{" DT: H5'", " C5'", NULL},
{" DT:1H2*", " C2*", NULL},
{" DT:2H2*", " C2*", NULL},
{" DT: H2*", " C2*", NULL},
{" DT: H2'", " C2'", NULL}, /* remediated names RMI 070722 */
{" DT:H2''", " C2'", NULL},
{" DT: H3T", " O3*: O3'", NULL},
{" DT: H5T", " O5*: O5'", NULL},
{" DT:3HO*", " O3*", NULL},
{" DT:*HO3", " O3*", NULL},
{" DT:HO3'", " O3'", NULL}, /* remediated names RMI 070718 */
{" DT:5HO*", " O5*", NULL},
{" DT:*HO5", " O5*", NULL},
{" DT:HO5'", " O5'", NULL}, /* remediated names RMI 070718 */
{" DT: HN3", " N3 ", NULL},
{" DT:1H5M", " C7 : C5A: CA5: C5M", NULL},
{" DT:2H5M", " C7 : C5A: CA5: C5M", NULL},
{" DT:3H5M", " C7 : C5A: CA5: C5M", NULL},
{" DT:1HM5", " C7 : C5A: CA5: C5M", NULL},
{" DT:2HM5", " C7 : C5A: CA5: C5M", NULL},
{" DT:3HM5", " C7 : C5A: CA5: C5M", NULL},
{" DT:H5M1", " C7 : C5A: CA5: C5M", NULL},
{" DT:H5M2", " C7 : C5A: CA5: C5M", NULL},
{" DT:H5M3", " C7 : C5A: CA5: C5M", NULL},
{" DT:1H5A", " C7 : C5A: CA5: C5M", NULL},
{" DT:2H5A", " C7 : C5A: CA5: C5M", NULL},
{" DT:3H5A", " C7 : C5A: CA5: C5M", NULL},
{" DT:1HA5", " C7 : C5A: CA5: C5M", NULL},
{" DT:2HA5", " C7 : C5A: CA5: C5M", NULL},
{" DT:3HA5", " C7 : C5A: CA5: C5M", NULL},
{" DT:H5A1", " C7 : C5A: CA5: C5M", NULL},
{" DT:H5A2", " C7 : C5A: CA5: C5M", NULL},
{" DT:H5A3", " C7 : C5A: CA5: C5M", NULL},
{" DT: H71", " C7 : C5A: CA5: C5M", NULL}, /* remediated names RMI 070718 */
{" DT: H72", " C7 : C5A: CA5: C5M", NULL}, /*  */
{" DT: H73", " C7 : C5A: CA5: C5M", NULL}, /*  */
{" DT: H6 ", " C6 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" DT: H1 ", " C1*: C1'", NULL},
{" DT: H4 ", " C4*: C4'", NULL},
{" DT:*H5 ", " C5*", NULL},
{" DT:H5' ", " C5'", NULL},
{" DT: H5 ", " C5': C5*: O5*: O5'", NULL},
{" DT:1H2 ", " C2*: C2'", NULL},
{" DT:2H2 ", " C2*: C2'", NULL},
{" DT: H2 ", " C2*: C2'", NULL},
{" DT:H2' ", " C2'", NULL},
{" DT:*HO ", " O3*: O5*", NULL},
{" DT:3HO ", " O3*: O3'", NULL},
{" DT:5HO ", " O5*: O5'", NULL},

{" DT: HN ", " N3 ", NULL},
{" DT:1H5 ", " C5M: C5A: C5*: C5'", NULL},
{" DT:2H5 ", " C5M: C5A: C5*: C5'", NULL},
{" DT:3H5 ", " C5M: C5A", NULL},
{" DT:1HM ", " C5M: CM5", NULL},
{" DT:2HM ", " C5M: CM5", NULL},
{" DT:3HM ", " C5M: CM5", NULL},
{" DT:1HA ", " C5A: CA5", NULL},
{" DT:2HA ", " C5A: CA5", NULL},
{" DT:3HA ", " C5A: CA5", NULL},
{" DT:H5M ", " C5M", NULL},
{" DT:H5A ", " C5A", NULL},

{" DT: H3 ", " N3 : C3*: C3': O3*: O3'", NULL},
#else
{" DT: H3 ", " N3 ", NULL},
#endif

{" DA: H1*", " C1*", NULL},
{" DA: H1'", " C1'", NULL},
{" DA: H3*", " C3*", NULL},
{" DA: H3'", " C3'", NULL},
{" DA: H4*", " C4*", NULL},
{" DA: H4'", " C4'", NULL},
{" DA:1H5*", " C5*", NULL},
{" DA:2H5*", " C5*", NULL},
{" DA:*H51", " C5*", NULL},
{" DA:*H52", " C5*", NULL},
{" DA:H5''", " C5'", NULL},
{" DA: H5'", " C5'", NULL},
{" DA:1H2*", " C2*", NULL},
{" DA:2H2*", " C2*", NULL},
{" DA: H2*", " C2*", NULL},
{" DA: H2'", " C2'", NULL}, /* remediated names RMI 070722 */
{" DA:H2''", " C2'", NULL},
{" DA: H3T", " O3*: O3'", NULL},
{" DA: H5T", " O5*: O5'", NULL},
{" DA:3HO*", " O3*", NULL},
{" DA:*HO3", " O3*", NULL},
{" DA:HO3'", " O3'", NULL}, /* remediated names RMI 070718 */
{" DA:5HO*", " O5*", NULL},
{" DA:*HO5", " O5*", NULL},
{" DA:HO5'", " O5'", NULL}, /* remediated names RMI 070718 */
{" DA:1H6 ", " N6 ", NULL},
{" DA:2H6 ", " N6 ", NULL},
{" DA:1HN6", " N6 ", NULL},
{" DA:2HN6", " N6 ", NULL},
{" DA: H61", " N6 ", NULL},
{" DA: H62", " N6 ", NULL},
{" DA: H8 ", " C8 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" DA: H1 ", " C1*: C1'", NULL},
{" DA: H3 ", " C3*: C3': O3*: O3'", NULL},
{" DA: H4 ", " C4*: C4'", NULL},
{" DA:1H5 ", " C5*: C5'", NULL},
{" DA:2H5 ", " C5*: C5'", NULL},
{" DA:*H5 ", " C5*", NULL},
{" DA:H5' ", " C5'", NULL},
{" DA: H5 ", " C5': C5*: O5*: O5'", NULL},
{" DA:1H2 ", " C2*: C2'", NULL},
{" DA:2H2 ", " C2*: C2'", NULL},
{" DA:H2' ", " C2'", NULL},
{" DA:*HO ", " O3*: O5*", NULL},
{" DA:3HO ", " O3*: O3'", NULL},
{" DA:5HO ", " O5*: O5'", NULL},

{" DA:1HN ", " N6 ", NULL},
{" DA:2HN ", " N6 ", NULL},
{" DA: H6 ", " N6 ", NULL},

{" DA: H2 ", " C2 : C2*: C2'", NULL},
#else
{" DA: H2 ", " C2 ", NULL},
#endif

{" DC: H1*", " C1*", NULL},
{" DC: H1'", " C1'", NULL},
{" DC: H3*", " C3*", NULL},
{" DC: H3'", " C3'", NULL},
{" DC: H4*", " C4*", NULL},
{" DC: H4'", " C4'", NULL},
{" DC:1H5*", " C5*", NULL},
{" DC:2H5*", " C5*", NULL},
{" DC:*H51", " C5*", NULL},
{" DC:*H52", " C5*", NULL},
{" DC:H5''", " C5'", NULL},
{" DC: H5'", " C5'", NULL},
{" DC:1H2*", " C2*", NULL},
{" DC:2H2*", " C2*", NULL},
{" DC: H2*", " C2*", NULL},
{" DC: H2'", " C2'", NULL}, /* remediated names RMI 070722 */
{" DC:H2''", " C2'", NULL},
{" DC: H3T", " O3*: O3'", NULL},
{" DC: H5T", " O5*: O5'", NULL},
{" DC:3HO*", " O3*", NULL},
{" DC:*HO3", " O3*", NULL},
{" DC:HO3'", " O3'", NULL}, /* remediated names RMI 070718 */
{" DC:5HO*", " O5*", NULL},
{" DC:*HO5", " O5*", NULL},
{" DC:HO5'", " O5'", NULL}, /* remediated names RMI 070718 */
{" DC:1H4 ", " N4 ", NULL},
{" DC:2H4 ", " N4 ", NULL},
{" DC:1HN4", " N4 ", NULL},
{" DC:2HN4", " N4 ", NULL},
{" DC: H41", " N4 ", NULL},
{" DC: H42", " N4 ", NULL},
{" DC: H6 ", " C6 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" DC: H1 ", " C1*: C1'", NULL},
{" DC: H3 ", " C3*: C3': O3*: O3'", NULL},
{" DC:1H5 ", " C5*: C5'", NULL},
{" DC:2H5 ", " C5*: C5'", NULL},
{" DC:*H5 ", " C5*", NULL},
{" DC:H5' ", " C5'", NULL},
{" DC:1H2 ", " C2*: C2'", NULL},
{" DC:2H2 ", " C2*: C2'", NULL},
{" DC: H2 ", " C2*: C2'", NULL},
{" DC:H2' ", " C2'", NULL},
{" DC:*HO ", " O3*: O5*", NULL},
{" DC:3HO ", " O3*: O3'", NULL},
{" DC:5HO ", " O5*: O5'", NULL},

{" DC:1HN ", " N4 ", NULL},
{" DC:2HN ", " N4 ", NULL},
{" DC: H4 ", " N4 : C4*: C4'", NULL},

{" DC: H5 ", " C5 : C5': C5*: O5*: O5'", NULL},
#else
{" DC: H5 ", " C5 ", NULL},
#endif

{" DG: H1*", " C1*", NULL},
{" DG: H1'", " C1'", NULL},
{" DG: H3*", " C3*", NULL},
{" DG: H3'", " C3'", NULL},
{" DG: H4*", " C4*", NULL},
{" DG: H4'", " C4'", NULL},
{" DG:1H5*", " C5*", NULL},
{" DG:2H5*", " C5*", NULL},
{" DG:*H51", " C5*", NULL},
{" DG:*H52", " C5*", NULL},
{" DG:H5''", " C5'", NULL},
{" DG: H5'", " C5'", NULL},
{" DG:1H2*", " C2*", NULL},
{" DG:2H2*", " C2*", NULL},
{" DG: H2*", " C2*", NULL},
{" DG: H2'", " C2'", NULL}, /* remediated names RMI 070722 */
{" DG:H2''", " C2'", NULL},
{" DG: H3T", " O3*: O3'", NULL},
{" DG: H5T", " O5*: O5'", NULL},
{" DG:3HO*", " O3*", NULL},
{" DG:*HO3", " O3*", NULL},
{" DG:HO3'", " O3'", NULL}, /* remediated names RMI 070718 */
{" DG:5HO*", " O5*", NULL},
{" DG:*HO5", " O5*", NULL},
{" DG:HO5'", " O5'", NULL}, /* remediated names RMI 070718 */
{" DG: HN1", " N1 ", NULL},
{" DG:1HN2", " N2 ", NULL},
{" DG:2HN2", " N2 ", NULL},
{" DG: H21", " N2 ", NULL},
{" DG: H22", " N2 ", NULL},
{" DG: H8 ", " C8 ", NULL},
#ifdef ALLOW_TRUNCATED_H_NAMES
{" DG: H3 ", " C3*: C3': O3*: O3'", NULL},
{" DG: H4 ", " C4*: C4'", NULL},
{" DG:1H5 ", " C5*: C5'", NULL},
{" DG:2H5 ", " C5*: C5'", NULL},
{" DG:*H5 ", " C5*", NULL},
{" DG:H5' ", " C5'", NULL},
{" DG: H5 ", " C5': C5*: O5*: O5'", NULL},
{" DG:H2' ", " C2'", NULL},
{" DG:*HO ", " O3*: O5*", NULL},
{" DG:3HO ", " O3*: O3'", NULL},
{" DG:5HO ", " O5*: O5'", NULL},

{" DG: HN ", " N1 ", NULL},
{" DG:1HN ", " N2 ", NULL},
{" DG:2HN ", " N2 ", NULL},
{" DG: H2 ", " N2 : C2*: C2'", NULL},

{" DG: H1 ", " N1 : C1*: C1'", NULL},
{" DG:1H2 ", " N2 : C2*: C2'", NULL},
{" DG:2H2 ", " N2 : C2*: C2'", NULL},
#else
{" DG: H1 ", " N1 ", NULL},
{" DG:1H2 ", " N2 ", NULL},
{" DG:2H2 ", " N2 ", NULL},
#endif

/*---------------------------------------*/
/* the following definitions are used to specify how       */
/* heavy atoms interact with other heavy atoms WITHIN a    */
/* standard amino acid residue or base.                    */
/* This is accessed only when the -STDBOND flag is on the  */
/* command line (UseStdBond == TRUE).                      */

{"ALA: N  ", " CA ",           NULL},
{"ALA: NT ", " CA ",           NULL},
{"ALA: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"ALA: O  ", " C  ",           NULL},
{"ALA: OXT", " C  ",           NULL},
{"ALA: OT1", " C  ",           NULL},
{"ALA: OT2", " C  ",           NULL},
{"ALA: CA ", " N  : NT : C  : CB ", NULL},
{"ALA: CB ", " CA ",           NULL},

{"CYS: N  ", " CA ",           NULL},
{"CYS: NT ", " CA ",           NULL},
{"CYS: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"CYS: O  ", " C  ",           NULL},
{"CYS: OXT", " C  ",           NULL},
{"CYS: OT1", " C  ",           NULL},
{"CYS: OT2", " C  ",           NULL},
{"CYS: CA ", " N  : NT : C  : CB ", NULL},
{"CYS: CB ", " CA : SG ",      NULL},
{"CYS: SG ", " CB ",           NULL},

{"ASP: N  ", " CA ",           NULL},
{"ASP: NT ", " CA ",           NULL},
{"ASP: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"ASP: O  ", " C  ",           NULL},
{"ASP: OXT", " C  ",           NULL},
{"ASP: OT1", " C  ",           NULL},
{"ASP: OT2", " C  ",           NULL},
{"ASP: CA ", " N  : NT : C  : CB ", NULL},
{"ASP: CB ", " CA : CG ",      NULL},
{"ASP: CG ", " CB : OD1: OD2", NULL},
{"ASP: OD1", " CG ",           NULL},
{"ASP: OD2", " CG ",           NULL},

{"GLU: N  ", " CA ",           NULL},
{"GLU: NT ", " CA ",           NULL},
{"GLU: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"GLU: O  ", " C  ",           NULL},
{"GLU: OXT", " C  ",           NULL},
{"GLU: OT1", " C  ",           NULL},
{"GLU: OT2", " C  ",           NULL},
{"GLU: CA ", " N  : NT : C  : CB ", NULL},
{"GLU: CB ", " CA : CG ",      NULL},
{"GLU: CG ", " CB : CD ",      NULL},
{"GLU: CD ", " CG : OE1: OE2", NULL},
{"GLU: OE1", " CD ",           NULL},
{"GLU: OE2", " CD ",           NULL},

{"PHE: N  ", " CA ",           NULL},
{"PHE: NT ", " CA ",           NULL},
{"PHE: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"PHE: O  ", " C  ",           NULL},
{"PHE: OXT", " C  ",           NULL},
{"PHE: OT1", " C  ",           NULL},
{"PHE: OT2", " C  ",           NULL},
{"PHE: CA ", " N  : NT : C  : CB ", NULL},
{"PHE: CB ", " CA : CG ",      NULL},
{"PHE: CG ", " CB : CD1: CD2", NULL},
{"PHE: CD1", " CG : CE1",      NULL},
{"PHE: CD2", " CG : CE2",      NULL},
{"PHE: CE1", " CD1: CZ ",      NULL},
{"PHE: CE2", " CD2: CZ ",      NULL},
{"PHE: CZ ", " CE1: CE2",      NULL},

{"GLY: N  ", " CA ",           NULL},
{"GLY: NT ", " CA ",           NULL},
{"GLY: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"GLY: O  ", " C  ",           NULL},
{"GLY: OXT", " C  ",           NULL},
{"GLY: OT1", " C  ",           NULL},
{"GLY: OT2", " C  ",           NULL},
{"GLY: CA ", " N  : NT : C  : CB ",      NULL},
/* note CA-CB for gly is to support mutations gracefully */
/* where a gly MC is attached to some other SC */

{"HIS: N  ", " CA ",           NULL},
{"HIS: NT ", " CA ",           NULL},
{"HIS: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"HIS: O  ", " C  ",           NULL},
{"HIS: OXT", " C  ",           NULL},
{"HIS: OT1", " C  ",           NULL},
{"HIS: OT2", " C  ",           NULL},
{"HIS: CA ", " N  : NT : C  : CB ", NULL},
{"HIS: CB ", " CA : CG ",      NULL},
{"HIS: CG ", " CB : ND1: CD2", NULL},
{"HIS: ND1", " CG : CE1",      NULL},
{"HIS: CD2", " CG : NE2",      NULL},
{"HIS: CE1", " ND1: NE2",      NULL},
{"HIS: NE2", " CD2: CE1",      NULL},

{"ILE: N  ", " CA ",           NULL},
{"ILE: NT ", " CA ",           NULL},
{"ILE: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"ILE: O  ", " C  ",           NULL},
{"ILE: OXT", " C  ",           NULL},
{"ILE: OT1", " C  ",           NULL},
{"ILE: OT2", " C  ",           NULL},
{"ILE: CA ", " N  : NT : C  : CB ", NULL},
{"ILE: CB ", " CA : CG1: CG2", NULL},
{"ILE: CG1", " CB : CD1",      NULL},
{"ILE: CG2", " CB ",           NULL},
{"ILE: CD1", " CG1",           NULL},

{"LYS: N  ", " CA ",           NULL},
{"LYS: NT ", " CA ",           NULL},
{"LYS: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"LYS: O  ", " C  ",           NULL},
{"LYS: OXT", " C  ",           NULL},
{"LYS: OT1", " C  ",           NULL},
{"LYS: OT2", " C  ",           NULL},
{"LYS: CA ", " N  : NT : C  : CB ", NULL},
{"LYS: CB ", " CA : CG ",      NULL},
{"LYS: CG ", " CB : CD ",      NULL},
{"LYS: CD ", " CG : CE ",      NULL},
{"LYS: CE ", " CD : NZ ",      NULL},
{"LYS: NZ ", " CE ",           NULL},

{"LEU: N  ", " CA ",           NULL},
{"LEU: NT ", " CA ",           NULL},
{"LEU: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"LEU: O  ", " C  ",           NULL},
{"LEU: OXT", " C  ",           NULL},
{"LEU: OT1", " C  ",           NULL},
{"LEU: OT2", " C  ",           NULL},
{"LEU: CA ", " N  : NT : C  : CB ", NULL},
{"LEU: CB ", " CA : CG ",      NULL},
{"LEU: CG ", " CB : CD1: CD2", NULL},
{"LEU: CD1", " CG ",           NULL},
{"LEU: CD2", " CG ",           NULL},

{"MET: N  ", " CA ",           NULL},
{"MET: NT ", " CA ",           NULL},
{"MET: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"MET: O  ", " C  ",           NULL},
{"MET: OXT", " C  ",           NULL},
{"MET: OT1", " C  ",           NULL},
{"MET: OT2", " C  ",           NULL},
{"MET: CA ", " N  : NT : C  : CB ", NULL},
{"MET: CB ", " CA : CG ",      NULL},
{"MET: CG ", " CB : SD ",      NULL},
{"MET: SD ", " CG : CE ",      NULL},
{"MET: CE ", " SD ",           NULL},

{"ASN: N  ", " CA ",           NULL},
{"ASN: NT ", " CA ",           NULL},
{"ASN: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"ASN: O  ", " C  ",           NULL},
{"ASN: OXT", " C  ",           NULL},
{"ASN: OT1", " C  ",           NULL},
{"ASN: OT2", " C  ",           NULL},
{"ASN: CA ", " N  : NT : C  : CB ", NULL},
{"ASN: CB ", " CA : CG ",      NULL},
{"ASN: CG ", " CB : OD1: ND2", NULL},
{"ASN: OD1", " CG ",           NULL},
{"ASN: ND2", " CG ",           NULL},

{"PRO: N  ", " CA : CD ",      NULL},
{"PRO: NT ", " CA : CD ",      NULL},
{"PRO: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"PRO: O  ", " C  ",           NULL},
{"PRO: OXT", " C  ",           NULL},
{"PRO: OT1", " C  ",           NULL},
{"PRO: OT2", " C  ",           NULL},
{"PRO: CA ", " N  : NT : C  : CB ", NULL},
{"PRO: CB ", " CA : CG ",      NULL},
{"PRO: CG ", " CB : CD ",      NULL},
{"PRO: CD ", " CG : N  : NT ", NULL},

{"GLN: N  ", " CA ",           NULL},
{"GLN: NT ", " CA ",           NULL},
{"GLN: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"GLN: O  ", " C  ",           NULL},
{"GLN: OXT", " C  ",           NULL},
{"GLN: OT1", " C  ",           NULL},
{"GLN: OT2", " C  ",           NULL},
{"GLN: CA ", " N  : NT : C  : CB ", NULL},
{"GLN: CB ", " CA : CG ",      NULL},
{"GLN: CG ", " CB : CD ",      NULL},
{"GLN: CD ", " CG : OE1: NE2", NULL},
{"GLN: OE1", " CD ",           NULL},
{"GLN: NE2", " CD ",           NULL},

{"ARG: N  ", " CA ",           NULL},
{"ARG: NT ", " CA ",           NULL},
{"ARG: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"ARG: O  ", " C  ",           NULL},
{"ARG: OXT", " C  ",           NULL},
{"ARG: OT1", " C  ",           NULL},
{"ARG: OT2", " C  ",           NULL},
{"ARG: CA ", " N  : NT : C  : CB ", NULL},
{"ARG: CB ", " CA : CG ",      NULL},
{"ARG: CG ", " CB : CD ",      NULL},
{"ARG: CD ", " CG : NE ",      NULL},
{"ARG: NE ", " CD : CZ ",      NULL},
{"ARG: CZ ", " NE : NH1: NH2", NULL},
{"ARG: NH1", " CZ ",           NULL},
{"ARG: NH2", " CZ ",           NULL},

{"SER: N  ", " CA ",           NULL},
{"SER: NT ", " CA ",           NULL},
{"SER: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"SER: O  ", " C  ",           NULL},
{"SER: OXT", " C  ",           NULL},
{"SER: OT1", " C  ",           NULL},
{"SER: OT2", " C  ",           NULL},
{"SER: CA ", " N  : NT : C  : CB ", NULL},
{"SER: CB ", " CA : OG ",      NULL},
{"SER: OG ", " CB ",           NULL},

{"THR: N  ", " CA ",           NULL},
{"THR: NT ", " CA ",           NULL},
{"THR: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"THR: O  ", " C  ",           NULL},
{"THR: OXT", " C  ",           NULL},
{"THR: OT1", " C  ",           NULL},
{"THR: OT2", " C  ",           NULL},
{"THR: CA ", " N  : NT : C  : CB ", NULL},
{"THR: CB ", " CA : OG1: CG2", NULL},
{"THR: OG1", " CB ",           NULL},
{"THR: CG2", " CB ",           NULL},

{"VAL: N  ", " CA ",           NULL},
{"VAL: NT ", " CA ",           NULL},
{"VAL: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"VAL: O  ", " C  ",           NULL},
{"VAL: OXT", " C  ",           NULL},
{"VAL: OT1", " C  ",           NULL},
{"VAL: OT2", " C  ",           NULL},
{"VAL: CA ", " N  : NT : C  : CB ", NULL},
{"VAL: CB ", " CA : CG1: CG2", NULL},
{"VAL: CG1", " CB ",           NULL},
{"VAL: CG2", " CB ",           NULL},

{"TRP: N  ", " CA ",           NULL},
{"TRP: NT ", " CA ",           NULL},
{"TRP: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"TRP: O  ", " C  ",           NULL},
{"TRP: OXT", " C  ",           NULL},
{"TRP: OT1", " C  ",           NULL},
{"TRP: OT2", " C  ",           NULL},
{"TRP: CA ", " N  : NT : C  : CB ", NULL},
{"TRP: CB ", " CA : CG ",      NULL},
{"TRP: CG ", " CB : CD1: CD2", NULL},
{"TRP: CD1", " CG : NE1",      NULL},
{"TRP: CD2", " CG : CE2: CE3", NULL},
{"TRP: NE1", " CD1: CE2",      NULL},
{"TRP: CE2", " NE1: CD2: CZ2", NULL},
{"TRP: CZ2", " CE2: CH2",      NULL},
{"TRP: CE3", " CD2: CZ3",      NULL},
{"TRP: CZ3", " CE3: CH2",      NULL},
{"TRP: CH2", " CZ2: CZ3",      NULL},

{"TYR: N  ", " CA ",           NULL},
{"TYR: NT ", " CA ",           NULL},
{"TYR: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"TYR: O  ", " C  ",           NULL},
{"TYR: OXT", " C  ",           NULL},
{"TYR: OT1", " C  ",           NULL},
{"TYR: OT2", " C  ",           NULL},
{"TYR: CA ", " N  : NT : C  : CB ", NULL},
{"TYR: CB ", " CA : CG ",      NULL},
{"TYR: CG ", " CB : CD1: CD2", NULL},
{"TYR: CD1", " CG : CE1",      NULL},
{"TYR: CD2", " CG : CE2",      NULL},
{"TYR: CE1", " CD1: CZ ",      NULL},
{"TYR: CE2", " CD2: CZ ",      NULL},
{"TYR: CZ ", " CE1: CE2: OH ", NULL},
{"TYR: OH ", " CZ ",           NULL},

{"ASX: N  ", " CA ",           NULL},
{"ASX: NT ", " CA ",           NULL},
{"ASX: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"ASX: O  ", " C  ",           NULL},
{"ASX: OXT", " C  ",           NULL},
{"ASX: OT1", " C  ",           NULL},
{"ASX: OT2", " C  ",           NULL},
{"ASX: CA ", " N  : NT : C  : CB ", NULL},
{"ASX: CB ", " CA : CG ",      NULL},
{"ASX: CG ", " CB : OD1: ND2: AD1: AD2: XD1: XD2", NULL},  /* remediated RMI 070718 */
{"ASX: OD1", " CG ",           NULL},
{"ASX: AD1", " CG ",           NULL},
{"ASX: XD1", " CG ",           NULL}, /*  */
{"ASX: ND2", " CG ",           NULL},
{"ASX: AD2", " CG ",           NULL},
{"ASX: XD2", " CG ",           NULL}, /*  */

{"GLX: N  ", " CA ",           NULL},
{"GLX: NT ", " CA ",           NULL},
{"GLX: C  ", " CA : O  : OXT: OT1: OT2", NULL},
{"GLX: O  ", " C  ",           NULL},
{"GLX: OXT", " C  ",           NULL},
{"GLX: OT1", " C  ",           NULL},
{"GLX: OT2", " C  ",           NULL},
{"GLX: CA ", " N  : NT : C  : CB ", NULL},
{"GLX: CB ", " CA : CG ",      NULL},
{"GLX: CG ", " CB : CD ",      NULL},
{"GLX: CD ", " CG : OE1: NE2: AE1: AE2: XE1: XE2", NULL},  /* remediated RMI 070718 */
{"GLX: OE1", " CD ",           NULL},
{"GLX: AE1", " CD ",           NULL},
{"GLX: XE1", " CD ",           NULL}, /*  */
{"GLX: NE2", " CD ",           NULL},
{"GLX: AE2", " CD ",           NULL},
{"GLX: XE2", " CD ",           NULL}, /*  */

{"  U: C1*", " N1 : C2*: O4*", NULL},
{"  U: C1'", " N1 : C2': O4'", NULL},
{"  U: C2*", " C3*: O2*: C1*", NULL},
{"  U: O2*", " C2*",           NULL},
{"  U: C3*", " C4*: O3*: C2*", NULL},
{"  U: O3*", " C3*",           NULL},
{"  U: C4*", " C5*: O4*: C3*", NULL},
{"  U: O4*", " C4*: C1*",      NULL},
{"  U: C5*", " C4*: O5*",      NULL},
{"  U: O5*", " C5*: P  : PA ", NULL},
{"  U: C2'", " C3': O2': C1'", NULL},
{"  U: O2'", " C2'",           NULL},
{"  U: C3'", " C4': O3': C2'", NULL},
{"  U: O3'", " C3'",           NULL},
{"  U: C4'", " C5': O4': C3'", NULL},
{"  U: O4'", " C4': C1'",      NULL},
{"  U: C5'", " C4': O5'",      NULL},
{"  U: O5'", " C5': P  : PA ", NULL},
{"  U: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{"  U: O1P", " P  ",           NULL}, /*  */
{"  U: O2P", " P  ",           NULL}, /*  */
{"  U: O3P", " P  ",           NULL}, /*  */
{"  U: OP1", " P  ",           NULL},
{"  U: OP2", " P  ",           NULL},
{"  U: OP3", " P  ",           NULL},
{"  U: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{"  U: O1A", " PA ",           NULL},
{"  U: O2A", " PA ",           NULL},
{"  U: O3A", " PA : PB ",      NULL},
{"  U: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{"  U: O1B", " PB ",           NULL},
{"  U: O2B", " PB ",           NULL},
{"  U: O3B", " PB : PG ",      NULL},
{"  U: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{"  U: O1G", " PG ",           NULL},
{"  U: O2G", " PG ",           NULL},
{"  U: O3G", " PG ",           NULL},
{"  U: S1G", " PG ",           NULL},
{"  U: N1 ", " C2 : C6 : C1*: C1'", NULL},
{"  U: C2 ", " N1 : N3 : O2 ", NULL},
{"  U: N3 ", " C2 : C4 ",      NULL},
{"  U: C4 ", " N3 : C5 : O4 ", NULL},
{"  U: C5 ", " C4 : C6 ",      NULL},
{"  U: C6 ", " C5 : N1 ",      NULL},
{"  U: O2 ", " C2 ",           NULL},
{"  U: O4 ", " C4 ",           NULL},

{"  T: C1*", " N1 : C2*: O4*", NULL},
{"  T: C1'", " N1 : C2': O4'", NULL},
{"  T: C2*", " C3*: O2*: C1*", NULL},
{"  T: O2*", " C2*",           NULL},
{"  T: C3*", " C4*: O3*: C2*", NULL},
{"  T: O3*", " C3*",           NULL},
{"  T: C4*", " C5*: O4*: C3*", NULL},
{"  T: O4*", " C4*: C1*",      NULL},
{"  T: C5*", " C4*: O5*",      NULL},
{"  T: O5*", " C5*: P  : PA ", NULL},
{"  T: C2'", " C3': O2': C1'", NULL},
{"  T: O2'", " C2'",           NULL},
{"  T: C3'", " C4': O3': C2'", NULL},
{"  T: O3'", " C3'",           NULL},
{"  T: C4'", " C5': O4': C3'", NULL},
{"  T: O4'", " C4': C1'",      NULL},
{"  T: C5'", " C4': O5'",      NULL},
{"  T: O5'", " C5': P  : PA ", NULL},
{"  T: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{"  T: O1P", " P  ",           NULL}, /*  */
{"  T: O2P", " P  ",           NULL}, /*  */
{"  T: O3P", " P  ",           NULL}, /*  */
{"  T: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{"  T: O1A", " PA ",           NULL},
{"  T: O2A", " PA ",           NULL},
{"  T: O3A", " PA : PB ",      NULL},
{"  T: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{"  T: O1B", " PB ",           NULL},
{"  T: O2B", " PB ",           NULL},
{"  T: O3B", " PB : PG ",      NULL},
{"  T: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{"  T: O1G", " PG ",           NULL},
{"  T: O2G", " PG ",           NULL},
{"  T: O3G", " PG ",           NULL},
{"  T: S1G", " PG ",           NULL},
{"  T: N1 ", " C2 : C6 : C1*: C1'", NULL},
{"  T: C2 ", " N1 : N3 : O2 ", NULL},
{"  T: N3 ", " C2 : C4 ",      NULL},
{"  T: C4 ", " N3 : C5 : O4 ", NULL},
{"  T: C5 ", " C4 : C6 : C5M: C5A: CM5: CA5: C7 ", NULL}, /* remediated RMI 070718 */
{"  T: C6 ", " C5 : N1 ",      NULL},
{"  T: O2 ", " C2 ",           NULL},
{"  T: O4 ", " C4 ",           NULL},
{"  T: C5M", " C5 ",           NULL},
{"  T: C5A", " C5 ",           NULL},
{"  T: CM5", " C5 ",           NULL},
{"  T: CA5", " C5 ",           NULL},
{"  T: C7 ", " C5 ",           NULL}, /*  */


{"  C: C1*", " N1 : C2*: O4*", NULL},
{"  C: C1'", " N1 : C2': O4'", NULL},
{"  C: C2*", " C3*: O2*: C1*", NULL},
{"  C: O2*", " C2*",           NULL},
{"  C: C3*", " C4*: O3*: C2*", NULL},
{"  C: O3*", " C3*",           NULL},
{"  C: C4*", " C5*: O4*: C3*", NULL},
{"  C: O4*", " C4*: C1*",      NULL},
{"  C: C5*", " C4*: O5*",      NULL},
{"  C: O5*", " C5*: P  : PA ", NULL},
{"  C: C2'", " C3': O2': C1'", NULL},
{"  C: O2'", " C2'",           NULL},
{"  C: C3'", " C4': O3': C2'", NULL},
{"  C: O3'", " C3'",           NULL},
{"  C: C4'", " C5': O4': C3'", NULL},
{"  C: O4'", " C4': C1'",      NULL},
{"  C: C5'", " C4': O5'",      NULL},
{"  C: O5'", " C5': P  : PA ", NULL},
{"  C: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{"  C: O1P", " P  ",           NULL}, /*  */
{"  C: O2P", " P  ",           NULL}, /*  */
{"  C: O3P", " P  ",           NULL}, /*  */
{"  C: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{"  C: O1A", " PA ",           NULL},
{"  C: O2A", " PA ",           NULL},
{"  C: O3A", " PA : PB ",      NULL},
{"  C: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{"  C: O1B", " PB ",           NULL},
{"  C: O2B", " PB ",           NULL},
{"  C: O3B", " PB : PG ",      NULL},
{"  C: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{"  C: O1G", " PG ",           NULL},
{"  C: O2G", " PG ",           NULL},
{"  C: O3G", " PG ",           NULL},
{"  C: S1G", " PG ",           NULL},
{"  C: N1 ", " C2 : C6 : C1*: C1'", NULL},
{"  C: C2 ", " N1 : N3 : O2 ", NULL},
{"  C: N3 ", " C2 : C4 ",      NULL},
{"  C: C4 ", " N3 : C5 : N4 ", NULL},
{"  C: C5 ", " C4 : C6 ",      NULL},
{"  C: C6 ", " C5 : N1 ",      NULL},
{"  C: O2 ", " C2 ",           NULL},
{"  C: N4 ", " C4 ",           NULL},

{"  A: C1*", " N9 : C2*: O4*", NULL},
{"  A: C1'", " N9 : C2': O4'", NULL},
{"  A: C2*", " C3*: O2*: C1*", NULL},
{"  A: O2*", " C2*",           NULL},
{"  A: C3*", " C4*: O3*: C2*", NULL},
{"  A: O3*", " C3*",           NULL},
{"  A: C4*", " C5*: O4*: C3*", NULL},
{"  A: O4*", " C4*: C1*",      NULL},
{"  A: C5*", " C4*: O5*",      NULL},
{"  A: O5*", " C5*: P  : PA ", NULL},
{"  A: C2'", " C3': O2': C1'", NULL},
{"  A: O2'", " C2'",           NULL},
{"  A: C3'", " C4': O3': C2'", NULL},
{"  A: O3'", " C3'",           NULL},
{"  A: C4'", " C5': O4': C3'", NULL},
{"  A: O4'", " C4': C1'",      NULL},
{"  A: C5'", " C4': O5'",      NULL},
{"  A: O5'", " C5': P  : PA ", NULL},
{"  A: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{"  A: O1P", " P  ",           NULL}, /*  */
{"  A: O2P", " P  ",           NULL}, /*  */
{"  A: O3P", " P  ",           NULL}, /*  */
{"  A: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{"  A: O1A", " PA ",           NULL},
{"  A: O2A", " PA ",           NULL},
{"  A: O3A", " PA : PB ",      NULL},
{"  A: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{"  A: O1B", " PB ",           NULL},
{"  A: O2B", " PB ",           NULL},
{"  A: O3B", " PB : PG ",      NULL},
{"  A: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{"  A: O1G", " PG ",           NULL},
{"  A: O2G", " PG ",           NULL},
{"  A: O3G", " PG ",           NULL},
{"  A: S1G", " PG ",           NULL},
{"  A: N9 ", " C8 : C4 : C1*: C1'", NULL},
{"  A: C8 ", " N9 : N7 ",      NULL},
{"  A: N7 ", " C8 : C5 ",      NULL},
{"  A: C5 ", " N7 : C6 : C4 ", NULL},
{"  A: C6 ", " C5 : N6 : N1 ", NULL},
{"  A: N6 ", " C6 ",           NULL},
{"  A: N1 ", " C6 : C2 ",      NULL},
{"  A: C2 ", " N1 : N3 ",      NULL},
{"  A: N3 ", " C2 : C4 ",      NULL},
{"  A: C4 ", " N9 : C5 : N3 ", NULL},

{"  G: C1*", " N9 : C2*: O4*", NULL},
{"  G: C1'", " N9 : C2': O4'", NULL},
{"  G: C2*", " C3*: O2*: C1*", NULL},
{"  G: O2*", " C2*",           NULL},
{"  G: C3*", " C4*: O3*: C2*", NULL},
{"  G: O3*", " C3*",           NULL},
{"  G: C4*", " C5*: O4*: C3*", NULL},
{"  G: O4*", " C4*: C1*",      NULL},
{"  G: C5*", " C4*: O5*",      NULL},
{"  G: O5*", " C5*: P  : PA ", NULL},
{"  G: C2'", " C3': O2': C1'", NULL},
{"  G: O2'", " C2'",           NULL},
{"  G: C3'", " C4': O3': C2'", NULL},
{"  G: O3'", " C3'",           NULL},
{"  G: C4'", " C5': O4': C3'", NULL},
{"  G: O4'", " C4': C1'",      NULL},
{"  G: C5'", " C4': O5'",      NULL},
{"  G: O5'", " C5': P  : PA ", NULL},
{"  G: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{"  G: O1P", " P  ",           NULL}, /*  */
{"  G: O2P", " P  ",           NULL}, /*  */
{"  G: O3P", " P  ",           NULL}, /*  */
{"  G: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{"  G: O1A", " PA ",           NULL},
{"  G: O2A", " PA ",           NULL},
{"  G: O3A", " PA : PB ",      NULL},
{"  G: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{"  G: O1B", " PB ",           NULL},
{"  G: O2B", " PB ",           NULL},
{"  G: O3B", " PB : PG ",      NULL},
{"  G: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{"  G: O1G", " PG ",           NULL},
{"  G: O2G", " PG ",           NULL},
{"  G: O3G", " PG ",           NULL},
{"  G: S1G", " PG ",           NULL},
{"  G: N9 ", " C8 : C4 : C1*: C1'", NULL},
{"  G: C8 ", " N9 : N7 ",      NULL},
{"  G: N7 ", " C8 : C5 ",      NULL},
{"  G: C5 ", " N7 : C6 : C4 ", NULL},
{"  G: C6 ", " C5 : O6 : N1 ", NULL},
{"  G: O6 ", " C6 ",           NULL},
{"  G: N1 ", " C6 : C2 ",      NULL},
{"  G: C2 ", " N1 : N2 : N3 ", NULL},
{"  G: N2 ", " C2 ",           NULL},
{"  G: N3 ", " C2 : C4 ",      NULL},
{"  G: C4 ", " N9 : C5 : N3 ", NULL},

/* DNA in the remediated system is DA, DT, DG, DC */

{" DT: C1*", " N1 : C2*: O4*", NULL},
{" DT: C1'", " N1 : C2': O4'", NULL},
{" DT: C2*", " C3*: C1*", NULL},
{" DT: C3*", " C4*: O3*: C2*", NULL},
{" DT: O3*", " C3*",           NULL},
{" DT: C4*", " C5*: O4*: C3*", NULL},
{" DT: O4*", " C4*: C1*",      NULL},
{" DT: C5*", " C4*: O5*",      NULL},
{" DT: O5*", " C5*: P  : PA ", NULL},
{" DT: C2'", " C3': C1'", NULL},
{" DT: C3'", " C4': O3': C2'", NULL},
{" DT: O3'", " C3'",           NULL},
{" DT: C4'", " C5': O4': C3'", NULL},
{" DT: O4'", " C4': C1'",      NULL},
{" DT: C5'", " C4': O5'",      NULL},
{" DT: O5'", " C5': P  : PA ", NULL},
{" DT: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" DT: O1P", " P  ",           NULL}, /*  */
{" DT: O2P", " P  ",           NULL}, /*  */
{" DT: O3P", " P  ",           NULL}, /*  */
{" DT: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" DT: O1A", " PA ",           NULL},
{" DT: O2A", " PA ",           NULL},
{" DT: O3A", " PA : PB ",      NULL},
{" DT: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" DT: O1B", " PB ",           NULL},
{" DT: O2B", " PB ",           NULL},
{" DT: O3B", " PB : PG ",      NULL},
{" DT: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" DT: O1G", " PG ",           NULL},
{" DT: O2G", " PG ",           NULL},
{" DT: O3G", " PG ",           NULL},
{" DT: S1G", " PG ",           NULL},
{" DT: N1 ", " C2 : C6 : C1*: C1'", NULL},
{" DT: C2 ", " N1 : N3 ", NULL},
{" DT: N3 ", " C2 : C4 ",      NULL},
{" DT: C4 ", " N3 : C5 : O4 ", NULL},
{" DT: C5 ", " C4 : C6 : C5M: C5A: CM5: CA5: C7 ", NULL}, /* remediated RMI 070718 */
{" DT: C6 ", " C5 : N1 ",      NULL},
{" DT: O4 ", " C4 ",           NULL},
{" DT: C5M", " C5 ",           NULL},
{" DT: C5A", " C5 ",           NULL},
{" DT: CM5", " C5 ",           NULL},
{" DT: CA5", " C5 ",           NULL},
{" DT: C7 ", " C5 ",           NULL}, /*  */


{" DC: C1*", " N1 : C2*: O4*", NULL},
{" DC: C1'", " N1 : C2': O4'", NULL},
{" DC: C2*", " C3*: C1*", NULL},
{" DC: C3*", " C4*: O3*: C2*", NULL},
{" DC: O3*", " C3*",           NULL},
{" DC: C4*", " C5*: O4*: C3*", NULL},
{" DC: O4*", " C4*: C1*",      NULL},
{" DC: C5*", " C4*: O5*",      NULL},
{" DC: O5*", " C5*: P  : PA ", NULL},
{" DC: C2'", " C3': C1'", NULL},
{" DC: C3'", " C4': O3': C2'", NULL},
{" DC: O3'", " C3'",           NULL},
{" DC: C4'", " C5': O4': C3'", NULL},
{" DC: O4'", " C4': C1'",      NULL},
{" DC: C5'", " C4': O5'",      NULL},
{" DC: O5'", " C5': P  : PA ", NULL},
{" DC: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" DC: O1P", " P  ",           NULL}, /*  */
{" DC: O2P", " P  ",           NULL}, /*  */
{" DC: O3P", " P  ",           NULL}, /*  */
{" DC: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" DC: O1A", " PA ",           NULL},
{" DC: O2A", " PA ",           NULL},
{" DC: O3A", " PA : PB ",      NULL},
{" DC: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" DC: O1B", " PB ",           NULL},
{" DC: O2B", " PB ",           NULL},
{" DC: O3B", " PB : PG ",      NULL},
{" DC: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" DC: O1G", " PG ",           NULL},
{" DC: O2G", " PG ",           NULL},
{" DC: O3G", " PG ",           NULL},
{" DC: S1G", " PG ",           NULL},
{" DC: N1 ", " C2 : C6 : C1*: C1'", NULL},
{" DC: C2 ", " N1 : N3 ", NULL},
{" DC: N3 ", " C2 : C4 ",      NULL},
{" DC: C4 ", " N3 : C5 : N4 ", NULL},
{" DC: C5 ", " C4 : C6 ",      NULL},
{" DC: C6 ", " C5 : N1 ",      NULL},
{" DC: N4 ", " C4 ",           NULL},

{" DA: C1*", " N9 : C2*: O4*", NULL},
{" DA: C1'", " N9 : C2': O4'", NULL},
{" DA: C2*", " C3*: C1*", NULL},
{" DA: C3*", " C4*: O3*: C2*", NULL},
{" DA: O3*", " C3*",           NULL},
{" DA: C4*", " C5*: O4*: C3*", NULL},
{" DA: O4*", " C4*: C1*",      NULL},
{" DA: C5*", " C4*: O5*",      NULL},
{" DA: O5*", " C5*: P  : PA ", NULL},
{" DA: C2'", " C3': C1'", NULL},
{" DA: C3'", " C4': O3': C2'", NULL},
{" DA: O3'", " C3'",           NULL},
{" DA: C4'", " C5': O4': C3'", NULL},
{" DA: O4'", " C4': C1'",      NULL},
{" DA: C5'", " C4': O5'",      NULL},
{" DA: O5'", " C5': P  : PA ", NULL},
{" DA: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" DA: O1P", " P  ",           NULL}, /*  */
{" DA: O2P", " P  ",           NULL}, /*  */
{" DA: O3P", " P  ",           NULL}, /*  */
{" DA: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" DA: O1A", " PA ",           NULL},
{" DA: O2A", " PA ",           NULL},
{" DA: O3A", " PA : PB ",      NULL},
{" DA: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" DA: O1B", " PB ",           NULL},
{" DA: O2B", " PB ",           NULL},
{" DA: O3B", " PB : PG ",      NULL},
{" DA: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" DA: O1G", " PG ",           NULL},
{" DA: O2G", " PG ",           NULL},
{" DA: O3G", " PG ",           NULL},
{" DA: S1G", " PG ",           NULL},
{" DA: N9 ", " C8 : C4 : C1*: C1'", NULL},
{" DA: C8 ", " N9 : N7 ",      NULL},
{" DA: N7 ", " C8 : C5 ",      NULL},
{" DA: C5 ", " N7 : C6 : C4 ", NULL},
{" DA: C6 ", " C5 : N6 : N1 ", NULL},
{" DA: N6 ", " C6 ",           NULL},
{" DA: N1 ", " C6 : C2 ",      NULL},
{" DA: C2 ", " N1 : N3 ",      NULL},
{" DA: N3 ", " C2 : C4 ",      NULL},
{" DA: C4 ", " N9 : C5 : N3 ", NULL},

{" DG: C1*", " N9 : C2*: O4*", NULL},
{" DG: C1'", " N9 : C2': O4'", NULL},
{" DG: C2*", " C3*: C1*", NULL},
{" DG: C3*", " C4*: O3*: C2*", NULL},
{" DG: O3*", " C3*",           NULL},
{" DG: C4*", " C5*: O4*: C3*", NULL},
{" DG: O4*", " C4*: C1*",      NULL},
{" DG: C5*", " C4*: O5*",      NULL},
{" DG: O5*", " C5*: P  : PA ", NULL},
{" DG: C2'", " C3': C1'", NULL},
{" DG: C3'", " C4': O3': C2'", NULL},
{" DG: O3'", " C3'",           NULL},
{" DG: C4'", " C5': O4': C3'", NULL},
{" DG: O4'", " C4': C1'",      NULL},
{" DG: C5'", " C4': O5'",      NULL},
{" DG: O5'", " C5': P  : PA ", NULL},
{" DG: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" DG: O1P", " P  ",           NULL}, /*  */
{" DG: O2P", " P  ",           NULL}, /*  */
{" DG: O3P", " P  ",           NULL}, /*  */
{" DG: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" DG: O1A", " PA ",           NULL},
{" DG: O2A", " PA ",           NULL},
{" DG: O3A", " PA : PB ",      NULL},
{" DG: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" DG: O1B", " PB ",           NULL},
{" DG: O2B", " PB ",           NULL},
{" DG: O3B", " PB : PG ",      NULL},
{" DG: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" DG: O1G", " PG ",           NULL},
{" DG: O2G", " PG ",           NULL},
{" DG: O3G", " PG ",           NULL},
{" DG: S1G", " PG ",           NULL},
{" DG: N9 ", " C8 : C4 : C1*: C1'", NULL},
{" DG: C8 ", " N9 : N7 ",      NULL},
{" DG: N7 ", " C8 : C5 ",      NULL},
{" DG: C5 ", " N7 : C6 : C4 ", NULL},
{" DG: C6 ", " C5 : O6 : N1 ", NULL},
{" DG: O6 ", " C6 ",           NULL},
{" DG: N1 ", " C6 : C2 ",      NULL},
{" DG: C2 ", " N1 : N2 : N3 ", NULL},
{" DG: N2 ", " C2 ",           NULL},
{" DG: N3 ", " C2 : C4 ",      NULL},
{" DG: C4 ", " N9 : C5 : N3 ", NULL},

/* RNA for Coot and CCP4 Ar, Ur, Tr, Cr, Gr */

{" UR: C1*", " N1 : C2*: O4*", NULL},
{" UR: C1'", " N1 : C2': O4'", NULL},
{" UR: C2*", " C3*: O2*: C1*", NULL},
{" UR: O2*", " C2*",           NULL},
{" UR: C3*", " C4*: O3*: C2*", NULL},
{" UR: O3*", " C3*",           NULL},
{" UR: C4*", " C5*: O4*: C3*", NULL},
{" UR: O4*", " C4*: C1*",      NULL},
{" UR: C5*", " C4*: O5*",      NULL},
{" UR: O5*", " C5*: P  : PA ", NULL},
{" UR: C2'", " C3': O2': C1'", NULL},
{" UR: O2'", " C2'",           NULL},
{" UR: C3'", " C4': O3': C2'", NULL},
{" UR: O3'", " C3'",           NULL},
{" UR: C4'", " C5': O4': C3'", NULL},
{" UR: O4'", " C4': C1'",      NULL},
{" UR: C5'", " C4': O5'",      NULL},
{" UR: O5'", " C5': P  : PA ", NULL},
{" UR: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" UR: O1P", " P  ",           NULL}, /*  */
{" UR: O2P", " P  ",           NULL}, /*  */
{" UR: O3P", " P  ",           NULL}, /*  */
{" UR: OP1", " P  ",           NULL},
{" UR: OP2", " P  ",           NULL},
{" UR: OP3", " P  ",           NULL},
{" UR: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" UR: O1A", " PA ",           NULL},
{" UR: O2A", " PA ",           NULL},
{" UR: O3A", " PA : PB ",      NULL},
{" UR: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" UR: O1B", " PB ",           NULL},
{" UR: O2B", " PB ",           NULL},
{" UR: O3B", " PB : PG ",      NULL},
{" UR: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" UR: O1G", " PG ",           NULL},
{" UR: O2G", " PG ",           NULL},
{" UR: O3G", " PG ",           NULL},
{" UR: S1G", " PG ",           NULL},
{" UR: N1 ", " C2 : C6 : C1*: C1'", NULL},
{" UR: C2 ", " N1 : N3 : O2 ", NULL},
{" UR: N3 ", " C2 : C4 ",      NULL},
{" UR: C4 ", " N3 : C5 : O4 ", NULL},
{" UR: C5 ", " C4 : C6 ",      NULL},
{" UR: C6 ", " C5 : N1 ",      NULL},
{" UR: O2 ", " C2 ",           NULL},
{" UR: O4 ", " C4 ",           NULL},

{" TR: C1*", " N1 : C2*: O4*", NULL},
{" TR: C1'", " N1 : C2': O4'", NULL},
{" TR: C2*", " C3*: O2*: C1*", NULL},
{" TR: O2*", " C2*",           NULL},
{" TR: C3*", " C4*: O3*: C2*", NULL},
{" TR: O3*", " C3*",           NULL},
{" TR: C4*", " C5*: O4*: C3*", NULL},
{" TR: O4*", " C4*: C1*",      NULL},
{" TR: C5*", " C4*: O5*",      NULL},
{" TR: O5*", " C5*: P  : PA ", NULL},
{" TR: C2'", " C3': O2': C1'", NULL},
{" TR: O2'", " C2'",           NULL},
{" TR: C3'", " C4': O3': C2'", NULL},
{" TR: O3'", " C3'",           NULL},
{" TR: C4'", " C5': O4': C3'", NULL},
{" TR: O4'", " C4': C1'",      NULL},
{" TR: C5'", " C4': O5'",      NULL},
{" TR: O5'", " C5': P  : PA ", NULL},
{" TR: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" TR: O1P", " P  ",           NULL}, /*  */
{" TR: O2P", " P  ",           NULL}, /*  */
{" TR: O3P", " P  ",           NULL}, /*  */
{" TR: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" TR: O1A", " PA ",           NULL},
{" TR: O2A", " PA ",           NULL},
{" TR: O3A", " PA : PB ",      NULL},
{" TR: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" TR: O1B", " PB ",           NULL},
{" TR: O2B", " PB ",           NULL},
{" TR: O3B", " PB : PG ",      NULL},
{" TR: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" TR: O1G", " PG ",           NULL},
{" TR: O2G", " PG ",           NULL},
{" TR: O3G", " PG ",           NULL},
{" TR: S1G", " PG ",           NULL},
{" TR: N1 ", " C2 : C6 : C1*: C1'", NULL},
{" TR: C2 ", " N1 : N3 : O2 ", NULL},
{" TR: N3 ", " C2 : C4 ",      NULL},
{" TR: C4 ", " N3 : C5 : O4 ", NULL},
{" TR: C5 ", " C4 : C6 : C5M: C5A: CM5: CA5: C7 ", NULL}, /* remediated RMI 070718 */
{" TR: C6 ", " C5 : N1 ",      NULL},
{" TR: O2 ", " C2 ",           NULL},
{" TR: O4 ", " C4 ",           NULL},
{" TR: C5M", " C5 ",           NULL},
{" TR: C5A", " C5 ",           NULL},
{" TR: CM5", " C5 ",           NULL},
{" TR: CA5", " C5 ",           NULL},
{" TR: C7 ", " C5 ",           NULL}, /*  */


{" CR: C1*", " N1 : C2*: O4*", NULL},
{" CR: C1'", " N1 : C2': O4'", NULL},
{" CR: C2*", " C3*: O2*: C1*", NULL},
{" CR: O2*", " C2*",           NULL},
{" CR: C3*", " C4*: O3*: C2*", NULL},
{" CR: O3*", " C3*",           NULL},
{" CR: C4*", " C5*: O4*: C3*", NULL},
{" CR: O4*", " C4*: C1*",      NULL},
{" CR: C5*", " C4*: O5*",      NULL},
{" CR: O5*", " C5*: P  : PA ", NULL},
{" CR: C2'", " C3': O2': C1'", NULL},
{" CR: O2'", " C2'",           NULL},
{" CR: C3'", " C4': O3': C2'", NULL},
{" CR: O3'", " C3'",           NULL},
{" CR: C4'", " C5': O4': C3'", NULL},
{" CR: O4'", " C4': C1'",      NULL},
{" CR: C5'", " C4': O5'",      NULL},
{" CR: O5'", " C5': P  : PA ", NULL},
{" CR: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" CR: O1P", " P  ",           NULL}, /*  */
{" CR: O2P", " P  ",           NULL}, /*  */
{" CR: O3P", " P  ",           NULL}, /*  */
{" CR: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" CR: O1A", " PA ",           NULL},
{" CR: O2A", " PA ",           NULL},
{" CR: O3A", " PA : PB ",      NULL},
{" CR: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" CR: O1B", " PB ",           NULL},
{" CR: O2B", " PB ",           NULL},
{" CR: O3B", " PB : PG ",      NULL},
{" CR: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" CR: O1G", " PG ",           NULL},
{" CR: O2G", " PG ",           NULL},
{" CR: O3G", " PG ",           NULL},
{" CR: S1G", " PG ",           NULL},
{" CR: N1 ", " C2 : C6 : C1*: C1'", NULL},
{" CR: C2 ", " N1 : N3 : O2 ", NULL},
{" CR: N3 ", " C2 : C4 ",      NULL},
{" CR: C4 ", " N3 : C5 : N4 ", NULL},
{" CR: C5 ", " C4 : C6 ",      NULL},
{" CR: C6 ", " C5 : N1 ",      NULL},
{" CR: O2 ", " C2 ",           NULL},
{" CR: N4 ", " C4 ",           NULL},

{" AR: C1*", " N9 : C2*: O4*", NULL},
{" AR: C1'", " N9 : C2': O4'", NULL},
{" AR: C2*", " C3*: O2*: C1*", NULL},
{" AR: O2*", " C2*",           NULL},
{" AR: C3*", " C4*: O3*: C2*", NULL},
{" AR: O3*", " C3*",           NULL},
{" AR: C4*", " C5*: O4*: C3*", NULL},
{" AR: O4*", " C4*: C1*",      NULL},
{" AR: C5*", " C4*: O5*",      NULL},
{" AR: O5*", " C5*: P  : PA ", NULL},
{" AR: C2'", " C3': O2': C1'", NULL},
{" AR: O2'", " C2'",           NULL},
{" AR: C3'", " C4': O3': C2'", NULL},
{" AR: O3'", " C3'",           NULL},
{" AR: C4'", " C5': O4': C3'", NULL},
{" AR: O4'", " C4': C1'",      NULL},
{" AR: C5'", " C4': O5'",      NULL},
{" AR: O5'", " C5': P  : PA ", NULL},
{" AR: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" AR: O1P", " P  ",           NULL}, /*  */
{" AR: O2P", " P  ",           NULL}, /*  */
{" AR: O3P", " P  ",           NULL}, /*  */
{" AR: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" AR: O1A", " PA ",           NULL},
{" AR: O2A", " PA ",           NULL},
{" AR: O3A", " PA : PB ",      NULL},
{" AR: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" AR: O1B", " PB ",           NULL},
{" AR: O2B", " PB ",           NULL},
{" AR: O3B", " PB : PG ",      NULL},
{" AR: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" AR: O1G", " PG ",           NULL},
{" AR: O2G", " PG ",           NULL},
{" AR: O3G", " PG ",           NULL},
{" AR: S1G", " PG ",           NULL},
{" AR: N9 ", " C8 : C4 : C1*: C1'", NULL},
{" AR: C8 ", " N9 : N7 ",      NULL},
{" AR: N7 ", " C8 : C5 ",      NULL},
{" AR: C5 ", " N7 : C6 : C4 ", NULL},
{" AR: C6 ", " C5 : N6 : N1 ", NULL},
{" AR: N6 ", " C6 ",           NULL},
{" AR: N1 ", " C6 : C2 ",      NULL},
{" AR: C2 ", " N1 : N3 ",      NULL},
{" AR: N3 ", " C2 : C4 ",      NULL},
{" AR: C4 ", " N9 : C5 : N3 ", NULL},

{" GR: C1*", " N9 : C2*: O4*", NULL},
{" GR: C1'", " N9 : C2': O4'", NULL},
{" GR: C2*", " C3*: O2*: C1*", NULL},
{" GR: O2*", " C2*",           NULL},
{" GR: C3*", " C4*: O3*: C2*", NULL},
{" GR: O3*", " C3*",           NULL},
{" GR: C4*", " C5*: O4*: C3*", NULL},
{" GR: O4*", " C4*: C1*",      NULL},
{" GR: C5*", " C4*: O5*",      NULL},
{" GR: O5*", " C5*: P  : PA ", NULL},
{" GR: C2'", " C3': O2': C1'", NULL},
{" GR: O2'", " C2'",           NULL},
{" GR: C3'", " C4': O3': C2'", NULL},
{" GR: O3'", " C3'",           NULL},
{" GR: C4'", " C5': O4': C3'", NULL},
{" GR: O4'", " C4': C1'",      NULL},
{" GR: C5'", " C4': O5'",      NULL},
{" GR: O5'", " C5': P  : PA ", NULL},
{" GR: P  ", " O5*: O5': O1P: O2P: O3P: OP1: OP2: OP3", NULL}, /* remediated RMI 070718 */
{" GR: O1P", " P  ",           NULL}, /*  */
{" GR: O2P", " P  ",           NULL}, /*  */
{" GR: O3P", " P  ",           NULL}, /*  */
{" GR: PA ", " O5*: O5': O1A: O2A: O3A", NULL},
{" GR: O1A", " PA ",           NULL},
{" GR: O2A", " PA ",           NULL},
{" GR: O3A", " PA : PB ",      NULL},
{" GR: PB ", " O5*: O5': O1B: O2B: O3B", NULL},
{" GR: O1B", " PB ",           NULL},
{" GR: O2B", " PB ",           NULL},
{" GR: O3B", " PB : PG ",      NULL},
{" GR: PG ", " O5*: O5': O1G: O2G: O3G : S1G", NULL},
{" GR: O1G", " PG ",           NULL},
{" GR: O2G", " PG ",           NULL},
{" GR: O3G", " PG ",           NULL},
{" GR: S1G", " PG ",           NULL},
{" GR: N9 ", " C8 : C4 : C1*: C1'", NULL},
{" GR: C8 ", " N9 : N7 ",      NULL},
{" GR: N7 ", " C8 : C5 ",      NULL},
{" GR: C5 ", " N7 : C6 : C4 ", NULL},
{" GR: C6 ", " C5 : O6 : N1 ", NULL},
{" GR: O6 ", " C6 ",           NULL},
{" GR: N1 ", " C6 : C2 ",      NULL},
{" GR: C2 ", " N1 : N2 : N3 ", NULL},
{" GR: N2 ", " C2 ",           NULL},
{" GR: N3 ", " C2 : C4 ",      NULL},
{" GR: C4 ", " N9 : C5 : N3 ", NULL},


{NULL, NULL, NULL}
};

/* initStdConnTable() - we have to do this once to setup the hash table */

void initStdConnTable() {
   StdResConnTableEntry_t *r;
   int i;

   for (i = 0; i < M; i++) { /* empty all the buckets */
      StdResTblBucket[i] = NULL;
   }

   for( r = StandardResAtomConnRec; r->key != NULL; r++) {
      InsertInStdResConnTable(r);
   }
}

/* searchForStdBondingPartner() - given a residue name and atom name find the bonded atom */
char * searchForStdBondingPartner(char *resname, char *atomname, int isAhydrogen) {
   char querystring[10], *p;

   /* a large part of this routine is devoted to converting residue names */
   /* and atom names into a search pattern. This pattern is 8 chars long  */
   /* and all uppercase. It has a colon as the forth character, separating*/
   /* residue name from atom name. Deuterum is converted to hydrogen.     */

   if (resname && atomname) {
      if (resname[0]) {
	 querystring[0] = toupper(resname[0]);
	 if (resname[1]) {
	    querystring[1] = toupper(resname[1]);
	    querystring[2] = resname[2] ? toupper(resname[2]) : ' ';
	 }
	 else { querystring[1] = querystring[2] = ' '; }
      }
      else { querystring[0] = querystring[1] = querystring[2] = ' '; }
      querystring[3] = ':';
      if (atomname[0]) {
	 querystring[4] = toupper(atomname[0]);
	 if (atomname[1]) {
	    querystring[5] = toupper(atomname[1]);
	    if (atomname[2]) {
	       querystring[6] = toupper(atomname[2]);
	       querystring[7] = atomname[3] ? toupper(atomname[3]) : ' ';
	    }
	    else { querystring[6] = querystring[7] = ' '; }
	 }
	 else { querystring[5] = querystring[6] = querystring[7] = ' '; }
      }
      else {
	 querystring[4] = querystring[5] = ' ';
	 querystring[6] = querystring[7] = ' ';
      }
      querystring[8] = '\0';

      if (isAhydrogen) {
	 /* because we know atomname is for a hydrogen atom we can */
	 /* make the following fixup for deuterium atoms           */
	 if ((querystring[4] == 'D') && (querystring[7] != ' ')) {
	    querystring[4] = 'H';
	 }
	 else if ((querystring[4] != 'H') && (querystring[5] == 'D')) {
	    querystring[5] = 'H';
	 }
      }

      p = SearchStdResConnTable(querystring);

      if (p == NULL) {
	 /* if failed because of a non-std nucleic acid residue name */
	 /* convert to "..X" and retry                               */

	 char *naALTlist = ":GUA:GTP:GDP:GMP:GSP:ADE:ATP:ADP:AMP:CYT:CTP:CDP:CMP:URA:UTP:UDP:UMP:THY:TTP:TDP:TMP: DA: DT: DC: DG: AR: UR: TR: CR: GR:";

	 querystring[3] = '\0'; /* temporarily look at just the first three chars */
	 if (strstr(naALTlist, querystring)) {
	    /* fixup common alternate nucleic acid residue names */
	         if(strstr(":GUA:GTP:GDP:GMP:GSP: DG: GR:", querystring)) {
	       querystring[2] = 'G';
	    }
	    else if(strstr(":ADE:ATP:ADP:AMP: DA: AR:",     querystring)) {
	       querystring[2] = 'A';
	    }
	    else if(strstr(":CYT:CTP:CDP:CMP: DC: CR:",     querystring)) {
	       querystring[2] = 'C';
	    }
	    else if(strstr(":URA:UTP:UDP:UMP: UR:",     querystring)) {
	       querystring[2] = 'U';
	    }
	    else if(strstr(":THY:TTP:TDP:TMP: DT: TR:",     querystring)) {
	       querystring[2] = 'T';
	    }
	    else { querystring[2] = ' '; }

	    querystring[0] = querystring[1] = ' '; /* i.e., blank-blank-[GACUT] */
				    /* put the string back together */
	    querystring[3] = ':';

	    p = SearchStdResConnTable(querystring); /* second time is a charm */
	 }
      }
   }
   else { p = NULL; }
   return p;
}

/* Hash function from PJ Weinberger's compiler as described          */
/* in Aho, Sethi & Ullman, Compilers, Principles, Techniques & Tools */
/* 1986, Addison Wesley, pg 436                                      */
int HashInStdResTbl(char *s) {
   unsigned h, g;

   for (h = 0; *s != '\0'; s++) {
      h = (h<<4) + *s;
      if ( (g = h & 0xf0000000) ) {
	 h ^= g >> 24;
	 h ^= g;
      }
   }
   return h % M;
}

/* InsertInStdResConnTable() - add an entry into the has table */
int InsertInStdResConnTable(StdResConnTableEntry_t *elem) {
   int r = 0, h = 0, rc = 1;
   StdResConnTableEntry_t *t = NULL;

   h = HashInStdResTbl(elem->key);
   t = StdResTblBucket[h];

   if (t == NULL) {           /* nothing else stored at this h value yet   */
           elem->next = NULL; /* the cdr of the new list is Z              */
      StdResTblBucket[h] = elem; /* store the element at the head of the list */
   }
   else if ((r = strcmp(elem->key, t->key)) > 0) { /* key(e) > key(t), look further down list */

      while (t->next && (r = strcmp(elem->key, t->next->key)) > 0) {
	 t = t->next;
      }

      if ((t->next == NULL) || (r != 0)) { /* hook elem in after t and before t->next */
	 elem->next = t->next;
	    t->next = elem;
      }
      else {
	 fprintf(stderr, "ERROR duplicate key: InsertInStdResConnTable(\"%s\", \"%s\") -- old(\"%s\" @%d)\n",
	       elem->key, elem->value, t->next->value, h);
	 rc = 0;
      }
   }
   else if (r < 0) {       /* the element comes ahead of the first list element */
           elem->next = t;    /* the cdr of the new list is the old list   */
      StdResTblBucket[h] = elem; /* store the element at the head of the list */
   }
   else { /* duplicate */
      fprintf(stderr, "ERROR duplicate key: InsertInStdResConnTable(\"%s\", \"%s\") - old(\"%s\" @%d)\n",
	    elem->key, elem->value, t->value, h);
      rc = 0;
   }
   return rc;
}

/* SearchStdResConnTable() - given a formatted search key, return the value in the hash table */
char *SearchStdResConnTable(char *key) {
   int r = 0;
   StdResConnTableEntry_t *t;

   t = StdResTblBucket[HashInStdResTbl(key)];

   while (t && (r = strcmp(key, t->key)) > 0) {
      t = t->next;
   }

   return (t ? (r ? NULL : t->value) : NULL);
}

/* dumpStdConnTable() - dump dictionary contents for debugging purposes */
void dumpStdConnTable(FILE * outf) {
   int h;
   StdResConnTableEntry_t *t;

   for (h = 0; h < M; h++) {
      fprintf(outf, "%5d] ", h);
      for (t = StdResTblBucket[h]; t; t = t->next) {
	 fprintf(outf, "(\"%s\" --> \"%s\") ", t->key, t->value);
      }
      fprintf(outf, "\n");
   }
}
