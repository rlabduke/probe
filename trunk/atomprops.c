/* name: atomprops.c                                     */
/* author: J. Michael Word     date written: 6/12/97     */
/* purpose: define atom properties                       */

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
#include <string.h>  /*060902 needs this for strlen() */
#include <ctype.h>
#define INIT_ATOM_TABLE 1
#include "atomprops.h"

static atomProp* AtomTblIndex[NUMATOMTYPES];
static float ExplMaxRad = 0.0;
static float ImplMaxRad = 0.0;

extern int NuclearRadii; /* defined in probe.c JJH */

void initalizeAtomTbl() {
   int i;
   atomProp *ap;

   ExplMaxRad = 0.0;
   ImplMaxRad = 0.0;
   for(i=0; i < NUMATOMTYPES; i++) { /* default to noAtom */
      AtomTblIndex[i] = &(AtomTbl[0]);
   }

   for(i=0; i < NUMATOMTYPES; i++) {
      ap = &(AtomTbl[i]);
      AtomTblIndex[ap->type] = ap;

      if (ap->iRad     > ImplMaxRad){ ImplMaxRad = ap->iRad;     }
      if (!NuclearRadii) {
        if (ap->eRad     > ExplMaxRad){ ExplMaxRad = ap->eRad;   }
      }
      else {
        if (ap->eRad_nuc > ExplMaxRad){ ExplMaxRad = ap->eRad_nuc; }
      }
   }
}

int   getAtno(int a)     { return AtomTblIndex[a]->atno;   }
char* getAtomName(int a) { return AtomTblIndex[a]->name;   }
float getExplRad(int a)  {
  if (NuclearRadii) {
    /*if ( (strcmp(AtomTblIndex[a]->name, "H")    == 0) ||
         (strcmp(AtomTblIndex[a]->name, "Har")  == 0) ||
         (strcmp(AtomTblIndex[a]->name, "Hpol") == 0) ||
         (strcmp(AtomTblIndex[a]->name, "HOd")  == 0) ) {
      return ( (AtomTblIndex[a]->eRad - 0.05) );
    }*/
    return AtomTblIndex[a]->eRad_nuc;
  }
  else {
    return AtomTblIndex[a]->eRad;
  }
}
float getImplRad(int a)  { return AtomTblIndex[a]->iRad;   }
float getCovRad(int a)   { return AtomTblIndex[a]->covRad; }
char* getColor(int a)    { return AtomTblIndex[a]->color;  }
float getMaxRadius(int isImpl){
   return (isImpl? ImplMaxRad : ExplMaxRad);
}
int   atomHasProp(int a, int f) {
   return AtomTblIndex[a]->flags & f;
}

int fixAtomName(const char* atomname, char resname[], int position) { /* no bool in C */
   char resn[6];
   char name[5] = "    ";
   int i;
   sprintf(resn, ":%-3.3s:", resname);
   for (i = 0; i < 4; i++) { /* uppercase the input */
      if (atomname[i] == '\0') { break; }
      name[i] = toupper(atomname[i]);
   }
   name[i] = '\0';
        switch(name[position]) {
           case 'E': if (strstr(HE_RESNAMES, resn) != NULL) { return 1; }
           case 'F': if (strstr(HF_RESNAMES, resn) != NULL) { return 1; }
           case 'G': if (strstr(HG_RESNAMES, resn) != NULL) { return 1; }
           case 'O': if (strstr(HO_RESNAMES, resn) != NULL) { return 1; }
           case 'S': if (strstr(HS_RESNAMES, resn) != NULL) { return 1; }
           default: break;
	}
   return 0;
}

int identifyAtom(char* name, char resname[], int Verbose) {  /*dcr041007 allow warning choice*/
   int n = -1, emitWarning = 0;

   switch(name[0]) {
   case '*': case '\'': case '"': case '`': case '_':
   case '+': case '-':  case ' ':
   case '0': case '1': case '2': case '3': case '4':
   case '5': case '6': case '7': case '8': case '9':
      switch(name[1]) {
      case 'A':
	 switch(name[3]) {
	 case '1': n = atomO; emitWarning = 1; break;
	 case '2': n = atomN; emitWarning = 1; break;
	 } break;
      case 'B': n = atomB; break;
      case 'C': n = atomC; break;
      case 'D': n = atomH; break;
      case 'F': n = atomF; break;
      case 'H':
        switch(name[2]) {
        case 'E': n = fixAtomName(name,resname,2) ? atomHe : atomH; break;
        case 'F': n = fixAtomName(name,resname,2) ? atomHf : atomH; break;
        case 'G': n = fixAtomName(name,resname,2) ? atomHg : atomH; break;
        case 'O': n = fixAtomName(name,resname,2) ? atomHo : atomH; break;
/*        case 'S': n = fixAtomName(name,resname,2) ? atomHs : atomH; break; */
        default : n = atomH; break;
        } break;
      case 'I': n = atomI; break;
      case 'K': n = atomK; break;
      case 'N': n = atomN; break;
      case 'O': n = atomO; break;
      case 'P': n = atomP; break;
    /*case 'S': n = atomS; break;*/
      case 'S':
         if(name[0] == ' ' && name[2] == 'E')
         {/*_SE likely refmac,cns missplaced Selenium atom name dcr041007*/
                n = atomSe; emitWarning = 1; break;
         }
         else  {n = atomS; break;} /*050121 NO warning for other S*/
      case 'U': n = atomU; break;
      case 'V': n = atomV; break;
      case 'W': n = atomW; break;
      case 'Y': n = atomY; break;
      } break;
   case 'A':
      switch(name[1]) {
      case 'C': n = atomC;  emitWarning = 1;break;/*nonstd!*/
      case 'G': n = atomAg; break;
      case 'H': n = atomH;  emitWarning = 1;break;
      case 'L': n = atomAl; break;
      case 'M': n = atomAm; break;
      case 'N': n = atomN;  emitWarning = 1;break;
      case 'O': n = atomO;  emitWarning = 1;break;
      case 'P': n = atomP;  emitWarning = 1;break;
      case 'R': n = atomAr; break;
      case 'S': n = atomAs; break;
      case 'T': n = atomAt; break;
      case 'U': n = atomAu; break;
      } break;
   case 'B':
      switch(name[1]) {
      case 'A': n = atomBa; break;
      case 'E': n = atomBe; break;
      case 'I': n = atomBi; break;
      case 'K': n = atomBk; break;
      case 'R': n = atomBr; break;
      } break;
   case 'C':
      switch(name[1]) {
      case 'A': n = atomCa; break;
      case 'C': n = atomC;  emitWarning = 1;break;
      case 'D': n = atomCd; break;
      case 'E': n = atomCe; break;
      case 'F': n = atomCf; break;
      case 'H': n = atomH;  emitWarning = 1;break;
      case 'L': n = atomCl; break;
      case 'M': n = atomCm; break;
      case 'N': n = atomN;  emitWarning = 1;break;
      case 'O': n = atomCo; emitWarning = 1;break;
      case 'P': n = atomP;  emitWarning = 1;break;
      case 'R': n = atomCr; break;
      case 'S': n = atomCs; break;
      case 'U': n = atomCu; break;
      default:  n = atomC;  emitWarning = 1;break;
      } break;
   case 'D':
      switch(name[1]) {
      case 'Y': n = atomDy; break;
      case 'C': n = atomC;  emitWarning = 1;break;
      case 'H': n = atomH;  emitWarning = 1;break;
      case 'N': n = atomN;  emitWarning = 1;break;
      case 'O': n = atomO;  emitWarning = 1;break;
      case 'P': n = atomP;  emitWarning = 1;break;
      default:  n = atomH;  emitWarning = 1;break;
      } break;
   case 'E':
      switch(name[1]) {
      case 'R': n = atomEr; break;
      case 'S': n = atomEs; break;
      case 'U': n = atomEu; break;
      case 'C': n = atomC;  emitWarning = 1;break;
      case 'H': n = atomH;  emitWarning = 1;break;
      case 'N': n = atomN;  emitWarning = 1;break;
      case 'O': n = atomO;  emitWarning = 1;break;
      case 'P': n = atomP;  emitWarning = 1;break;
      } break;
   case 'F':
      switch(name[1]) {
      case 'E': n = atomFe; break;
      case 'M': n = atomFm; break;
      case 'R': n = atomFr; break;
      case 'C': n = atomC;  emitWarning = 1;break;
      case 'H': n = atomH;  emitWarning = 1;break;
      case 'N': n = atomN;  emitWarning = 1;break;
      case 'O': n = atomO;  emitWarning = 1;break;
      case 'P': n = atomP;  emitWarning = 1;break;
      } break;
   case 'G':
      switch(name[1]) {
      case 'A': n = atomGa; break;
      case 'D': n = atomGd; break;
      case 'E': n = atomGe; break;
      case 'C': n = atomC;  emitWarning = 1;break;
      case 'H': n = atomH;  emitWarning = 1;break;
      case 'N': n = atomN;  emitWarning = 1;break;
      case 'O': n = atomO;  emitWarning = 1;break;
      case 'P': n = atomP;  emitWarning = 1;break;
      } break;
   case 'H':
      switch(name[1]) {
      case 'E': n = fixAtomName(name,resname,1) ? atomHe : atomH; break;
      case 'F': n = fixAtomName(name,resname,1) ? atomHf : atomH; break;
      case 'G': n = fixAtomName(name,resname,1) ? atomHg : atomH; break;
      case 'O': n = fixAtomName(name,resname,1) ? atomHo : atomH; break;
/*      case 'S': n = fixAtomName(name,resname,1) ? atomHs : atomH; break; */
      default : n = atomH; break;
      } break;

/*    case 'E': */
/*  if (isdigit(name[2])) { n = atomH; emitWarning = 1;} */ /* Hepsilon?? */
/*  else {n = atomHe; emitWarning = 1;} */
/*  break; */
/*      case 'F': n = atomHf; emitWarning = 1;break; */
/*      case 'G': */
/*  if (isdigit(name[2])) { n = atomH; emitWarning = 1;} */ /* Hgamma?? */
/*  else {n = atomHg; emitWarning = 1;} */
/*  break; */
/*      case 'O': n = atomHo; emitWarning = 1;break; */
/*      default:  n = atomH;  emitWarning = 1;break; */
/*      } break; */

   case 'I':
      switch(name[1]) {
      case 'N': n = atomIn; break;
      case 'R': n = atomIr; break;
      } break;
   case 'K':
      if (name[1] == 'R') n = atomKr; break;
   case 'L':
      switch(name[1]) {
      case 'A': n = atomLa; break;
      case 'I': n = atomLi; break;
      case 'U': n = atomLu; break;
      } break;
   case 'M':
      switch(name[1]) {
      case 'D': n = atomMd; break;
      case 'G': n = atomMg; break;
      case 'N': n = atomMn; break;
      case 'O': n = atomMo; break;
      } break;
   case 'N':
      switch(name[1]) {
      case 'A': n = atomNa; emitWarning = 1;break;
      case 'B': n = atomNb; emitWarning = 1;break;
      case 'C': n = atomC;  emitWarning = 1;break;
      case 'D': n = atomNd; emitWarning = 1;break;
      case 'E': n = atomNe; emitWarning = 1;break;
      case 'H': n = atomH;  emitWarning = 1;break;
      case 'I': n = atomNi; break;
      case 'N': n = atomN;  emitWarning = 1;break;
      case 'O': n = atomO;  emitWarning = 1;break;/*nonstd!*/
      case 'P': n = atomP;  emitWarning = 1;break;/*nonstd!*/
      case 'S': n = atomS;  emitWarning = 1;break;
      default:  n = atomN;  emitWarning = 1;break;
      } break;
   case 'O':
      switch(name[1]) {
      case 'S': n = atomOs; break;
      default:  n = atomO;  emitWarning = 1;break;
      } break;
   case 'P':
      switch(name[1]) {
      case 'A': n = atomPa; emitWarning = 1;break;
      case 'B': n = atomPb; emitWarning = 1;break;
      case 'D': n = atomPd; emitWarning = 1;break;
      case 'M': n = atomPm; break;
      case 'O': n = atomPo; break;
      case 'R': n = atomPr; break;
      case 'T': n = atomPt; break;
      case 'U': n = atomPu; break;
      default:  n = atomP;  emitWarning = 1;break;
      } break;
   case 'R':
      switch(name[1]) {
      case 'A': n = atomRa; break;
      case 'B': n = atomRb; break;
      case 'E': n = atomRe; break;
      case 'H': n = atomRh; break;
      case 'N': n = atomRn; break;
      case 'U': n = atomRu; break;
      } break;
   case 'S':
      switch(name[1]) {
      case 'B': n = atomSb; emitWarning = 1;break;
      case 'C': n = atomSc; break;
      case 'E': n = atomSe; emitWarning = 1;break;
      case 'I': n = atomSi; break;
      case 'M': n = atomSm; break;
      case 'N': n = atomSn; break;
      case 'R': n = atomSr; break;
      default:  n = atomS;  emitWarning = 1;break;
      } break;
   case 'T':
      switch(name[1]) {
      case 'A': n = atomTa; break;
      case 'B': n = atomTb; break;
      case 'C': n = atomTc; break;
      case 'E': n = atomTe; break;
      case 'H': n = atomTh; break;
      case 'I': n = atomTi; break;
      case 'L': n = atomTl; break;
      case 'M': n = atomTm; break;
      } break;
   case 'X':
      if (name[1] == 'E') n = atomXe; break;
   case 'Y':
      if (name[1] == 'B') n = atomYb; break;
   case 'Z':
      switch(name[1]) {
      case 'N': n = atomZn; break;
      case 'R': n = atomZr; break;
      } break;
   default: break;
   }

   if (n < 0) { emitWarning = 1;
      n = atomC;

      switch(name[1]) {
      case 'H': case 'D': n = atomH; break;
      case 'C': n = atomC; break;
      case 'N': n = atomN; break;
      case 'O': n = atomO; break;
      case 'P': n = atomP; break;
      case 'S': n = atomS; break;

      case 'I': n = atomI; break;
      case 'K': n = atomK; break;
      case 'V': n = atomV; break;
      case 'W': n = atomW; break;
      case 'U': n = atomU; break;

      case 'A':
	 switch(name[2]) {
	 case 'G': n = atomAg; break;
	 case 'L': n = atomAl; break;
	 case 'S': n = atomAs; break;
	 case 'U': n = atomAu; break;
	 } break;
      case 'F':
	 if (name[2] == 'E') n = atomFe; break;
      case 'G':
	 if (name[2] == 'D') n = atomGd; break;
      case 'L':
	 if (name[2] == 'I') n = atomLi; break;
      case 'M':
	 switch(name[2]) {
	 case 'G': n = atomMg; break;
	 case 'N': n = atomMn; break;
	 case 'O': n = atomMo; break;
	 } break;
      case 'Z':
	 if (name[2] == 'N') n = atomZn; break;
/* --------- default if we fall through... ---------------------*/
      default: n = atomC; break;
      }
   }
   if (emitWarning && Verbose) {/*Verbose controlled dcr041007*/
      char warnstr[5]={'_','_','_','_','\0'}; /*dcr041007*/
      char warnresn[5]={'_','_','_','_','\0'}; /*rmi070719*/
      char* atstr;
      int j=0;
      while(name[j]!='\0'&&warnstr[j]!='\0')
      {/*explicitly show blanks so can understand atom name problem*/
         if(name[j]==' ') {warnstr[j]='_';}
         else             {warnstr[j]=name[j];}
         j++;
      }
      j=0;
      while(resname[j]!='\0'&&warnresn[j]!='\0')
      {/*explicitly show blanks so can understand residue/atom name problem*/
         if(name[j]==' ') {warnresn[j]='_';}
         else             {warnresn[j]=resname[j];}
         j++;
      }
      atstr = getAtomName(n);
      if(strlen(atstr)==1)
           {fprintf(stderr,"WARNING: atom %s from resn %s will be treated as  %s\n",warnstr,warnresn,atstr);}
      else {fprintf(stderr,"WARNING: atom %s from resn %s will be treated as %s\n",warnstr,warnresn,atstr);}
         /*name, getAtomName(n));*/
   }

   return n;
}

