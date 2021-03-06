/* name: autobondrot.c                                          */
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utility.h"
#include "autobondrot.h"
#include "readPDBrecs.h"

/*describeXformDB() writes header-comments to the .map file! */
void describeXformDB(FILE *outf, xformDatabase* xdb, char *cch);
void fillinTransformationRec(transformData* xform, char* rec);
xformDatabase* initXformDB();
transformData* newTransformation(char* rec);
transformData* deleteTransformation(transformData* xform,
            abrAtomListProc delAtomProc, void *otherstuff);
void appendTransformation(xformDatabase* xdb, transformData* xform);
void appendAtomRec(transformData* xform, char* rec, abrMkAtomProc mkAtom, void *atomstuff);
xformAtomRecords * newXformAtomRecord(char *rec, abrMkAtomProc mkAtom, void *otherstuff);
xformAtomRecords* deleteXformAtomRecord(xformAtomRecords *xr,
            abrAtomListProc delAtomProc, void *deletestuff);
void indexLevels(xformDatabase* xdb);
void buildXDBsortedAtomRecs(xformDatabase* xdb);
void rebuildTmat(transformData* xform);
atom* orientedAtoms(xformDatabase* xdb);
void setOrientation(xformDatabase* xdb, goToRec* g);
void cleanupListExtras(atom *alst, xformDatabase* xdb,
            abrAtomListProc delAtomProc, void *deletestuff);
void initializeXFStack(xformDatabase* xdb);
int pushXFStack(xformDatabase* xdb, matrix4* ctm);
void popXFStack(xformDatabase* xdb);
matrix4* topXFStack(xformDatabase* xdb);
double runnameAndTorsion(xformDatabase* xdb, char* buf, int max);
void addGoToRecord(xformDatabase* xdb, char *rec);
void XFdescribeAtomList(FILE *outf, atom* a);
void appendBiasFunction(transformData* xform, char *rec);
int compareXDBAtomRecs(const void *a, const void *b);
void scanXDBinputText(xformDatabase* xdb, FILE *inf, FILE *outf,
	    abrMkAtomProc mkAtom, void *atomstuff, char *cch);
void doListCleanupProcessing(xformDatabase* xdb, abrAtomListProc inputListProc, void *liststuff);

#define MAX_XF_REC_LEN  300
static char XFworkBuf[MAX_XF_REC_LEN + 1]; 
/*e.g. angle values held in XFworkBuf for output to .map file*/

/* read in records and interpret them to fillout a database */
xformDatabase* readTransformationDatabase(FILE *inf, FILE *outf,
	    abrMkAtomProc mkAtom, void *atomstuff,
	    abrAtomListProc inputListProc, void *liststuff,
	    char *cch) 
{
/*
        mabis.inf                == FILE *inf
        outf                     == FILE *outf
        newMovingAtom            == abrMkAtomProc mkAtom
        &mabis                   == void *atomstuff
        movingAtomListProcessing == abrAtomListProc inputListProc
        NULL                     == void *liststuff
        RAW_HEADER_COMMENT       == char *cch
*/

   xformDatabase* xdb = NULL;

   xdb = initXformDB();
   if (xdb) 
   {
      scanXDBinputText(xdb, inf, outf, mkAtom, atomstuff, cch);

      indexLevels(xdb);

      buildXDBsortedAtomRecs(xdb);

      doListCleanupProcessing(xdb, inputListProc, liststuff);

      /*describeXformDB() writes header-comments to the .map file! */
      describeXformDB(outf, xdb, cch);
   }
   return xdb;
}

void scanXDBinputText(xformDatabase* xdb, FILE *inf, FILE *outf,
	    abrMkAtomProc mkAtom, void *atomstuff, char *cch) {
   transformData* xform = NULL;
   char *s = NULL;
   char saynull[] = "NULL:identity transformation:";

   while(readRecord(inf, XFworkBuf, MAX_XF_REC_LEN) >= 0) {
      if (isAtom(XFworkBuf) || isHet(XFworkBuf)) {
	 if (! xform) { /* make someplace to hang these atoms */
	    xform = newTransformation(saynull);
	    appendTransformation(xdb, xform);
	 }
	 if (! xform) {
	    errmsg("could not create Identity transformation - Atom skipped");
	 }
	 else if (isPseudoAtom(XFworkBuf)) { /* ignore pseudo atoms */ }
	 else {
	    appendAtomRec(xform, XFworkBuf, mkAtom, atomstuff);
	 }
      }
      else {
	 s = XFworkBuf + strspn(XFworkBuf, " \t\n\r");
	 if (!strncasecmp(s, "BONDROT",7)
	  || !strncasecmp(s, "ROT",  3)
	  || !strncasecmp(s, "TRANS",5)
	  || !strncasecmp(s, "SAVE", 4)   || (s[0] == '(')
	  || !strncasecmp(s, "RESTORE",7) || (s[0] == ')') ) {
	    xform = newTransformation(s);
	    appendTransformation(xdb, xform);
	 }
	 else if (!strncasecmp(s, "GO", 2)) {
	    addGoToRecord(xdb, s);
	 }
	 else if (!strncasecmp(s, "COS", 3)
	       || !strncasecmp(s, "POLY",   4)
	       || !strncasecmp(s, "CONST",  5) ) {
	    if (xform) {
	       appendBiasFunction(xform, s);
	    }
	    else {
	       warn("no transformation to attach bias function:");
	       warn(s);
	    }
	 }
	 else if (s[0] == '#' || s[0] == '\0') {
	       /* comments are ok - we do nothing */
	 }
	 else if (s[0] == '@') { /* include file */
	    FILE *if2 = NULL;
	    if2 = fopen(s+1, "rb");
	    if (if2) { /* recursive call */
	       fprintf(outf, "%s%s\n", cch, s);

	       scanXDBinputText(xdb, if2, outf, mkAtom, atomstuff, cch);
	       fclose(if2);
	    }
	    else {
	       warn("could not open include file:");
	       warn(s+1);
	    }
	 }
	 else {
	    warn("unknown record type:");
	    warn(s);
	 }
      }
   }
}

void discardXformDB(xformDatabase* xdb,
            abrAtomListProc delAtomProc, void *deletestuff) {
   int i = 0;
   transformData *p = NULL, *next = NULL;
   goToRec *g = NULL, *nextg = NULL;

   if (xdb) {
      for(p = xdb->xforms; p; p = next) {
         next = p->next;
         deleteTransformation(p, delAtomProc, deletestuff);
      }
      xdb->xforms = NULL;
      xdb->last   = NULL;

      for(i = 0; i < MAX_XFORM_LEVELS; i++) {
         xdb->level[i] = NULL;
      }
      xdb->nlevs  = 0;
      xdb->nstack = 0;
      free(xdb->ctmstack[0]); /* all stack elements allocated in a block */
      for(g = xdb->golist; g; g = nextg) {
         nextg = g->next;
         free(g);
      }
      if (xdb->sortedByRes) { free(xdb->sortedByRes); }
      xdb->sortedByRes = NULL;

      free(xdb);
   }
}

void autobondrot(FILE *outf, xformDatabase* xdb,
	       abrProbeProc    probeProc,   void *probestuff,
	       abrAtomListProc delAtomProc, void *deletestuff,
	       int dumpGoToAtoms) 
{/*autobondrot()*/
/*where:  (see probe/autobondrot(...); )
   outf  == stderr
   xdb   == xdb
   probeProc == movingDoCommand
   probestuff == mcis
   delAtomProc == deleteMovingAtom
   deletestuff == mabis
   dumpGoToAtoms == Verbose
*/

   int overflow = FALSE, cursor = 0;
   transformData* xform = NULL;
   goToRec* g = NULL;
   atom* alst = NULL;
   double torsionScore = 0.0;

   if (xdb)
   {
      if (xdb->golist)  /*working from a list of GO statements */
      {/*if using GO statements: specific orientations requested */
         for(g = xdb->golist; g; g = g->next) 
         {
            setOrientation(xdb, g);

            alst = orientedAtoms(xdb); /* atom* alst == atom *allMovingAtoms */
	    torsionScore = runnameAndTorsion(xdb, XFworkBuf, MAX_XF_REC_LEN);

            /*torsionScore is NOT a score, holds angle-bias to modify score*/

	    if (dumpGoToAtoms) 
            {
	       fprintf(outf, "REMARK Atoms Rotated To %s\n", XFworkBuf);
	       XFdescribeAtomList(outf, alst);
	    }

            probeProc(XFworkBuf, torsionScore, alst, probestuff); 
               /* process the oriented atoms */
               /*probeProc == movingDoCommand() which calls doCommand()*/
               /*where: (see probe.c/movingDoCommand() )
                    XFworkBuf == char* orientationName  holds angle values
                    torsionScore == double scoreBias
                    alst == atom *allMovingAtoms existance implies autobondrot
                    probestuff == void *userdata == mcis
               */

            cleanupListExtras(alst, xdb, delAtomProc, deletestuff);
         }
         return; /* can quit when we are done with the 'go to's */
      }/*if using GO statements: specific orientations requested */

      /*else: working from intermixed atoms and transformations */

      /* initialize transformation values */
      for(cursor=0; cursor < xdb->nlevs; cursor++) {
         xform = xdb->level[cursor];
         xform->currVal = xform->startVal;
         rebuildTmat(xform);
      }

#define ABR_COMP_EPSILON 0.00001

      overflow = FALSE;
      cursor = xdb->nlevs - 1;
      while(cursor >= 0) 
      {
         while(! overflow) 
         {

            alst = orientedAtoms(xdb); /* atom* alst == atom *allMovingAtoms */
	    torsionScore = runnameAndTorsion(xdb, XFworkBuf, MAX_XF_REC_LEN);

            /*torsionScore is NOT a score, holds angle-bias to modify score*/

            probeProc(XFworkBuf, torsionScore, alst, probestuff); 
               /* process the oriented atoms */
               /*probeProc == movingDoCommand() which calls doCommand()*/
               /*where: (see probe.c/movingDoCommand() )
                    XFworkBuf == char* orientationName  holds angle values
                    torsionScore == double scoreBias
                    alst == atom *allMovingAtoms existance implies autobondrot
                    probestuff == void *userdata == mcis
               */
 
            cleanupListExtras(alst, xdb, delAtomProc, deletestuff);

            cursor = xdb->nlevs - 1;
            xform = xdb->level[cursor];
            xform->currVal += xform->stepVal;
            rebuildTmat(xform);
            overflow = (xform->stepVal > 0)
                       ?  (xform->currVal > (xform->endVal + ABR_COMP_EPSILON))
                       :  (xform->currVal < (xform->endVal - ABR_COMP_EPSILON));
         }

         xform = xdb->level[cursor];
         xform->currVal = xform->startVal;
         rebuildTmat(xform);

         overflow = FALSE;

         if (--cursor >= 0) 
         {
            xform = xdb->level[cursor];
            xform->currVal += xform->stepVal;
            rebuildTmat(xform);
            overflow = (xform->stepVal > 0)
                      ?  (xform->currVal > (xform->endVal + ABR_COMP_EPSILON))
                      :  (xform->currVal < (xform->endVal - ABR_COMP_EPSILON));
         }
      }
   }
}/*autobondrot()*/

   /*describeXformDB() writes header-comments to the .map file! */
void describeXformDB(FILE *outf, xformDatabase* xdb, char *cch) 
{
   int n = 0, slev = 0;
   transformData *p = NULL;
   xformAtomRecords *xa = NULL;
   biasFunction *f = NULL;
   goToRec *g = NULL;

   if (cch == NULL) { cch = "#"; }

   if (xdb) 
   {/*xdb exists, write comments to .map file*/
      slev = 0;

      fprintf(outf, "%sAUTOBONDROT with %d variables\n", cch, xdb->nlevs);

      for(p = xdb->xforms; p; p = p->next) {
	 if (p->type == NULLxform) {
#ifdef DUMP_EXTRA_XDB_DESCR
	    fprintf(outf, "%sNULL transformation\n", cch);
#else
	    fprintf(outf, "%s\n", cch);
#endif
	 }
	 else if (p->type == ROTxform || p->type == TRANSxform) {
	    if (p->type == ROTxform) {
	       fprintf(outf, "%sROT ", cch);
	    }
	    else if (p->type == TRANSxform) {
	       fprintf(outf, "%sTRANS ", cch);
	    }
	    fprintf(outf, "%s ", ((p->name)?(p->name):""));
#ifdef DUMP_EXTRA_XDB_DESCR
	    fprintf(outf, "orig %g scan from %g to %g by %g (curr %g)\n",
			   p->phase, p->startVal, p->endVal,
			   p->stepVal, p->currVal);
#else
	    fprintf(outf, "orig %g scan from %g to %g by %g\n",
			   p->phase, p->startVal, p->endVal,
			   p->stepVal);
#endif
	    fprintf(outf, "%s        AXIS(%g, %g, %g -> %g, %g, %g)\n", cch,
			   p->a1.x, p->a1.y, p->a1.z,
			   p->a2.x, p->a2.y, p->a2.z);
#ifdef DUMP_EXTRA_XDB_DESCR
	    fprintf(outf, "%s        %12g %12g %12g %12g\n", cch,
			   p->tmat.element[0][0],
			   p->tmat.element[0][1],
			   p->tmat.element[0][2],
			   p->tmat.element[0][3]);
	    fprintf(outf, "%s        %12g %12g %12g %12g\n", cch,
			   p->tmat.element[1][0],
			   p->tmat.element[1][1],
			   p->tmat.element[1][2],
			   p->tmat.element[1][3]);
	    fprintf(outf, "%s        %12g %12g %12g %12g\n", cch,
			   p->tmat.element[2][0],
			   p->tmat.element[2][1],
			   p->tmat.element[2][2],
			   p->tmat.element[2][3]);
	    fprintf(outf, "%s        %12g %12g %12g %12g\n", cch,
			   p->tmat.element[3][0],
			   p->tmat.element[3][1],
			   p->tmat.element[3][2],
			   p->tmat.element[3][3]);
#endif
	 }
	 else if (p->type == SAVExform) {
	    fprintf(outf, "%sSAVE transformation (level %d)\n", cch, ++slev);
	 }
	 else if (p->type == RESTORExform) {
	    fprintf(outf, "%sRESTORE transformation (level %d)\n", cch, slev);
	    if (slev > 0) { slev--; }
	 }
	 else {
	    fprintf(outf, "%sUNKNOWN transformation (type %d)\n", cch, p->type);
	 }
	 for(f = p->funcs; f; f = f->next) {
	    if (f->type == CONSTbiasfunc) {
	       fprintf(outf, "%s        CONST bias function %g\n", cch, f->v);
	    }
	    else if (f->type == POLYbiasfunc) {
	       fprintf(outf, "%s        POLYNOMIAL bias function %g*(x - %g)^%g\n", cch,
					  f->v, f->ph, f->freq);
	    }
	    else if (f->type == COSbiasfunc) {
	       fprintf(outf, "%s        COSINE bias function %g*(%g - cos(%g*(x - %g)))/2\n", cch,
		  f->v, f->shift, f->freq, f->ph);
	    }
	    else {
	       fprintf(outf, "%s        UNKNOWN bias function type: %d\n", cch, f->type);
	    }
	 }

	 for(xa = p->recs; xa; xa = xa->next) {
	    fprintf(outf, "%s    %4.4s%c%3.3s%2s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %s\n", cch,
	       xa->a->atomname, xa->a->altConf,
	       xa->a->r->resname, xa->a->r->chain,
	       xa->a->r->resid, xa->a->r->resInsCode,
	       xa->loc.x, xa->loc.y, xa->loc.z,
	       xa->a->occ, xa->a->bval, xa->a->r->segid);
	 }
      }
#ifdef DUMP_EXTRA_XDB_DESCR
      for(i = 0; i < MAX_XFORM_LEVELS; i++) {
         if (xdb->level[i]) {
	    fprintf(outf, "%sLevel%d (%s %d)\n", cch, i,
	       ((xdb->level[i]->name)?(xdb->level[i]->name):""),
	       xdb->level[i]->num);
	 }
      }
#endif
      n = 0;
      for(g = xdb->golist; g; g = g->next) {
	 ++n;
#ifdef DUMP_EXTRA_XDB_DESCR
	 fprintf(outf, "%sGO %d", cch, n);
	 for(i = 0; i < xdb->nlevs; i++) {
	    fprintf(outf, ": %g ", g->level[i]);
	 }
	 fprintf(outf, "\n");
#endif
      }
      fprintf(outf, "%sread %d GO statements\n", cch, n);
   }/*xdb exists, write comments to .map file*/
   else 
   {
      fprintf(outf, "NULL xformDatabase\n");
   }
}

void XFdescribeAtomList(FILE *outf, atom* a) {
   int i = 0;
   for(; a; a = a->next) {
      fprintf(outf, "ATOM   %4d %4.4s%c%3.3s%2s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f      auto\n",
	       ++i, a->atomname, a->altConf,
	       a->r->resname, a->r->chain,
	       a->r->resid, a->r->resInsCode,
	       a->loc.x, a->loc.y, a->loc.z,
	       a->occ, a->bval);
   }
}

xformDatabase* initXformDB() {
   xformDatabase *xdb = NULL;
   matrix4 *mp = NULL;
   int i = 0;

   xdb = (xformDatabase *)malloc(sizeof(xformDatabase));
   if (!xdb) {
      warn("could not create transformation database");
      return NULL;
   }

   mp = (matrix4*)malloc(sizeof(matrix4)*MAX_XFORM_STACK_DEPTH);
   if (!mp) {
      warn("could not create transformation stack");
      if (xdb) { free(xdb); }
      return NULL;
   }

   xdb->xforms = xdb->last = NULL;

   xdb->nlevs = 0;
   for(i = 0; i < MAX_XFORM_LEVELS; i++) {
      xdb->level[i] = NULL;
   }

   xdb->nstack = 0;
   for(i = 0; i < MAX_XFORM_STACK_DEPTH; i++) {
      xdb->ctmstack[i] = &(mp[i]);
   }
   initializeXFStack(xdb);

   xdb->golist = xdb->glast = NULL;

   xdb->sortedByRes = NULL;

   return xdb;
}

transformData* newTransformation(char* rec) {
   transformData* xform = NULL;
   if (rec) {
      xform = (transformData *)malloc(sizeof(transformData));
      if (!xform) {
         warn("could not create transformation");
         return NULL;
      }
      xform->next    = NULL;
      xform->recs    = NULL;
      xform->lastrec = NULL;

      xform->funcs   = NULL;

      fillinTransformationRec(xform, rec);
   }
   return xform;
}

void fillinTransformationRec(transformData* xform, char* rec) {
   char *p = NULL, *s = NULL;

   xform->name    = NULL;
   xform->num     = -1; /* filled in during indexing */
   xform->type    =  0;

   xform->phase   = 0.0;
   xform->startVal= 0.0;
   xform->endVal  = 0.0;
   xform->stepVal = 0.0;
   xform->currVal = 0.0; /* varied during processing */
   xform->a1.x    = 0.0;
   xform->a1.y    = 0.0;
   xform->a1.z    = 0.0;
   xform->a2.x    = 0.0;
   xform->a2.y    = 0.0;
   xform->a2.z    = 0.0;
   v3identityMat(&(xform->tmat));  /* varied during processing */

   p = strtok(rec, XFORMREC_DELIM);
   if (p) {
      s = p + strspn(rec, " \t");

      if (!strncasecmp(s, "NULL", 4)) {
         xform->type = NULLxform;
      }
      else if (!strncasecmp(s, "BONDROT", 7) || !strncasecmp(s, "ROT", 3)) {
         xform->type = ROTxform;
      }
      else if (!strncasecmp(s, "TRANS", 5)) {
         xform->type = TRANSxform;
      }
      else if (!strncasecmp(s, "SAVE", 4) || (s[0] == '(')) {
         xform->type = SAVExform;
      }
      else if (!strncasecmp(s, "RESTORE", 7)  || (s[0] == ')')) {
         xform->type = RESTORExform;
      }

      if (xform->type == ROTxform || xform->type == TRANSxform) {
         p = strtok(NULL, XFORMREC_DELIM);
         xform->name = strdup(p?p:"/NULL/");
         if (p) {
            p = strtok(NULL, XFORMREC_DELIM);
            if (p) {
               xform->phase = atof(p);
	       p = strtok(NULL, XFORMREC_DELIM);
               if (p) {
                  xform->startVal = atof(p);
                  p = strtok(NULL, XFORMREC_DELIM);
                  if (p) {
                     xform->endVal = atof(p);
                     p = strtok(NULL, XFORMREC_DELIM);
                     if (p) {
                        xform->stepVal = atof(p);
                        p = strtok(NULL, XFORMREC_DELIM);
                        if (p) {
                           xform->a1.x = atof(p);
                           p = strtok(NULL, XFORMREC_DELIM);
                           if (p) {
                              xform->a1.y = atof(p);
                              p = strtok(NULL, XFORMREC_DELIM);
                              if (p) {
                                 xform->a1.z = atof(p);
                                 p = strtok(NULL, XFORMREC_DELIM);
                                 if (p) {
                                    xform->a2.x = atof(p);
                                    p = strtok(NULL, XFORMREC_DELIM);
                                    if (p) {
                                       xform->a2.y = atof(p);
                                       p = strtok(NULL, XFORMREC_DELIM);
                                       if (p) {
                                          xform->a2.z = atof(p);
                                       } /* a2z */
                                    } /* a2y */
                                 } /* a2x */
                              } /* a1z */
                           } /* a1y */
                        } /* a1x */
                     } /* stepVal */
                  } /* endVal */
	       } /* startVal */
	    } /* phase */
         } /* name */
	 
	 /* insure agreement between range and step direction */
         if (xform->startVal < xform->endVal) {
	    xform->stepVal = fabs(xform->stepVal);
	 }
         else if (xform->startVal > xform->endVal) {
	    xform->stepVal = -fabs(xform->stepVal);
	 }
	 else { /* must be equal */
	    warn("TRANSFORMATION where start == end (step reset)");
	    xform->stepVal = 9999; /* error */
	 }
      } /* ROT || TRANS */
   } /* type */

   if (xform->name == NULL) { xform->name = strdup(""); }
}

transformData* deleteTransformation(transformData* xform,
            abrAtomListProc delAtomProc, void *deletestuff) {
   transformData* following = NULL;
   xformAtomRecords *p = NULL, *next = NULL;
   biasFunction *f = NULL, *nextf = NULL;

   if (xform) {
      for (p = xform->recs; p; p = next) {
         next = p->next;
         deleteXformAtomRecord(p, delAtomProc, deletestuff);
      }
      xform->recs = xform->lastrec = NULL;

      if (xform->name != NULL) {
         free(xform->name);
         xform->name = NULL;
      }

      for (f = xform->funcs; f; f = nextf) {
	 nextf = f->next;
	 f->next = NULL;
	 free(f);
      }
      xform->funcs = NULL;

      following = xform->next;
      xform->next = NULL;
      free(xform);
   }
   return following;
}

void appendTransformation(xformDatabase* xdb, transformData* xform) {
   if (xdb && xform) {
      if (xdb->xforms == NULL) {
         xdb->xforms = xform;
         xdb->last   = NULL;
      }
      if (xdb->last != NULL) {
         xdb->last->next = xform;
      }
      xdb->last = xform;
   }
}

void appendAtomRec(transformData* xform, char* rec, abrMkAtomProc mkAtom, void *atomstuff) {
   xformAtomRecords * xr = NULL;

   if (xform && rec) {
      xr = newXformAtomRecord(rec, mkAtom, atomstuff);

      if (xr) {
	 if (xform->recs == NULL) {
	    xform->recs    = xr;
	    xform->lastrec = NULL;
	 }
	 if (xform->lastrec != NULL) {
	    xform->lastrec->next = xr;
	 }
	 xform->lastrec = xr;
      }
   }
}

xformAtomRecords * newXformAtomRecord(char *rec, abrMkAtomProc mkAtom, void *atomstuff) {
   xformAtomRecords *xr = NULL;
   int getAtno(int);

   if (rec) {
      xr = (xformAtomRecords *)malloc(sizeof(xformAtomRecords));
      if (!xr) {
         warn("could not create Xform atom record");
         return NULL;
      }
      xr->next = NULL;
      xr->a = mkAtom(rec, atomstuff);

      if (xr->a) {
	                       /* record the original position on the atom */
	 xr->loc = xr->a->loc; /* so that we can continually return to it  */
      }
      else {
	 free(xr); /* may be an H when implicitH processing is being used */
	 xr = NULL;
      }
   }
   return xr;
}

void doListCleanupProcessing(xformDatabase* xdb, abrAtomListProc inputListProc, void *liststuff) {
   atom* alst = NULL, *a = NULL;
   int i = 0;

   if (xdb) {
      alst = NULL;
      if (xdb->sortedByRes) { /* (re-)compose the atom list from the sort index */
	 for(i=0; xdb->sortedByRes[i]; i++) {
	    a = xdb->sortedByRes[i]->a;
	    a->next = alst;
	    alst = a;
	 }
      }
      if (alst && inputListProc) {
	 inputListProc(alst, liststuff);
      }
   }
}

xformAtomRecords* deleteXformAtomRecord(xformAtomRecords *xr,
            abrAtomListProc delAtomProc, void *deletestuff) {
   xformAtomRecords* following = NULL;
   if (xr) {
      if (xr->a != NULL) {
         delAtomProc(xr->a, deletestuff);
         xr->a = NULL;
      }

      following = xr->next;
      xr->next  = NULL;
      free(xr);
   }
   return following;
}

void indexLevels(xformDatabase* xdb) {
   transformData* xform = NULL;
   int n = 0;
   
   if (xdb) {
      n = 0;
      for(xform = xdb->xforms; xform; xform = xform->next) {
         if (xform->type == ROTxform || xform->type == TRANSxform) {
            /* the NULL transformation does not get indexed (it is not a scan dimension) */
	    xform->num = n;
            xdb->level[n++] = xform;
         }
      }
      xdb->nlevs = n;
   }
}

void buildXDBsortedAtomRecs(xformDatabase* xdb) {
   transformData *xform = NULL;
   xformAtomRecords *xa = NULL;
   xformAtomRecords **srtvec = NULL;
   int n = 0;

   if (xdb) {
      n = 0; /* first count the number of atoms */
      for (xform = xdb->xforms; xform; xform = xform->next) {
	 for(xa = xform->recs; xa; xa = xa->next) {
	    n++;
	 }
      }

      /* build an array with extra room for a terminating null pointer */
      srtvec = (xformAtomRecords**)malloc((n + 1)*sizeof(xformAtomRecords*));
      if (!srtvec) {
	 errmsg("could not create atom sorted-order array");
	 xdb->sortedByRes = NULL;
	 return;
      }

      n = 0; /* now we store pointers to the AtomRecs in the array */
      for (xform = xdb->xforms; xform; xform = xform->next) {
	 for(xa = xform->recs; xa; xa = xa->next) {
	    srtvec[n] = xa;
	    n++;
	 }
      }
      srtvec[n] = NULL; /* array ends with a NULL */

      qsort(&(srtvec[0]), n, sizeof(xformAtomRecords*), compareXDBAtomRecs);

      xdb->sortedByRes = srtvec; /* stash away the sorted array */
   }
}

/* compare function used in sorting the transformation atom records */

int compareXDBAtomRecs(const void *avptr, const void* bvptr) {
   xformAtomRecords **ar = (xformAtomRecords **)avptr;
   xformAtomRecords **br = (xformAtomRecords **)bvptr;
   atom *a = NULL, *b = NULL;
   int cmpval = 0, t = 0;

   a = ar[0]->a;
   b = br[0]->a;
   
   /* Divide into separate residues. Atom sorting within a residue is less important. */

   if      (strcmp(a->r->chain, b->r->chain) > 0) { cmpval =  1; }
   else if (strcmp(a->r->chain, b->r->chain) < 0) { cmpval = -1; }
   else {
      if      (a->r->resid > b->r->resid) { cmpval =  1; }
      else if (a->r->resid < b->r->resid) { cmpval = -1; }
      else {
	 if      (a->r->resInsCode > b->r->resInsCode) { cmpval =  1; }
	 else if (a->r->resInsCode < b->r->resInsCode) { cmpval = -1; }
	 else {
	    if      (a->atomname[2] > b->atomname[2]) { cmpval =  1; }
	    else if (a->atomname[2] < b->atomname[2]) { cmpval = -1; }
	    else {
	       t = strcmp(a->atomname, b->atomname);
	       if (t != 0) { cmpval = t; }
	       else {
		  if      (a->altConf > b->altConf) { cmpval =  1; }
		  else if (a->altConf < b->altConf) { cmpval = -1; }
		  else {
		     if      (a->loc.x - 0.001 > b->loc.x) { cmpval =  1; }
		     else if (a->loc.x + 0.001 < b->loc.x) { cmpval = -1; }
		     else {
			if      (a->loc.y - 0.001 > b->loc.y) { cmpval =  1; }
			else if (a->loc.y + 0.001 < b->loc.y) { cmpval = -1; }
			else {
			   if      (a->loc.z - 0.001 > b->loc.z) { cmpval =  1; }
			   else if (a->loc.z + 0.001 < b->loc.z) { cmpval = -1; }
			   else {
			      /* records seem to have the same name and position */

			      note("atom duplicated in autobondrot input");

			      cmpval = 0;
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   return -cmpval; /* invert sort order because the list will be built "CONSed" (appended) */
}

void rebuildTmat(transformData* xform) {
   vector3d v;

   if (xform) {
      if (xform->type == ROTxform) {
	 vector3d axis;
	 matrix4  m1, m2, m3, m4;
	 double twistangle = 0.0;

	 v = xform->a2;
	 v3translationMat(&v, &m1);

	 twistangle = xform->currVal - xform->phase;
	 v3makeVec(&(xform->a2), &(xform->a1), &axis);
	 v3rotationMat(&axis, twistangle, &m2);

	 v3negate(&v);
	 v3translationMat(&v, &m4);

	 v3matMul(&m2, &m1, &m3);
	 v3matMul(&m4, &m3, &(xform->tmat));
      }
      else if (xform->type == TRANSxform) {

	 v3sub(&(xform->a2), &(xform->a1), &v);
	 if (v3squaredLength(&v) > 0.000001) {
	    v3scale(&v, xform->currVal - xform->phase);
	    v3translationMat(&v, &(xform->tmat));
	 }
	 else { v3identityMat(&(xform->tmat)); }
      }
   }
}

atom* orientedAtoms(xformDatabase* xdb) {
   transformData* xform = NULL;
   xformAtomRecords* xa = NULL;
   matrix4 ctm, rsltMat;
   atom* alst = NULL, *a = NULL;

   if (xdb) {
      initializeXFStack(xdb);

      v3identityMat(&ctm);

      for (xform = xdb->xforms; xform; xform = xform->next) {
	 if (xform->type == SAVExform) {
	       pushXFStack(xdb, &ctm);
	 }
	 else if (xform->type == RESTORExform) {
	    ctm = * topXFStack(xdb);
	    popXFStack(xdb);
	 }
	 else if (xform != NULLxform) {
	    v3matMul(&(xform->tmat), &ctm, &rsltMat);
	    ctm = rsltMat;
	 }

	 for(xa = xform->recs; xa; xa = xa->next) {
	    a = xa->a;

	    a->loc = xa->loc;                /* apply ctm to the original loc of each atom  */
	    v3mulPointMat(&(a->loc), &ctm); /* and use to set the new position of the atom */
	 }
      }

      alst = NULL;
      if (xdb->sortedByRes) { /* (re-)compose the atom list from the sort index */
	 int i = 0;
	 for(i=0; xdb->sortedByRes[i]; i++) {
	    a = xdb->sortedByRes[i]->a;
	    a->next = alst;
	    alst = a;
	 }
      }
   }

   return alst; /* returns the list of (pre-)sorted atoms */
}

double runnameAndTorsion(xformDatabase* xdb, char* buf, int max) {
   char *p = NULL;
   int i = 0, n = 0;
   double bias = 0.0, phi = 0.0;
   biasFunction *f = NULL;

   if (xdb && buf) {
      bias = 0.0;
      p = buf; /*dcr?: (dcr guesses) that p now points to beginning of buf */
      for(i = 0; i < xdb->nlevs; i++) {
	 if (max > 15) { /* add value to the output buffer */
	    n = sprintf(p, "%g ", xdb->level[i]->currVal);
	            /*dcr?: that sprintf wrote n char into buf starting at p*/
	    p += n; 
                    /*dcr?: p now points n nchar further along in buf*/
                    /*note the space after the %g conversion of double currVal*/
                    /*so when next char written they will have a space between*/
	    max -= n; 

            /*dcr?: max tracks available space in buf, now decremented by n */
            /*buf holds the angle values in order of their use*/
            /*it would be nice to reorder those values for output to .map file*/
            /*e.g. so plotting of phi,tau,psi could be as phi,psi,tau ! */
            /*probably this should be done at the output stage anyway*/
            /*autobondrot output is from probe.c/rawEnumerate() where */
            /*angle values are called char *rawname  which seems to == buf */
	 }

	 /* sum up any bias functions */

	 for(f = xdb->level[i]->funcs; f; f = f->next) {
	    if (f->type == CONSTbiasfunc) {
	       bias += f->v;
	    }
	    else if (f->type == POLYbiasfunc) {
	       bias += f->v * pow(xdb->level[i]->currVal - f->ph, f->freq);
	    }
	    else if (f->type == COSbiasfunc) {
	       phi = f->freq * DEG2RAD * (xdb->level[i]->currVal - f->ph);
	       bias += f->v * 0.5 * (f->shift - cos(phi));
	    }
	 }
      }
   }
   return bias;
}

void addGoToRecord(xformDatabase* xdb, char *rec) {
   char *p = NULL;
   int i = 0;
   goToRec* g = NULL;

   if (xdb && rec) {
      g = (goToRec *)malloc(sizeof(goToRec));
      if (!g) {
         warn("could not create 'go to' record");
         return;
      }

      /* initialize */
      g->next = NULL;
      for(i = 0; i < MAX_XFORM_LEVELS; i++) {
         g->level[i] = 0.0;
      }

      if (! xdb->golist) { xdb->golist = g; }
      if (xdb->glast) { xdb->glast->next = g; }
      xdb->glast = g;

      p = strtok(rec, XFORMREC_DELIM); /* type id */
      if (p) {
	 p = strtok(NULL, XFORMREC_DELIM);
	 for(i=0; p && (i < MAX_XFORM_LEVELS); i++) {
	    g->level[i] = atof(p);
	    p = strtok(NULL, XFORMREC_DELIM);
	 }
      }
   }
}

/* transformation described in the database based on the 'go to' entry */
void setOrientation(xformDatabase* xdb, goToRec* g) {
   int i = 0;

   if (xdb && g) {
      for(i = 0; i < xdb->nlevs; i++) {
         xdb->level[i]->currVal = g->level[i];
         rebuildTmat(xdb->level[i]);
      }
   }
}

void appendBiasFunction(transformData* xform, char *rec) {
   biasFunction *f = NULL;
   char *p = NULL, *s = NULL;

   if (xform && rec) {
      f = (biasFunction *)malloc(sizeof(biasFunction));
      if (!f) {
         warn("could not create 'basis function' record");
         return;
      }

      /* append to the transformation */
      f->next = xform->funcs;
      xform->funcs = f;

      /* initialize */
      f->type  = CONSTbiasfunc;
      f->v     =  0.0;
      f->ph    =  0.0;
      f->freq  =  0;
      f->shift =  0.0;

      p = strtok(rec, XFORMREC_DELIM);
      if (p) {
	 s = p + strspn(rec, " \t");

	 if (!strncasecmp(s, "CONST", 5)) {
	    f->type = CONSTbiasfunc;
	 }
	 else if (!strncasecmp(s, "POLY", 4)) {
	    f->type = POLYbiasfunc;
	    f->ph    =  0.0;
	    f->freq  =  2; /* default is a quadratic */
	 }
	 else if (!strncasecmp(s, "COS", 3)) {
	    f->type = COSbiasfunc;
	    f->ph    =  60.0; /* default cos parameters (still need a scale value) */
	    f->freq  =  3;
	    f->shift =  1.0;
	 }

	 p = strtok(NULL, XFORMREC_DELIM);
	 if (p) {
	    if (nonblankstr(p)) { f->v = atof(p); }
	    p = strtok(NULL, XFORMREC_DELIM);
	    if (p) {
	       if (nonblankstr(p)) { f->ph = atof(p); }
	       p = strtok(NULL, XFORMREC_DELIM);
	       if (p) {
		  if (nonblankstr(p)) { f->freq = atof(p); }
		  p = strtok(NULL, XFORMREC_DELIM);
		  if (p) {
		     if (nonblankstr(p)) { f->shift = atof(p); }
		  }
	       }
	    }
	 }
      }
   }
}

/* cull out any new atoms added in processing */

void cleanupListExtras(atom *alst, xformDatabase* xdb,
            abrAtomListProc delAtomProc, void *deletestuff) {

   atom *a = NULL, *prev= NULL, *next = NULL;
   transformData *xform = NULL;
   xformAtomRecords *xa = NULL;

   if (alst && xdb) {
      /* mark everything in the atom list */
      for(a = alst; a; a = a->next) {
	 a->mark = 1;
      }

      /* unmark atoms reachable from the transformation database */
      for(xform = xdb->xforms; xform; xform = xform->next) {
	 for(xa = xform->recs; xa; xa = xa->next) {
	    xa->a->mark = 0;
	 }
      }

      /* if still marked - must be extra - delete it */
      prev = NULL;
      for(a = alst; a; a = next) {
	 next = a->next;
	 if (a->mark) {
	    delAtomProc(a, deletestuff);
	    if (prev) { prev->next = next; }
	 }
	 else { prev = a; }
      }
   }
}

/* some stack manipulation functions */
void initializeXFStack(xformDatabase* xdb) {
   xdb->nstack = 1;
   v3identityMat(xdb->ctmstack[0]);
}

int pushXFStack(xformDatabase* xdb, matrix4* ctm) {
   if (xdb->nstack > 0 && xdb->nstack < MAX_XFORM_STACK_DEPTH) {
      xdb->nstack++;
      *(xdb->ctmstack[xdb->nstack - 1]) = *ctm;
      return TRUE;
   }
   else {
      return FALSE;
   }
}

void popXFStack(xformDatabase* xdb) {
   if (xdb->nstack > 1) { xdb->nstack--; }
}

matrix4* topXFStack(xformDatabase* xdb) {
   return (xdb->nstack > 0) ? xdb->ctmstack[xdb->nstack - 1] : NULL;
}
