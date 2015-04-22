/* name: probe.h                         */
/* author: J. Michael Word               */
/* date written: 2/26/96                 */
/* modified: 10/18/96, 11/6/96           */
/* purpose: compute intersecton surfaces */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/

#ifndef PROBE_H
#define PROBE_H 1

#include <stdio.h>
#include "utility.h"
#include "dots.h"
#include "abin.h"
#include "readPDBrecs.h"
#include "select.h"
#include "atomprops.h"
#include "stdconntable.h"
#include "autobondrot.h"

/* number of dot categories */
/*0 wide contact,1 close contact,2 small overlap,3 bad overlap,4 H-bonds - original, changed by SJ see below*/
/*0 wide contact,1 close contact,2 weak H bonds, 3 small overlap,4 bad overlap (0.4-0.5),5 worse overlap (>= 0.5), 6 H-bonds*/ /*04/08/2015 SJ moved weak H bonds catergory with the contacts, and separated bad overlap into bad and worse*/
/*#define NODEWIDTH 6*/ /*20111215dcr change 5 to 6 for optional weak H bonds*/
#define NODEWIDTH 7 /*04/08/2015 SJ adding category for worse overlap*/
/* the new numbers are: 0 wide contact, 1 close contact, 2 weak H-bond, 3 small overlap, 4 bad overlap, 5 worse overlap, 6 H-bonds - weak H bonds and worse overlaps only separated when LweakHbonds and LworseOverlaps is true. This is false by default*/

/* selection identifiers */
#define SET1 1
#define SET2 2
#define IGNORE_FLAG 4

/* mode types */
#define EXTERNALSURFACE   0
#define INTERSECTONCE     1
#define INTERSECTBOTHWAYS 2
#define SELFINTERSECT     3
#define DUMPATOMCOUNT     4

#define MCMCCHR   'M'   /*dcr041017*/
#define SCSCCHR   'S'   /*dcr041017*/
#define MCSCCHR   'P'   /*dcr041017*/
#define OTHERCHR  'O'   /*dcr041017*/

/* before output, dots are stored in dotNodes */
typedef struct dotNode_t {
    struct dotNode_t *next; /* link to next dot */
    atom *a;       /* dot's owner  */
    atom *t;       /* dot's cause  */
    point3d loc;   /* dot position */
    point3d spike; /* end of spike */
    int type;      /* -1 bump, 0 touch, +1 H bond */
    float gap;     /* vdw-vdw distance */
    char ptmaster; /* point master, kinemage output, (M for mc) dcr041009*/
    int dotCount;  /*added by SJ -10/07/2011 for keeping count of number of dots in the interaction, used only in writeRaw and Condense functions*/
    float angle;   /*dcr20120120 angle: dot's: parent, self, cause, esp for H */
} dotNode;

/* data structure used to determine if N and O are at chain ends */
typedef struct {
   atom *ambigN[8];    /* first residue's ambiguous N in this chain [0-3]*/
                       /* along with the last residue's Ns [4-7]         */
   atom *ambigO[8]; /*  last residue's ambiguous Os in this chain [0-7]  */
   int res_mc_oxy_cnt;       /* how many mainchain Oxygens on last res?  */
   int first, last;          /* where do the chain ids break?            */
   int Ntmarkers, Ctmarkers; /* where were N and C term markers observed? */
   int hisNHcount[4];
} chainEndData_t;

typedef struct { /* buffer to pass data to newMovingAtom */
   int filenum;
   FILE *inf;
   int close; /* do we need to close the file */
   pattern *srcPat;
   pattern *targPat;
   residue **reslstptr;
   residue *scratchRes;
} movingAtomBuildInfo;

typedef struct { /* buffer to pass data to movingDoCommand */
   int firstPass;
   int keepUnselected;
   FILE *outf;
   int method;
   atom *allMainAtoms;
   atom *waterClones;
   atomBins *abins;
   pointSet *dots;
   float probeRad;
   float density;
   float spikelen;
   int countDots;
   int rawOutput;
   int drawSpike;
   int sayGroup;
   char* groupLabel;
   int argc;
   char **argv;
   char *message;
} movingCommandInfo;

int mainProbeProc(int argc, char **argv, FILE *outf);
void doCommand(FILE *outf, int method,
   atom *allMainAtoms, atomBins *abins,
   atom *allMovingAtoms, atomBins *bbins,
   pointSet dots[], float probeRad, float density, float spikelen,
   int countDots, int rawOutput, int conFlag, char* rawname, double scoreBias,
   int drawSpike, int sayGroup, char* groupLabel,
   int argc, char **argv, char message[]);//conFlag added by SJ 10/07/2011
void descrCommand(FILE *fp, char* hdr1, char* hdr2, int argc, char **argv);
void loadDotSpheres(pointSet dots[], float density);
void unloadDotSpheres(pointSet dots[]);
atom* processCommandline(int argc, char **argv, int *method, region *bboxA,
            float *density, float *probeRad,
            int *drawSpike, float *spikelen, int *countDots,
            int *keepUnselected,
            char **srcArg, char **targArg, char **extraArg, char **ignoreArg,
            char **groupLabel, int *rawOutput, int * conFlag, int *sayGroup,
            int *addKinToFile, movingAtomBuildInfo *mabip,
            residue **reslstptr);//conFlag added by SJ - 01/07/2011
atom* loadAtoms(FILE *fp, atom *atomlist, region *boundingBox, int file,
            residue **resDataLst);
atomBins* binAtoms(atom* allAtoms, region *boundingBox, char serialNum,
            float probeRad, int keepUnselected, int selflags);
float getRadius(int at, int useCOScale);
atom * newAtom(char *rec, int file, int model, residue* resDataBlk);
int atomsClose(atom *a, atom *b, float probeRad);
int inRange(point3d *p, point3d *q, float lim);
float gapSize(point3d *p, point3d *q, float qrad);
void selectSource(atom *allAtoms, pattern *sp, int srcFlag,
	pattern *tp, int objFlag, pattern *ignorePat);
atom* findTouchingAtoms(atom *src, atom *head, atomBins *bins, float probeRad, int flag,int *ok);

void saveDot(atom *src, atom *targ, int type, point3d *loc, point3d *spike,
  dotNode *results[][NODEWIDTH], int ovrlaptype, float mingap, char ptmaster, float XHTangle);/*dcr041009,XHTangle dcr20120120 */
  dotNode * newDot(atom *src,atom *targ, point3d *loc, point3d *spike, int ovrlaptype, float gap, char ptmaster, float angle);/*dcr041009*/  /*XHTangle dcr20120120*/

void examineOneDotEach(atom *src, int type, atom *scratch,
        pointSet dots[], float probeRad, float spikelen,
        int objFlag, dotNode *results[][NODEWIDTH], atom *allMainAtoms);
                /*allMainAtoms20120120*/
void examineDots(atom *src, int type, atom *scratch,
		pointSet dots[], float probeRad, float spikelen,
		int objFlag, dotNode *results[][NODEWIDTH]);
void markBonds(atom *src, atom *neighbors, int distcount, int max);
int dotType(atom *src, atom *atomList, int recalcOnly);
void genDotIntersect(atom *allMainAtoms, atomBins *abins,
			atom *allMovingAtoms, atomBins *bbins,
			pointSet dots[],
			float probeRad, float spikelen,
			int srcFlag, int targFlg, dotNode *results[][NODEWIDTH]);
void genDotSurface(atom *allMainAtoms, atomBins *abins,
			atom *allMovingAtoms, atomBins *bbins,
			pointSet dots[],
			float probeRad, float spikelen, int srcFlag,
			dotNode *results[][NODEWIDTH]);
void surfDots(atom *src, int type, atom *scratch,
	pointSet dots[], float probeRad,
	float spikelen, dotNode *results[][NODEWIDTH]);
void initResults(dotNode *results[][NODEWIDTH]);
void freeResults(dotNode *results[][NODEWIDTH]);
void writeOutput(FILE *outf, char *groupname, dotNode *results[][NODEWIDTH],
                 int drawSpike, int method, char *extramastername, float probeRad);
                /*041020 method for better kinemage keywords*/
                /*060129 extra master name controls original vs fitted dots*/ 
                /*probeRad added 20111220dcr*/
void writeAltFmtO(FILE *outf, int showBegin, int showEnd,
      char* groupname, dotNode *results[][NODEWIDTH], int drawSpike);
void writeAltFmtXV(FILE *outf, int showBegin, int showEnd,
      char* groupname, dotNode *results[][NODEWIDTH], int drawSpike);

void writeRaw(FILE *outf, char *groupname, dotNode *results[][NODEWIDTH],
               float rp, char*s, float,int conFlag); 
dotNode * Condense(dotNode * head,int conFlag); //added 10/04/11 - SJ for condensing the rawOutput
void enumerate(FILE *outf, char* groupname, dotNode *results[][NODEWIDTH],
               float probeRad, int method,
	       int nsel, int spike, int outdots, int numSkinDots,
	       float density);
void rawEnumerate(FILE *outf, char* groupname, dotNode *results[][NODEWIDTH],
               int method, int nsel, int spike, int outdots, int numSkinDots,
	       float density, char *namestring, char *rawname, double scoreBias);
int countSelected(atom *allAtoms, int srcFlag);
int enumDotSkin(atom *allMainAtoms, atomBins *abins,
	    atom *allMovingAtoms, atomBins *bbins,
	    pointSet dots[], int srcFlag);
int countSkin(atom *src, atom *scratch, pointSet dots[]);
atom* updateHydrogenInfo(FILE *outf, atom *allMainAtoms,   atomBins *abins,
				    atom *allMovingAtoms, atomBins *bbins,
				    int selectedFlag, int mustSaveMainWater);
int dotClassIndex(int t, float mingap);
void fixupLongBondChains(atom *src, atom *neighbors, int cutoff);
float dot2bullsEye(point3d *dot, atom *src, atom *targ);
float dot2srcCenter(point3d *dot, atom *src, atom *targ);
float kissEdge2bullsEye(float ra, float rb, float rp);
char* assignGapColorForKin(float gap, int class);
char* assignGapColorForO(float gap, int class);
char* assignGapColorForXV(float gap, int class);
char* convertKinColorToO(char* incolor);
char* convertKinColorToXV(char* incolor);

atom * newMovingAtom(char *rec, void* userdata);
void deleteMovingAtom(atom *a, void*);
void movingAtomListProcessing(atom *atomlst, void *userdata);
void movingDoCommand(char* name, double scoreBias, atom *allMovingAtoms, void* userdata);

void initEndData(chainEndData_t *ed);
void hisNHcheck(chainEndData_t *ed, atom *atomlist, int rescnt);
void resCheck(chainEndData_t *ed, atom *atomlist, int rescnt);
void CtermCheck(chainEndData_t *ed, int rescnt, int isChainEnd);
void noticedNt(chainEndData_t *ed, int rescnt);
void noticedCt(chainEndData_t *ed, int rescnt);
void NtermCheck(chainEndData_t *ed, int rescnt, int isChainEnd);
void ProcessResInfo(chainEndData_t *ed, atom *a);


ringInfo * newRingInfo(ringInfo **head, point3d* ctr, point3d* norm);
void deleteRingInfoList(ringInfo *ri);
residue * newResidueData();
void deleteResidueData(residue *r);
void disposeListOfResidues(residue *theRes);
void deleteAtom(atom *a);
void disposeListOfAtoms(atom *theAtom);
int resDiffersFromPrev(residue *r1, residue *r2);
void dumpRes(residue *theRes);
void dump_changes(FILE *outf);
void countsummary(FILE *outf, char* modestr, int Lines, int Pass); /*dcr041101*/
#endif
