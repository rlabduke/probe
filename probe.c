
/*{{{--probe.c general comments                                              */

/* name: probe.c                         */
/* author: J. Michael Word               */
/* date written: 2/26/96                 */
/* purpose: compute intersecton surfaces */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1996-2013 J. Michael Word                       */
/*****************************************************************/

/*Updated to work with the remediated PDB names rmi 070727 */
/* Essentially ennumerated additional possible atom names and added */
/* the residue name info to the routine that assigns element type */

/*Jan 2006 DCR compromise comments to enable jEdit type folding {{{ }}}  */
/* without messing up vi type beginning-end bracket identification shift % */
/* So jEdit folds have to stay within code regions, esp subroutines 060129*/

/*Oct 2004 DCR modifications to version 2.10.031014  starting with dcr041007 */
/*dcr: sorry Mike, but changing to Berkeley-Altman brackets really helps me!*/
/*i.e. some routines changed: opening and closing brackets vertically in line*/
/*also substituted spaces for tabs which are not robust across text editors*/
/*also tend to shorten lines to <=80 col, standard card-image window size*/
/*which many text editors use as a default width*/

/*050111:  in trying to track behavior of autobondrot mode, the lack of */
/*global mode state is a severe limitation.  */
/*probe trys to be so object oriented as to not have global modal states */
/*but given the existance of a whole slew of probe.c global parameters,  */
/*a few more would only help clear up dependencies and program flow! */
/*hence logical Lautobondrot as a probe.c global*/


/* modifications: see dump_changes() */
/*}}}--probe.c general comments */

/*{{{--includes and declarations   (globals for probe.c)                     */
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "probe.h"

#define INLINE_FOR_SPEED 1

static char *versionString = "probe: version 2.16.160404, Copyright 1996-2016, J. Michael Word";
/*static char *versionString = "probe: version 2.16.130520, Copyright 1996-2013, J. Michael Word";*/
/*"probe: version 2.15.130427, merged probeVector, Copyright 1996-2013, J. Michael Word";*/
/*"probe: version 2.14.130116, Copyright 1996-2013, J. Michael Word";*/
/*"probe: version 2.13.120907, Copyright 1996-2012, J. Michael Word";*/
/*"probe: version 2.14.120109, Copyright 1996-2011, J. Michael Word";probeVector version*/
/*"probe: version 2.13.110909, Copyright 1996-2011, J. Michael Word";*/
/*"probe: version 2.13.110830, Copyright 1996-2011, J. Michael Word";*/
/*"probe: version 2.12.110413, Copyright 1996-2007, J. Michael Word";*/
/*"probe: version 2.12.070821, Copyright 1996-2007, J. Michael Word";*/
/*"probe: version 2.11.061018, Copyright 1996-2006, J. Michael Word";*/
/*"probe: version 2.11.060831, Copyright 1996-2006, J. Michael Word";*/
/*minor work: 2.11.060212*/
/*"probe: version 2.11.060129, Copyright 1996-2006, J. Michael Word";*/
/*"probe: version 2.11.050121, Copyright 1996-2005, J. Michael Word";*/
/*"probe: version 2.11.041112, Copyright 1996-2004, J. Michael Word";*/
/*"probe: version 2.10.031014dcr041101, Copyright 1996-2004, J. Michael Word";*/
/*"probe: version 2.10  10/14/2003, Copyright 1996-2003, J. Michael Word";*/
   /*jmw & dcr agreement on version name and maintenance by dcr 041110*/
static char *shortVersionStr = "probe.2.16.160404";
/*static char *shortVersionStr = "probe.2.16.130520";*/
/*static char *shortVersionStr = "probe.2.15.130427";*/
/*static char *shortVersionStr = "probe.2.14.130116";*/
/*static char *shortVersionStr = "probe.2.13.120907"; gjk changed %s for OUTCOLor to make -OUT work again*/
/*static char *shortVersionStr = "probe.2.14.120109";probeVector version*/
/*static char *shortVersionStr = "probe.2.13.110909";*/
/*static char *shortVersionStr = "probe.2.13.110830";*/
/*static char *shortVersionStr = "probe.2.12.110413";*/
/*static char *shortVersionStr = "probe.2.11.061018";*/
/*static char *shortVersionStr = "probe.2.11.060831";*/
/*static char *shortVersionStr = "Probe V2.11.060129";*/
/*static char *shortVersionStr = "Probe V2.11.050121";*/
/*static char *shortVersionStr = "Probe V2.11.041112";*/ /*041112 version change*/
/*static char *shortVersionStr = "Probe V2.10.031014dcr041101";*/
/*static char *shortVersionStr = "Probe V2.9 ( 10/14/2003)";*/
static char *referenceString = "Word, et. al. (1999) J. Mol. Biol. 285, 1711-1733.";
static char *electronicReference = "http://kinemage.biochem.duke.edu/";

static int LweakHbonds = FALSE; /*111215dcr global Logical to group weak H-bonds*/
static int LworseOverlap = FALSE; /* 04/09/2015 SJ - to separate the bad and worse overlaps if set to true*/
static int LOneDotEach = FALSE; /*111013 global One Dot for Each neighbor */
static int Lnearest = FALSE;    /*111021dcr selection only on nearest neighbor*/
static int Lfilternearest = FALSE;/*111022dcr final filter nearest by extraPat*/
static int Lfilterocc1 = FALSE; /*20120120dcr final filter use only OCC==1   */
static int Lfilternoalt = FALSE; /*20120120dcr final filter exclude any alts */
static int Lgapbins = FALSE;    /*111020dcr gapbin ptmasters*/
static int LXHvector = FALSE;   /*20120120dcr for H: vector parent--H--Target angle*/
static int Ldotdump    = FALSE; /* global flag for Ldotdump messages in examineDots*/
static int Lhonlyocclude = FALSE; /*130427dcr global flag: only Hatoms can occlude*/
static pattern *extraPat = NULL;/*111022dcr finally filter nearest by extraPat*/
static pattern *causePat = NULL;/*111022dcr initial select nearest by targPat*/
static pattern *sourcePat= NULL;/*111021dcr select src bases on srcPat */

static int LMasterName = FALSE; /*global extra master={name} on lists 060129*/
static int LMergeContacts = TRUE; /*global combine wide & close contacts 060129*/
/*static int Lautobondrot = FALSE;*/ /* global flag for AUTOBONDROT mode 050111*/
static int ImplicitH = FALSE; /* global controling dot radii     */
static int UsePolarH = TRUE;  /* global controling VDW radii of polar Hs */
static int Verbose   = TRUE;  /* global flag for output messages */
static int DoMcMc    = FALSE; /* global controling mainchain->mainchain dots */
static int DoHet     = TRUE;  /* global controling non-Water HET dots   */
static int DoH2O     = TRUE;  /* global controling Water dots           */
static int DoWatWat  = FALSE; /* global controling water->water dots */
static int LimitDots = TRUE;  /* global controling chopping around bumps*/
static int LensDots  = FALSE; /* global controling lens keyword in kin file*/
static int ColorGap  = TRUE;  /* global controling coloration of Gaps   */
static int ByNABaseColor= FALSE;  /* global controling coloration of dots */
static int DumpNewHO = FALSE; /*global controling dump of water newH atom data*/
static int Maxbonded = 4;     /* global controling exclusion bond count */
static int UseStdBond= FALSE; /* global flag: if TRUE std AA and NA bonding assumed */
static int UseHParent= TRUE;  /* global flag: bonding based on parent names on std H atoms */
static int ShowTicks = TRUE;  /* global controling display of residue name ticker */
static int OldStyleU = FALSE; /* global controling -u output (true gives kiss edge stuff) */
static int SplitBySegID=FALSE;/* global controling the splitting of residues (-SEGID) */
static int HB2aromFace=TRUE;  /* global controling aromatic face Hbonding */
static int AtomMasters=FALSE; /* global controling use of masters */
static int PermitCHXHB=FALSE; /* global controling use CHX (CHO) hbonds */
static int DebugLevel= 0;     /* global controling debugging information */
int NuclearRadii = FALSE; /* global controling Hydrogen vdW radii, cannot be static */
#ifdef JACKnANDREA
static int writeHbonds = FALSE; /* global flag: if true, write out Hbonds every time they are detected */
#endif

static int OutputFmtType = 0;    /* global selecting output format */
                                 /*        0 ==> kinemage format   */
                                 /*        1 ==> O format          */
                                 /*        2 ==> XtalView format   */
                                 /*        3 ==> oneline multi:count:dcr041101*/
static int OutputHBs     = TRUE; /* global controling display of contact categories */
static int OutputClashes = TRUE; /* global controling display of contact categories */
static int OutputVDWs    = TRUE; /* global controling display of contact categories */
static int OnlyBadOut = FALSE; /* global controling display of contact categories dcr041010*/
static int ContactSUMMARY = FALSE; /* global controling output of contact summaries dcr041101*/

static float Min_regular_hb_cutoff=0.6; /* globals controling hbond cutoff */
static float Min_charged_hb_cutoff=0.8; /* defaults set in processCommandline() */
static float RadScaleFactor       =1.0; /* global VDW radius scale Factor r*f */
static float RadScaleOffset       =0.0; /* global VDW radius scale Offset r+o */
static float CORadScale   =(1.65/1.75); /* global VDW radius scale Factor for C=O */
static float GAPweight            =0.25;/* global raw GAP score weight */
static float BUMPweight           =10.0;/* global raw BUMP score scale Factor */
static float HBweight             = 4.0;/* global raw HBond score scale Factor */
static float CHOHBfactor          = 0.5;/* global CH..O HBond score scale Factor */

static float LowGoodCut           =-0.4;/* global cutoff for good bin */
static float HighGoodCut          = 0.25;/* global cutoff for good bin */
static float WorseBumpCut         =-0.5; /* SJ 04/10/2015 - global cutoff for worse overlap*/

static float OccupancyCutoff      =0.02; /* global occupancy below which atom has presence but no clashes */
   /*unfortunately, one needs intervening atoms to avoid clashes between atoms*/
   /*that are separated by <= Maxbonded, which thus needs to pass through any*/
   /*atom no matter what it's own occupancy happens to be */
   /*050118 atom->occ < OccupancyCutoff neither clash nor transfer bonded info*/
   /* so flanking atoms show horrible spurious clashes (e.g. 1BKR lys 5 lys 78*/
   /*mechanism seems to be atom->flags | IGNORE_FLAG */
   /*050119 atoms irrespective of occ now put into neighbor bins for bonding*/
   /* but atom->occ < OccupancyCutoff are not allowed to show clashes*/

static int OccupancyCriteria = TRUE; /*global: use occupancy criteria 060831*/
   /*explicit default: unset by -noocc flag, 060831*/

static char *OutPointColor = "gray";

static pointSet COdots; /* adjustible radius dots for C=O carbon */

/*Old: dcr041020 NODEWIDTH defined 5 probe.h, hard code 6 to match initialization*/
/*Old: 0 wide contact,1 close contact,2 small overlap,3 bad overlap,4 H-bonds*/
/*Old: 5 running sum which can be tested for existance of any of these spikes*/
/* New 04/15/2015 SJ: 0 wide contact,1 close contact,2 weak H bonds, 3 small overlap, 4 bad overlap, 5 worse overlap, 6 H-bonds*/
/* changing the size of the array to 8*/
static long mcmccont[2][8] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; /*dcr041017,dcr041020, SJ 04/15/2015*/
static long scsccont[2][8] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; /*dcr041017,dcr041020, SJ 04/15/2015*/
static long mcsccont[2][8] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; /*dcr041017,dcr041020, SJ 04/15/2015*/
static long othrcont[2][8] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; /*dcr041017,dcr041020, SJ 04/15/2015*/
static long summcont[2][8] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; /*dcr041017,dcr041020, SJ 04/15/2015*/

/* comment: the above number 8 should not be hardcoded, maybe initialize in the main function? - SJ 04/09/2015*/

static char inputfilename[256]; /*dcr041023 global, subroutines can report it*/
static int  modelNumber[100]; /*global record of incoming model # 041114*/
static int  modelCount = 0;   /*global count of such records, all files 041114*/
static int  modelLimit = 0;   /*global loop's flag for multiple models 041114*/
  /*model counting simple-minded only good for one input file cases like SELF*/
static int  modelSrc  = 0; /*global model specified for source 041114*/
static int  modelTarg = 0; /*global model specified for target 041114*/
static int  modelToProcess = 0; /*global model specified for processing 041114*/

#define ATOMS_IN_THE_SAME_RES(aa, bb) IS_THE_SAME_RES((aa)->r, (bb)->r)
/*041112 IS_THE_SAME_RES() cannot avoid modeltest, i.e. implied Lmodeltest==1*/
#define IS_THE_SAME_RES(ra, rb)  \
 (  ((ra) == (rb))                \
  ||(   (ra)->resid      == (rb)->resid      \
     && (ra)->resInsCode == (rb)->resInsCode \
     && strcmp((ra)->chain,  (rb)->chain) == 0     \
     && (ra)->model      == (rb)->model )    \
  && ((SplitBySegID == FALSE) || ( strcmp((ra)->segid, (rb)->segid) == 0) ) )

#define TRACEATOM(where, txt, atptr, endtxt) \
(fprintf((where), "%s{%s%c %s%d%c}[%d]%s", (txt),\
         (atptr)->atomname, (atptr)->altConf,\
         (atptr)->r->resname, (atptr)->r->resid,\
         (atptr)->r->resInsCode, (atptr)->r->rescnt, (endtxt)))
/*}}}--includes and declarations */

/*{{{main()***************************** MAIN ********************************/
int main(int argc, char **argv) {
   int rc = 0;
   FILE *outf = stdout;

   initalizeAtomTbl(); /*atomprops.c*/
   initStdConnTable(); /*stdconntable.c*/

   inputfilename[0] = '\0'; /*global dcr041023*/
   rc = mainProbeProc(argc, argv, outf);
   return rc;
}
/*}}}main()__________________________________________________________________*/

/*{{{mainProbeProc()******** called from main() ******************************/
int mainProbeProc(int argc, char **argv, FILE *outf)
{/*mainProbeProc()*/
   int method;
   int keepUnselected;
   region boundingBoxA;
   atom *allMainAtoms = NULL;
   atomBins *abins = NULL;
   pointSet dots[NUMATOMTYPES];
   char *srcArg=NULL, *targArg=NULL, *extraArg=NULL, *ignoreArg=NULL; /* extraArg dcr 111022*/
   char *groupLabel=NULL;
   pattern *srcPat = NULL, *targPat = NULL, *ignorePat = NULL;
   float density, probeRad, spikelen;
   int drawSpike, countDots, rawOutput, sayGroup, addKinToFile, conFlag; //conFlag by SJ - 10/07/2011
   char message[200];
   movingAtomBuildInfo mabis;
   residue *resLst = NULL;

   /*mabis moving atom build info structure */
   /*buffer to pass data to newMovingAtom for autobondrot*/

   allMainAtoms = processCommandline(argc, argv, &method,
        &boundingBoxA, &density, &probeRad,
        &drawSpike, &spikelen, &countDots, &keepUnselected,
        &srcArg, &targArg, &extraArg, &ignoreArg, &groupLabel, &rawOutput, &conFlag,
        &sayGroup, &addKinToFile, &mabis, &resLst);//conFlag by SJ -10/07/2011

        /*called with address of mabis, i.e. mabip moving atom build info ptr */
        /*autobondrot mode: not read in any mobile atoms yet*/
        /*though may have static set of atoms already in. */

/* dumpRes(resLst); */

   if ((allMainAtoms == NULL) && (mabis.inf == NULL))
   {
      /*static atoms read into allMainAtoms linked list, */
      /* and/or autobondrot file name for mobile atoms held in mabis.inf */

      sprintf(message, "No atom data in input.");
      note(message);
   }
   else
   {/*there are atoms for probe to work on...*/
      if (Verbose)
      {/*Verbose: lots of stuff to screen...*/

	 note(versionString);
	 sprintf(message, "density: %g, probe radius: %g, maxBonded: %d",
			density, probeRad, Maxbonded);
	 note(message);
	 sprintf(message, "regular HB cutoff: %g, charged HB cutoff: %g",
			Min_regular_hb_cutoff, Min_charged_hb_cutoff);
	 note(message);
	 sprintf(message,
	 "Dot gap bins: low to %g, %g to 0, 0 to %g, %g to high, & Hbonds",
			LowGoodCut, LowGoodCut, HighGoodCut, HighGoodCut);
	 note(message);
	 if (RadScaleOffset < -0.0001 || RadScaleOffset > 0.0001) {
	    sprintf(message, "Add %g to VDW Rradius",
			RadScaleOffset);
	    note(message);
	 }
	 if (RadScaleFactor < 0.9999 || RadScaleFactor > 1.0001) {
	    sprintf(message, "Scale VDW Rradius by %g",
			RadScaleFactor);
	    note(message);
	 }
	 sprintf(message, "C=O carbon VDW scaled by %.3f to a radius of %g A",
			   CORadScale, CORadScale*1.75);
	 note(message);
	 if (rawOutput || countDots) {
	    sprintf(message, "Score Weights: gapWt=%g, bumpWt=%g, HBWt=%g",
			GAPweight, BUMPweight, HBweight);
	    note(message);
     }
	 if (drawSpike) {
	    sprintf(message, "draw spike, len: %g", spikelen);
	    note(message);
	 }
	 if (ImplicitH) {
	    note("implicit hydrogens");
	 }
	 if (PermitCHXHB) {
	    sprintf(message, "CH..O type hbonds recognized (scale factor %g)", CHOHBfactor);
	    note(message);
	 }
	 if (!DoMcMc) {
	    note("Excluding mainchain-mainchain interactions");
	 }

	 if (!DoHet && !DoH2O) {
	    note("Excluding both HET groups and waters");
	 }
	 else {
	    if (!DoHet) {
	       note("Excluding non-water HET groups");
	    }
	    if (!DoH2O) {
	       note("Excluding waters");
	    }
	    else if (!DoWatWat) {
	       note("Excluding water-water interactions");
	    }
	 }
	 if (LimitDots) {
	    note("Limiting dots from bumps to max non-bump cap");
	 }

	 if (rawOutput) { /* nothing required for raw format */
	 }
	 else if (OutputFmtType == 1) {
	    note("Formatting for display in O");
	 }
	 else if (OutputFmtType == 2) {
	    note("Formatting for display in XtalView");
	 }
	 else if (OutputFmtType == 3) { /*dcr041101*/
	    note("outputing one line colon:deliminated:counts:by:severity:type:");
	 }

	 if (ByNABaseColor) {
	    if (ColorGap) {
	       note("Coloring dots by gap (& list color by nucleic acid base)");
	    }
	    else {
	       note("Coloring dots by nucleic acid base");
	    }
	 }
	 else if (ColorGap) {
	    note("Coloring dots by gap (& list color by atom)");
	 }
	 if (!UsePolarH) {
	    note("Using old (long) radii for polar H atoms");
	 }
	 if (!OutputHBs) {
	    note("Attention: ***Output not generated for H-Bond contacts***");
	 }
	 if (!OutputClashes) {
	    note("Attention: ***Output not generated for clashes***");
	 }
	 if (!OutputVDWs) {
	    note("Attention: ***Output not generated for van der Waals contacts***");
	 }
	 if (OnlyBadOut) {
	    note("Attention: ***Output only generated for Bad Clashes***");/*dcr041010*/
     }
	 if (!keepUnselected) {
	    note("Attention: ***Unselected atoms dropped (ignored)***");
	 }
    }/*Verbose: lots of stuff to screen...*/

    /*now do some work...*/
    if(   addKinToFile  && (OutputFmtType == 0)
       && (! countDots) && (! rawOutput) ) {
       fprintf(outf, "@kinemage 1\n");
    }

    srcPat = getPat(srcArg,   "pattern 1", Verbose);
    sourcePat = getPat(srcArg,   "pattern 1", Verbose);/*111021dcr*/
    modelSrc = modelInPattern(srcPat); /*parse.c 041114*/
    if(Verbose && modelSrc > 0)
    {
       sprintf(message, "specified src  model==%d", modelSrc);
       note(message);
    }
    targPat = getPat(targArg,  "pattern 2", Verbose);
    causePat = getPat(targArg,  "pattern 2", Verbose); /*111021dcr nearest*/
    modelTarg = modelInPattern(targPat); /*parse.c 041114*/
    if(Verbose && modelSrc > 0)
    {
       sprintf(message, "specified targ model==%d", modelTarg);
       note(message);
    }
    extraPat = getPat(extraArg,"extra", Verbose); /*111023dcr*/
    ignorePat = getPat(ignoreArg,"ignore",    Verbose);

    loadDotSpheres(dots, density);

      if (allMainAtoms) /*atoms usually read in from processCommandline()*/
      {  /*build the bins for all static atoms*/
         selectSource(allMainAtoms, srcPat, SET1, targPat, SET2, ignorePat);
         abins = binAtoms(allMainAtoms, &boundingBoxA,
                          'a', probeRad, keepUnselected, SET1|SET2);
      }
      else /*but autobondrot mode may not yet have read in any atoms*/
      {/*minimal boundingBoxes setup*/
         boundingBoxA.min.x = boundingBoxA.min.y = boundingBoxA.min.z = 0.0;
         boundingBoxA.max.x = boundingBoxA.max.y = boundingBoxA.max.z = 0.1;
         abins = initBins('a', &boundingBoxA, 1); /* dummy bbox */
      }

      if (mabis.inf == NULL) /*pseudo modal flag for autobondrot*/
      {/*there is NOT an autobondrot atom info input file */
         if (! ImplicitH)
         { /*preliminary step acts on all atoms*/
            updateHydrogenInfo(outf, allMainAtoms, abins,
                               NULL, NULL, SET1|SET2, FALSE);
            /*(returned param used on separate call for autobondrot processing)*/
            /*e.g. creates: newH->elem = atomHOd , dummy Hydrogen to make an H-Bond*/
         }
         /*now do the real work: */
	 doCommand(outf, method,
	    allMainAtoms, abins, NULL, NULL,
	    dots, probeRad, density, spikelen,
	    countDots, rawOutput, conFlag, "", 0.0, drawSpike,
	    sayGroup, groupLabel, argc, argv, message);
      }
      else
      {/*setup to call autobondrot*/
	 xformDatabase* xdb = NULL;
	 movingCommandInfo mcis;

	 mcis.firstPass  = TRUE;
	 mcis.keepUnselected = keepUnselected;
	 mcis.outf       = outf;
	 mcis.method     = method;
	 mcis.allMainAtoms = allMainAtoms;
	 mcis.waterClones= NULL;
	 mcis.abins      = abins;
	 mcis.dots       = dots;
	 mcis.probeRad   = probeRad;
	 mcis.density    = density;
	 mcis.spikelen   = spikelen;
	 mcis.countDots  = countDots;
	 mcis.rawOutput  = rawOutput;
	 mcis.drawSpike  = drawSpike;
	 mcis.sayGroup   = sayGroup;
	 mcis.groupLabel = groupLabel;
	 mcis.argc       = argc;
	 mcis.argv       = argv;
	 mcis.message    = message;

	 mabis.srcPat    = srcPat;
	 mabis.targPat   = targPat;

#define RAW_HEADER_COMMENT "#"

        descrCommand(outf, RAW_HEADER_COMMENT, RAW_HEADER_COMMENT, argc, argv);

        /*now actually read in both the mobile atoms and transformation info*/

	xdb = readTransformationDatabase(mabis.inf, outf, newMovingAtom, &mabis,
                                       movingAtomListProcessing, NULL,
                                       RAW_HEADER_COMMENT);

/*        mabis.inf                == FILE *inf                      */
/*        outf                     == FILE *outf                     */
/*        newMovingAtom            == abrMkAtomProc mkAtom           */
/*        &mabis                   == void *atomstuff                */
/*        movingAtomListProcessing == abrAtomListProc inputListProc  */
/*        NULL                     == void *liststuff                */
/*        RAW_HEADER_COMMENT       == char *cch                      */

         /*autobondrot/readTransformationDatabase()             */
         /* also calls autobondrot/describeXformDB()            */
         /*  which writes the header-comments to the .map file! */

	 if (mabis.close) { fclose(mabis.inf); }

	 autobondrot(stderr, xdb, movingDoCommand, &mcis,
		  deleteMovingAtom, &mabis, Verbose);

         /*autobondrot.c/autobondrot() is the call to do the autobondrot stuff*/
         /*movingDoCommand is the name of the probeProc() called from there */
         /*probe.c/movingDoCommand() in turn calls probe.c/doCommand()*/

	 discardXformDB(xdb, deleteMovingAtom, &mabis);
      }/*setup to call autobondrot*/

      /*after finish, then drop through to here to release memory, etc.*/

      disposeBins(abins);         abins = NULL;
      freeDotSphere(&COdots);
      unloadDotSpheres(dots);
      freePattern(ignorePat); ignorePat = NULL;
      freePattern(targPat);     targPat = NULL;
      freePattern(srcPat);       srcPat = NULL;
      freePattern(sourcePat); sourcePat = NULL; /*111023dcr*/
      freePattern(causePat);   causePat = NULL; /*111023dcr*/
      freePattern(extraPat);   extraPat = NULL; /*111023dcr*/

      if (Verbose)
      {
	 note("If you publish work which uses probe, please cite:");
	 note(referenceString);
	 sprintf(message, "For more information see %s", electronicReference);
	 note(message);
      }
   }/*there are atoms for probe to work on...*/

   disposeListOfResidues(resLst);          resLst = NULL;
   disposeListOfAtoms(allMainAtoms); allMainAtoms = NULL;

   return 0; /*to probe.c/main() */
}/*mainProbeProc()*/
/*}}}mainProbeProc()_________________________________________________________*/

/*{{{newMovingAtom()**********************************************************/
atom * newMovingAtom(char *rec, void *userdata)
{
   atom* a = NULL;
   movingAtomBuildInfo *m = (movingAtomBuildInfo *)userdata;

   if (rec == NULL) { return NULL; }

   if (m->scratchRes == NULL) {
      m->scratchRes = newResidueData();
   }

   /*atom * newAtom(char *rec, int file, int model, residue * resDataBuf)*/
   a = newAtom(rec, m->filenum, 0, m->scratchRes);
   if (a) {
      a->next = NULL;
      if (ImplicitH && isHatom(a->elem)) {
	 deleteAtom(a);  /* filter out implicit hydrogens */
	 return NULL;
      }

      selectSource(a, m->srcPat, SET1, m->targPat, SET2, NULL);
             /* note: new moving atom can not be ignored using -IGNORE flag */

      if (*(m->reslstptr) == NULL) {/* first residue goes on residue list */
	 *(m->reslstptr) = m->scratchRes;
	 m->scratchRes = NULL;
      }
      else { /* otherwise we have to compare residue data blocks */
	 if (resDiffersFromPrev(*(m->reslstptr), m->scratchRes)) {
	    m->scratchRes->nextRes = *(m->reslstptr);
	    *(m->reslstptr) = m->scratchRes;
	    m->scratchRes = NULL;

	 }
	 else {
	    a->r = *(m->reslstptr); /* same, so point to prior block */
	    a->r->a = a; /* makes this atom the head of atom list for this res */
	 }
      }
   }
   return a;
}
/*}}}newMovingAtom()_________________________________________________________*/

/*{{{deleteMovingAtom()*******************************************************/
void deleteMovingAtom(atom *a, void *userdata)
{
   movingAtomBuildInfo *m = (movingAtomBuildInfo *)userdata;
   if (m && m->scratchRes) {
      deleteResidueData(m->scratchRes); /* clean up extra memory */
      m->scratchRes = NULL;
   }
   deleteAtom(a);
}
/*}}}deleteMovingAtom()______________________________________________________*/

/*{{{movingDoCommand()********************************************************/
void movingDoCommand(char* orientationName, double scoreBias,
		  atom *allMovingAtoms, void *userdata)
{/*movingDoCommand() for autobondrot*/
   /*orientationName holds the angle values as a character string*/

/* NOTE: allMainAtoms/abins & allMovingAtoms/bbins must be disjoint */
/*       sets of atoms (none in common) or allMovingAtoms/bbins can */
/*       be NULL.                                                   */
/*       allMovingAtoms refers to autobondrot set of atoms          */

   movingCommandInfo *m = (movingCommandInfo *)userdata;
   region boundingBoxB, nonMovingBB;
   atom *a = NULL;
   atomBins *bbins = NULL;

   /*autobondrot.c: alst == atom *allMovingAtoms existance implies autobondrot*/
   if (allMovingAtoms) {
      boundingBoxB.min = allMovingAtoms->loc;
      boundingBoxB.max = allMovingAtoms->loc;
      for(a = allMovingAtoms; a; a = a->next) {
	 updateBoundingBox(&(a->loc), &boundingBoxB);
      }

      /* Expand bounding box because of extra water hydrogens */
      /*                        added in updateHydrogenInfo() */
      nonMovingBB.min = m->abins->min;
      nonMovingBB.max = m->abins->max;
      addBBox2BBox(&nonMovingBB, &boundingBoxB);

      bbins = binAtoms(allMovingAtoms, &boundingBoxB, 'b', m->probeRad,
			m->keepUnselected, SET1|SET2);

      if (m->abins == NULL || bbins == NULL) {
	 halt("no atoms for autobondrot processing");
      }

      if (! ImplicitH) {
	 if (m->firstPass) {
	    m->waterClones = updateHydrogenInfo(m->outf, m->allMainAtoms,
                             m->abins, allMovingAtoms, bbins, SET1|SET2, TRUE);
	    m->firstPass = FALSE;
	 }
	 else {
	    /*the DONOR prop tells updateHydrogenInfo() to build H? atoms*/
	    for(a = m->waterClones; a; a = a->next) {
	       if ((a->props & AMBIGWATER_PROP) && (a->elem == atomO)) {
		  a->props |= DONOR_PROP; /* turn water donor property back on*/
	       }
	    } /* since each pass through updateHydrogenInfo turns it off */

	    for(a = allMovingAtoms; a; a = a->next) {
	       if ((a->props & AMBIGWATER_PROP) && (a->elem == atomO)) {
		  a->props |= DONOR_PROP; /* turn water donor property back on*/
	       }
	    } /* since each pass through updateHydrogenInfo turns it off */

	    updateHydrogenInfo(m->outf, m->waterClones, m->abins,
	       allMovingAtoms, bbins, SET1|SET2, FALSE);
	 }
      }

      /*autobondrot at this stage has m->method == INTERSECTONCE == 1 */
      /* or m->method == SELFINTERSECT == 3 (seems as given in probe command)*/
      /* and m->countDots==1, m->rawOutput==1, m->drawSpike==1 */
      /* orientationName==rawname== char* holding angle values */

      doCommand(m->outf, m->method,
	 m->allMainAtoms, m->abins,
	 allMovingAtoms, bbins,
	 m->dots, m->probeRad, m->density, m->spikelen,
	 m->countDots, m->rawOutput, FALSE, orientationName, scoreBias,
	 m->drawSpike, m->sayGroup, m->groupLabel,
	 m->argc, m->argv, m->message); // SJ -10/07/2011 by default no condensed output for this

         /*allMovingAtoms is closest thing to a flag for autobondrot mode*/

      disposeBins(bbins);
   }
}/*movingDoCommand() for autobondrot*/
/*}}}movingDoCommand()_______________________________________________________*/

/*{{{doCommand()******* called from mainProbeProc() ***calls writeOutput()****/
void doCommand(FILE *outf, int method,
   atom *allMainAtoms, atomBins *abins,
   atom *allMovingAtoms, atomBins *bbins,
   pointSet dots[], float probeRad, float density, float spikelen,
   int countDots, int rawOutput, int conFlag, char* rawname, double scoreBias,
   int drawSpike, int sayGroup, char* groupLabel,
   int argc, char **argv, char message[])
{

/* NOTE: allMainAtoms/abins & allMovingAtoms/bbins must be disjoint */
/*       sets of atoms (none in common) or allMovingAtoms/bbins can */
/*       be NULL.                                                   */
/*       allMovingAtoms refers to autobondrot set of atoms          */

   /*doCommand is called from mainProbeProc() with rawname=="" */
   /* and called from movingDoCommand() with rawname==orientationName */
   /* which holds the string of the current autobondrot cycle's angle values*/
   /*for autobondrot: method is set by probe command, e.g. */
   /* INTERSECTONCE==1 as used in example of mobile sidechain in a protein */
   /* SELFINTERSECT==3 as used by alaphitaupsi ramachandran calculation*/
   /* autobondrot option defines: countDots==1, rawOutput==1*/
   /* autobondrot seems to get here with: drawSpike==1*/

   dotNode *rslts[NUMATOMTYPES][NODEWIDTH];
   int nsel = 0, numSkinDots = 0;
   int usesMovingAtoms = FALSE;

   int  j=0;

   char extrastr[32]; /*060129 for extra master to control orig vs fitted dots*/

   /*allMovingAtoms is closest thing to a flag for autobondrot mode*/
   /*so if it exists and bbins exist (dcr?: mobile atoms near enough?) */
   /*then usesMovingAtoms becomes the flag for autobondrot mode */
   /*in effect probe is modal but too object oriented to admit such globalness*/

   usesMovingAtoms = ((allMovingAtoms != NULL) && (bbins != NULL));
   if (! usesMovingAtoms)
   {
	allMovingAtoms = NULL;
	bbins = NULL;
   }

   initResults(rslts);
   if (!countDots && rawOutput && Verbose)
   { /*write rawOutput col headers*/
        if (OldStyleU)
        {
        	if(conFlag)
        		note(">>name:pat:type:srcAtom:targAtom:dot-count:min-gap:gap:"
        		"kissEdge2BullsEye:dot2BE:dot2SC:spike:score:"
        		"stype:ttype:x:y:z:sBval:tBval");
        	else
        		note(">>name:pat:type:srcAtom:targAtom:min-gap:gap:"
        		"kissEdge2BullsEye:dot2BE:dot2SC:spike:score:"
        		"stype:ttype:x:y:z:sBval:tBval");
       }
       else
       {
       	       if(conFlag)
       	      		 note(">>name:pat:type:srcAtom:targAtom:dot-count:min-gap:gap:"
       	       	       "spX:spY:spZ:spikeLen:score:"
              	       "stype:ttype:x:y:z:sBval:tBval");
              else
              	   note(">>name:pat:type:srcAtom:targAtom:min-gap:gap:"
       	       	       "spX:spY:spZ:spikeLen:score:"
              	       "stype:ttype:x:y:z:sBval:tBval");
       }
   }/*write rawOutput col headers*/
   if (method == SELFINTERSECT)
   {/*{{{(method == SELFINTERSECT)___expected default method___****/
   	   if (Verbose && !(countDots && rawOutput))
   		note("SelfIntersect");

    /*SELFINTERSECT case, using one input file (and where srcPat == targPat)*/
    /*is where unwashed NMR files with multiple models can be processed */
    /*when there is not a single model specified in the pattern.*/
    /*This uses a loop through all the model numbers actually in the file.*/

    	  if(modelSrc == 0 && modelTarg == 0 && modelCount > 1) /*041114*/ /*multiple models in file but no model specified in pattern*/
    	        modelLimit = modelCount;
          else
   		modelLimit = 1;

   	  for(j = 1; j <= modelLimit; j++) /*041114*/
   	  {/*loop over models*/

   	  	  if(modelLimit > 1)
   	  	  {/*defines the multiple model case*/
   	  	  	  modelToProcess = modelNumber[j]; /*models can be in any order 041114*/
   	  	  	  if(Verbose)
   	  	  	    	  fprintf(stderr,"processing modelNumber== %d\n",modelNumber[j]);
   	  	  }
   	  	  if(j > 1)
   	  	  {
   	  	  	  freeResults(rslts);
   	  	  	  initResults(rslts);
   	  	  }
   	  	  genDotIntersect(allMainAtoms, abins, allMovingAtoms, bbins,dots, probeRad, spikelen,SET1, SET1, rslts);
   	  	  /*does this for all atoms... calls examineDots()*/

          
   	  	  if (countDots)
   	  	  {/*autobondrot sets this, commandline can set this*/
   	  	  	  if (!rawOutput) {descrCommand(outf, "program:", "command:", argc, argv);}
   	  	  	  numSkinDots = enumDotSkin(allMainAtoms, abins, allMovingAtoms, bbins, dots, SET1);
   	  	  	  /*numSkinDots used to normalize output score*/
   	  	  	  if (!rawOutput)
   	  	  	  {
   	  	  	  	  fprintf(outf, "selection: self\nname: %s\n", groupLabel?groupLabel:"dots");
   	  	  	  	  fprintf(outf, "density: %.1f dots per A^2\nprobeRad: %.3f A\nVDWrad: (r * %.3f) + %.3f A\n",
   	  	  	  	  	  density, probeRad, RadScaleFactor, RadScaleOffset);
   	  	  	  	  fprintf(outf, "score weights: gapWt=%g, bumpWt=%g, HBWt=%g\n",
   	  	  	  	  	  GAPweight, BUMPweight, HBweight);
   	  	  	  }
   	  	  	  if (usesMovingAtoms)
   	  	  	    	  nsel = countSelected(allMainAtoms,   SET1) + countSelected(allMovingAtoms, SET1);
   	  	  	  else
   	  	  	    	  nsel = countSelected(allMainAtoms, SET1);

   	  	  	  if (rawOutput)
   	  	  	  {/*autobondrot sets this*/
   	  	  	  	  rawEnumerate(outf, "", rslts, method,
   	  	  	  	  	  nsel, drawSpike, FALSE, numSkinDots, density,
   	  	  	  	  	  groupLabel?groupLabel:"", rawname, scoreBias);
   	  	  	  }
   	  	  	  else
   	  	  	  {
   	  	  	  	  enumerate(outf, "self dots", rslts, probeRad, method, nsel,
   	  	  	  	  	  drawSpike, FALSE, numSkinDots, density);
   	  	  	  }
   	  	  }/*countDots*/
   	  	  else
   	  	  {
   	  	  	  if (rawOutput)
				{
   	  	  	    	  writeRaw(outf, "1->1", rslts, probeRad,groupLabel?groupLabel:"", density,conFlag);
				}
   	  	  	  else if (OutputFmtType == 1)
   	  	  	   	  writeAltFmtO(outf, TRUE, TRUE, "self_dots", rslts, drawSpike);
   	  	  	  else if (OutputFmtType == 2)
   	  	  	  {
   	  	  	  	  descrCommand(outf, "# software:", "# command:", argc, argv);
   	  	  	  	  writeAltFmtXV(outf, TRUE, TRUE, "self_dots", rslts, drawSpike);
   	  	  	  }
   	  	  	  else if (OutputFmtType == 3)
   	  	  	  { /*dcr041101 ONELINE :count:summaries:*/
   	  	  	  	  countsummary(outf,"SelfIntersect", 1, 1); /*dcr041101 Lines,Pass*/
   	  	  	  }
   	  	  	  else
   	  	  	  {/*kinemage output*/
   	  	  	  	  if(ContactSUMMARY) /*dcr041101 Lines,Pass*/
   	  	  	  	    	countsummary(outf,"SelfIntersect", 9, 1);

   	  	  	  	  descrCommand(outf, "@caption", " command:", argc, argv);

   	  	  	  	  if (sayGroup)
   	  	  	  	  {
   	  	  	  	  	  if(modelLimit > 1)
   	  	  	  	  	  {/*doing jth of multiple models of an ensemble*/
   	  	  	  	  	  	  fprintf(outf, "@group dominant {%s M%d} animate\n",groupLabel?groupLabel:"dots",j);
   	  	  	  	  	  }
   	  	  	  	  	  else
   	  	  	  	  	    	  fprintf(outf, "@group dominant {%s}\n",groupLabel?groupLabel:"dots");
   	  	  	  	  }
   	  	  	  	  sprintf(extrastr,"%s",groupLabel?groupLabel:"dots"); /*060129*/
                  writeOutput(outf, "self dots", rslts, drawSpike, method, extrastr, probeRad);
                           /*add probeRad 20111220dcr*/
   	  	  	  }/*kinemage output*/
   	  	  }
   	  }/*loop over models*/
   }/*}}}(method == SELFINTERSECT)________________________________*/

   else if (method == INTERSECTONCE)
   {/*{{{(method == INTERSECTONCE)*********************************/
      if (Verbose && !(countDots && rawOutput))
      	      note("IntersectOnce");

      genDotIntersect(allMainAtoms, abins, allMovingAtoms, bbins,dots, probeRad, spikelen, SET1, SET2, rslts);

      if (countDots)
      {/*autobondrot sets this*/
      	      if (!rawOutput)
      	      	      descrCommand(outf, "program:", "command:", argc, argv);

      	      numSkinDots = enumDotSkin(allMainAtoms, abins, allMovingAtoms, bbins, dots, SET1);
         /*numSkinDots used to normalize output score*/
              if (!rawOutput)
	      {
	 	     fprintf(outf, "selection: once\nname: %s\n", groupLabel?groupLabel:"dots");
	 	     fprintf(outf, "density: %.1f dots per A^2\nprobeRad: %.3f A\nVDWrad: (r * %.3f) + %.3f A\n",
	 	     	     density, probeRad, RadScaleFactor, RadScaleOffset);
		     fprintf(outf, "score weights: gapWt=%g, bumpWt=%g, HBWt=%g\n",
		     	     GAPweight, BUMPweight, HBweight);
	      }
	      if (usesMovingAtoms)
	      	      nsel = countSelected(allMainAtoms,   SET1) + countSelected(allMovingAtoms, SET1);
              else
              	      nsel = countSelected(allMainAtoms, SET1);

              if (rawOutput)
              {/*autobondrot sets this*/
              	      rawEnumerate(outf,"", rslts, method, nsel,
              	      	      drawSpike, FALSE, numSkinDots, density,
              	      	      groupLabel?groupLabel:"", rawname, scoreBias);
              }
              else
              {
              	      enumerate(outf,"once dots", rslts, probeRad, method, nsel,
              	      	      drawSpike, FALSE, numSkinDots, density);
              }
      }
      else
      {
      	      if (rawOutput)
		{
      	      	    writeRaw(outf, "1->2", rslts, probeRad,groupLabel?groupLabel:"", density,conFlag);

		}
		else if (OutputFmtType == 1)
	      	      writeAltFmtO(outf, TRUE, TRUE, "once_dots", rslts, drawSpike);
	      else if (OutputFmtType == 2)
	      {
	      	      descrCommand(outf, "# software:", "# command:", argc, argv);
	      	      writeAltFmtXV(outf, TRUE, TRUE, "once_dots", rslts, drawSpike);
	      }
	      else if (OutputFmtType == 3)
	      { /*dcr041101 ONELINE :count:summaries:*/
	      	      countsummary(outf,"IntersectOnce", 1, 1); /*dcr041101 Lines,Pass*/
	      }
	      else
	      {/*write kinemage*/
	      	      if(ContactSUMMARY) /*dcr041101 Lines,Pass*/
               			countsummary(outf,"IntersectOnce", 9, 1);

               	      descrCommand(outf, "@caption", " command:", argc, argv);
               	      if (sayGroup)
               	      	      fprintf(outf, "@group dominant {%s}\n", groupLabel?groupLabel:"dots");

               	      sprintf(extrastr,"%s",groupLabel?groupLabel:"dots"); /*060129*/
                         writeOutput(outf, "once dots", rslts, drawSpike, method, extrastr, probeRad);
                         /*add probeRad 20111220dcr*/
              }
      }
   }/*}}}(method == INTERSECTONCE)________________________________*/

   else if (method == INTERSECTBOTHWAYS)
   {/*{{{(method == INTERSECTBOTHWAYS)*****************************/
      if (Verbose && !(countDots && rawOutput))
      		note("IntersectBothWays");

      if (countDots)
      {
	 if (!rawOutput)
	 {
	    descrCommand(outf, "program:", "command:", argc, argv);
	    fprintf(outf, "selection: both\nname: %s\n", groupLabel?groupLabel:"dots");
	    fprintf(outf, "density: %.1f dots per A^2\nprobeRad: %.3f A\nVDWrad: (r * %.3f) + %.3f A\n",
				 density, probeRad, RadScaleFactor, RadScaleOffset);
	    fprintf(outf, "score weights: gapWt=%g, bumpWt=%g, HBWt=%g\n",
				 GAPweight, BUMPweight, HBweight);
	 }
      }
      else
      {
      	      if (rawOutput) {}
		   /* do nothing on purpose */
	      else if (OutputFmtType == 1) {}
	      else if (OutputFmtType == 2)
	      	      descrCommand(outf, "# software:", "# command:", argc, argv);
	      else
	      {/*kinemage: keywords before double pass*/
	      	      descrCommand(outf, "@caption", " command:", argc, argv);
	      	      if (sayGroup)
	      	      		fprintf(outf, "@group {%s}\n",groupLabel?groupLabel:"dots");
	    }
      }

      genDotIntersect(allMainAtoms, abins, allMovingAtoms, bbins,
		  dots, probeRad, spikelen, SET1, SET2, rslts);

      if (countDots)
      {
	 numSkinDots = enumDotSkin(allMainAtoms, abins, allMovingAtoms, bbins, dots, SET1);
         /*numSkinDots used to normalize output score*/
	 if (usesMovingAtoms)
	     nsel = countSelected(allMainAtoms,   SET1) + countSelected(allMovingAtoms, SET1);
	 else
             nsel = countSelected(allMainAtoms, SET1);

	 if (rawOutput)
	 {
	    rawEnumerate(outf, "1->2", rslts, method, nsel,
			drawSpike, FALSE, numSkinDots, density,
			groupLabel?groupLabel:"", rawname, scoreBias);
	 }
	 else
	 {
	    enumerate(outf, "1->2", rslts, probeRad, method, nsel,
			drawSpike, FALSE, numSkinDots, density);
	 }
      }
      else
      {
	 if (rawOutput)
	{
	    writeRaw(outf, "1->2", rslts, probeRad, groupLabel?groupLabel:"", density,conFlag);

	}
	else if (OutputFmtType == 1)
	    writeAltFmtO(outf, TRUE, !sayGroup, "1->2", rslts, drawSpike);
	 else if (OutputFmtType == 2)
	    writeAltFmtXV(outf, TRUE, !sayGroup, "1->2", rslts, drawSpike);
	 else if (OutputFmtType == 3)
	 { /*dcr041101 ONELINE :count:summaries:*/
            countsummary(outf,"IntersectBothWays 1->2", 1, 0);
            /*dcr041101 Lines,Pass==0 no output for this pass*/
	 }
	 else
	 {
            sprintf(extrastr,"%s",groupLabel?groupLabel:"dots"); /*060129*/
        writeOutput(outf, "1->2", rslts, drawSpike, method, extrastr,probeRad);
        /*add probeRad 20111220dcr*/

            if(ContactSUMMARY) /*dcr041101 Lines,Pass*/
             	countsummary(outf,"IntersectBothWays 1->2", 9, 0);
               /*dcr041101 Lines,Pass==0 no output for this pass*/
	 }
      }
      freeResults(rslts);
      initResults(rslts);

      genDotIntersect(allMainAtoms, abins, allMovingAtoms, bbins,
			dots, probeRad, spikelen, SET2, SET1, rslts);

      if (countDots)
      {
	 numSkinDots = enumDotSkin(allMainAtoms, abins, allMovingAtoms, bbins, dots, SET2);
         /*numSkinDots used to normalize output score*/
	 if (usesMovingAtoms)
	 	 nsel = countSelected(allMainAtoms,   SET2) + countSelected(allMovingAtoms, SET2);
	 else
	 	 nsel = countSelected(allMainAtoms, SET2);

	 if (rawOutput)
	 {
	    rawEnumerate(outf, "2->1", rslts, method, nsel,
			drawSpike, FALSE, numSkinDots, density,
			groupLabel?groupLabel:"", rawname, scoreBias);
	 }
	 else
	 {
	    enumerate(outf, "2->1", rslts, probeRad, method, nsel,
			drawSpike, FALSE, numSkinDots, density);
	 }
      }
      else
      {
	 if (rawOutput)
	 {
         writeRaw(outf, "2->1", rslts, probeRad,
			groupLabel?groupLabel:"", density,conFlag);

	 }
	 else if (OutputFmtType == 1)
	 {
	    writeAltFmtO(outf, !sayGroup, TRUE, "2->1", rslts, drawSpike);
	 }
	 else if (OutputFmtType == 2)
	 {
	    writeAltFmtXV(outf, !sayGroup, TRUE, "2->1", rslts, drawSpike);
	 }
	 else if (OutputFmtType == 3)
	 { /*dcr041101 ONELINE :count:summaries:*/
            countsummary(outf,"IntersectBothWays 2->1", 1, 2);
            /*dcr041101 Lines,Pass==2 output sum of 2 passes*/
	 }
	 else
	 {
            sprintf(extrastr,"%s",groupLabel?groupLabel:"dots"); /*060129*/
        writeOutput(outf, "2->1", rslts, drawSpike, method, extrastr, probeRad);
        /*add probeRad 20111220dcr*/

            if(ContactSUMMARY) /*dcr041101 Lines,Pass*/
               countsummary(outf,"IntersectBothWays 2->1", 9, 2);
               /*dcr041101 Lines,Pass==2 output sum of 2 passes*/
	 }
      }
   }/*}}}(method == INTERSECTBOTHWAYS)____________________________*/

   else if (method == EXTERNALSURFACE)
   {/*{{{(method == EXTERNALSURFACE)*******************************/
      if (Verbose && !(countDots && rawOutput))
      	      note("ExternalSurface");

      genDotSurface(allMainAtoms, abins, allMovingAtoms, bbins, dots, probeRad, spikelen, SET1, rslts);

      if (countDots)
      {
	 if (!rawOutput)
	    descrCommand(outf, "program:", "command:", argc, argv);

	 numSkinDots = enumDotSkin(allMainAtoms, abins, allMovingAtoms, bbins, dots, SET1);
         /*numSkinDots used to normalize output score*/
	 if (!rawOutput)
	 {
	    fprintf(outf, "selection: external\nname: %s\n",groupLabel?groupLabel:"dots");
	    fprintf(outf, "density: %.1f dots per A^2\nprobeRad: %.3f A\nVDWrad: (r * %.3f) + %.3f A\n",
		     density, probeRad, RadScaleFactor, RadScaleOffset);
	 }
	 if (usesMovingAtoms)
	 	 nsel = countSelected(allMainAtoms,   SET1) + countSelected(allMovingAtoms, SET1);
	 else
	 	 nsel = countSelected(allMainAtoms, SET1);

	 if (rawOutput)
	 {
	    rawEnumerate(outf, "", rslts, method, nsel,
		     FALSE, TRUE, numSkinDots, density,
		     groupLabel?groupLabel:"", rawname, scoreBias);
	 }
	 else
	 {
	    enumerate(outf, "extern dots", rslts, probeRad, method, nsel,
		     FALSE, TRUE, numSkinDots, density);
	 }
      }
      else
      {
	 if (rawOutput)
	 {
	    writeRaw(outf, "1->none", rslts, probeRad,
		  groupLabel?groupLabel:"", density,conFlag);

	 }
	 else if (OutputFmtType == 1)
	 {
	    writeAltFmtO(outf, TRUE, TRUE, "extern_dots", rslts, FALSE);
	 }
	 else if (OutputFmtType == 2)
	 {
	    descrCommand(outf, "# software:", "# command:", argc, argv);
	    writeAltFmtXV(outf, TRUE, TRUE, "extern_dots", rslts, FALSE);
	 }
	 else
	 {
	    descrCommand(outf, "@caption", " command:", argc, argv);
	    if (sayGroup)
	       	fprintf(outf, "@group dominant {%s}\n", groupLabel?groupLabel:"dots");

            sprintf(extrastr,"%s",groupLabel?groupLabel:"dots"); /*060129*/
        writeOutput(outf, "extern dots", rslts, FALSE, method, extrastr, probeRad);
        /*add probeRad 20111220dcr*/
	 }
      }
   }/*}}}(method == EXTERNALSURFACE)______________________________*/

   else if (method == DUMPATOMCOUNT)
   {/*{{{(method == DUMPATOMCOUNT)*********************************/
      if (Verbose && !rawOutput)
      {
	 note("dumpAtomInfo");
	 descrCommand(outf, "program:", "command:", argc, argv);
	 fprintf(outf, "selection: self\nname: %s\n", groupLabel?groupLabel:"dots");
	 fprintf(outf, "density: %.1f dots per A^2\nprobeRad: %.3f A\nVDWrad: (r * %.3f) + %.3f A\n",
				 density, probeRad, RadScaleFactor, RadScaleOffset);
	 fprintf(outf, "score weights: gapWt=%g, bumpWt=%g, HBWt=%g\n",
				 GAPweight, BUMPweight, HBweight);
      }
      if (usesMovingAtoms)
	 nsel = countSelected(allMainAtoms,   SET1) + countSelected(allMovingAtoms, SET1);
      else
	 nsel = countSelected(allMainAtoms, SET1);

      if (rawOutput)
      {
	 if (groupLabel)
	 {
	    fprintf(outf, "%d %s %s%s\n", nsel, rawname,
	       RAW_HEADER_COMMENT, groupLabel);
	 }
	 else
	 {
	    fprintf(outf, "%d %s\n", nsel, rawname);
	 }
      }
      else
      {
	 fprintf(outf, "atoms selected: %d\n", nsel);
      }
   }/*}}}(method == DUMPATOMCOUNT)_______________________________*/

   freeResults(rslts);

}/*doCommand()*/
/*}}}doCommand()_____________________________________________________________*/

/*{{{descrCommand()***********************************************************/
void descrCommand(FILE *outf, char *hdr1, char *hdr2, int argc, char **argv)
{
   int i = 0;
   time_t t;
   char *ts = NULL;

   t = time(NULL);
   ts = asctime(localtime(&t)); /* system dependent time string */

   fprintf(outf, "%s %s, run %s", hdr1, shortVersionStr, ts);
   fprintf(outf, "%s %s", hdr2, argv[0]);
   for (i = 1; i < argc; i++) {
      int dontQuote = (strpbrk(argv[i], " \t\n\b\r\v\f\a*\'") == NULL);
      fprintf(outf,(dontQuote?" %s":" \"%s\""), argv[i]);
   }
   fprintf(outf, "\n");
}
/*}}}descrCommand()__________________________________________________________*/

/*{{{loadDotSpheres()*********************************************************/
void loadDotSpheres(pointSet dots[], float density)
{
   int i = 0;

   for (i = 0; i < NUMATOMTYPES; i++) {
      dotSphere(&(dots[i]), getRadius(i, 0), density);
   }
   dotSphere(&COdots, getRadius(atomC, 1), density);
}
/*}}}loadDotSpheres()________________________________________________________*/

/*{{{unloadDotSpheres()*******************************************************/
void unloadDotSpheres(pointSet dots[])
{
   int i = 0;

   for (i = 0; i < NUMATOMTYPES; i++) {
      freeDotSphere(&(dots[i]));
   }
}
/*}}}unloadDotSpheres()______________________________________________________*/

/*{{{probehelp()**************************************************************/
void probehelp(int longlist)
{
   fprintf(stderr, "\nSyntax: probe input.pdb >> out.kin\n");
   fprintf(stderr, "    or: probe [flags] \"src pattern\" [\"target pattern\"] pdbfiles... >> out.kin\n\n");
   fprintf(stderr, "Flags:\n");
   fprintf(stderr, "  -SElf  self intersection:   src  -> src (default)\n");
   fprintf(stderr, "  -Both  intersect both ways: src <=> targ\n");
   fprintf(stderr, "  -ONce  single intersection: src  -> targ\n");
   fprintf(stderr, "  -OUt   external van der Waals surface of src (solvent contact surface)\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "  -AUTObondrot filename    read and process an autobondrot file\n");
if (! longlist) {
   fprintf(stderr, "\n  shortcuts: <<NO FLAGS>>, -SCSurface, -EXPOsed, -ASurface, -ACCESS, -SCAN0, -SCAN1\n");
   fprintf(stderr, "\n most simple dotkin: probe -self all -kinemage input.pdb > out.kin\n");
}
else { /*longlist option*/
   fprintf(stderr, "\n  shortcuts:\n");
   fprintf(stderr, " <<NO FLAGS>>same as: -4H -mc -het -self \"altA ogt33\"\n");
                            /*changed from 3 : dcr041017*/
   fprintf(stderr, "  -DEFAULTs  same as: <<NO FLAGS>>, but allows some other flags\n"); /*dcr041101*/
   fprintf(stderr, "  -SCSurface same as: -drop -rad1.4         -out \"not water\"\n");
   fprintf(stderr, "  -EXPOsed   same as: -drop -rad1.4         -out      (note: user supplies pattern)\n");
   fprintf(stderr, "  -ASurface  same as: -drop -rad0.0 -add1.4 -out \"not water\"\n");
   fprintf(stderr, "  -ACCESS    same as: -drop -rad0.0 -add1.4 -out      (note: user supplies pattern)\n");
   fprintf(stderr, "  -SCAN0 same as: -4H -mc -self \"alta blt40 ogt33\"\n");
                        /*changed from 3 : dcr041017*/
   fprintf(stderr, "  -SCAN1 same as: -4H -once \"sc alta blt40 ogt33\" \"alta blt40 ogt65,(not water ogt33)\"\n");
                        /*changed from 3 : dcr041017*/
   fprintf(stderr, "\n");
   fprintf(stderr, "  -DUMPAtominfo   count the atoms in the selection: src\n\n");
   fprintf(stderr, "  (note that BOTH and ONCE require two patterns while\n");
   fprintf(stderr, "   OUT, SELF and DUMPATOMINFO require just one pattern)\n\n");
   fprintf(stderr, "  -Implicit    implicit hydrogens\n");
   fprintf(stderr, "  -Explicit    explicit hydrogens (default)\n");
   fprintf(stderr, "  -DEnsity#    set dot density (default 16 dots/sq A)\n");
   fprintf(stderr, "  -Radius#.#   set probe radius (default 0.25 A)\n");
   fprintf(stderr, "  -ADDvdw#.#   offset added to Van der Waals radii (default 0.0)\n");
   fprintf(stderr, "  -SCALEvdw#.# scale factor for Van der Waals radii (default 1.0)\n");
   fprintf(stderr, "  -COSCale#.#  scale C=O carbon Van der Waals radii (default 0.94)\n");
   fprintf(stderr, "  -SPike       draw spike instead of dots (default)\n");
   fprintf(stderr, "  -SPike#.#    set spike scale (default=0.5)\n");
   fprintf(stderr, "  -NOSpike     draw only dots\n");
   fprintf(stderr, "  -HBRegular#.# max overlap for regular Hbonds(default=%.1f)\n", Min_regular_hb_cutoff);
   fprintf(stderr, "  -HBCharged#.# max overlap for charged Hbonds(default=%.1f)\n", Min_charged_hb_cutoff);
   fprintf(stderr, "  -Keep        keep nonselected atoms (default)\n");
   fprintf(stderr, "  -DRop        drop nonselected atoms\n");
   fprintf(stderr, "  -LIMit       limit bump dots to max dist when kissing (default)\n");
   fprintf(stderr, "  -NOLIMit     do not limit bump dots\n");
   fprintf(stderr, "  -LENs        add lens keyword to kin file\n");
   fprintf(stderr, "  -NOLENs      do not add lens keyword to kin file (default)\n");
   fprintf(stderr, "  -MC          include mainchain->mainchain interactions\n");
   fprintf(stderr, "  -HETs        include dots to non-water HET groups (default)\n");
   fprintf(stderr, "  -NOHETs      exclude dots to non-water HET groups\n");
   fprintf(stderr, "  -WATers      include dots to water (default)\n");
   fprintf(stderr, "  -NOWATers    exclude dots to water\n");
   fprintf(stderr, "  -WAT2wat     show dots between waters\n");
   fprintf(stderr, "  -DUMPH2O     include water H? vectorlist in output\n");
   fprintf(stderr, "  -4H          extend bond chain dot removal to 4 for H (default)\n");
   fprintf(stderr, "  -3           limit bond chain dot removal to 3\n");
   fprintf(stderr, "  -2           limit bond chain dot removal to 2\n");
   fprintf(stderr, "  -1           limit bond chain dot removal to 1\n");
   fprintf(stderr, "  -IGNORE \"pattern\" explicit drop: ignore atoms selected by pattern\n");
   fprintf(stderr, "  -DOCHO       recognize CH..O Hbonds\n");
   fprintf(stderr, "  -CHO#.#      scale factor for CH..O Hbond score (default=%.1f)\n", CHOHBfactor);
   fprintf(stderr, "  -PolarH      use short radii of polar hydrogens (default)\n");
   fprintf(stderr, "  -NOPolarH    do not shorten radii of polar hydrogens\n");
   fprintf(stderr, "  -NUClear     use nuclear position vdW radii\n"); // 130114 JJH
   fprintf(stderr, "               (default is electron cloud positions)\n"); // 130114 JJH
   fprintf(stderr, "  -NOFACEhbond do not identify HBonds to aromatic faces\n\n");

   fprintf(stderr, "  -Name \"name\" specify the group name (default \"dots\")\n");
   fprintf(stderr, "  -DOTMASTER  group name used as extra master={name} on lists\n"); /*060129*/
   fprintf(stderr, "  -NOGroup     do not generate @group statement in .kin format output\n");
   fprintf(stderr, "  -KINemage    add @kinemage 1 statement to top of .kin format output\n\n");

   fprintf(stderr, "  -Countdots   produce a count of dots-not a dotlist\n");
   fprintf(stderr, "  -Unformated  output raw dot info\n");
   fprintf(stderr, "     name:pat:type:srcAtom:targAtom:mingap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval:\n");
   fprintf(stderr, "  -CONdense    raw output in condensed format, i.e. one line per source->target atom. Also gives the dot count for that interaction. Works only with -Unformated flag\n");
   fprintf(stderr, "     name:pat:type:srcAtom:targAtom:dot-count:mingap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval:\n");
   fprintf(stderr, "  -OFORMAT     output dot info formatted for display in O\n");
   fprintf(stderr, "  -XVFORMAT    output dot info formatted for display in XtalView\n");
   fprintf(stderr, "  -ONELINE     output one line :contacts:by:severity:type:\n"); /*dcr041101*/
   fprintf(stderr, "  -ONEDOTeach  output one dot for each src-to-neighbor\n"); /*dcr111013*/
   fprintf(stderr, "               With -Unformated flag: for H gives angle parent-src-cause,\n");
   fprintf(stderr, "     name:pat:type:srcAtom:targAtom:XHTangle:mingap:gap:spX:spY:spZ:spikeLen:score:stype:ttype:x:y:z:sBval:tBval:\n");
   fprintf(stderr, "  -GAPBINs     2nd pointmaster character, -.5...0...+.5\n"); /*dcr111020*/
   fprintf(stderr, "  -NEAREST     apply selection only on nearest neighbor\n"); /*dcr111020*/
   fprintf(stderr, "  -FILTERnearest apply 3rd selection on neighbor\n"); /*dcr111020*/
   fprintf(stderr, " e.g. -nearest -filternearest patt1 all patt3\n"); /*dcr111020*/
   fprintf(stderr, "  -OCC1filter  final filter exclude any w occ < 1\n");  /*dcr111020*/
   fprintf(stderr, "  -NOALTfilter final filter exclude any alts\n");  /*dcr111020*/
   fprintf(stderr, "  -GAPcolor    color dots by gap amount (default)\n");
   fprintf(stderr, "  -ATOMcolor   color dots by atom type\n");
   fprintf(stderr, "  -BASEcolor   color dots by nucleic acid base type\n");
   fprintf(stderr, "  -COLORBase   color dots by gap and nucleic acid base type\n");
   fprintf(stderr, "  -OUTCOLor \"name\" specify the point color for -OUT (default \"gray\")\n");
   fprintf(stderr, "  -GAPWeight#  set weight for scoring gaps (default 0.25)\n");
   fprintf(stderr, "  -BUMPWeight# set relative scale for scoring bumps (default 10.0)\n");
   fprintf(stderr, "  -HBWeight#   set relative scale for scoring Hbonds (default 4.0)\n");
   fprintf(stderr, "  -DIVLow#.#   Division for Bump categories    (default -0.4)\n");
   fprintf(stderr, "  -DIVHigh#.#  Division for Contact categories (default  0.25)\n");
   fprintf(stderr, "  -DIVWorse#.#  Division for regarding a clash as a worse clash (default  -0.5)\n"); /* 04/21/2015 SJ*/
   fprintf(stderr, "  -SEPWORSE    To separate bad overlaps and worse overlaps (default false, if true, deafult value to separate the worse clashes is -0.5)\n");
   fprintf(stderr, "  -MINOCCupancy#.#  Occupancy below this is same as zero (default  0.02)\n");
   fprintf(stderr, "  -ELEMent     add master buttons for different elements in kin output\n");
   fprintf(stderr, "  -NOHBOUT     do not output contacts for HBonds\n");
   fprintf(stderr, "  -NOCLASHOUT  do not output contacts for clashes\n");
   fprintf(stderr, "  -NOVDWOUT    do not output contacts for van der Waals interactions\n");
   fprintf(stderr, "  -ONLYBADOUT  onlybadout output bad clashes (severe overlap contacts)\n"); /*dcr041010*/
   fprintf(stderr, "  -SUMMARY     output summary list of contacts and clashes\n"); /*dcr041101*/
   fprintf(stderr, "  -ONELINE     output summary list on oneline\n"); /*dcr041101*/
   fprintf(stderr, "  -NOTICKs     do not display the residue name ticker during processing\n");
   fprintf(stderr, "  -STDBONDs    assume only standard bonding patterns in standard residues\n");
   fprintf(stderr, "  -NOPARENT    do not bond hydrogens based on table of parent heavy atoms\n\n");
   fprintf(stderr, "  -SEGID       use the PDB SegID field to descriminate between residues\n");
   fprintf(stderr, "  -OLDU        generate old style -u output: kissEdge2BullsEye, etc\n");
   fprintf(stderr, "  -VErbose     verbose mode (default)\n");
   fprintf(stderr, "  -DOTDUMP     dump dot info while doing examineDots()\n");
   fprintf(stderr, "  -REFerence   display reference string\n");
   fprintf(stderr, "  -CHANGEs     display a list of program changes\n");
   fprintf(stderr, "  -Quiet       quiet mode\n");
}/*longlist option*/
   fprintf(stderr, "\n  -Help  show expanded help notice (includes other flags)\n");
   fprintf(stderr, "\n  -VERSION   one line version to stdout\n\n");/*dcr041009*/
   fprintf(stderr, "Pattern elements:  (should be put in quotes on the command line)\n");
if (longlist) { /*longlist option continues*/
   fprintf(stderr, "   FILE#     within file #\n");
   fprintf(stderr, "   MODEL#    within model #\n");
   fprintf(stderr, "   CHAINaa   within chain a\n");
   fprintf(stderr, "   SEGaaaa   segment identifier aaaa (where _ represents blank)\n");
   fprintf(stderr, "   ALTa      alternate conformation a\n");
   fprintf(stderr, "   ATOMaaaa  atom name aaaa (where _ represents blank)\n");
   fprintf(stderr, "             (all 4 characters are used so H would be ATOM_H__)\n");
   fprintf(stderr, "   RESaaa    residue aaa\n");
   fprintf(stderr, "   #         residue #\n");
   fprintf(stderr, "   #a        residue #, insert a\n");
   fprintf(stderr, "   #-#       residue range # (insert codes ignored)\n");
   fprintf(stderr, "   a         residue type by one letter codes   (eg. y)\n");
   fprintf(stderr, "   aaa       residue type by three letter codes (eg. tyr)\n");
   fprintf(stderr, "   ALL,PROTEIN,MC,SC,BASE,ALPHA,BETA,NITROGEN,CARBON,OXYGEN,\n");
   fprintf(stderr, "   SULFUR,PHOSPHORUS,HYDROGEN,METAL,POLAR,NONPOLAR,CHARGED,\n");
   fprintf(stderr, "   DONOR,ACCEPTOR,AROMATIC,METHYLENE,METHYL,HET,WATER,DNA,RNA\n");
   fprintf(stderr, "             all or a subset of the atoms\n");
   fprintf(stderr, "   OLT#      Occupancy less than # (integer percent)\n");
   fprintf(stderr, "   OGT#      Occupancy greater than # (integer percent)\n");
   fprintf(stderr, "   BLT#      B-value less than # (integer)\n");
   fprintf(stderr, "   BGT#      B-value greater than # (integer)\n");
   fprintf(stderr, "   INSa      Insert code a (where _ represents blank)\n");
   fprintf(stderr, "   \n");
   fprintf(stderr, "   WITHIN #.# OF #.#, #.#, #.#   atoms within distance from point\n");
}/*longlist option continues*/
else { /*NOT longlist*/
   fprintf(stderr, "   #         residue number\n");
   fprintf(stderr, "   #a        residue #, insert a\n");
   fprintf(stderr, "   #-#       residue number range\n");
   fprintf(stderr, "   a OR aaa  residue type by one or three letter codes\n");
   fprintf(stderr, "   FILE#,MODEL#,CHAINa,SEGaaaa,ALTa,ATOMaaaa,RESaaa,\n");
   fprintf(stderr, "   ALL,PROTEIN,MC,SC,BASE,ALPHA,BETA,NITROGEN,CARBON,OXYGEN,\n");
   fprintf(stderr, "   SULFUR,PHOSPHORUS,HYDROGEN,METAL,POLAR,NONPOLAR,CHARGED,\n");
   fprintf(stderr, "   DONOR,ACCEPTOR,AROMATIC,METHYLENE,METHYL,HET,WATER,DNA,RNA,\n");
   fprintf(stderr, "   OLT#, OGT#, BLT#, BGT#, INSa, WITHIN #.# OF #.#, #.#, #.#\n");
}
   fprintf(stderr, "   \n");
   fprintf(stderr, "   Patterns can be combined into comma separated lists\n");
   fprintf(stderr, "   such as \"trp,phe,tyr\" meaning TRP or PHE or TYR.\n");
   fprintf(stderr, "   \n");
   fprintf(stderr, "   Patterns that are sepatated by blanks must all be true\n");
   fprintf(stderr, "   such as \"chainb 1-5\" meaning residues 1 to 5 in chain B.\n");
   fprintf(stderr, "   \n");
   fprintf(stderr, "   You can also group patterns with parenthesis, separate\n");
   fprintf(stderr, "   multiple patterns with | meaning 'or' and choose the\n");
   fprintf(stderr, "   complement with NOT as in \"not file1\" meaning not in file 1.\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "   An autobondrot file is similar to other PDB input files\n");
   fprintf(stderr, "   but it includes information identifying atoms subject to rotations\n");
   fprintf(stderr, "   and other transformations.\n");
if (longlist) { /*more longlist*/
   fprintf(stderr, "\n   Example autobondrot file fragment showing Calpha-Cbeta bond rotation\n");
   fprintf(stderr, "   and a periodic torsion penalty function for this rotation\n");
   fprintf(stderr, "     ATOM      1  CB  TYR    61      34.219  17.937   4.659  1.00  0.00\n");
   fprintf(stderr, "     bondrot:chi1:78.7:  0:359:5:33.138:18.517: 5.531:34.219:17.937: 4.659\n");
   fprintf(stderr, "     cos:-3:60:3:\n");
   fprintf(stderr, "     ATOM      1 1HB  TYR    61      34.766  18.777   4.206  1.00  0.00\n");
   fprintf(stderr, "     ATOM      1 2HB  TYR    61      34.927  17.409   5.315  1.00  0.00\n");
   fprintf(stderr, "     ATOM      1  CG  TYR    61      33.836  16.989   3.546  1.00  0.00\n");
   fprintf(stderr, "     ...\n");

   fprintf(stderr, "   Autobondrot commands use colons to separate values\n");
   fprintf(stderr, "   Transformations: BONDROT:id:currAng:start:end:stepSz:x1:y1:z1:x2:y2:z2\n");
   fprintf(stderr, "                    TRANS:  id:currpos:start:end:stepSz:x1:y1:z1:x2:y2:z2\n");
   fprintf(stderr, "                    NULL  # dummy\n");
   fprintf(stderr, "   Bias functions:  COS:scale:phaseOffset:frequency\n");
   fprintf(stderr, "                    POLY:scale:offset:polynomialDegree\n");
   fprintf(stderr, "                    CONST:value\n");
   fprintf(stderr, "   Branching:       SAVE and RESTORE or \"(\" and \")\"\n");
   fprintf(stderr, "           (e.g. to rotate each Chi and the methyls for isoleucine the\n");
   fprintf(stderr, "            sequence is: rotChi1/SAVE/rotChi2/rotCD1/RESTORE/rotCG2)\n");
   fprintf(stderr, "   Set orientation: GO:angle1:angle2:...\n");
   fprintf(stderr, "   Include files:   @filename\n");
   fprintf(stderr, "   Comments:        # comment text\n");
}/*more longlist*/
   fprintf(stderr, "\n%s\n", versionString);
   exit(0);
}/*probehelp()*/
/*}}}probehelp()_____________________________________________________________*/

/*{{{probeversion()***********************************************************/
void probeversion(void) /*VERSION  dcr041009*/
{/*stdout one liner version for, e.g., MolProbity  dcr041009*/
   /*beware: orig defaults get superseded by shortcuts, like raw "probe"*/
   printf("%s\n", shortVersionStr);
   exit(0);
}
/*}}}probeversion()__________________________________________________________*/

/*{{{processCommandline()**** called from mainProbeProc() ********************/
atom* processCommandline(int argc, char **argv, int *method, region *bboxA,
			float *density, float *probeRad,
			int *drawSpike, float *spikelen, int *countDots,
			int *keepUnselected,
            char **srcArg, char **targArg, char **extraArg, char **ignoreArg,
            char **groupLabel, int *rawOutput, int *conFlag, int *sayGroup,
			int *addKinToFile, movingAtomBuildInfo *mabip,
			residue **reslstptr) //conFlag added SJ - 10/07/2011, extraArg DCR
{/*processCommandline()*/
   /*gets running conditions from the commandline as well as */
   /*loads atomlist from in file, returns atomlist which becomes allMainAtoms*/
   /* atomlist = probe.c/loadAtoms() */

   /*called with address of mabis, i.e. mabip moving atom build info ptr */

   int file = 0, nargs = 0, argcnt = 0, i = 0, n = 0;
   char *p, message[200];
   FILE *inf = stdin;
   atom *atomlist = NULL;

   *method = SELFINTERSECT;
   nargs = 1;

   *density = 16.0;
   *probeRad = 0.25;
   *drawSpike = TRUE;
   *spikelen = 0.50;
   *countDots = FALSE;
   *keepUnselected = TRUE;
   *rawOutput = FALSE;
   *conFlag = FALSE; // added by SJ - 10/07/2011
   *sayGroup = TRUE;
   *addKinToFile = FALSE;

    *srcArg = NULL;
   *targArg = NULL;
   *ignoreArg=NULL;
   *groupLabel = NULL;

   /*&mabis of mainProbeProc()== *mabip here */
   if (mabip) {
      mabip->filenum  = 0;
      mabip->inf      = NULL;
      mabip->close    = FALSE;
      mabip->srcPat   = NULL;
      mabip->targPat  = NULL;
      mabip->reslstptr= reslstptr;
      mabip->scratchRes= NULL;
   }

   file = 0;
   argcnt = 0;
   for (i = 1; i < argc; i++) {
      p = argv[i];
      if (p[0] == '-') { /* -flag item */
	 if ((argcnt >= nargs) && (p[1] == '\0')) {
	    file++;
	    atomlist = loadAtoms(stdin,atomlist,bboxA,file, reslstptr);
	 }
	 else if(compArgStr(p+1, "HELP", 1)){
	    probehelp(1);
	 }
	 else if(compArgStr(p+1, "VERSION", 7)){/*before VERBOSE,  dcr041009*/
	    probeversion();
	 }
	 else if(compArgStr(p+1, "VERBOSE", 2)){
	    Verbose = TRUE;
	 }
     else if(compArgStr(p+1, "DOTDUMP", 7)){
        Ldotdump = TRUE; /*111202dcr Ldotdump in examineOneDotEach() */
     }
	 else if(compArgStr(p+1, "QUIET", 1)){
	    Verbose = FALSE;
	 }
	 else if(compArgStr(p+1, "EXPLICIT", 1)){
	    ImplicitH = FALSE;
	 }
	 else if(compArgStr(p+1, "IMPLICIT", 1)){
	    ImplicitH = TRUE;
	 }
	 else if(compArgStr(p+1, "DROP", 2)){
	    *keepUnselected = FALSE;
	 }
	 else if(compArgStr(p+1, "KEEP", 2)){
	    *keepUnselected = TRUE;
	 }
	 else if((n = compArgStr(p+1, "DENSITY", 2))){
	    *density = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "RADIUS", 1))){
	    *probeRad   = parseReal(p, n+1, 10);

		  /* parameters correlated with radius */
	    if ((*probeRad) > GAPweight) {
	       GAPweight   = *probeRad;
	    }
	    if ((*probeRad) > HighGoodCut) {
	       HighGoodCut = *probeRad;
	    }
	 }
	 else if((n = compArgStr(p+1, "SCALEVDW", 5))   /* new name,*/
	      || (n = compArgStr(p+1, "VDWSCALE", 3))){ /* old name */
	    RadScaleFactor = parseReal(p, n+1, 10);
	    if (RadScaleFactor < 0.0001
	    ||  RadScaleFactor > 1000.0) {
	       sprintf(message, "invalid scale factor: -scalevdw%g",
			RadScaleFactor);
	       halt(message);
	    }
	 }
	 else if((n = compArgStr(p+1, "ADDVDW", 3))){
	    RadScaleOffset = parseReal(p, n+1, 10);
	    if (RadScaleOffset < -10.0
	    ||  RadScaleOffset > 1000.0) {
	       sprintf(message, "invalid vdw offset: -addvdw%g",
			RadScaleOffset);
	       halt(message);
	    }
	 }
	 else if((n = compArgStr(p+1, "SPIKE", 2))){
	    *drawSpike = TRUE;
	    if (p[n+1]) {
	       *spikelen = parseReal(p, n+1, 10);
	    }
	 }
	 else if(compArgStr(p+1, "NOSPIKE", 3)){
	    *drawSpike = FALSE;
	 }
	 else if(compArgStr(p+1, "KINEMAGE", 3)){
	    *addKinToFile = TRUE;
	    *countDots = FALSE; /* also forces kin format */
	    *rawOutput = FALSE;
	    *conFlag = FALSE; // by SJ 10/07/2011
	    OutputFmtType = 0;
	 }
	 else if(compArgStr(p+1, "NOGROUP", 3)){
	    *sayGroup = FALSE;
	 }
	 else if(compArgStr(p+1, "ELEMENT", 4)){
	    AtomMasters = TRUE;
	 }
	 else if(compArgStr(p+1, "COUNTDOTS", 1)){
	    *countDots = TRUE;
	 }
	 else if(compArgStr(p+1, "UNFORMATED", 1)){
	    *rawOutput = TRUE;
	 }
	 else if(compArgStr(p+1, "CONDENSE", 3)){
         *conFlag = TRUE;// added by SJ 10/07/2011 - to parse the flag for condened raw output
	 }
    else if(compArgStr(p+1, "WEAKHbonds", 5)){
         LweakHbonds = TRUE; /*111215dcr*/
    }
    else if(compArgStr(p+1,"SEPWORSE",8)){
        LworseOverlap = TRUE; /*04/09/2015 SJ - for separating the worse overlaps(>=0.5A) into their own list*/
    }
    else if(compArgStr(p+1, "XHvector", 2)){
        LXHvector = TRUE; /*20120120dcr*/
    }
    else if(compArgStr(p+1, "ONEDOTeach", 6)){
        LOneDotEach = TRUE; /*111013dcr*/
    }
    else if(compArgStr(p+1, "GAPBINs", 6)){
        Lgapbins = TRUE; /*111020dcr*/
    }
    else if(compArgStr(p+1, "NEAREST", 7)){
        Lnearest = TRUE; /*111021dcr*/
    }
    else if(compArgStr(p+1, "FILTERnearest", 6)){
        Lfilternearest = TRUE; /*111022dcr*/
        Lnearest = TRUE; /*to make the logic explicit*/
    }
    else if(compArgStr(p+1, "honlyocclude", 5)){
        Lhonlyocclude = TRUE; /*130427dcr re ala-dipept gamma turn vs C occlusion*/
    }
    else if(compArgStr(p+1, "OCC1filter", 4)){
        Lfilterocc1 = TRUE; /*20120120dcr*/
    }
    else if(compArgStr(p+1, "NOALTfilter", 5)){
        Lfilternoalt = TRUE; /*20120120dcr*/
    }
	 else if(compArgStr(p+1, "OUTSIDE", 2)){
	    *method = EXTERNALSURFACE;
	    nargs = 1;
	 }
	 else if(compArgStr(p+1, "ONCE", 2)){
	    *method = INTERSECTONCE;
	    nargs = 2;
	 }
	 else if(compArgStr(p+1, "BOTH", 1)){
	    *method = INTERSECTBOTHWAYS;
	    nargs = 2;
	 }
	 else if(compArgStr(p+1, "SELF", 2)){
	    *method = SELFINTERSECT;
	    nargs = 1;
	 }
	 else if(compArgStr(p+1, "DUMPATOMINFO", 5)){
	    *method = DUMPATOMCOUNT;
	    nargs = 1;
	 }
	 else if(compArgStr(p+1, "SCSURFACE", 3)){
	    *method = EXTERNALSURFACE;
	    nargs = 1;
	    argcnt = 1;
	    * srcArg = "not water";
	    *keepUnselected = FALSE;
	    *probeRad      = 1.4;
	    *groupLabel = "SCS";
	 }
	 else if(compArgStr(p+1, "EXPOSED", 4)){
	    *method = EXTERNALSURFACE;
	    nargs = 1;
	    *keepUnselected = FALSE;
	    *probeRad      = 1.4;
	    *groupLabel = "SCS";
	 }
	 else if(compArgStr(p+1, "ASURFACE", 2)){
	    *method = EXTERNALSURFACE;
	    nargs = 1;
	    argcnt = 1;
	    * srcArg = "not water";
	    *keepUnselected = FALSE;
	    RadScaleOffset = 1.4;
	    *probeRad      = 0.0;
	    *groupLabel = "AS";
	 }
	 else if(compArgStr(p+1, "ACCESS", 6)){
	    *method = EXTERNALSURFACE;
	    nargs = 1;
	    *keepUnselected = FALSE;
	    RadScaleOffset = 1.4;
	    *probeRad      = 0.0;
	    *groupLabel = "AS";
	 }
	 else if((n = compArgStr(p+1, "SCAN", 4))){
	    int scanType = parseInteger(p, n+1, 10);
	    switch(scanType) {
	    case 0:
	       *method = SELFINTERSECT;
	       nargs = 1;
	       argcnt = 1;
	       * srcArg = "alta blt40 ogt33";
	       Maxbonded = 4; /*changed from 3 : dcr041017*/
	       DoMcMc = TRUE;
	       break;
	    case 1:
	       *method = INTERSECTONCE;
	       nargs = 2;
	       argcnt = 2;
	       * srcArg = "sc alta blt40 ogt33";
	       *targArg = "alta blt40 ogt65,(not water ogt33)";
	       Maxbonded = 4; /*changed from 3 : dcr041017*/
	       break;
	    default:
	       sprintf(message,
		  "invalid scan type %d, for more info use -help",
		     scanType);
	       halt(message);
	       break;
	    }
	 }
	 else if(compArgStr(p+1, "SEGID", 5)){
	    SplitBySegID = TRUE;
	 }
	 else if(compArgStr(p+1, "LIMIT", 3)){
	    LimitDots = TRUE;
	 }
	 else if(compArgStr(p+1, "NOLIMIT", 5)){
	    LimitDots = FALSE;
	 }
	 else if(compArgStr(p+1, "LENS", 3)){
	    LensDots = TRUE;
	 }
	 else if(compArgStr(p+1, "NOLENS", 5)){
	    LensDots = FALSE;
	 }
	 else if(compArgStr(p+1, "MC", 2)){
	    DoMcMc = TRUE;
	 }
	 else if(compArgStr(p+1, "DUMPH2O", 5)){
	    DumpNewHO = TRUE;
	 }
	 else if(compArgStr(p+1, "HETS", 3)){
	    DoHet = TRUE;
	 }
	 else if(compArgStr(p+1, "NOHETS", 5)){
	    DoHet = FALSE;
	 }
	 else if(compArgStr(p+1, "WAT2WAT", 4)){
	    DoWatWat = TRUE; /* should be before WATERS */
	    DoH2O = TRUE;
	 }
	 else if(compArgStr(p+1, "WATERS", 3)){
	    DoH2O = TRUE;
	 }
	 else if(compArgStr(p+1, "NOWATERS", 5)){
	    DoH2O = FALSE;
	 }
	 else if(compArgStr(p+1, "4H", 1)){
	    Maxbonded = 4;
	 }
	 else if(compArgStr(p+1, "3", 1)){
	    Maxbonded = 3;
	 }
	 else if(compArgStr(p+1, "2", 1)){
	    Maxbonded = 2;
	 }
	 else if(compArgStr(p+1, "1", 1)){
	    Maxbonded = 1;
	 }
	 else if(compArgStr(p+1, "POLARH", 1)){
	    UsePolarH = TRUE;
	 }
	 else if(compArgStr(p+1, "NOPOLARH", 3)){
	    UsePolarH = FALSE;
	 }
	 else if(compArgStr(p+1, "NUCLEAR", 3)){
	   NuclearRadii = TRUE;
	 }
	 else if(compArgStr(p+1, "NOFACEHBOND", 6)){
	    HB2aromFace = FALSE;
	 }
	 else if(compArgStr(p+1, "NAME", 1)){
	    if (++i < argc) {
	       *groupLabel = argv[i];
	    }
	    else {
	       halt("no group name after -Name flag");
	    }
	 }
	 else if(compArgStr(p+1, "DOTMASTER", 9)){
                 LMasterName = TRUE; /*extra master={name} on lists 060129*/
	 }
	 else if(compArgStr(p+1, "IGNORE", 6)){
	    if (++i < argc) {
	       *ignoreArg = argv[i];
	    }
	    else {
	       halt("no pattern after -IGNORE flag");
	    }
	 }
	 else if((n = compArgStr(p+1, "COSCALE", 4))){
	    CORadScale = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "HBREGULAR", 3))){
	    Min_regular_hb_cutoff = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "HBCHARGED", 3))){
	    Min_charged_hb_cutoff = parseReal(p, n+1, 10);
	 }
	 else if(compArgStr(p+1, "DOCHO", 5)){
	    PermitCHXHB = TRUE;
	 }
	 else if((n = compArgStr(p+1, "CHO", 3))){
	    CHOHBfactor = parseReal(p, n+1, 10);
	 }
	 else if(compArgStr(p+1, "ATOMCOLOR", 4)){
	    ColorGap = FALSE;
	 }
	 else if(compArgStr(p+1, "GAPCOLOR", 3)){
	    ColorGap = TRUE;
	 }
	 else if(compArgStr(p+1, "BASECOLOR", 4)){
	    ByNABaseColor = TRUE;
	    ColorGap = FALSE; /* forces base color to be only type of display */
	 }
	 else if(compArgStr(p+1, "COLORBASE", 6)){
	    ByNABaseColor = TRUE;
	 }
	 else if((n = compArgStr(p+1, "GAPWEIGHT", 4))){
	    GAPweight = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "BUMPWEIGHT", 5))){
	    BUMPweight = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "HBWEIGHT", 3))){
	    HBweight = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "DIVLow", 4))){
	    LowGoodCut = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "DIVHigh", 4))){
	    HighGoodCut = parseReal(p, n+1, 10);
	 }
     else if((n = compArgStr(p+1, "DIVWorse", 4))){// 04/21/2015 SJ to change the default of the worse overlap
        WorseBumpCut = parseReal(p, n+1, 10);
     }
	 else if((n = compArgStr(p+1, "MINOCCupancy", 6))){
	    OccupancyCutoff = parseReal(p, n+1, 10);
	 }
	 else if((n = compArgStr(p+1, "NOOCCupancy", 5))){
	    OccupancyCriteria = FALSE; /*060831*/
	 }
	 else if(compArgStr(p+1, "OFORMAT", 7)){
	    OutputFmtType = 1;
	 }
	 else if(compArgStr(p+1, "XVFORMAT", 8)){
	    OutputFmtType = 2;
	 }
	 else if(compArgStr(p+1, "ONELINE", 7)){ /*oneline dcr041101*/
            ContactSUMMARY = TRUE;
	    OutputFmtType = 3;
            argcnt = 1; /*recognize next naked item as an input file name*/
	 }
	 else if(compArgStr(p+1, "SUMMARY", 7)){ /*summary dcr041101*/
            ContactSUMMARY = TRUE;
            argcnt = 1; /*recognize next naked item as an input file name*/
	 }
	 else if(compArgStr(p+1, "DEFAULTS", 7)){ /*defaults dcr041101*/
	    *srcArg = "altA ogt33";
	    *method = SELFINTERSECT;
	    Maxbonded = 4; /*changed from 3 : dcr041017*/
	    DoMcMc = TRUE;
	    DoHet  = TRUE;
            argcnt = 1; /*recognize next naked item as an input file name*/
            /*like naked defaults but allow flags like -summary or -oneline*/
	 }
	 else if(compArgStr(p+1, "NOHBOUT", 7)){
	    OutputHBs = FALSE;
	 }
	 else if(compArgStr(p+1, "NOCLASHOUT", 10)){
	    OutputClashes = FALSE;
	 }
	 else if(compArgStr(p+1, "NOVDWOUT", 8)){
	    OutputVDWs = FALSE;
	 }
	 else if(compArgStr(p+1, "ONLYBADOUT", 10)){ /*dcr041010*/
	    OnlyBadOut = TRUE;
	 }
	 else if(compArgStr(p+1, "NOTICKS", 6)){
	    ShowTicks = FALSE;
	 }
	 else if(compArgStr(p+1, "STDBONDS", 7)){
	    UseStdBond = TRUE;
	 }
	 else if(compArgStr(p+1, "OLDU", 4)){
	    OldStyleU = TRUE;
	    *rawOutput = TRUE;
	 }
#ifdef JACKnANDREA
         else if(compArgStr(p+1, "WRITEHBONDS", 6)){
            writeHbonds = TRUE;
         }
#endif
	 else if(compArgStr(p+1, "NOPARENT", 8)){
	    UseHParent = FALSE;
	 }
	 else if(compArgStr(p+1, "OUTCOLOR", 6)){
	    if (++i < argc) {
	       OutPointColor = argv[i];
	    }
	    else {
	       halt("no color name after -OUTCOLolor flag");
	    }
	 }
	 else if(compArgStr(p+1, "AUTObondrot", 4)){
	    if (++i < argc) {
	       p = argv[i];
	       if (p[0] == '-') {
		  if (p[1] == '\0') {
		     file++;
		     if (mabip) {
			mabip->filenum  = file;
			mabip->inf      = stdin; /*input of autobondrot atoms*/
			mabip->close    = FALSE;
		     }
		  }
		  else {
		     sprintf(message, "%s flag instead of filename after -AUTObondrot", p);
		     halt(message);
		  }
	       }
	       else {
		  file++;
		  inf = fopen(p, "rb"); /*p holds autobondrot input file name*/
		  if (inf) {
		     if (mabip) {
			mabip->filenum  = file;
			mabip->inf      = inf;
			mabip->close    = TRUE;
		     }
		  }
		  else {
                    sprintf(message, "could not open -AUTOBONDROT file: %s", p);
                    halt(message);
		  }
	       }
               /*Lautobondrot = TRUE;*/ /*050111 probe.c global logical */
	       *countDots = TRUE;
	       *rawOutput = TRUE;
               /*autobondrot mode does NOT read from input file at this time*/
	    }
	    else {
	       halt("no filename after -AUTObondrot flag");
	    }
	 }
	 else if(compArgStr(p+1, "REFerence", 3)){
	    sprintf(message, "Please cite: %s", referenceString);
	    note(message);
	    sprintf(message, "For more information see %s", electronicReference);
	    note(message);
	    exit(0);
	 }
	 else if((n = compArgStr(p+1, "CHANGEs", 6))){
	    dump_changes(stderr);
	 }
	 else if((n = compArgStr(p+1, "DEBUG", 5))){
	    DebugLevel = parseInteger(p, n+1, 10);
	 }
	 else { /*either naked -  or unrecognized character string*/
	    sprintf(message, "unrecognized flag, %s", p);
	    halt(message);
	 }
      }/* -flag item */
    else if (argcnt <= 0)
    { /*lonely naked item: store as a source argument*/
       *srcArg = p; /*if nothing else comes in, use as input file name*/
       argcnt = 1;
    }
    else if (argcnt == 1 && nargs == 2)
    {
       if((p[0] == '=') && (p[1] == '\0'))
       {
          *targArg = *srcArg;
       }
       else
       {
          *targArg = p;
       }
       argcnt = 2;
    }
    else if(Lfilternearest && argcnt == 2)
    {
       *extraArg = p; /*111021dcr*/
       /* -onedotonly -nearest "srcPat" "targPat" "extraPat" output filter*/
       argcnt = 3;
    }
    else
    {
       file++;
       inf = fopen(p, "rb");
       if (inf)
       {
          strcpy(inputfilename,p); /*dcr041023*/
            /*read atoms from input file now*/
          atomlist = loadAtoms(inf, atomlist, bboxA, file, reslstptr);
          fclose(inf);
       }
       else
       {
          sprintf(message, "could not open file %d: %s", file, p);
          halt(message);
       }
     }
   }
   if (file < 1)
   {
      if (argc == 2 && argcnt == 1 && *srcArg)
      { /*naked item stored as srcArg*/
         p = *srcArg; /*interpret as an input file name, then presume defaults*/
         note("*SIMPLE SELF DOTS (for more info use -help flag and look for <<NO FLAGS>>)");
         *srcArg = "altA ogt33";
         *method = SELFINTERSECT;
         Maxbonded = 4; /*changed from 3 : dcr041017*/
         DoMcMc = TRUE;
         DoHet  = TRUE;
         file++;
         inf = fopen(p, "rb");
         if (inf)
         {
            strcpy(inputfilename,p); /*dcr041023*/
            /*read atoms from input file now*/
            atomlist = loadAtoms(inf, atomlist, bboxA, file, reslstptr);
            fclose(inf);
         }
         else
         {
            sprintf(message, "could not open file %d: %s", file, p);
            halt(message);
         }
      }
      else
      {
         errmsg("no input files");
         probehelp(0);
      }
   }

   return atomlist; /*loaded under some conditions, but not if autobondrot!*/
}/*processCommandline()*/
/*}}}processCommandline()____________________________________________________*/

/*{{{initEndData()************************************************************/
void initEndData(chainEndData_t *ed)
{
   int i = 0;
   for (i = 0; i < 8; i++) { ed->ambigN[i] = NULL; ed->ambigO[i] = NULL; }
   ed->res_mc_oxy_cnt = 0;
   ed->first = 0;
   ed->last  = 0;
   ed->Ntmarkers = 0;
   ed->Ctmarkers = 0;
   for (i = 0; i < 4; i++) { ed->hisNHcount[i] = 0; }
}
/*}}}initEndData()___________________________________________________________*/

/*{{{hisNHcheck()*************************************************************/
/* His ring atoms start out with a positive charge.                       */
/* Here we remove the charge for His residues without two ring NH protons.*/
void hisNHcheck(chainEndData_t *ed, atom *atomlist, int rescnt)
{
   int i = 0;
   atom *rat = NULL;

   if (ed->hisNHcount[0]) { /* the zeroth counter flags: Is this a HIS? */
      /*is a histidine*/
      for(rat = atomlist; rat; rat = rat->next) {
	 if ((rat->r->rescnt == rescnt) && (rat->props & POSITIVE_PROP)) {
	    char altc = toupper(rat->altConf);
	    int ialt = 0;

	    if ((altc==' ')||(altc=='A')||(altc=='1')) { ialt = 1; }
	    else if ((altc=='B')||(altc=='2'))         { ialt = 2; }
	    else if ((altc=='C')||(altc=='3'))         { ialt = 3; }

	    if (ialt > 0) {
	       if (ed->hisNHcount[ialt] < 2) { /* his not positive */
		  rat->props &= ~POSITIVE_PROP;
	       }
	       else {  /* his positive but not an acceptor */
		  rat->props &= ~ACCEPTOR_PROP;
	       }
	    }
	 }
      }

      for (i = 0; i < 4; i++) { ed->hisNHcount[i] = 0; }
   }/*is a histidine*/
}
/*}}}hisNHcheck()____________________________________________________________*/

/*{{{resCheck()***************************************************************/
/* called after each residue is read */
void resCheck(chainEndData_t *ed, atom *atomlist, int rescnt)
{
    hisNHcheck(ed, atomlist, rescnt);
}
/*}}}resCheck()______________________________________________________________*/

/*{{{CtermCheck()*************************************************************/
/* called when a new residue begins and at the end of all processing */
/* to process the previous residue                                   */
void CtermCheck(chainEndData_t *ed, int rescnt, int isChainEnd)
{
   int i = 0;
   /* avoid processing the first time 'cause there ain't no chain pending */

   if (isChainEnd && rescnt) {
      ed->last = rescnt; /* last residue */

      /* see if we can put the pieces together to determine end charge */

      for (i = 0; i < 4; i++) { /* only array[0-3] contains first Ns */
	 if (ed->ambigN[i]) {
	    if (ed->ambigN[i]->r->rescnt == ed->first
	    &&  ed->Ntmarkers         == ed->first) {
	       ed->ambigN[i]->props |= POSITIVE_PROP;
	    }
	 }
	 else { break; }
      }
      for (i = 0; i < 8; i++) { /* array[0-7] contains last Os */
	 if (ed->ambigO[i]) {
	    if (ed->ambigO[i]->r->rescnt == ed->last
	    &&  ed->Ctmarkers         == ed->last) {
	       ed->ambigO[i]->props |= NEGATIVE_PROP;
	    }
	 }
	 else { break; }
      }

      initEndData(ed); /* reset the end data record */
   }
}
/*}}}CtermCheck()____________________________________________________________*/

/*{{{noticedNt()**************************************************************/
/* called when we have an indicator that residue is at N terminus */
/* to look for ambigN and mark charged */
void noticedNt(chainEndData_t *ed, int rescnt)
{
   int i = 0;

   for (i = 4; i < 8; i++) { /* only array[4-7] contains last Ns */

      if (ed->ambigN[i]) {
	 if (ed->ambigN[i]->r->rescnt == rescnt) {
	    ed->ambigN[i]->props |= POSITIVE_PROP;
	 }
      }
      else { break; }
   }
}
/*}}}noticedNt()_____________________________________________________________*/

/*{{{noticedCt()**************************************************************/
/* called when we have an indicator that residue is at C terminus */
/* to look for ambigO and mark charged */
void noticedCt(chainEndData_t *ed, int rescnt)
{
   int i = 0;

   for (i = 0; i < 8; i++) { /* array[0-7] contains last Os */
      if (ed->ambigO[i]) {
	 if (ed->ambigO[i]->r->rescnt == rescnt) {
	    ed->ambigO[i]->props |= NEGATIVE_PROP;
	 }
      }
      else { break; }
   }
}
/*}}}noticedCt()_____________________________________________________________*/

/*{{{NtermCheck()*************************************************************/
/* called whenever a new residue begins */
void NtermCheck(chainEndData_t *ed, int rescnt, int isChainEnd)
{
   int i = 0;

   ed->res_mc_oxy_cnt = 0; /* reset Oxygen data for each residue */
   for (i = 0; i < 8; i++) { ed->ambigO[i] = NULL; }

   /* we have reserved the top 4 nitrogen slots for most recent res */
   /* so we have to clear this here                                 */
   for (i = 4; i < 8; i++) { ed->ambigN[i] = NULL; }

   if (isChainEnd) {
      ed->first = rescnt; /* first residue */
   }
}
/*}}}NtermCheck()____________________________________________________________*/

/*{{{ProcessResInfo()*** called from loadAtoms() & movingAtomListProcessing()*/
/* called for each atom read in */
void ProcessResInfo(chainEndData_t *ed, atom *a)
{
   int i = 0;

   if (strstr("HIS:his", a->r->resname)) {/*is a histidine*/
      ed->hisNHcount[0] = 1; /* indicates this is a HIS */
      if ( strstr(" HD1: hd1: HE2: he2", a->atomname)) {
	 char altc = toupper(a->altConf);
	 if ((altc==' ')||(altc=='A')||(altc=='1')) {
	    ed->hisNHcount[1]++;
	 }
	 else if ((altc=='B')||(altc=='2')) {
	    ed->hisNHcount[2]++;
	 }
	 else if ((altc=='C')||(altc=='3')) {
	    ed->hisNHcount[3]++;
	 }
      }
   }/*is a histidine*/

   /* look for indicators that we are at the end of a chain */

   if (a->props & CHECK_ENDS_PROP) {
      if (isHatom(a->elem)) {
	 if (! ed->Ntmarkers) { ed->Ntmarkers = a->r->rescnt; }
	 noticedNt(ed, a->r->rescnt);
      }
      else if ( (a->elem == atomO) &&
	       ((a->altConf==' ')||(a->altConf=='A')||(a->altConf=='1')) ){
	 ed->res_mc_oxy_cnt++;
	 if (ed->res_mc_oxy_cnt == 2) {
	    ed->Ctmarkers = a->r->rescnt;
	    noticedCt(ed, a->r->rescnt);
	 }
      }
      else if (strstr(a->atomname, "OXT")) { /* check SUBSET of name */
	 ed->Ctmarkers = a->r->rescnt;
	 noticedCt(ed, a->r->rescnt);
      }
   }

   /* save pointers to the first ambigNs and the last ambigOs and Ns */
   /* we save multiple atoms to both find both possible mc oxygens   */
   /* and to handle multiple conformations (we can do 4 max)         */
   /* All atoms must have the same residue counter                   */

   if (a->props & MAYBECHG_PROP) {
      if (!strcmp(a->atomname, " N  ")) {

         if (ed->first == a->r->rescnt) { /* mark the Nterm -- jmw 20011001 */
            a->props |= CHECK_ENDS_PROP;
         }

	 for (i = 0; i < 4; i++) {/* first Ns[0-3] (cleared at chain ends)*/
	    if (ed->ambigN[i] == NULL) {
	       if ((i == 0) || (ed->ambigN[i-1]->r == a->r)) {
		  ed->ambigN[i] = a;
	       }
	       break;
	    }
	 }
	 for (i = 4; i < 8; i++) {   /* last Ns[4-7] (cleared each res)*/
	    if (ed->ambigN[i] == NULL) {
	       if ((i == 4) || (ed->ambigN[i-1]->r == a->r)) {
		  ed->ambigN[i] = a;
	       }
	       break;
	    }
	 }
      }
      if (!strcmp(a->atomname, " O  ")) { /* last Os[0-7]  (cleared each res)*/
	 for (i = 0; i < 8; i++) {
	    if (ed->ambigO[i] == NULL) {
	       if ((i == 0)|| (ed->ambigN[i-1]->r == a->r)) {
		  ed->ambigO[i] = a;
	       }
	       break;
	    }
	 }
      }
   }

}
/*}}}ProcessResInfo()________________________________________________________*/

/*{{{movingAtomListProcessing()***********************************************/
void movingAtomListProcessing(atom *initAtomLst, void *userdata)
{
   int previd = 0, rescnt = 0;
   char prevInsCode = ' ';
   char prevChain[5];
   atom *a = NULL, *prevAtomLst = NULL, *nexta = NULL;
   chainEndData_t endData;

   initEndData(&endData);

   prevChain[0] = prevChain[1] = '?';
   prevChain[2] = '\0';

   for(a=initAtomLst; a; a = nexta) {
      nexta = a->next;

      a->next = prevAtomLst; /* remove a from initial atom list */
      prevAtomLst = a;       /* and build up a reversed list of */
                             /* only atoms previously seen      */

      if (a->r->resid != previd || a->r->resInsCode != prevInsCode
                             || strcmp(a->r->chain, prevChain) != 0) {

	 if (rescnt){
	    resCheck(&endData, prevAtomLst, rescnt);
	 }
	 CtermCheck(&endData, rescnt, (strcmp(a->r->chain, prevChain) != 0));

	 ++rescnt;

	 NtermCheck(&endData, rescnt, (strcmp(a->r->chain, prevChain) != 0));
	 previd = a->r->resid;
	 prevInsCode = a->r->resInsCode;
	 strcpy(prevChain, a->r->chain);
         //prevChain[2] = '\0';
      }
      a->r->rescnt = rescnt;

      ProcessResInfo(&endData, a);
   }

   if (rescnt) { resCheck(&endData, prevAtomLst, rescnt); }
   CtermCheck(&endData, rescnt, TRUE);
}
/*}}}movingAtomListProcessing()______________________________________________*/

/*{{{loadAtoms()************* called from processCommandline() ***************/
atom* loadAtoms(FILE *fp, atom *atomlist, region *boundingBox, int file,
		  residue **reslstptr)
{
   int previd = 0, rescnt = 0, model = 0;
   char *rec, prevInsCode = ' ';
   char prevChain[5];
   atom *a = NULL;
   residue *scratchRes;
   chainEndData_t endData;

   prevChain[0] = prevChain[1] = '?';
   prevChain[2] = '\0';

   scratchRes = newResidueData();

   initEndData(&endData);

   while((rec = getPDBrecord(fp))) { /*loop over all records in pdb file*/
      if (isTer(rec)) {
         prevChain[0] = prevChain[1] = '?';
         prevChain[2] = '\0';
      }
      if ((isAtom(rec) || isHet(rec)) && ! isPseudoAtom(rec)) {/*atom or htatm*/
	 a = newAtom(rec, file, model, scratchRes);
	 if (!a) { break;}

	 if (a->r->resid != previd || a->r->resInsCode != prevInsCode
	                        || strcmp(a->r->chain, prevChain) != 0) {

	    if (rescnt){
	       resCheck(&endData, atomlist, rescnt);
	    }
	    CtermCheck(&endData, rescnt, (strcmp(a->r->chain, prevChain) != 0));

	    ++rescnt;

	    NtermCheck(&endData, rescnt, (strcmp(a->r->chain, prevChain) != 0));
	    previd = a->r->resid;
	    prevInsCode = a->r->resInsCode;
	    strcpy(prevChain, a->r->chain);
            //prevChain[2] = '\0';
	 }
	 a->r->rescnt = rescnt;

	 ProcessResInfo(&endData, a);

	 if (ImplicitH && isHatom(a->elem)) {
	    deleteAtom(a); /* filter out implicit hydrogens */
	    continue;
	 }
	 if (isDummyAtom(*a)) {
	    deleteAtom(a);  /* filter out dummy Xplor atoms */
	    continue;
	 }

	 if (!atomlist) {
	    boundingBox->min = a->loc;
	    boundingBox->max = a->loc;
	 }
	 a->next = atomlist;
	 atomlist = a;

	 updateBoundingBox(&(a->loc), boundingBox);

	 if (*reslstptr == NULL) {/* first residue goes on residue list */
	    *reslstptr = scratchRes;
	    scratchRes = newResidueData();
	 }
	 else { /* otherwise we have to compare residue data blocks */
	    if (resDiffersFromPrev(*reslstptr, scratchRes)) {
	       scratchRes->nextRes = *reslstptr;
	       *reslstptr = scratchRes;
	       scratchRes = newResidueData();

	    }
	    else {
	       a->r = *reslstptr; /* same, so point to prior block */
	       a->r->a = a; /* makes this atom the head of atom list for this res */
	    }
	 }

      }/*atom or htatm*/
      if (isModel(rec)) {
	 model = parseModel(rec);
         modelNumber[++modelCount] = model; /*041114*/
if(Verbose)
{
   fprintf(stderr,"modelNumber[%d]==%d\n",modelCount,modelNumber[modelCount]);
}
      }
   }/*loop over all records in pdb file*/
   if (rescnt) { resCheck(&endData, atomlist, rescnt); }
   CtermCheck(&endData, rescnt, TRUE);

   deleteResidueData(scratchRes); /* clean up any extra residue */

   return atomlist;
}
/*}}}loadAtoms()_____________________________________________________________*/

/*{{{binAtoms()***************************************************************/
atomBins* binAtoms(atom *theAtoms, region *boundingBox, char serialNum,
		  float probeRad, int keepUnselected, int selflags)
{
   atom *a = NULL;
   atomBins *bins = NULL;

   bins = initBins(serialNum, boundingBox,
	    2.0*(getMaxRadius(ImplicitH)+probeRad)+0.2);
   if (bins) {
      for(a = theAtoms; a; a = a->next) {
	 if (keepUnselected || (a->flags & selflags)) {
	    if (! (a->flags & IGNORE_FLAG)) {
	       addNeighbor(a, bins);
	    }
            /*050118 this seems the only place that uses IGNORE_FLAG*/
	 }
      }
   }
   /* note: addNeighbor() also called in updateHydrogenInfo() */
   return bins;
}
/*}}}binAtoms()______________________________________________________________*/

/*{{{getRadius()**************************************************************/
float getRadius(int at, int useCOScale)
{
   float rad = 0.0, sf = RadScaleFactor;

   if (useCOScale){ sf *= CORadScale; }

   if (ImplicitH) { rad = getImplRad(at); }
   else           { rad = getExplRad(at); }

   return (rad*sf) + RadScaleOffset;
}
/*}}}getRadius()_____________________________________________________________*/

/*{{{newRingInfo()************************************************************/
ringInfo * newRingInfo(ringInfo **head, point3d* ctr, point3d* norm)
{
   ringInfo *ri = NULL;
   ri = (ringInfo *)malloc(sizeof(ringInfo));
   if (ri) {
      ri->nextRing = NULL;
      ri->ringCenter = *ctr;
      ri->ringNormal = *norm;

      if (head) { /* if list passed, link new ring as new head */
	 ri->nextRing = *head;
	 *head = ri;
      }
   }
   return ri;
}
/*}}}newRingInfo()___________________________________________________________*/

/*{{{deleteRingInfoList()*****************************************************/
void deleteRingInfoList(ringInfo *ri)
{ /* kill entire list */
   ringInfo *p = NULL, *nxt = NULL;
   p = ri;
   while(p) {
      nxt = p->nextRing;
      free(p);
      p = nxt;
   }
}
/*}}}deleteRingInfoList()____________________________________________________*/

/*{{{newResidueData()*********************************************************/
residue * newResidueData()
{ /* new blank residue */
   residue *r = NULL;
   r = (residue *)malloc(sizeof(residue));
   if (r) {
      r->nextRes = NULL;
      r->a = NULL;
      r->file = 0;
      r->model = 0;
      //r->chain[0] = ' ';   /*RMI I don't know what to do here*/
      r->chain[0] = ' ';
      r->chain[1] = ' ';
      r->chain[2] = '\0';
      r->resid = ' ';
      r->resInsCode = ' ';
      r->rescnt = 0;
      r->segid[0] = '\0';

      r->ring = NULL;
   }
   return r;
}
/*}}}newResidueData()________________________________________________________*/

/*{{{resDiffersFromPrev()*****************************************************/
int resDiffersFromPrev(residue *r1, residue *r2)
{
   if (r1 == NULL || r2 == NULL) { return 1; }
   return (! IS_THE_SAME_RES(r1, r2));
}
/*}}}resDiffersFromPrev()____________________________________________________*/

/*{{{deleteResidueData()******************************************************/
void deleteResidueData(residue *r)
{
   if (r) {
      deleteRingInfoList(r->ring);
      free(r);
   }
}
/*}}}deleteResidueData()_____________________________________________________*/

/*{{{disposeListOfResidues()**************************************************/
void disposeListOfResidues(residue *theRes)
{
   residue *r = NULL, *nextr = NULL;

   for(r=theRes; r; r = nextr) {
      nextr = r->nextRes;
      deleteResidueData(r);
   }
}
/*}}}disposeListOfResidues()_________________________________________________*/

/*{{{dumpRes()****************************************************************/
/* dumpRes() used for debugging */
void dumpRes(residue *theRes)
{
   residue *r = NULL;
   atom *a = NULL;
   for(r = theRes; r; r = r->nextRes) {
      fprintf(stderr, "RES{%s%d%c}[%d] ring(%p)\n",
		r->resname, r->resid, r->resInsCode, r->rescnt,
		r->ring);
      fprintf(stderr, "  ");
      for(a = r->a; a && (a->r == r); a = a->next) {
	 fprintf(stderr, "%s ", a->atomname);
      }
      fprintf(stderr, "\n");
   }
}
/*}}}dumpRes()_______________________________________________________________*/

/*{{{deleteAtom()*************************************************************/
void deleteAtom(atom *a)
{
   if (a) {
      free(a);
   }
}
/*}}}deleteAtom()____________________________________________________________*/

/*{{{disposeListOfAtoms()*****************************************************/
void disposeListOfAtoms(atom *theAtoms)
{
   atom *a = NULL, *nexta = NULL;

   for(a=theAtoms; a; a = nexta) {
      nexta = a->next;
      deleteAtom(a);
   }
}
/*}}}disposeListOfAtoms()____________________________________________________*/

/*{{{newAtom() <--loadAtoms(),newMovingAtom(); -->various property tests *****/
atom * newAtom(char *rec, int file, int model, residue * resDataBuf)
{/*newAtom()*/
   atom *a = NULL;
   char msg[100];

   a = (atom *)malloc(sizeof(atom));
   if (a) {
      a->r = resDataBuf;

      a->next     = NULL;
      a->nextInBin= NULL;
      a->scratch  = NULL;
      a->mark = 0;    /* used to mark bonds */
      a->flags = 0;
      a->props = 0;
      parseResidueName(rec, a->r->resname);
      parseAtomName(rec, a->atomname);
      a->r->a = a; /* residue points back to this atom (can change for res) */
      a->r->file = file;
      a->r->model = model;
      parseChain(rec, a->r->chain);
     /* a->r->chain = parseChain(rec);*/
      a->r->resid = parseResidueNumber(rec);
      parseResidueHy36Num(rec, a->r->Hy36resno);
      a->r->resInsCode = parseResidueInsertionCode(rec);
      a->r->rescnt = 0;
      a->altConf = parseAltLocCode(rec);
      a->binSerialNum = '?'; /* set when we bin */
      if (strstr(":ASX:GLX:ASN:GLN:", a->r->resname) /* special case treats undecided */
      && (a->atomname[1] == 'A')) {             /* as an oxygen */
	 a->elem = identifyAtom(" O  ", a->r->resname, Verbose); /*dcr041007 allow warning  add resname to call */
	 sprintf(msg, "atom %s will be treated as oxygen", a->atomname);
	 warn(msg);
      }
      else { /* normal case */
	 a->elem = identifyAtom(a->atomname, a->r->resname, Verbose);/*dcr041007 allow warning  add resname to call */
      }

      /*next section seems to be the only place where atom->bondedto is set.*/
      /*select.c/setHydrogenParentName() and select.c/setMainchainBonding()*/
      /* both merely call stdconntable.c/searchForStdBondingPartner() */
      /* which constructs a search string from (a->r->resname, a->atomname) */
      /* and calls stdconntable.c/SearchStdResConnTable()  */
      /* which returns the string found in the StdResTblBucket[] hash table.*/
      /*This hash table is a marvelous tour de force listing of names of*/
      /*residue/atoms and what atoms they can be bonded to */

      a->bondedto = NULL;
      if (isHatom(a->elem)) {
	 a->atomclass = -1; /* to be specified later */
	 if (UseHParent) {
	    a->bondedto = setHydrogenParentName(a->r->resname, a->atomname);
	 }
      }
      else {
	 a->atomclass = a->elem;
	 if (UseStdBond) {
	    a->bondedto = setMainchainBonding(a->r->resname, a->atomname);

            /*setMainchainBonding() does not seem to have a distance limit*/
	 }
      }

      setProperties(a, isHet(rec), HB2aromFace, PermitCHXHB); /*select.c*/
        /*This is where (e.g.) both HET_PROP and DNA_PROP are set for the atom*/
      if (ByNABaseColor) {  /* forces coloring by base rather than atom type */
	 a->atomclass = naBaseCategory(a);
      }

      parseXYZ(rec, &(a->loc));
      a->ix = a->iy = a->iz = 0;

   /* note: H DONOR_PROP assignment done later in updateHydrogenInfo() */

      a->occ  = parseOccupancy(rec);
      a->bval = parseTempFactor(rec);
      parseSegID(rec, a->r->segid);

      //use SEGID instead of chainID if it exists and chainID is blank
      if ( (strncmp(a->r->segid, "    ",4) != 0) &&
           (strncmp(a->r->chain, "  ",2)   == 0)){
        //parseSegID(rec, a->r->chain);
        a->r->chain[0] = a->r->segid[0];
        a->r->chain[1] = a->r->segid[1];
        a->r->chain[2] = a->r->segid[2];
        a->r->chain[3] = a->r->segid[3];
        a->r->chain[4] = '\0';
      }

      a->radius = getRadius(a->elem, isCarbonylAtom(a));
      a->covRad = getCovRad(a->elem);
   }
   else {
      warn("could not allocate space for new atom");
   }
   return a;
}/*newAtom()*/
/*}}}newAtom()_______________________________________________________________*/

/*{{{selectSource() <--mainProbeProc(),newMovingAtom(); -->select/matchPat() */
void selectSource(atom *theAtoms, pattern *srcPat, int srcFlag,
		  pattern *targPat, int targFlg, pattern *ignorePat)
{
   atom *src = NULL;

   for(src = theAtoms; src; src = src->next)
   {
      if(   (DoHet || ! (src->props &   HET_PROP) )
         && (DoH2O || ! (src->props & WATER_PROP) ) )
      {
         /*060212 hets marked as prot,dna,rna could be excluded by these tests*/
         /* if src||target was specified as not-prot,dna,rna,... */
         /* seemingly because atoms of hets are also assigned those types */
	 if (srcPat && matchPat(src,  srcPat)) { src->flags |= srcFlag; }
	 if(targPat && matchPat(src, targPat)) { src->flags |= targFlg; }

	 /*if ((FABS(src->occ) <= OccupancyCutoff)        */
	 /*   || (ignorePat && matchPat(src, ignorePat))) */
         /*050119 not use occ cutoff here, so 0-occ atom will be put into bins*/
	 if (ignorePat && matchPat(src, ignorePat))
         {
	    src->flags = IGNORE_FLAG; /* overwrite any other settings */
            /*050118 this seems the only place that sets IGNORE_FLAG*/
	 }
      }
   }
}
/*}}}selectSource()__________________________________________________________*/

/*{{{atomsClose()*************************************************************/
/*atomsClose() only called from findTouchingAtoms() */
int atomsClose(atom *a, atom *b, float probeRad)
{
   int nearpt = FALSE;
   float lim = 0.0, dsq = 0.0;

   lim = a->radius + b->radius + probeRad + probeRad;

   dsq = v3distanceSq(&(a->loc), &(b->loc));

   /* if too close they must be the same atom...dummy atom?? */

/*removing the (dsq > 0.001) test actually removes one side of a clash!!*/
/*20111207 dcr changed test to .003*.003 = .00001 i.e. .001 in each coord */
/* to minimize false same-atom-ness but still allow filtering of dummy atoms*/
   if ((dsq > 0.00001) && (dsq <= (lim*lim)))
   {
      nearpt = TRUE;
   }
   return nearpt;
}
/*}}}atomsClose()____________________________________________________________*/

/*{{{inRange()****************************************************************/
int inRange(point3d *p, point3d *q, float lim)
{
   return v3distanceSq(p, q) <= (lim*lim);
}
/*}}}inRange()_______________________________________________________________*/

/*{{{perpendicular()**********************************************************/
void perpendicular(point3d *p, point3d *q, point3d *r, point3d *s)
{  /*cribbed from prekin perpendicular, 111018 */
   /*construct perpendicular from point to line*/
   /*perpendicular from 3rd point to line of points 1 and 2 */
   /* fourth point is on line at perpendicular to third point*/

   float  A1X,A1Y,A1Z,A2X,A2Y,A2Z,B1X,B1Y,B1Z,ADX,ADY,ADZ,B2X,B2Y,B2Z,ART;

        A1X = p->x; A1Y = p->y; A1Z = p->z;
        A2X = q->x; A2Y = q->y; A2Z = q->z;
        B1X = r->x; B1Y = r->y; B1Z = r->z;
        ADX = A2X - A1X;ADY = A2Y - A1Y;ADZ = A2Z - A1Z;
/*arbitrary point on line A
   B2X = A1X + ADX*ART;B2Y = A1Y + ADY*ART;B2Z = A1Z + ADZ*ART;

   direction components of line B1 --- B2
   BDX = A1X-B1X+ADX*ART; BDY = A1Y-B1Y+ADY*ART; BDZ = A1Z-B1Z+ADZ*ART;
   for this to be perpendicular to line A, sum of product of components ==0
   ADX*(A1X-B1X+ADX*ART)+ADY*(A1Y-B1Y+ADY*ART)+ADZ*(A1Z-B1Z+ADZ*ART)=0;
   rearrange
   ADX*A1X-ADX*B1X+ADX*ADX*ART
   +ADY*A1Y-ADY*B1Y+ADY*ADY*ART
   +ADZ*A1Z-ADZ*B1Z+ADZ*ADZ*ART=0;
   regroup
   ART*(ADX*ADX++ADY*ADY+ADZ*ADZ)
         =ADX*B1X-ADX*A1X+ADY*B1Y-ADY*A1Y+ADZ*B1Z-ADZ*A1Z;
*/
         ART =  (ADX*B1X-ADX*A1X+ADY*B1Y-ADY*A1Y+ADZ*B1Z-ADZ*A1Z)
               /(ADX*ADX+ADY*ADY+ADZ*ADZ);

         B2X = A1X + ADX*ART;
         B2Y = A1Y + ADY*ART;
         B2Z = A1Z + ADZ*ART;
         s->x = B2X;
         s->y = B2Y;
         s->z = B2Z;
}
/*}}}perpendicular()_________________________________________________________*/

/*{{{gapSize()****************************************************************/
float gapSize(point3d *p, point3d *q, float qrad)
{
   return v3distance(p, q) - qrad;
}
/*}}}gapSize()_______________________________________________________________*/

/*{{{findTouchingAtoms()******************************************************/
/* findTouchingAtoms() - Note: resets bond marks! Returns list of close atoms */
/*that are in the set of bins here called abins */
/*in the form of the first atom that is close to the src-atom where */
/*a linked list of all close atoms is made by the scratch member of each atom*/
/*which indicates another neighbor of the src-atom */

atom* findTouchingAtoms(atom *src, atom *head, atomBins *abins,
			   float probeRad, int targFlg, int *ok)
{/*findTouchingAtoms()*/
   atom *a = NULL;
   int i = 0, j = 0, k = 0, nearpt = 0, targcount = 0;
   int jx = (src->ix), jy = (src->iy), jz = (src->iz);
   int imin = 0, jmin = 0, kmin = 0;
   int imax = 0, jmax = 0, kmax = 0;

   if (src->binSerialNum != abins->binSerialNum)
   {/*in foreign bins */
      /*dcr?: why do they have to be in foreign==different bins??*/
      /*this requirement seems not to affect loss of close clashes*/

      /* must look up the positions */

      if(  (src->loc.x < (abins->min.x - abins->delta))
	|| (src->loc.y < (abins->min.y - abins->delta))
	|| (src->loc.z < (abins->min.z - abins->delta))
	|| (src->loc.x > (abins->max.x + abins->delta))
	|| (src->loc.y > (abins->max.y + abins->delta))
	|| (src->loc.z > (abins->max.z + abins->delta)) )
      {
	 return NULL; /* nothing anywhere nearby */
      }
      binLoc(&(src->loc), &jx, &jy, &jz, abins);

#ifdef DEBUG_A2B
fprintf(stderr, "f(%c != %c) %s%d %s [%d, %d, %d] <1..%d, 1..%d, 1..%d>\n",
   src->binSerialNum, abins->binSerialNum,
   src->r->resname, src->r->resid, src->atomname, jx, jy, jz, abins->nx-2, abins->ny-2, abins->nz-2);
#endif

   }/*in foreign bins */

   imin = jx-1; imax = jx+1; /* bin ranges */
   jmin = jy-1; jmax = jy+1;
   kmin = jz-1; kmax = jz+1;

   /* trim any excess edges */
   if (imin < 1)               { imin = 1; }
   if (jmin < 1)               { jmin = 1; }
   if (kmin < 1)               { kmin = 1; }
   if (imax > (abins->nx - 2)) { imax = (abins->nx - 2); }
   if (jmax > (abins->ny - 2)) { jmax = (abins->ny - 2); }
   if (kmax > (abins->nz - 2)) { kmax = (abins->nz - 2); }

   for(i = imin; i <= imax; i++) {
      for(j = jmin; j <= jmax; j++) {
         for(k = kmin; k <= kmax; k++) {

	    for(a = abins->list[i][j][k]; a; a = a->nextInBin) {

               /*Lmodeltest for interaction of diff. models 041112*/
               /*050121 Lmodeltest superseded, mage sends model# */
               /*...within subroutine findTouchingAtoms()... */
	       if( (src->r->model != a->r->model) ) { continue; }
               if( modelLimit > 1 ) /*041114*/
               {/*looping through multiple models*/
                  if(src->r->model != modelToProcess) { continue; }
               }

               /*if (( src->r  == a->r)            && -- alternate      */
                                                    /* alternate      */
	       if (( src->altConf != a->altConf) && /* conformations  */
	           ( src->altConf != ' ')        && /* don't interact */
	           (   a->altConf != ' ')     ) { continue; }
            
    
/*nearpt seems to be essential to get any dots or spikes */
/*...seems not to be where too-close atoms fail to get any dots or spikes*/
	       nearpt = atomsClose(src, a, probeRad);/*only call to atomsClose*/
	       if (nearpt)
               {/*add this atom to the linked-list of those touching src atom*/
		  a->mark = 0;       /* clear the bonding marker.   */
		  a->scratch = head; /* link using scratch pointers */
		  head = a;          /* with <a> at head of list    */

		  if (a->flags & targFlg) { targcount++; }
           }
        }
     }
    }
   }
   src->mark = 0;       /* clear the src atom bonding marker.   */

   *ok = !(targFlg && (targcount < 1));

   return head;
}/*findTouchingAtoms*/
/*}}}findTouchingAtoms()_____________________________________________________*/

/*{{{dotClassIndex()**********************************************************/
/* dotClassIndex() - maps dot type to class index   */
/*   OLD             t == 0  --> 0,1 : contact dot  */
/*                   t <  0  --> 2,3 : bump         */
/*                   t >  0  --> 4   : hbond        */
/*   NEW - SJ        t == 0  --> 0,1   : contact dot    */
/*  04/10/2015       t <  0  --> 3,4,5 : small, bad, worse bump*/
/*                   t >  0  --> 2,6   : weak Hbonds, Hbonds*/
int dotClassIndex(int t, float mingap)
{
   int idx = 0;
    
   // Old code for old definition
   //if      (t == 0) { idx = (mingap > HighGoodCut) ? 0 : 1; } /* contact */
   //else if (t  < 0) { idx = (mingap >  LowGoodCut) ? 2 : 3; } /* clash   */
   //else if (LweakHbonds)
   //{
   //        if (t==1) {idx = 4;                     } /*good hbond*/
   //        else      {idx = 5;                     } /*weak hbond 20111215dcr*/
   //}
   //else             { idx = 4;                              } /* hbonds  */
   
   /* new code SJ 04/10/2105 to implement the new definiion above*/
    if (t == 0) { /* contact */
        if (mingap > HighGoodCut)
            idx = 0; /* wide contact*/
        else
            idx = 1; /* close contact*/
    }
    else if(t < 0){ /* clash*/
        if (mingap > LowGoodCut)
            idx = 3; /* small overlap*/
        else if(LworseOverlap) /* divide bad and worse overlap only if this is true*/
        {
            if(mingap > WorseBumpCut)
                idx = 4; /* bad overlap */
            else
                idx = 5; /* worse overlap */
        }
        else
            idx = 4; /* bad overlap*/
    }
    else{ /* H bond */
        if(LweakHbonds){/* divide H bonds into weak only if this is true.*/
            if(mingap > 0)
                idx = 2; /* weak H bond*/
            else
                idx = 6; /* H bonds*/
        }
        else
            idx = 6; /* H bonds*/
    }
    return idx;
}
/*}}}dotClassIndex()_________________________________________________________*/

/*{{{saveDot() called from examineDots() and from SurfDots() *****************/
void saveDot(atom *src, atom *targ, int type, point3d *loc, point3d *spike,
   dotNode *results[][NODEWIDTH], int ovrlaptype, float mingap, char ptmaster, float angle)
{/*saveDot()*/
   /*ptmaster dcr041009*/
   dotNode* dot = NULL;
   int which = 0, idx = 0;

   /*overlaptype:  -1 bump, 0 touch, +1 H bond */
   /*entering from SurfDots needs no further Logical filtering*/
   /*move this decision making back up into examineDots()*/
   /*where we now will accummulate count numbers also*/
#ifdef OLDCODE  /*041020 saveDot() decisions moved up into examineDots*/
   if ((OutputHBs     && ovrlaptype  > 0)
   ||  (OutputClashes && ovrlaptype  < 0)
   ||  (OutputVDWs    && ovrlaptype == 0)) {
#endif

      idx = dotClassIndex(ovrlaptype, mingap);
      /*OLD: 0 wide contact,1 close contact,2 small overlap,3 bad overlap,4 H-bonds*/
      /*OLD: LweakHbonds: 4 good H-bonds, 5 weak (distant) H-bonds 20111215dcr*/
      /* NEW 04/10/2015 SJ: 0 wide contact, 1 close contact, 2 weak H-bonds, 3 small overlap, 4 bad overlap, 5 worse overlap, 6 H-bonds*/
      /*idx=1 for EXTERNALSURFACE ovrlaptype==0, mingap==0.0 */

         dot = newDot(src, targ, loc, spike, ovrlaptype, mingap, ptmaster, angle);
         /*ptmaster dcr041009*/ /*XHTangle dcr20120120*/
         if (dot)
         {/*add new dot to head of lists arrayed by type and severity */
             which = type;
             dot->next = results[which][idx];
             results[which][idx] = dot;
         }
#ifdef OLDCODE
   } /* otherwise we forget about the dot */
#endif
}/*saveDot()*/
/*}}}saveDot_________________________________________________________________*/

/*{{{newDot()*****************************************************************/
dotNode* newDot(atom *src, atom* targ, point3d *loc, point3d *spike,
          int ovrlaptype, float gap, char masterchr, float angle)
{ /*ptmaster dcr041009, XHTangle dcr20120120*/
   dotNode* d = NULL;

   d = (dotNode *)malloc(sizeof(dotNode));
   if (d) {
      d->next  = NULL;
      d->a     = src;
      d->t     = targ;
      d->loc   = *loc;
      d->spike = *spike;
      d->type  = ovrlaptype; /* -1 bump, 0 touch, +1 H bond */
      d->gap   = gap;
      d->ptmaster   = masterchr; /*dcr041009*/
      d->angle = angle; /* XHTangle dcr20120120*/
   }
   else warn("could not allocate space for new dot");

   return d;
}
/*}}}newDot()________________________________________________________________*/

/*{{{INRANGEINLINE()**********************************************************/
#define INRANGEINLINE(a, b, lim) \
(((((a).x-(b).x)*((a).x-(b).x))+\
  (((a).y-(b).y)*((a).y-(b).y))+\
  (((a).z-(b).z)*((a).z-(b).z))) <= ((lim)*(lim)))

  /*returns true when distance-sq between a and b <= limit-sq */
/*}}}INRANGEINLINE()_________________________________________________________*/

/*{{{examineOneDotEach()****** called from genDotIntersect() LOneDotEach *****/
void examineOneDotEach(atom *src, int type, atom *scratch,
      pointSet dots[], float probeRad,
      float spikelen, int targFlg, dotNode *results[][NODEWIDTH], atom *allMainAtoms)
                /*allMainAtoms20120120*/
{
   /*called from genDotIntersect() with same names except scratch==atomList2 */
   /*cloned and modified from examineDots to do only the OneDotEach option 111205*/
   /*111013 switch to dump details... 111202 now input param: -dotdump */
   float neighborhoodsq,neighbor2dotsq;
   atom *holdsrc;
   atom *holdcause;
   int holdtype;
   point3d holddotloc;
   point3d holdspikeloc;
   dotNode *holdresults;
   int holdovrlaptype;
   float holdmingap = 999.9;
   char holdptmaster;
   int Lholdingdot = FALSE; /*111021dcr src and targ can have diff patterns*/
   int Lbetweenatoms = FALSE; /*111202dcr*/

   atom *targ = NULL, *cause = NULL;
   atom *extratarg = NULL; /*extra outer loop over neighbors, 111013*/
   int i = 0, nearpt = 0, within = 0, ok = 0, maskHB = 0;
   int ndots = 1; /*for only One Dot for Each neighbor LOneDotEach 111008*/
   int skipthisMCMC = 0, isaHB = 0, tooCloseHB = 0, ovrlaptype = 0;
   float gap = 0.0, mingap = 0.0, sl = 0.0, hbondbumpgap = 0.0;
   float srcCauseDist = 0.0; /*check for disaster overlap condition*/
   float hbcutoff = 0.0, targVDW = 0.0;
   point3d targetloc, dotloc, dotvect, exploc, spikeloc;
   point3d gapvect, gaploc;
   point3d perploc; /*LOneDotEach occluder perpendicular to src-targ 111018*/
   point3d distalloc,distalvec; /*point across gap from source dot 111024dcr*/
   pointSet *srcDots = NULL;
   char ptmaster = ' '; /*pointmaster character  dcr041009*/
   int idx = 0; /*dcr041020*/
   atom *parent = NULL; /*dcr20120120*/
   /*int LXHvector = TRUE;*/  /*dcr20120120  eventually commandline param*/
   float XHdistance = 99.9;  /*dcr20120120*/
   float XHTangle = 0.0;  /*dcr20120120*/

   if (src->elem == ignoreAtom)
      { return; }
   /*050119 but there seems no way to ever set an atom->elem = ignoreAtom !?*/

   if (OccupancyCriteria && src->occ < OccupancyCutoff) /*Occ.Crit. 060831*/
      {/*050119*/ return; }

   cause = NULL; /* who is responsible for the dot */

   srcDots = &(dots[src->elem]);

   if (src->elem == atomC && isCarbonylAtom(src))
      {srcDots = &COdots;}

   maskHB = 0;
   if(src->props & DONOR_PROP)
   {
      maskHB = ACCEPTOR_PROP;
      if(src->props & ACCEPTOR_PROP)
      { /* ambig */
         maskHB |= DONOR_PROP;
      }
   }
   else if (src->props & ACCEPTOR_PROP)
   {
      maskHB = DONOR_PROP;
   }

   /*111013 original mode: for each src atom (each call to this subroutine) */
   /*   place dot ball at src, */
   /*   loop over all dots, for each dot loop over neighbors*/

   /*OneDotEach mode: for each src atom (each call to this subroutine) */
   /*   Outer loop over all neighbors, each in turn taken to cause a dot*/
   /*   for each neighbor, calculate dot at src surface in line with neighbor */
   /*   pass dot to Inner loop for this single dot,*/
   /*   then loop over neighbors to see if any occlude the line-of-sight src--cause */
   /*   I.E. wrap orig "Inner" dot loop with extra "Outer" neighbor loop, */
   /*   alternate modes do only one pass at one or the other loop. */


   if(Ldotdump)
       fprintf(stderr,"\n\n******* SOURCE Atom at src:%7.3f,%7.3f,%7.3f, name %s %s %s SOURCE\n",src->loc.x,src->loc.y,src->loc.z,src->atomname,src->r->resname,src->r->Hy36resno);
/* */
   /*{{{OUTER LOOP over target atoms*/
   for(extratarg = scratch; extratarg; extratarg = extratarg->scratch) /*OUTER LOOP*/
   {/* OUTER LOOP over target atoms */
      /* scan over neighbors and calc in-sight-line dot position for each */

       if(Ldotdump)
           fprintf(stderr,"[src %s %s %s] TARGET:%7.3f,%7.3f,%7.3f, name %s %s %s TARGET\n",src->atomname,src->r->resname,src->r->Hy36resno,extratarg->loc.x,extratarg->loc.y,extratarg->loc.z,extratarg->atomname,extratarg->r->resname,extratarg->r->Hy36resno);
/* */

      /*operate on a single dot from outer loop */
      /*first get src to target vector, call this working dot vector*/
      v3sub(&(extratarg->loc), &(src->loc), &dotvect);
      /*then scale to src vdW to get src ctr to close point vector*/
      v3scale(&dotvect, src->radius); /*the actual dot vector*/
      /*add this to src position to get close point position in 3 space*/
      v3add(&(src->loc), &dotvect, &dotloc);/*dotlocation on line to target*/
      /*save actual dotvector, scaled later to get dot on other sde of gap*/
      gapvect.x = dotvect.x;
      gapvect.y = dotvect.y;
      gapvect.z = dotvect.z;
      /*scale this out by probeRadius to get exploc */
      v3scale(&dotvect, src->radius + probeRad);
      v3add(&(src->loc), &dotvect, &exploc); /*probeRad out from dotloc*/

       if(Ldotdump)
           fprintf(stderr,"[src %s %s %s] DOT [%s %s %s] dotloc:%7.3f,%7.3f,%7.3f DOT\n",src->atomname,src->r->resname,src->r->Hy36resno,extratarg->atomname,extratarg->r->resname,extratarg->r->Hy36resno,dotloc.x,dotloc.y,dotloc.z);
/* */

      mingap = 999.9;
      isaHB = tooCloseHB = FALSE;
      hbondbumpgap = 0.0;
      ok = FALSE; /*expect dot to be rejected*/

      /*{{{attributes of dot wrt causal target: sets ok */

         /*LOneDotEach only considers target atom that caused this dot*/
         /* to assign attributes that might filter out this dot */
         /*so at the end of this section, we will break out and then */
         /*if the dot is of the right sort so far, we will scan over */
         /*neighbors to see if this dot is occluded. */
         targ = extratarg; /*override: evaluate causual target for this dot*/

         hbcutoff = Min_regular_hb_cutoff; /* redefined later on */
         targetloc = targ->loc;
         targVDW = targ->radius;
#ifdef INLINE_FOR_SPEED
         nearpt = INRANGEINLINE(exploc, targetloc, probeRad+targVDW);
#else
         nearpt = inRange(&exploc, &targetloc, probeRad+targVDW);
#endif
         /*returns true when distance-sq between a and b <= limit-sq */
         /*exploc is location of probe ctr with probe surface at dotloc*/
         /*Fails erratically when atoms superimpose so distance back from a */
         /*probe distance beyond surface of src*/
         /*is a full probe radius + vdW to the center of target atom */
         /*So full overlap disaster condition needs to be tested*/
         /*before filter on nearpt value. */


      if(Ldotdump)
      {
          fprintf(stderr,"nearpt=1  if within reach of probe\n");
          fprintf(stderr,"         targ->mark>=1 if target within BONDED set of src\n");
          fprintf(stderr,"nearpt %d, targ->mark %d [src: %s %s %s] %8.3f,%8.3f,%8.3f ; [targ %s %s %s] %8.3f,%8.3f,%8.3f ; exp:%8.3f,%8.3f,%8.3f\n",nearpt,targ->mark,src->atomname,src->r->resname,src->r->Hy36resno,src->loc.x,src->loc.y,src->loc.z,targ->atomname,targ->r->resname,targ->r->Hy36resno,targetloc.x,targetloc.y,targetloc.z,exploc.x,exploc.y,exploc.z);
      }
      
      skipthisMCMC = 0;

         /*just checking if this target atom acceptable for current src dot*/

         if(    DoMcMc == FALSE
             && ( src->props & MC_PROP)
             && (targ->props & MC_PROP) )
         {/*potential MC-MC for us to skip this target atom*/
            if (! (   ( src->props & HET_PROP)
                   || (targ->props & HET_PROP) ) )
            {/* neither can be a HET */
               if (strcmp(src->r->chain, targ->r->chain) == 0)
               {/* must be the same chain*/ /*problem for NMR models ????*/
                  skipthisMCMC = 1;
               }
            }
         }

         /*now enter if...if else... chain of decisions...*/
         /*drop through defaults to ok=0, dot is REJECTED*/

         if (targ->elem == ignoreAtom) {/*NOP: condition not contribute to decision*/}

         else if (OccupancyCriteria && targ->occ < OccupancyCutoff)
            {/*050119 will not contribute */}   /*Occ.Crit. 060831*/

         else if (skipthisMCMC) {/*NOP: condition not contribute to decision*/}

         else if(    DoWatWat == FALSE
                  && ( src->props & WATER_PROP)
                  && (targ->props & WATER_PROP))
         {
            ; /* water will not contribute to decision when not being considered*/
         }

         if(targ->mark == 0)
         {/*not flagged as a bonded atom*/
                 /*20111207dcr: unfortunately, can't distinguish between atoms*/
                 /* that by element type might be bonding and those in separate*/
                 /* groups that really are too close*/
           /*overlap disaster: when src and target(==cause) atoms on top of each other*/
           /*tests done from the edge of src to ctr of target are ill behaved here*/
           /*double v3distance(point3d* a, point3d* b) */
           /*{return(sqrt(v3distanceSq(a, b)));}*/
           srcCauseDist = v3distance(&(extratarg->loc), &(src->loc));
           if( (srcCauseDist <= (src->radius))||(srcCauseDist <= (extratarg->radius)))
           {
              nearpt = 2; /*disasterously near*/
              ok = TRUE; /*really want to see this clash*/
              cause = targ;
              /*targ gets reused, use cause as identifier of extratarg atom*/
               if(Ldotdump) fprintf(stderr,"nearpt==2, atoms overlap!\n");
           }
         }
         /*do next as plain if, rather than else if 20111207 dcr*/
         if((nearpt==1) && (targ->mark == 0)) /*filter on nearpt value*/
         {/* close but nonbonded, so might be an OK dot...*/
             if(Ldotdump) fprintf(stderr,"close but nonbonded: calc gap...\n");
             gap = gapSize(&dotloc, &targetloc, targVDW);
/*111013*/
             if(Ldotdump)
                 fprintf(stderr,"[src %s %s %s]  gap== %8.3f,mingap== %8.3f  [etarg %s %s %s]\n",src->atomname,src->r->resname,src->r->Hy36resno,gap,mingap,targ->atomname,targ->r->resname,targ->r->Hy36resno);

            /*identify target with minimum gap*/
            /*set properties*/
            if (gap < mingap) /*otherwise ok wll remain FALSE*/
            {/*if: gap < mingap*/

               int bothCharged = FALSE;
               int chargeComplement = FALSE;

               if(  ( src->props & (NEGATIVE_PROP|POSITIVE_PROP))
                  &&(targ->props & (NEGATIVE_PROP|POSITIVE_PROP)) )
               {
                  bothCharged = TRUE;
                  chargeComplement =
                     (   (   ( src->props & POSITIVE_PROP)
                          && (targ->props & NEGATIVE_PROP))
                      || (   ( src->props & NEGATIVE_PROP)
                          && (targ->props & POSITIVE_PROP)) );
               }

               if(    (targ->props & maskHB) /*test for H-BOND*/
                   && ( (! bothCharged) || chargeComplement) )
               {
                  hbcutoff = bothCharged ?
                        Min_charged_hb_cutoff :
                        Min_regular_hb_cutoff;

                  if (gap < -hbcutoff)
                  {
                     tooCloseHB = TRUE; /* treat like a bump */
                     hbondbumpgap = gap;
                     isaHB = TRUE;
                  }
                  else
                  { /* hbond or contact */
                     isaHB = TRUE;
                     tooCloseHB = FALSE;
                  }
               }
               else
               {/* clash or contact, not an H-BOND */
                  /*NOTE: atomHOd : hb-only-dummy, phantom H atom*/
                  if (src->elem == atomHOd)  { continue; } /*with for: loop*/
                  if (targ->elem == atomHOd) { continue; } /*with for: loop*/
                  isaHB = FALSE;
               }
               mingap = gap; /* update minimum gap */
               ok = TRUE;
               cause = targ; /*targ is loop counter, about to go out of scope*/
            }/*if: gap < mingap*/
         }/*else if: close but nonbonded*/

         if(nearpt == 2) /*disasterously near, i.e. on top of each other*/
         {
            ok = TRUE;
            if(Ldotdump) fprintf(stderr,"OK by courtesy, nearpt==2, keep horrible clash\n");
         }

         if(ok)
         {
             if(Ldotdump) fprintf(stderr,"OK so far ...  (targ->mark== %d)\n",targ->mark);
         }
         else
         {
             if(Ldotdump) fprintf(stderr,"NOT OK: REJECT (targ->mark== %d)\n",targ->mark);
         }

      /*}}}attributes of dot wrt causal target: sets ok */

       if(Ldotdump) fprintf(stderr,"Mark1 ok== %d, nearpt== %d\n",ok,nearpt);

        /*cause is atom identified as extratarg which invokes the dot on src atom*/
        /*targ is current atom being considered, overloaded name: reused below */
        /*as putative occluding 3rd atom from the surrounding neighbor atoms*/
        /*Does OneDotEach always identify cause atom as a targFlg atom, need to check??*/
        /*BUT probe has several ways to define sets of atoms, some as src, some as targ*/

      if (ok && !(cause->flags & targFlg))
      { /* drop non-target atom dot */
         ok = FALSE;
      }

      if (ok && (nearpt==1))
      {/*ok so far!: check for occluding atoms */
          if(Ldotdump) fprintf(stderr,"DOT OK so far, check for OCCLUDING neighbors\n");

         /*{{{scan over putative occluding neighbors*/

         /*targ is an atom*, scratch is an atom*,  atom* defined in abin.h */
         /* atom* scratch is a member of atom*: which has member scratch...*/
         /* so for-loop moves through all these atoms while an atom has a scratch*/

         for(targ = scratch; targ; targ = targ->scratch)
         {/*for: scan over target atoms*/

            if(LOneDotEach  && (targ == extratarg))continue;/*to interate loop*/
            /*Looking from src atom toward extratarg, check neighbors for occlusion*/
/*111013*/
            if(Ldotdump)
                 fprintf(stderr,"----------OCCLUSION test: [src %s %s %s]dot vs [neighbor %s %s %s elem==%d]\n",src->atomname,src->r->resname,src->r->Hy36resno,targ->atomname,targ->r->resname, targ->r->Hy36resno, targ->elem); /*elem 130427dcr*/

            /* "targ" atoms here are neighbors potentially occluding*/

            if(Lhonlyocclude && !isHatom(targ->elem) )continue;/*to interate loop 130427dcr*/

            if (targ->mark && (targ->elem != atomHOd) &&!LOneDotEach)/*111122 ODE test*/
            /*if (targ->mark && (targ->elem != atomHOd))*/ /*111122 ODE test*/
            {/*bonded atom*/
               /*NOTE: atomHOd : hb-only-dummy, phantom H atom*/
#ifdef INLINE_FOR_SPEED
               within = INRANGEINLINE(dotloc, (targ->loc), targ->radius);
#else            /*returns true when distance-sq between a and b <= limit-sq */
               within = inRange(&dotloc, &(targ->loc), targ->radius);
#endif
               if (within)
               {
if(Ldotdump)
fprintf(stderr," NEIGHBOR [%s %s %s] OCCLUDS    src---targ interaction\n",targ->atomname, targ->r->resname,targ->r->Hy36resno);
                  ok = FALSE;
                  break;
               }
               else
               {
if(Ldotdump)
fprintf(stderr," NEIGHBOR [%s %s %s] BENIGN for src---targ interaction\n",targ->atomname, targ->r->resname,targ->r->Hy36resno);
               }
            }/*bonded atom*/
            else if ((src->elem == atomHOd) && !(targ->props & ACCEPTOR_PROP))
            {/*dummy HYDROGEN*/   /* atomHOd : hb-only-dummy, phantom H atom*/
               /*eliminate if H? within non-Acceptor atom! WIERD*/
               /*111122dcr maybe this is for dummy H on water: only make H-bonds????*/
#ifdef INLINE_FOR_SPEED
               within = INRANGEINLINE(dotloc, (targ->loc), targ->radius);
#else            /*returns true when distance-sq between a and b <= limit-sq */
               within = inRange(&dotloc, &(targ->loc), targ->radius);
#endif
               if (within)
               {
                   if(Ldotdump)
                       fprintf(stderr,"WIERD:HYDROGEN?[%s %s %s] OCCLUDS src---target\n",targ->atomname,targ->r->resname,targ->r->Hy36resno);
                  ok = FALSE;
                  break;
               }
               else
               {
                   if(Ldotdump)
                       fprintf(stderr,"WIERD:HYDROGEN?[%s %s %s] BENIGN  src---target\n",targ->atomname,targ->r->resname,targ->r->Hy36resno);
               }
            }/*dummy HYDROGEN*/   /* atomHOd : hb-only-dummy, phantom H atom*/
            else if(LOneDotEach)
            {/*LOneDotEach*/
              {/*project putative occluding atom onto src--target line*/
                 /*drop perpendicular from putative occluder to src--target line*/
                 /*NOTE: one atom is src, the other is extratarg == target causes dot*/
                 /*putative occluding atom == "targ", i.e. scan all targ for each src*/
                  if(Ldotdump)
                      fprintf(stderr,"check PERPENDICULAR distance of occluding atom to src---targ line\n");
                 /*get position on line defined by atoms where 3rd atom is perpendicular*/
                 perpendicular(&(src->loc),&(extratarg->loc),&(targ->loc),&perploc);

                 /*now check distance of putative occluder (targ) to perploc*/
#ifdef INLINE_FOR_SPEED
                 within = INRANGEINLINE(perploc, (targ->loc), targ->radius);
#else
                 within = inRange(&perploc, &(targ->loc), targ->radius);
#endif
                 if(within)
                 {/*occluder within vdW of line defined by atoms (src and extratarg)*/
                     if(Ldotdump)
                         fprintf(stderr,"occluding atom within it's own vdW of src---targ line\n");

                   /*TWO conditions can cause occlusion: */
                   /* 1:block gap */
                    /*close enough to occlude interaction if perploc within the gap*/
                    /*src loc -- perp loc -- target-causing-dot loc */
                   /* 2:eat dots */
                    /*3rd atom overlaps either one of the gap defining points*/
                    /* which can happen even if perploc not between atoms*/

                     /*TWO possible occluding situations: block gap or eat dots*/
                     /*need to know location of gap dot on target atom side*/
                     v3scale(&gapvect, src->radius + gap);
                     v3add(&(src->loc), &gapvect, &gaploc);

                     if(  (  v3distance(&gaploc,&perploc) < fabs(gap) )
                        &&(  v3distance(&dotloc,&perploc) < fabs(gap) ) )
                     {/*perp hits within gap so does occlude line of sight*/
                         if(Ldotdump)
                             fprintf(stderr,"atom [%s %s %s] OCCLUDES src---targ line at %7.3f %7.3f %7.3f, DOT REJECTED\n" ,targ->atomname,targ->r->resname,targ->r->Hy36resno,perploc.x,perploc.y,perploc.z);
                       ok = FALSE;
                       break;
                     }/*perp hits within gap so does occlude line of sight*/

                     /*check if neighbor close enough to occlude gapdots*/
                     /*need to know if perploc actually between atom centers*/
                     /*check that perpendicular drop point is between the two atoms*/
                     if(  (v3distance(&perploc,&(src->loc))
                           < v3distance(&(src->loc),&(extratarg->loc)))
                        &&(v3distance(&perploc,&(extratarg->loc))
                           < v3distance(&(src->loc),&(extratarg->loc))))
                     {/*perp hits line inbetween atoms*/
                        Lbetweenatoms = TRUE;
                     }/*perp hits line inbetween atoms*/
                     else {Lbetweenatoms = FALSE;}
                     /*now see if 3rd atom eats dot*/
                     if(  (  (v3distance(&dotloc,&(targ->loc)) < targ->radius)
                                ||(v3distance(&gaploc,&(targ->loc)) < targ->radius) )
                             &&Lbetweenatoms)
                     {/*neighbor vdW engulfs gap dot(s),i.e. OCCLUDES*/
                         if(Ldotdump)
                             fprintf(stderr,"atom [%s %s %s] OCCLUDES gap dots,distance: %7.3f & %7.3f,vs vdW %7.3f: DOT REJECTED\n" ,targ->atomname,targ->r->resname,targ->r->Hy36resno, v3distance(&dotloc,&(targ->loc)),v3distance(&gaploc,&(targ->loc)),targ->radius );
                          ok = FALSE;
                          break;
                     }/*neighbor vdW engulfs gap dot(s),i.e. OCCLUDES*/
                     else
                     {
                         if(Ldotdump)
                             fprintf(stderr,"atom [%s %s %s] NOT occlude dots,distance: %7.3f & %7.3f,vs vdW %7.3f: DOT OK\n" ,targ->atomname,targ->r->resname,targ->r->Hy36resno, v3distance(&dotloc,&(targ->loc)),v3distance(&gaploc,&(targ->loc)),targ->radius );
                     }
                     /*check if neighbor close enough to occlude gapdots*/

                 }/*occluder within vdW of line defined by atoms (src and extratarg)*/
                 else
                 {
                     if(Ldotdump)
                         fprintf(stderr,"this neighbor atom too far away to occlude src---targ line\n");
                 }
              }/*project putative occluding atom onto src--target line*/
            }/*LOneDotEach*/

         }/*for: scan over target atoms*/

         /*}}}scan over putative occluding neighbors*/

      }/*ok so far?: check for occluding atoms */

      if(Ldotdump) fprintf(stderr,"OK2 == %d, DOT to be outputed if ok\n",ok);

      if (ok)
      {/*dot is ok, atom that this dot hits is called cause */ /* Contact */

         /* ovrlaptype : -1 bump, 0 touch, +1 H bond */
         /* mingap originally minimum distance of src dot to cause atom*/

         if(nearpt == 2) /*disasterously near, i.e. on top of each other*/
         {
            sl = -((src->radius) +  (extratarg->radius))/2.0;
            mingap = 2.0*LowGoodCut; /*force arbitrary bad clash*/
            ovrlaptype =-1; /*disaster overlap == clash!! */
            /*spikeloc calc from src->radius + sl, sl=spikelen*mingap,*/
            /*spikelen = 0.50 OR input param*/
            /*usually: sl = spikelen*mingap*/
         }
         else if (mingap > 0.0)
         {
            sl = 0.0;
            if(LweakHbonds && isaHB) {ovrlaptype = 2;} /*20111215dcr*/
             /* TODO SJ: Why is this 2*/
            else {ovrlaptype = 0;}
         }
         else if (isaHB && tooCloseHB)
         {
            mingap = hbondbumpgap + hbcutoff;
            sl = spikelen*mingap;
            ovrlaptype =-1; /* Hbond which is too close is a clash */
         }
         else if (isaHB)                           /* Hbond */
         {
            sl = spikelen*mingap;
            ovrlaptype =+1;

            /* must test angle for Hbond in some cases */
            if (        src->props & TEST_ACCEPT_ANGLE_PROP)
            {
               /* cause->loc vs src->aromaticRing; ovrlaptype =-1;?? */
            }
            else if ( cause->props & TEST_ACCEPT_ANGLE_PROP)
            {
               /* src->loc vs cause->aromaticRing; ovrlaptype =-1;?? */
            }

            if(   src->props & CH_DONOR_PROP
               || cause->props & CH_DONOR_PROP)
            {/* CH..O type Hbond */
               sl *= CHOHBfactor; /*down scale score for these types of Hbonds*/
            }

#ifdef JACKnANDREA
            if (writeHbonds)
              printf("{%s %s %s %d}P %f %f %f {%s %s %s %d} %f %f %f\n",
                     src->atomname, src->r->resname, src->r->chain, src->r->resid, src->loc.x, src->loc.y, src->loc.z,
                     cause->atomname, cause->r->resname, cause->r->chain, cause->r->resid, cause->loc.x, cause->loc.y, cause->loc.z);
/*this makes a monster list of all possible H-bond dots as heavy-atom--H vecs*/
/*which must then be pared down to the unique atom-pair vectors*/
/*also includes all the atomOHd dummy phantoms which junks up the display*/
#endif

         }
         else                                    /* Clash */
         {
            sl = spikelen*mingap;
            ovrlaptype =-1;
         }


          if(Ldotdump)
              fprintf(stderr,"[src %s %s %s] OK2 passed, ovrlaptype == %d********************\n",src->atomname,src->r->resname,src->r->Hy36resno,ovrlaptype);
/* */
         /*NOTE: atomHOd : hb-only-dummy, phantom H atom*/
         /* ovrlaptype : -1 bump, 0 touch, +1 H bond */
         if((ovrlaptype == 1)||(src->elem != atomHOd && cause->elem != atomHOd))
         { /* TODO SJ: Figure out what will happen is weak H bonds is true*/
           v3scale(&dotvect, src->radius + sl);
           v3add(&(src->loc), &dotvect, &spikeloc);

           /* possibly limit contact dots... */
           /* ovrlaptype : -1 bump, 0 touch, +1 H bond */

      /*kissEdge2bullsEye seems not responsible for loss of very close clashes*/
           if (   ovrlaptype != 0
               || LimitDots == FALSE
               || dot2srcCenter(&dotloc, src, cause) <=
                  kissEdge2bullsEye(src->radius, cause->radius, probeRad))
           {/*could use this dot/spike, */
             /*saveDot decisions moved here to examineDots 041020 */
             /*overlaptype:  -1 bump, 0 touch, +1 H bond */
             idx = dotClassIndex(ovrlaptype, mingap);
               /* OLD: idx: 0 wide contact,1 close contact,2 small overlap,3 bad overlap,4 H-bonds*/
               /* New SJ 04/21/2015: idx: 0 wide contact,1 close contact,2 weak H-bonds, 3 small overlap,4 bad overlap, 5 worse overlaps, 6 H-bonds*/
             /*OnlyBadOut option dcr041010*/

            /* ovrlaptype : -1 bump, 0 touch, +1 H bond */
            if(  (!OnlyBadOut && OutputHBs && ovrlaptype > 0)
               ||(!OnlyBadOut && OutputClashes && ovrlaptype  < 0)
               ||( OnlyBadOut && OutputClashes && /*idx==3*/ (idx ==4 || idx == 5))/*Bad Clash dcr041020*/ /* SJ 04/21/2015 changed idx according to new values*/
               ||(!OnlyBadOut && OutputVDWs && ovrlaptype == 0))
             { /*logicals allow the overlaptype of this dot*/

               /*now derive ptmaster*/
               if( (src->props & MC_PROP)  && (cause->props & MC_PROP) )
               {/*potential MC-MC to flag*/
                  if( !(( src->props & HET_PROP) || (cause->props & HET_PROP)) )
                  {/*neither can be a HET, what wierd condition is avoided?*/
                     ptmaster = 'M'; /*McMc interaction*/
                     mcmccont[0][idx]++;
                     //mcmccont[0][5]++;
                     mcmccont[0][7]++; /* SJ 04/21/2015 changing the second dimension to 7 to match new initialization*/
                  }
               }
               else if( (src->props & SC_PROP)  && (cause->props & SC_PROP) )
               {/*potential SC-SC to flag*/
                  if( !(( src->props & HET_PROP) || (cause->props & HET_PROP)) )
                  {/*neither can be a HET, what wierd condition is avoided?*/
                     ptmaster = 'S'; /*ScSc interaction*/
                     scsccont[0][idx]++;
                     //scsccont[0][5]++;
                     scsccont[0][7]++;/* SJ 04/21/2015 changing the second dimension to 7 to match new initialization*/
                  }
               }
               else if(  ((src->props & SC_PROP)  && (cause->props & MC_PROP))
                       ||((src->props & MC_PROP)  && (cause->props & SC_PROP)))
               {/*potential MC-SC or SC-MC to flag*/
                  if( !(( src->props & HET_PROP) || (cause->props & HET_PROP)) )
                  {/*neither can be a HET, what wierd condition is avoided?*/
                     ptmaster = 'P'; /*McSc or ScMc interaction*/
                     mcsccont[0][idx]++;
                     //mcsccont[0][5]++;
                     mcsccont[0][7]++;/* SJ 04/15/2015 changing the second dimension to 7 to match new initialization*/
                  }
               }
               else
               {
                     ptmaster = 'O'; /*Oh for Other, general default*/
                     othrcont[0][idx]++;
                     //othrcont[0][5]++;
                     othrcont[0][7]++;/* SJ 04/15/2015 changing the second dimension to 7 to match new initialization*/
               }
/*111013*/
                 if(Ldotdump)
                     fprintf(stderr,"[src %s %s %s] :idx==%d: [neigh %s %s %s]*** saveDot:: gap %8.3f, mingap %8.3f, holdmingap %8.3f, dotloc: %7.3f,%7.3f,%7.3f; spikeloc: %7.3f,%7.3f,%7.3f \n",src->atomname,src->r->resname,src->r->Hy36resno,idx,cause->atomname, cause->r->resname,cause->r->Hy36resno,gap,mingap,holdmingap, dotloc.x,dotloc.y,dotloc.z,spikeloc.x,spikeloc.y,spikeloc.z);
               if(Lnearest)
               {
                 if(mingap < holdmingap)
                 {
                    holdsrc = src;
                    holdcause = cause;
                    holdtype = type;
                    holddotloc.x = dotloc.x;
                    holddotloc.y = dotloc.y;
                    holddotloc.z = dotloc.z;
                    holdspikeloc.x = spikeloc.x;
                    holdspikeloc.y = spikeloc.y;
                    holdspikeloc.z = spikeloc.z;
                    holdovrlaptype = ovrlaptype;
                    holdmingap = mingap;
                    holdptmaster = ptmaster;
                    Lholdingdot = TRUE;
                    if(Ldotdump)
                         fprintf(stderr,"reset holdings: holdmingap %8.3f\n",holdmingap);
                 }
               }
               else
               {
                   if(Ldotdump)
                       fprintf(stderr,"saveDot CALL+++++++++++++++++++++++\n");
                  saveDot(src, cause, type, &dotloc, &spikeloc,
                          results, ovrlaptype, mingap, ptmaster,XHTangle); /*dcr041009*/
               }
             }/*logicals allow the overlaptype of this dot*/
           }/*could use this dot/spike, */
         }
      }/*dot is ok, atom that this dot hits is called cause */
      /*operate on a single dot from outer loop */

/*111013*/
       if(Ldotdump)
           fprintf(stderr,"+++++++++++END of pass: loop on extratarg, finished with %s %s %s\n\n",extratarg->atomname,extratarg->r->resname,extratarg->r->Hy36resno);
/* */
   }/* OUTER LOOP over target atoms */
   /*}}}OUTER LOOP over target atoms*/

   if(Lnearest && Lholdingdot)
   {/*Lnearest mode and has an OK dot which is nearest that is allowed by targPat*/
      ok = TRUE;
      if(Lfilternearest) /*behavior if extraPat==NULL ???? */
      {
         if(!matchPat(holdcause, extraPat)) /*use 3rd Pattern 111022dcr*/
         {/*src deemed OK before call examineDots, unless cause OK, ignore atom*/
            ok = FALSE;
         }
      }
      if(ok && Lfilterocc1) /*final filter exclude any with occ < 1  20120120dcr*/
      {
         if(holdsrc->occ < 1 || holdcause->occ < 1)
         {
            ok = FALSE;
            fprintf(stderr,"EXCLUDE[src %s OCC:%4.2f, %s %s] :: [cause %s OCC:%4.2f, %s %s]\n",holdsrc->atomname,holdsrc->occ,holdsrc->r->resname,holdsrc->r->Hy36resno
                    ,holdcause->atomname,holdcause->occ,holdcause->r->resname,holdcause->r->Hy36resno);
         }
      }
      if(ok && Lfilternoalt) /*final filter exclude any alts  20120120dcr*/
      {
         if(holdsrc->altConf != ' ' || holdcause->altConf != ' ')
         {
            ok = FALSE;
            fprintf(stderr,"EXCLUDE[src %s ALT:%c, %s %s] :: [cause %s ALT:%c, %s %s]\n"
                    ,holdsrc->atomname,holdsrc->altConf,holdsrc->r->resname,holdsrc->r->Hy36resno
                    ,holdcause->atomname,holdcause->altConf,holdcause->r->resname,holdcause->r->Hy36resno);
         }
      }
      if(ok)
      {
          if(Ldotdump)
              fprintf(stderr,"[src %s %s %s] :: [neigh %s %s %s]*YES**about to saveDot:: mingap %8.3f HOLD\n",holdsrc->atomname,holdsrc->r->resname,holdsrc->r->Hy36resno,holdcause->atomname,holdcause->r->resname,holdcause->r->Hy36resno,holdmingap);

         if(LXHvector && holdsrc->elem == atomH)
         {
            for(parent = allMainAtoms; parent; parent = parent->next) /*dcr20120120*/
            {/*for: loop over all main atoms taking each in turn as the parent-atom*/
               XHdistance = v3distance(&(parent->loc), &(src->loc));
               if(XHdistance < 1.2 && XHdistance > 0.1 )
               {/*parent atom for this Hydrogen*/
                  XHTangle = v3angle(&(parent->loc), &(src->loc), &(holdcause->loc));
                  XHTangle = XHTangle*360/(2*3.14159265);

                   if(Ldotdump)
                       fprintf(stderr,"[src %s %s %s] :: [parent %s %s %s]*YES**about to saveDot:: mingap %8.3f, XHdistance: %8.3f, XHTangle: %8.3f HOLD\n",holdsrc->atomname,holdsrc->r->resname,holdsrc->r->Hy36resno,parent->atomname,parent->r->resname,parent->r->Hy36resno,holdmingap,XHdistance,XHTangle);

                 /*save the XHvector and XHTangle construction*/
                /* saveDot(holdsrc, parent, holdtype, &parent->loc, &holdsrc->loc,
                    results, -1, holdmingap, holdptmaster, XHTangle);*/
               }
            }
         }
         else {XHTangle = 0.0;}
         /*now save the real dot*/

         saveDot(holdsrc, holdcause, holdtype, &holddotloc, &holdspikeloc,
            results, holdovrlaptype, holdmingap, holdptmaster,XHTangle);

      }
   }/*Lnearest mode and has an OK dot*/
}/*examineOneDotEach()*/
/*}}}examineOneDotEach()_____________________________________________________*/

/*{{{examineDots()*********** called from genDotIntersect() ******************/
void examineDots(atom *src, int type, atom *scratch,
		pointSet dots[], float probeRad,
		float spikelen, int targFlg, dotNode *results[][NODEWIDTH])
{
   /*called from genDotIntersect() with same names except scratch==atomList2 */

   atom *targ = NULL, *cause = NULL;
   int i = 0, nearpt = 0, within = 0, ok = 0, maskHB = 0;
   int skipthisMCMC = 0, isaHB = 0, tooCloseHB = 0, ovrlaptype = 0;
   float gap = 0.0, mingap = 0.0, sl = 0.0, hbondbumpgap = 0.0;
   float hbcutoff = 0.0, targVDW = 0.0;
   point3d targetloc, dotloc, dotvect, exploc, spikeloc;
   pointSet *srcDots = NULL;
   char ptmaster = ' '; /*pointmaster character  dcr041009*/
   int idx = 0; /*dcr041020*/

   if (src->elem == ignoreAtom) { return; }
   /*050119 but there seems no way to ever set an atom->elem = ignoreAtom !?*/

   if (OccupancyCriteria && src->occ < OccupancyCutoff) /*Occ.Crit. 060831*/
      {/*050119*/ return; }

   cause = NULL; /* who is responsible for the dot */

   srcDots = &(dots[src->elem]);

   if (src->elem == atomC && isCarbonylAtom(src)) {
      srcDots = &COdots;
   }

   maskHB = 0;
   if(src->props & DONOR_PROP) {
      maskHB = ACCEPTOR_PROP;
      if(src->props & ACCEPTOR_PROP) { /* ambig */
            maskHB |= DONOR_PROP;
      }
   }
   else if (src->props & ACCEPTOR_PROP) {
      maskHB = DONOR_PROP;
   }

   for(i=0; i < srcDots->n; i++)
   {/*loop over dots on src dot-ball */
      dotvect = srcDots->p[i];

      v3add(&(src->loc), &dotvect, &dotloc);

      v3scale(&dotvect, src->radius + probeRad);
      v3add(&(src->loc), &dotvect, &exploc);

      ok = FALSE;
      mingap = 999.9;
      isaHB = tooCloseHB = FALSE;
      hbondbumpgap = 0.0;

      /*targ is an atom*, scratch is an atom*,  atom* defined in abin.h */
      /* atom* scratch is a member of atom*: which has member scratch...*/
      /* so for-loop moves through all these atoms while an atom has a scratch*/

      for(targ = scratch; targ; targ = targ->scratch)
      {/*for: loop over target atoms: set ok=1 for some subset*/
         hbcutoff = Min_regular_hb_cutoff; /* redefined later on */

         targetloc = targ->loc;
         targVDW = targ->radius;

#ifdef INLINE_FOR_SPEED
         nearpt = INRANGEINLINE(exploc, targetloc, probeRad+targVDW);
#else
         nearpt = inRange(&exploc, &targetloc, probeRad+targVDW);
#endif
         /*returns true when distance-sq between a and b <= limit-sq */
         /*...in subroutine examineDots()...*/

         skipthisMCMC = 0;

         /*just checking if this target atom acceptable for current src dot*/
         /*(premature to assign pointmasters when target not yet determined*/
         if(    DoMcMc == FALSE
             && ( src->props & MC_PROP)
             && (targ->props & MC_PROP) )
         {
            /*potential MC-MC for us to skip this target atom*/
            if (! (   ( src->props & HET_PROP)
                   || (targ->props & HET_PROP) ) )
            {/* neither can be a HET */
               if (strcmp(src->r->chain, targ->r->chain) == 0)
               {/* must be the same chain*/ /*problem for NMR models ????*/
        	          skipthisMCMC = 1;
               }
            }
         }

         /*now enter if...if else... chain of decisions...*/

         if (targ->elem == ignoreAtom) {/* will not contribute */}

         else if (OccupancyCriteria && targ->occ < OccupancyCutoff)
            {/*050119 will not contribute */}   /*Occ.Crit. 060831*/

         else if (skipthisMCMC) {/* will not contribute */}

         else if(    DoWatWat == FALSE
                  && ( src->props & WATER_PROP)
                  && (targ->props & WATER_PROP))
         {
            ; /* water will not contribute when not being considered*/
         }

         /*Lautobondrot 050111 investigation point */

         else if (nearpt && (targ->mark == 0))
         {/*else if: close but nonbonded*/
            gap = gapSize(&dotloc, &targetloc, targVDW);

            /*identify target with minimum gap*/

            if (gap < mingap)
            {/*if: gap < mingap*/

               int bothCharged = FALSE;
               int chargeComplement = FALSE;

               if(  ( src->props & (NEGATIVE_PROP|POSITIVE_PROP))
                  &&(targ->props & (NEGATIVE_PROP|POSITIVE_PROP)) )
               {
                  bothCharged = TRUE;
                  chargeComplement =
                     (   (   ( src->props & POSITIVE_PROP)
                          && (targ->props & NEGATIVE_PROP))
                      || (   ( src->props & NEGATIVE_PROP)
                          && (targ->props & POSITIVE_PROP)) );
               }

               if(    (targ->props & maskHB)
                   && ( (! bothCharged) || chargeComplement) )
               {
                  hbcutoff = bothCharged ?
                        Min_charged_hb_cutoff :
                        Min_regular_hb_cutoff;

                  if (gap < -hbcutoff)
                  {
                     tooCloseHB = TRUE; /* treat like a bump */
                     hbondbumpgap = gap;
                     isaHB = TRUE;
                  }
                  else
                  { /* hbond or contact */
                     isaHB = TRUE;
                     tooCloseHB = FALSE;
                  }
               }
               else
               { /* clash or contact */
                  /*NOTE: atomHOd : hb-only-dummy, phantom H atom*/
                  if (src->elem == atomHOd)  { continue; } /*with for: loop*/
                  if (targ->elem == atomHOd) { continue; } /*with for: loop*/
                  isaHB = FALSE;
               }
               mingap = gap; /* update minimum gap */
               ok = TRUE;
               cause = targ;
            }/*if: gap < mingap*/
         }/*else if: close but nonbonded*/
      }/*for: loop over target atoms: set ok=1 for some subset*/
      if (ok && !(cause->flags & targFlg)) { /* drop non-target atom dot */
         ok = FALSE;
      }
      if (ok)
      {/*ok: apply bonded atom dot filters */

      /*targ is an atom*, scratch is an atom*,  atom* defined in abin.h */
      /* atom* scratch is a member of atom*: which has member scratch...*/
      /* so for-loop moves through all these atoms while an atom has a scratch*/

         for(targ = scratch; targ; targ = targ->scratch)
         {/*for: scan over target atoms*/
            /* eliminate dot if within bonded atom! */
            if (targ->mark && (targ->elem != atomHOd))
            {  /*NOTE: atomHOd : hb-only-dummy, phantom H atom*/
#ifdef INLINE_FOR_SPEED
               within = INRANGEINLINE(dotloc, (targ->loc), targ->radius);
#else
               within = inRange(&dotloc, &(targ->loc), targ->radius);
#endif
               /*returns true when distance-sq between a and b <= limit-sq */
               /*...in subroutine examineDots()...*/

               if (within) { ok = FALSE; break; }
             }
             else if ((src->elem == atomHOd) && !(targ->props & ACCEPTOR_PROP))
             {   /*NOTE: atomHOd : hb-only-dummy, phantom H atom*/
                /* eliminate if H? within non-Acceptor atom! */
#ifdef INLINE_FOR_SPEED
               within = INRANGEINLINE(dotloc, (targ->loc), targ->radius);
#else
               within = inRange(&dotloc, &(targ->loc), targ->radius);
#endif
               /*returns true when distance-sq between a and b <= limit-sq */
               /*...in subroutine examineDots()...*/

               if (within) { ok = FALSE; break; }
            }
         }/*for: scan over target atoms*/
      }/*ok: apply bonded atom dot filters */
      if (ok)
      {/*dot is ok, atom that this dot hits is called cause */ /* Contact */
         /* ovrlaptype : -1 bump, 0 touch, +1 H bond */
         /*dcr?: mingap seems to do what ???? 050111*/
          /* 04/10/2015 - SJ, including weak H bonds in overlaptype +1, if LweakHbonds is true*/
         if (mingap > 0.0)
         {
            sl = 0.0; /* TODO: SJ not sure waht this sl is. Find out and change accordingly*/
            if(LweakHbonds && isaHB) /* 04/10/2015 ovrlaptype is H bond if detecting weak H bonds*/
                ovrlaptype = 1;
            else
                ovrlaptype = 0;
         }
         else if (isaHB && tooCloseHB)
         {
            mingap = hbondbumpgap + hbcutoff;
            sl = spikelen*mingap;
            ovrlaptype =-1; /* Hbond which is too close is a clash */
         }
         else if (isaHB)                           /* Hbond */
         {
            sl = spikelen*mingap;
            ovrlaptype =+1;

            /* must test angle for Hbond in some cases */
            if (        src->props & TEST_ACCEPT_ANGLE_PROP)
            {
               /* cause->loc vs src->aromaticRing; ovrlaptype =-1;?? */
            }
            else if ( cause->props & TEST_ACCEPT_ANGLE_PROP)
            {
               /* src->loc vs cause->aromaticRing; ovrlaptype =-1;?? */
            }

            if(   src->props & CH_DONOR_PROP
               || cause->props & CH_DONOR_PROP)
            {/* CH..O type Hbond */
               sl *= CHOHBfactor; /*down scale score for these types of Hbonds*/
            }

#ifdef JACKnANDREA
            if (writeHbonds)
              printf("{%s %s %s %d}P %f %f %f {%s %s %s %d} %f %f %f\n",
                     src->atomname, src->r->resname, src->r->chain, src->r->resid, src->loc.x, src->loc.y, src->loc.z,
                     cause->atomname, cause->r->resname, cause->r->chain, cause->r->resid, cause->loc.x, cause->loc.y, cause->loc.z);
            /*this makes a monster list of all possible H-bond dots as heavy-atom--H vecs*/
            /*which must then be pared down to the unique atom-pair vectors*/
            /*also includes all the atomOHd dummy phantoms which junks up the display*/
#endif

         }
         else                                    /* Clash */
         {
            sl = spikelen*mingap;
            ovrlaptype =-1;
         }

         /*NOTE: atomHOd : hb-only-dummy, phantom H atom*/
         /* ovrlaptype : -1 bump, 0 touch, +1 H bond */
         if(  (ovrlaptype == 1)/* TODO: SJ Have to figue out what will happen here is LweakHbonds is true*/
            ||(src->elem != atomHOd && cause->elem != atomHOd) )
         {
           v3scale(&dotvect, src->radius + sl);
           v3add(&(src->loc), &dotvect, &spikeloc);

           /* possibly limit contact dots... */
           /* ovrlaptype : -1 bump, 0 touch, +1 H bond */

           /*kissEdge2bullsEye seems not responsible for loss of very close clashes*/
           if (   ovrlaptype != 0
               || LimitDots == FALSE
               || dot2srcCenter(&dotloc, src, cause) <=
                  kissEdge2bullsEye(src->radius, cause->radius, probeRad))
           {/*could use this dot/spike, */
             /*saveDot decisions moved here to examineDots 041020 */
             /*overlaptype:  -1 bump, 0 touch, +1 H bond */
             idx = dotClassIndex(ovrlaptype, mingap);
             /*OLD idx: 0 wide contact,1 close contact,2 small overlap,3 bad overlap,4 H-bonds*/
               /* new 04/10/2015 SJ: 0 wide contact, 1 close contact, 2 weak H bonds, 3 small overlap, 4 bad overlap, 5 worse overlap, 6 H-bonds*/
             /*OnlyBadOut option dcr041010*/

            /* ovrlaptype : -1 bump, 0 touch, +1 H bond */
                if(  (!OnlyBadOut && OutputHBs && ovrlaptype > 0)
                     ||(!OnlyBadOut && OutputClashes && ovrlaptype  < 0)
                     ||( OnlyBadOut && OutputClashes && /*idx==3 SJ*/ (idx == 4 || idx == 5))/*Bad Clash dcr041020*/ /* idx value changes according to new above, SJ*/
                     ||(!OnlyBadOut && OutputVDWs && ovrlaptype == 0))
                { /*logicals allow the overlaptype of this dot*/

                    /*now derive ptmaster*/
                     if( (src->props & MC_PROP)  && (cause->props & MC_PROP) )
                     {/*potential MC-MC to flag*/
                        if( !(( src->props & HET_PROP) || (cause->props & HET_PROP)) )
                        {/*neither can be a HET, what wierd condition is avoided?*/
                           ptmaster = 'M'; /*McMc interaction*/
                           mcmccont[0][idx]++;
                           //mcmccont[0][5]++;
                           mcmccont[0][7]++; /* SJ 04/15/2015 changing the second dimension to 7 to match new initialization*/
                        }
                     }
                     else if( (src->props & SC_PROP)  && (cause->props & SC_PROP) )
                     {/*potential SC-SC to flag*/
                        if( !(( src->props & HET_PROP) || (cause->props & HET_PROP)) )
                        {/*neither can be a HET, what wierd condition is avoided?*/
                           ptmaster = 'S'; /*ScSc interaction*/
                           scsccont[0][idx]++;
                           //scsccont[0][5]++;
                           scsccont[0][7]++;/* SJ 04/15/2015 changing the second dimension to 7 to match new initialization*/
                        }
                     }
                     else if(  ((src->props & SC_PROP)  && (cause->props & MC_PROP))
                               ||((src->props & MC_PROP)  && (cause->props & SC_PROP)))
                     {/*potential MC-SC or SC-MC to flag*/
                        if( !(( src->props & HET_PROP) || (cause->props & HET_PROP)) )
                        {/*neither can be a HET, what wierd condition is avoided?*/
                           ptmaster = 'P'; /*McSc or ScMc interaction*/
                           mcsccont[0][idx]++;
                           //mcsccont[0][5]++;
                           mcsccont[0][7]++;/* SJ 04/15/2015 changing the second dimension to 7 to match new initialization*/
                        }
                     }
                     else
                     {
                        ptmaster = 'O'; /*Oh for Other, general default*/
                        othrcont[0][idx]++;
                        //othrcont[0][5]++;
                        othrcont[0][7]++;/* SJ 04/15/2015 changing the second dimension to 7 to match new initialization*/
                     }

                     saveDot(src, cause, type, &dotloc, &spikeloc,
        	               results, ovrlaptype, mingap, ptmaster,0); /* dcr041009 */
                   /*XHTangle dummy dcr20120120*/
               }/*logicals allow the overlaptype of this dot*/
            }/*could use this dot/spike, */
         }
      }/*dot is ok, atom that this dot hits is called cause */
   }/*loop over dots on src dot-ball */
}/*examineDots()*/
/*}}}examineDots()___________________________________________________________*/

/*{{{markBonds()**************************************************************/
/*markBonds() decides which atoms are bonded to each other and thus */
/*decides under what conditions close atoms will NOT be considered to clash*/
/*this seems to be the only subroutine that uses atom->bondedto information*/
/*for a given src atom: mark all atoms that are close enough to bond and are */
/*allowed to bond and extend this to thier neighbors through Maxbonded bonds*/
/*including marking the src atom itself if such a Maxbonded-loop returns to it*/

void markBonds(atom *src, atom *neighbors, int distcount, int max)
{/*markBonds  reformated 041112*/
   atom *targ = NULL, *theHatom = NULL, *theOtherAtom = NULL;
   int nearpt, tooclose, isAanH, isBanH, letItBond = FALSE;

   /*when making surface dots, max==1 */
   /*when looking for contacts and clashes, max==Maxbonded */

   for(targ = neighbors; targ; targ = targ->scratch)
   {/*for(loop over neighbors)*/

      /*first pass look at neighbors of external calling atom: targ->mark == 0*/
      /*later look at neighbors of neighbors: targ->mark > max */
      /* until Maxbonded neighbor when targ->mark == max */

      if (targ->mark == 0 || targ->mark > distcount)
      {/*if(targ->mark == 0 || targ->mark > distcount) */

	 /* note: when we are relying on the StdBond table   */
	 /*       to figure out which bonds are allowed      */
	 /*       we can be more generous with the distance  */
	 /*       cuttofs to permit more oddball distortions */
         /*dcr: so it seems that with UseStdBond 2*COVRADFUDGE is added */
         /*atomprops.h/defines COVRADFUDGE 0.2        050111 */

#ifdef INLINE_FOR_SPEED
	 nearpt = INRANGEINLINE((src->loc), (targ->loc),
		     src->covRad + targ->covRad + COVRADFUDGE
		     + (UseStdBond ? COVRADFUDGE : 0.0) );

	 tooclose = INRANGEINLINE((src->loc), (targ->loc),
			   src->covRad + targ->covRad - 0.6
			   - (UseStdBond ? COVRADFUDGE : 0.0) );
#else
	 nearpt = inRange(&(src->loc), &(targ->loc),
		     src->covRad + targ->covRad + COVRADFUDGE
		     + (UseStdBond ? COVRADFUDGE : 0.0) );

	 tooclose = inRange(&(src->loc), &(targ->loc),
			   src->covRad + targ->covRad - 0.6
			   - (UseStdBond ? COVRADFUDGE : 0.0) );
#endif
            /*returns true when distance-sq between a and b <= limit-sq */

/*050111 tooclose=   does NOT seem to be the call where */
/*atoms that move too close together suddenly stop showing bad clashes ! */
/*i.e. force tooclose = 0  here does not affect very close clash cut-off*/
/* and indeed atoms not bonded because of tooclose should show clashes !? */
/*but actual close distance of cut-off seems to be dependent on UseStdBond */
/*so since this is the only place where COVRADFUDGE is used, one suspects */
/*that this markBonds() is responsible for showing as bonded those atoms */
/*that fail to show clashes as they get near to each other*/
/*However, although this markBonds seems to check each atom pair many times, */
/*it does look as if it declares the letItBond flag correctly even when close!*/
/*050111 BUT this seems where all bonding info is calc, so maybe letItBond is */
/*not interpreted correctly by this analysis...*/
/*050117 specific case of ala_dipeptide at phi==0,psi==0, tau<= 102.25 seems */
/*the N of res 2 is just within bonding criteria of C of res 0 and with*/
/*UseStdBond N -- C  bonds are allowed between atoms of different residues*/
/*thus the covalent neighbors of that N and that C, e.g. NH of 2 and O of 0 */
/*are also considered to NOT clash! */
/*For autobondrot and mage/probe strict succession would cure this...*/
/*BUT one does need to allow for deletions in residue order!*/
/*SO this is a FEATURE, NOT a BUG */
/*and the work-around is to restrict the tau deviation to what is more */
/*reasonable for tau anyway of  -8, thus -5 to +10 around ideal 111.1 */
/*might be the best range for getting a feel for allowable phi,psi regions*/

	 /* conditions for allowing a bond to be formed... */

	 if (nearpt && (! tooclose) )
         {/*close enough for a bond*/
	    isAanH = isHatom(src->elem);
	    isBanH = isHatom(targ->elem);

#ifdef OLD_BONDING_CALC
	    if(  (   ! isAanH && ! isBanH ) /* neither atom is a hydrogen  */
	       ||(  (! isAanH || ! isBanH) /*or one is an H */
	          && ATOMS_IN_THE_SAME_RES(src, targ) /* and both are in same*/
                                                      /* res with parent name*/
	          &&(  (theHatom->bondedto == NULL) /*and either unknown */
                     or as expected*/
	             || strstr(theHatom->bondedto, theOtherAtom->atomname)
                                                                           ) ) )
            {
	       targ->mark = distcount;
	       if (distcount < max) {
		  markBonds(targ, neighbors, distcount+1, max);
	       }
	    }
#else /*NOT OLD_BONDING_CALC*/
	    letItBond = FALSE; /* we start out skeptical */
	    if ((! isAanH) && (! isBanH))
            {/* neither atom is a hydrogen - both are heavy atoms */

	       if (UseStdBond)
               {/* if we are allowing only "standard" bonding */

		  if (! ATOMS_IN_THE_SAME_RES(src, targ))
                  { /* and are in diff res */

		     if ((src->props & MC_PROP) && (targ->props & MC_PROP))
                     {/* both are mc in diff res - */
                        /*only specific bonding patterns allowed */

                        /*strcmp(a,b) returns 0 when strings == (a==b) */
                        /*strstr returns 1 when any str in colon delineated*/
                        /* list matches reference str */

			letItBond =
                           (    !strcmp(" N  ", src->atomname)
                             && !strcmp(" C  ",targ->atomname) )
                        || (    !strcmp(" N  ",targ->atomname)
                             && !strcmp(" C  ", src->atomname) )
			|| (    !strcmp(" P  ", src->atomname)
                             &&  strstr(" O3*: O3'", targ->atomname) )
			|| (    !strcmp(" P  ",targ->atomname)
                             &&  strstr(" O3*: O3'",  src->atomname) );
		     }
		     else if(   !(src->props & MC_PROP)
                             && !(targ->props & MC_PROP) )
                     {
			/* both are sc in diff res they must be CYS-CYS SS */
			letItBond = (   !strcmp(src->r->resname,  "CYS")
                                     && !strcmp(src->atomname,   " SG ")
				     && !strcmp(targ->r->resname, "CYS")
                                     && !strcmp(targ->atomname,  " SG ") );
		     }

                     /*050121 remove dcr041112 code insertion */
                        /*inserted test to allow mage/probe fitting of a*/
                        /*particular sidechain in context of an NMR ensemble*/
                     /*Lmodeltest not needed, mage sends model # when needed*/

                  }/* and are in diff res */
		  else
                  {/* both heavy atoms in same res */
                        /*strstr returns 1 when any str in colon delineated*/
                        /* list matches reference str */
#ifdef STRICT_CONN_TABLE /*050111 seems to be not defined*/
		     /* strict - heavy atom must be as exactly expected */
		     letItBond = (   (src->bondedto != NULL)
                                  && strstr(src->bondedto, targ->atomname) );
#else /*NOT STRICT_CONN_TABLE*/
		     /* heavy atom either unknown or as expected */
		     letItBond = (   (src->bondedto == NULL)
                                  || strstr(src->bondedto, targ->atomname) );
#endif /*def? STRICT_CONN_TABLE*/
		  }
               }/* if we are allowing only "standard" bonding */
	       else
               {/* not using std bond table - */
                  /*so we let heavy atoms bond when they are close enough */
		  letItBond = TRUE;
	       }
            }/* neither atom is a hydrogen - both are heavy atoms */
	    else if(   (! isAanH || ! isBanH) /* only one atom is a hydrogen */
	            && ATOMS_IN_THE_SAME_RES(src, targ) ) /* and both atoms*/
                                                       /* are in the same res*/
            {
	       if (isAanH) { theHatom = src;  theOtherAtom = targ; }
	       else        { theHatom = targ; theOtherAtom = src;  }

	       if(  (theHatom->bondedto == NULL) /* heavy atom either unknown*/
                                                 /*  or as expected */
	          || strstr(theHatom->bondedto, theOtherAtom->atomname)
	          || (UseHParent == FALSE) )
               {
		  letItBond = TRUE;
	       }
	    }
	    else
            {/*either both are hydrogens in the same residue- no bond possible*/
	       /*  or one is a H and they are in diff res - no bond possible */
	       letItBond = FALSE;
	    }

	    if (letItBond) /*check connectivity through Maxbonded neighbors*/
            {
	       targ->mark = distcount; /*initially == 1 as called externally*/
	       if (distcount < max)  /*max == Maxbonded  as input or defaulted*/
               {/*reentrant: no clashes between atoms within Maxbonded bonds*/
		  markBonds(targ, neighbors, distcount+1, max);
	       }
               /*else targ->mark == max > discount and loop will end*/
#ifdef DUMP_DEBUG_EXTRA
	       if (DebugLevel > 7)
               {
		  fprintf(stdout, "{%4.4s%c%3.3s%s%4.4s%c}P %8.3f%8.3f%8.3f\n",
			   src->atomname, src->altConf,
			   src->r->resname, src->r->chain,
			   src->r->Hy36resno, src->r->resInsCode,
			   src->loc.x, src->loc.y, src->loc.z);
		  fprintf(stdout, "{%4.4s%c%3.3s%s%4.4s%c}L %8.3f%8.3f%8.3f\n",
			   targ->atomname, targ->altConf,
			   targ->r->resname, targ->r->chain,
			   targ->r->Hy36resno, targ->r->resInsCode,
			   targ->loc.x, targ->loc.y, targ->loc.z);
	       }
#endif /*DUMP_DEBUG_EXTRA*/
	    }
#endif /*def? OLD_BONDING_CALC not defined as of 050111*/
         }/*close enough for a bond*/
      }/*if(targ->mark == 0 || targ->mark > distcount) */

   }/*for(loop over neighbors)*/
}/*markBonds*/
/*}}}markBonds()_____________________________________________________________*/

/*{{{fixupLongBondChains()****************************************************/
/* fixupLongBondChains() - unmark and remove from the bonded set    */
/*                         atoms which are more than cutoff bonds   */
/*                         from the source when neither src or targ */
/*                         is a hydrogen                            */
void fixupLongBondChains(atom *src, atom *neighbors, int cutoff)
{
   atom *targ = NULL;

   if (! isHatom(src->elem)) {
      for(targ = neighbors; targ; targ = targ->scratch) {
	 if (targ->mark > cutoff && ! isHatom(targ->elem)) {
	    targ->mark = 0;
	 }
      }
   }
}
/*}}}fixupLongBondChains()___________________________________________________*/

/*{{{dotType()****************************************************************/
int dotType(atom *src, atom *atomList, int recalcOnly)
{
   atom *a = NULL;
   int rslt = src->atomclass;

   /* When ByNABaseColor is TRUE, the atomclass is not the real element type */
   /* but instead is the nucleic acid base type. In this case, the times we  */
   /* need the actual element (or parent element for hydrogens) recalcOnly is*/
   /* TRUE, we determine the element & do not molest the atomclass assignment*/

   if (src->atomclass < 0 || recalcOnly) {
      rslt = src->elem;

      if (isHatom(src->elem)) {
	 for(a = atomList; a; a = a->scratch) {
	    if (a->mark == 1 && ! isHatom(a->elem)
	                     && src->r == a->r) {
	       rslt = a->elem; break;
	    }
	 }
      }
      if (!recalcOnly) { /* usually reassign atomclass, but not for recalc */
         src->atomclass = rslt;
      }
   }

   return rslt;
}
/*}}}dotType()_______________________________________________________________*/

/*{{{debugBondingLists()******************************************************/
void debugBondingLists(atom *src, atom *neighbors)
{
   atom *targ = NULL;
   int i = 0;

   fprintf(stdout, "%4.4s%c%3.3s%s%4.4s%c(%c%d)",
			   src->atomname, src->altConf,
			   src->r->resname, src->r->chain,
			   src->r->Hy36resno, src->r->resInsCode,
			   src->binSerialNum, src->mark);
   for(targ = neighbors; targ; targ = targ->scratch) {
      fprintf(stdout, ", %d: %4.4s%c%3.3s%s%4.4s%c(%c%d)", ++i,
			   targ->atomname, targ->altConf,
			   targ->r->resname, targ->r->chain,
			   targ->r->Hy36resno, targ->r->resInsCode,
			   targ->binSerialNum, targ->mark);
   }
   fprintf(stdout, "\n");
}
/*}}}debugBondingLists()_____________________________________________________*/

/*{{{genDotIntersect()*** called from doCommand() *** calls examineDots() ****/
/*genDotIntersect() called from doCommand during processing of modes: */
/*once with srcFlag==1, targFlag==2 in modes SELFINTERSECT,INTERSECTONCE, */
/*and twice from INTERSECTBOTHWAYS first with srcFlag,targFlag == 1,2 then 2,1*/
/*autobondrot static atoms are allMainAtoms with abin bins*/
/* mobile atoms are allMovingAtoms with bbins bins*/

/* NOTE: allMainAtoms/abins & allMovingAtoms/bbins must be disjoint */
/*       sets of atoms (none in common) or allMovingAtoms/bbins can */
/*       be NULL.                                                   */
/*       allMovingAtoms refers to autobondrot set of atoms          */

void genDotIntersect(atom *allMainAtoms, atomBins *abins,
			atom *allMovingAtoms, atomBins *bbins,
			pointSet dots[],
			float probeRad, float spikelen,
			int srcFlag, int targFlg, dotNode *results[][NODEWIDTH])
{/*genDotIntersect()*/

   atom *src = NULL, *atomList = NULL, *atomList2 = NULL;
   int type = 0, usesMovingAtoms = FALSE;
   int oktargsA = TRUE, oktargsB = TRUE;

    int i=0,j=0;/* SJ - 05/06/2015 moved the variables here from below for compatibility with windows C compiler*/ 
    atom *source = NULL, *target=NULL;
    dotNode *head = NULL, *curNode = NULL, *prevNode = NULL;

   if(Ldotdump) {
    fprintf(stderr,"entered genDotIntersect()\n");
   }

   usesMovingAtoms = ((allMovingAtoms != NULL) && (bbins != NULL));

   for(src = allMainAtoms; src; src = src->next) /*main==autobondrotstatic*/
   {/*for: loop over all main atoms taking each in turn as the src-atom*/

      if (src->flags & srcFlag)  /*dcr?: seems srcFlag always either 1 or 2 ??*/
      {/*for each src atom*/

          oktargsA = TRUE; oktargsB = TRUE;
          atomList = findTouchingAtoms(src, NULL,abins, probeRad, targFlg, &oktargsA);

           /*20111206 findTouchingAtoms() sets ~exactly overlapping atoms to be NOT oktargsA*/
           /*this is where too close things fail to show dots !! ?? */

          if (usesMovingAtoms) /*autobondrot*/
          {
	            atomList2 = findTouchingAtoms(src, atomList, bbins, probeRad, targFlg, &oktargsB);
	        }
	        else
          {
               atomList2 = atomList;
          }

          if(Ldotdump)
              fprintf(stderr,"findTouchingAtoms() sets oktargsA== %d \n",oktargsA);
          if (atomList2 && (oktargsA || (usesMovingAtoms && oktargsB) ))
          {
	           markBonds(src, atomList2, 1, Maxbonded); /*in genDotIntersect()*/
              /*markBonds identifies bonds between atoms - */
              /*seems not to have lower distance limit*/
              /*but this is where autobondrot fails to spike very close atoms*/

              if (Maxbonded > 3)
                   fixupLongBondChains(src, atomList2, 3);

	           type = dotType(src, atomList2, FALSE);

              if(LOneDotEach)
              {
                   if(Ldotdump)
                     fprintf(stderr,"genDotIntersect() calls examineOneDotEach()\n");
                
                   examineOneDotEach(src,type,atomList2,dots,probeRad,spikelen,targFlg,results,allMainAtoms);
                  /*allMainAtoms20120120*/
              }
              else
              {
                  examineDots(src, type, atomList2, dots, probeRad, spikelen, targFlg, results);
              }
              if(Verbose && ShowTicks)
              {
		            fprintf(stderr, "%s%d   \r",
			          src->r->resname, src->r->resid);
              }
          }
      }/*for each src atom*/
   }/*for: loop over all static atoms taking each in turn as the src-atom*/

   if(usesMovingAtoms)
   {/*if: usesMovingAtoms==autobondrot*/

      for(src = allMovingAtoms; src; src = src->next)
      {/*for: loop over allMovingAtoms taking each in turn as the src-atom*/

	       if (src->flags & srcFlag)  /*dcr?: what is being flagged here??*/
         {

	           oktargsA = TRUE;
              oktargsB = TRUE;
	           atomList  = findTouchingAtoms(src, NULL,     abins, probeRad, targFlg, &oktargsA);
	           atomList2 = findTouchingAtoms(src, atomList, bbins, probeRad, targFlg, &oktargsB);
	           if (atomList2 && (oktargsA || oktargsB))
              {
	               markBonds(src, atomList2, 1, Maxbonded);
               /*markBonds identifies bonds between atoms - */
               /*seems not to have lower distance limit*/
              /*but is where autobondrot fails to spike very close atoms*/

	               if (Maxbonded > 3)
                       fixupLongBondChains(src, atomList2, 3);

                 type = dotType(src, atomList2, FALSE);

                  /*usesMovingAtoms==autobondrot does not use LOneDotEach option 111205*/
                 examineDots(src,type,atomList2,dots,probeRad,spikelen,targFlg,results);

	               if(Verbose && ShowTicks)
                 {
		                 fprintf(stderr, "%s%d   \r",
			              src->r->resname, src->r->resid);
                 }
	          }
	      }
      }/*for: loop over allMovingAtoms taking each in turn as the src-atom*/
   }/*if: usesMovingAtoms==autobondrot*/
    
    // SJ - 10/03/2014 - Code to delete contact dots (except H bonds) between atom pairs when atleast one of them has an altConf " " and the sum of their occpancies is <= 1
    // UNKNOWN BUGS: Not sure what this code will do if the results is genearted using autobondrot
    
    for(i=0;i<NUMATOMTYPES;i++)
    {/*loop over atom types*/
       // for(j=0;j<4;j++) // Old: first j:0-3 correspond to wc, cc, so, bo. Not looping over H-bond dotNode, last one is H bonds
        for(j=0; j< NODEWIDTH - 1; j++) // changed on 04/09/2015 SJ to accomodate for extra list for weak H bonds and worse overlaps
            /* see NODEWIDTH definition in probe.h - SJ*/
        {/*loop over wc, cc, weak H bonds (wh, if LweakHbonds is true), so, bo, and worse overlap (wo, if LworseOverlap is true)*/
            
            head = results[i][j];
            prevNode=NULL;
            curNode=head;
            while(curNode) // if the list is empty, it will not enter this loop
            {
                source = curNode->a;
                target = curNode->t;
                if((source->altConf == ' ' || target->altConf == ' ') && source->occ + target->occ <= 1)
                { //this node has to be deleted
                    if(curNode==head)
                    {
                        head=curNode->next;
                        free(curNode);
                        prevNode=NULL;
                        curNode=head;
                    }
                    else
                    {
                        prevNode->next=curNode->next;
                        free(curNode);
                        curNode=prevNode->next;
                    }
                }
                else
                {//move forward when no deletion
                    prevNode=curNode;
                    curNode=curNode->next;
                }
            }/*while loop*/
            results[i][j]=head;
        }/*loop over wc, cc, weak H bonds (wh, if LweakHbonds is true), so, bo, and worse overlap (wo, if LworseOverlap is true)*/
    }/*for loop over atomtypes*/
    //results now is filtered
    
}/*genDotIntersect()*/
/*}}}genDotIntersect()_______________________________________________________*/

/*{{{genDotSurface() only called when method == EXTERNALSURFACE **************/

/* NOTE: allMainAtoms/abins & allMovingAtoms/bbins must be disjoint */
/*       sets of atoms (none in common) or allMovingAtoms/bbins can */
/*       be NULL.                                                   */
/*       allMovingAtoms refers to autobondrot set of atoms          */

void genDotSurface(atom *allMainAtoms, atomBins *abins,
			atom *allMovingAtoms, atomBins *bbins,
			pointSet dots[],
			float probeRad, float spikelen, int srcFlag,
			dotNode *results[][NODEWIDTH])
{
   atom *src = NULL, *atomList = NULL, *atomList2 = NULL;
   int type = 0, usesMovingAtoms = FALSE;
   int dummy = TRUE;

   usesMovingAtoms = ((allMovingAtoms != NULL) && (bbins != NULL));

   for(src = allMainAtoms; src; src = src->next) {
      if (src->flags & srcFlag) {

	 atomList  = findTouchingAtoms(src, NULL,     abins, probeRad, 0, &dummy);
	 if (usesMovingAtoms) {
	    atomList2 = findTouchingAtoms(src, atomList, bbins, probeRad, 0, &dummy);
	 }
	 else { atomList2 = atomList; }

	 if (atomList2) {
	    markBonds(src, atomList2, 1, 1); /*in genDotSurface()*/
	 }

	 type = dotType(src, atomList2, FALSE);

	 surfDots(src, type, atomList2, dots,
		     probeRad, spikelen, results);

	 if(Verbose && ShowTicks) {
	    fprintf(stderr, "%s%d   \r",
			src->r->resname, src->r->resid);
	 }
      }
   }

   if (usesMovingAtoms) {
      for(src = allMovingAtoms; src; src = src->next) {
	 if (src->flags & srcFlag) {

	    atomList  = findTouchingAtoms(src, NULL,     abins, probeRad, 0, &dummy);
	    atomList2 = findTouchingAtoms(src, atomList, bbins, probeRad, 0, &dummy);
	    if (atomList2) {
	       markBonds(src, atomList2, 1, 1); /*in genDotSurface()*/
	    }

	    type = dotType(src, atomList2, FALSE);

	    surfDots(src, type, atomList2, dots,
		     probeRad, spikelen, results);

	    if(Verbose && ShowTicks) {
	       fprintf(stderr, "%s%d   \r",
			src->r->resname, src->r->resid);
	    }
	 }
      }
   }
}
/*}}}genDotIntersect()_______________________________________________________*/

/*{{{surfDots() only called from genDotSurface(), method == EXTERNALSURFACE **/
void surfDots(atom *src, int type, atom *scratch, pointSet dots[],
		float probeRad, float spikelen, dotNode *results[][NODEWIDTH])
{
   atom *targ;
   int i, nearpt, ok;
   point3d dotloc, dotvect, exploc, spikeloc;
   pointSet *srcDots;
   char ptmaster = ' '; /*pointmaster character  dcr041009*/

   if (src->elem == ignoreAtom) { return; }

   srcDots = &(dots[src->elem]);
   if (src->elem == atomC && isCarbonylAtom(src)) {
      srcDots = &COdots;
#ifdef DUMP_DEBUG_EXTRA
      if (DebugLevel>8) {
	 fprintf(stdout, "DEBUG surfDots(%s%d %s radius = %.3f)\n",
	       src->r->resname,src->r->resid,src->atomname,src->radius);
      }
#endif
   }

   for(i=0; i < srcDots->n; i++) {
      dotvect = srcDots->p[i];

      v3add(&(src->loc), &dotvect, &dotloc);

      v3scale(&dotvect, src->radius + probeRad);
      v3add(&(src->loc), &dotvect, &exploc);

      v3scale(&dotvect, src->radius + 0.0);
      v3add(&(src->loc), &dotvect, &spikeloc);

      ok = TRUE;

      /*targ is an atom*, scratch is an atom*,  atom* defined in abin.h */
      /* atom* scratch is a member of atom*: which has member scratch...*/
      /* so for-loop moves through all these atoms while an atom has a scratch*/

      for(targ = scratch; targ; targ = targ->scratch) {
	 if ( (targ->elem == ignoreAtom)
	   || ((!DoWatWat) && (src->props &  WATER_PROP)
	                   && (targ->props & WATER_PROP)
	          && (src->r != targ->r)) /* we do see HOH hydrogens */
	   || ((!DoHet)    && (targ->props &   HET_PROP))
	   || ((!DoH2O)    && (targ->props & WATER_PROP)) ) {
	       /* ignore */
	 }
	 else {
#ifdef INLINE_FOR_SPEED
	    nearpt = INRANGEINLINE(exploc, (targ->loc),
			probeRad+targ->radius);
#else
	    nearpt = inRange(&exploc, &(targ->loc),
			probeRad+targ->radius);
#endif

	    if (nearpt) { ok = FALSE; break; }
	 }
      }
      if (ok) {
	 saveDot(src, NULL, type, &dotloc, &spikeloc, results, 0, 0.0,ptmaster,0);
         /*XHTangle dummy dcr20120120*/
         /*ptmaster not used to distinguish anykind of surface dots  dcr041009*/
         /*Note:  here only when method == EXTERNALSURFACE */
         /*and ovrlaptype==0, mingap==0.0, ptmaster==' ' */
      }
   }
}/*surfDots()*/
/*}}}genDotIntersect()_______________________________________________________*/

/*{{{initResults()************************************************************/
void initResults(dotNode *results[][NODEWIDTH])
{
   int i, j;

   for (i = 0; i < NUMATOMTYPES; i++) {
      for (j = 0; j < NODEWIDTH; j++) {
	 results[i][j] = NULL;
      }
   }
}
/*}}}initResults()___________________________________________________________*/

/*{{{freeResults()************************************************************/
void freeResults(dotNode *results[][NODEWIDTH])
{
   int i, j;
   dotNode *node, *next;
   atom * a;

   for (i = 0; i < NUMATOMTYPES; i++) {
      for (j = 0; j < NODEWIDTH; j++) {
	 for (node = results[i][j]; node; node = next) {
	    next = node->next;
	    free(node);
	 }
	 results[i][j] = NULL;
      }
   }
}
/*}}}freeResults()___________________________________________________________*/

/*{{{assignGapColorForKin()***************************************************/
char* assignGapColorForKin(float gap, int class)
{
   char *colorValue = "";
   //if (class == 4)     { colorValue = "greentint "; } /* hbond */ /* 04/16/2015 SJ changed 4 to 6, as H-bond are now 6 in NODEWIDTH, see probe.h*/
   if (class == 6)     { colorValue = "greentint "; }
   else if (gap > 0.35){ colorValue = "blue ";      }
   else if (gap > 0.25){ colorValue = "sky ";       }
   else if (gap > 0.15){ colorValue = "sea ";       }
   else if (gap > 0.0) { colorValue = "green ";     }
   else if (gap >-0.1) { colorValue = "yellowtint ";}
   else if (gap >-0.2) { colorValue = "yellow ";    }
   else if (gap >-0.3) { colorValue = "orange ";    }
   else if (gap >-0.4) { colorValue = "red ";       }
   else                { colorValue = "hotpink ";   }
   return colorValue;
}
/*}}}assignGapColorForKin()__________________________________________________*/

/*{{{assignGapColorForO()*****************************************************/
char* assignGapColorForO(float gap, int class)
{
   char *colorValue = "";
   //if (class == 4)     { colorValue = "pale_green ";      } /* hbond *//* 04/16/2015 SJ changed 4 to 6, as H-bond are now 6 in NODEWIDTH, see probe.h*/
   if (class == 6)     { colorValue = "pale_green ";      }
   else if (gap > 0.35){ colorValue = "cornflower_blue "; }
   else if (gap > 0.25){ colorValue = "sky_blue ";        }
   else if (gap > 0.15){ colorValue = "aquamarine ";      }
   else if (gap > 0.0) { colorValue = "green ";           }
   else if (gap >-0.1) { colorValue = "yellow_green ";    }
   else if (gap >-0.2) { colorValue = "yellow ";          }
   else if (gap >-0.3) { colorValue = "orange ";          }
   else if (gap >-0.4) { colorValue = "red ";             }
   else                { colorValue = "orange_red ";      }
   return colorValue;
}
/*}}}assignGapColorForO()____________________________________________________*/

/*{{{assignGapColorForXV()****************************************************/
char* assignGapColorForXV(float gap, int class)
{
   return assignGapColorForKin(gap, class);
}
/*}}}assignGapColorForXV()___________________________________________________*/

/*{{{convertKinColorToO()*****************************************************/
char* convertKinColorToO(char* incolor)
{
   char *outcolor = "light_gray";
        if (strcmp(incolor, "red"       ) == 0) { outcolor = "red";             }
   else if (strcmp(incolor, "green"     ) == 0) { outcolor = "green";           }
   else if (strcmp(incolor, "blue"      ) == 0) { outcolor = "cornflower_blue"; }
   else if (strcmp(incolor, "cyan"      ) == 0) { outcolor = "cyan";            }
   else if (strcmp(incolor, "yellow"    ) == 0) { outcolor = "yellow";          }
   else if (strcmp(incolor, "magenta"   ) == 0) { outcolor = "magenta";         }
   else if (strcmp(incolor, "white"     ) == 0) { outcolor = "white";           }
   else if (strcmp(incolor, "pink"      ) == 0) { outcolor = "salmon";          }
   else if (strcmp(incolor, "orange"    ) == 0) { outcolor = "orange";          }
   else if (strcmp(incolor, "purple"    ) == 0) { outcolor = "purple";          }
   else if (strcmp(incolor, "sky"       ) == 0) { outcolor = "sky_blue";        }
   else if (strcmp(incolor, "brown"     ) == 0) { outcolor = "brown";           }
   else if (strcmp(incolor, "gray"      ) == 0) { outcolor = "light_gray";      }
   else if (strcmp(incolor, "black"     ) == 0) { outcolor = "black";           }
   else if (strcmp(incolor, "gold"      ) == 0) { outcolor = "gold";            }
   else if (strcmp(incolor, "yellowtint") == 0) { outcolor = "yellow_green";    }
   else if (strcmp(incolor, "sea"       ) == 0) { outcolor = "aquamarine";      }
   else if (strcmp(incolor, "pinktint"  ) == 0) { outcolor = "pink";            }
   else if (strcmp(incolor, "bluetint"  ) == 0) { outcolor = "light_blue";      }
   else if (strcmp(incolor, "greentint" ) == 0) { outcolor = "pale_green";      }
   else if (strcmp(incolor, "hotpink"   ) == 0) { outcolor = "orange_red";      }
   else if (strcmp(incolor, "invisible" ) == 0) { outcolor = "black";           }
   else                                         { outcolor = "light_gray";      }
   return outcolor;
}
/*}}}convertKinColorToO()____________________________________________________*/

/*{{{convertKinColorToXV()****************************************************/
char* convertKinColorToXV(char* incolor)
{
   return incolor;
}
/*}}}convertKinColorToXV()___________________________________________________*/

/*{{{writeOutput()************************************************************/
void writeOutput(FILE *outf, char* groupname, dotNode *results[][NODEWIDTH], int spike, int method, char* extrastr, float probeRad)
{ /*writeOutput for kinemage*/
   /*add probeRad 20111220dcr*/
   /*dcr041020 need method for better kinemage keywords*/
   /* 060129 extrastr for extra master to control orig vs fitted dots */
   int i, j;
   int k; /*counter for gapbins 111020dcr*/
   char ptm[2] = {' ','\0'}; /*20120120dcr ptmaster for bval*/
   dotNode *node;
   atom *a; /*owner of dot*/
   atom *t; /*cause of dot, output for OneDotEach mode 110113*/
   char *color = "";
  // OLD: char *mast[NODEWIDTH] = {"wide  contact", "close contact", "small overlap",
                           // "bad overlap", "H-bonds", "H-weaks"};
                            /*111215dcr weakHbonds option*/
   char *mast[NODEWIDTH] = {"wide contact", "close contact", "weak H-bonds", "small overlap", "bad overlap",
                            "worse overlap", "H-bonds"}; /* 04/13/2015 SJ changed to match new NODEWIDTH, see probe.h*/
   char *surfacestr = "surface"; /*041020*/
   char *contactstr = "vdw contact"; /*060129*/
   char extraMstr[32]; /*fill below with extra master name 060129*/
   char pointid[100], lastpointid[100];
   char ptmast[7]={'\0','\0','\0','\0','\0','\0','\0'}; /*111020dcr 6 to 7 */
   /*dcr041009 string to hold 'M' (McMc), 'S' (ScSc), 'P' (McSc), 'O' (other)*/
   /*dcr041009 ptmaster part of each dotNode member of results[][]*/
   char masterchr=' '; /*dcr041017*/
/*111020dcr gapbin ptmasters*/
   int ngapbins = 22; /*(std probe diameter .5 +clash overlap .5=1.0/0.05interval) +2*/
   /*+1 for 0th fence post, another +1 for X excess*/
   char gapnamestr[22] =
{'z','y','x','w','v','u','t','g','r','q','f','F','Q','R','G','T','U','V','W','X','Y','Z'};
/*0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21*/
   /*std gapbins scope at least -.5 to +.5, wider if probeRad > 0.25 standard*/
   int  gaplimit = ((2*(probeRad>0.25?probeRad:0.25) +0.5)/0.05)+2;  /*20111220dcr*/
   long *gapcounts = (long *) malloc(gaplimit * sizeof(long));
   long maxgapcounts = 0;
   int  Lgotgapbin = FALSE;
   char strcName[6]; /*first nendstr char of inputfilename 20120120dcr*/
   int  nendstr = 6;

   if(Lgapbins) /*111020dcr gapbins pointmasters, 111215dcr count bin entries*/
   {
      for(k=0; k<gaplimit;k++)
      {
         gapcounts[k] = 0;
      }
   }

   if(LMergeContacts)  /*060129*/
   {
      mast[0] = contactstr;
      mast[1] = contactstr;
   }
   if(LMasterName)  /*060129*/
   {
      sprintf(extraMstr," master={%s}",extrastr);
   }
   else{extraMstr[0] = '\0';} /*occupies no space in output*/

   fprintf(outf, "@subgroup dominant {%s}\n", groupname);
   /*list masters for contact types dcr041020*/
   if(method == EXTERNALSURFACE) /*dcr041020*/
   {/*a bit crude, but better than calling all surfaces a close contact*/
      mast[1] = surfacestr; /*the only master button invoked dcr041020*/
      fprintf(outf, "@master {%s}\n", mast[1]); /*was close contact*/
   }
   else
   {/*fine control over masters names   dcr041020*/
      if (OutputVDWs && !OnlyBadOut) { /*dcr041010*/
         fprintf(outf, "@master {%s}\n", mast[0]);
        if(!LMergeContacts)
         fprintf(outf, "@master {%s}\n", mast[1]);
      }
      if (OutputClashes || OnlyBadOut) { /*dcr041010*/
         if(!OnlyBadOut) {
             //fprintf(outf, "@master {%s}\n", mast[2]);
             fprintf(outf, "@master {%s}\n", mast[3]); /* 04/13/2015 SJ changed to match new number for small overlaps*/
         }
         //fprintf(outf, "@master {%s}\n", mast[3]); /*the Bad clashes*/
         fprintf(outf, "@master {%s}\n", mast[4]); /* 04/13/2015 SJ changed to match new number for bad overlaps*/
          if (LworseOverlap) {/* 04/13/2015 SJ if the worse overlaps are separated*/
              fprintf(outf, "@master {%s}\n", mast[5]);
          }
      }
      if (OutputHBs && !OnlyBadOut)
      { /*dcr041010*/
         //fprintf(outf, "@master {%s}\n", mast[4]);/* 04/13/2015 SJ changed to match new number for H bonds*/
         fprintf(outf, "@master {%s}\n", mast[6]);
         if(LweakHbonds) /*111215dcr*/
         {
            //fprintf(outf, "@master {%s}\n", mast[5]);
             fprintf(outf, "@master {%s}\n", mast[2]);/* 04/13/2015 SJ changed to match new number for weak H bonds*/
         }
      }
   }/*fine control over masters names*/ /* 04/15/2015 SJ changing all 5s to 7s to reflect new initialization*/
   //if(mcmccont[0][5] > 0)
   if(mcmccont[0][7] > 0)
   {
      masterchr = MCMCCHR; /*single char*/
      fprintf(outf, "@pointmaster '%c' {McMc contacts}\n",masterchr);
   }
   //if(scsccont[0][5] > 0)
   if(scsccont[0][7] > 0)
   {
      masterchr = SCSCCHR; /*single char*/
      fprintf(outf, "@pointmaster '%c' {ScSc contacts}\n",masterchr);
   }
   //if(mcsccont[0][5] > 0)
   if(mcsccont[0][7] > 0)
   {
      masterchr = MCSCCHR; /*single char*/
      fprintf(outf, "@pointmaster '%c' {McSc contacts}\n",masterchr);
   }
   //if(othrcont[0][5] > 0)
   if(othrcont[0][7] > 0)
   {
      masterchr = OTHERCHR; /*single char*/
      fprintf(outf, "@pointmaster '%c' {Hets contacts}\n",masterchr);
   }
   if(Lgapbins) /*111020dcr*/
   {
      for(i = 0; i < ngapbins; i++)
      {
        fprintf(outf, "@pointmaster '%c' {gap %3.2f}\n"
                ,gapnamestr[i],(((float)i-11.0)/20.0)+0.05); /*20111220dcr +0.05*/
      }
   }
   for (i = 0; i < NUMATOMTYPES; i++)
   {/*i: Hatom varients + rest of periodic table + base types, see atomprops.h*/
      for (j = 0; j < NODEWIDTH; j++)
      {/* OLD j: wide contact, close contact, small overlap, bad overlap, H-bonds*/
              /*OLD optionally: H-weaks if LweakHbonds  111215dcr*/
          /* NEW 04/13/2015 - SJ j: wide contact, close contact, weak H bonds, small overlap, bad overlap, worse overlap, H-bonds*/
         node = results[i][j]; /*head of list of dot nodes, i.e. first node*/
         if (node)
         {/*node exists, write list headers*/
           if (j == 0 || j == 1) /*wide and close contacts*/
           {/* contact */
             if (AtomMasters)
             {
                if(LOneDotEach)
                {
                   fprintf(outf,
                   "@balllist {x} radius=0.2 color=%s master={%s dots} master={%s}%s%s\n",
                   getColor(i), getAtomName(i), mast[j],extraMstr,
                  LensDots ? " lens" : "");
                }
                else
                {
                   fprintf(outf,
                   "@dotlist {x} color=%s master={%s dots} master={%s}%s%s\n",
                   getColor(i), getAtomName(i), mast[j],extraMstr,
                   LensDots ? " lens" : "");
                }
             }
             else
             {
                if(LOneDotEach)
                {
                   fprintf(outf,
                     "@balllist {x} radius=0.2 color=%s master={%s}%s%s\n",
                   getColor(i), mast[j],extraMstr,
                   LensDots ? " lens" : "");
                }
                else
                {
                    fprintf(outf,
                      "@dotlist {x} color=%s master={%s}%s%s\n",
                    getColor(i), mast[j],extraMstr,
                    LensDots ? " lens" : "");
                }
             }
           }/* contact */
          // else if ((j == 2 || j == 3) && spike)
           else if ((j == 3 || j == 4 || j == 5) && spike) /* 04/13/2015 SJ changed to reflect new numbers for different categories*/
           {/* bump w/ spike */
			  if (AtomMasters) 
              {
                 if(LOneDotEach)
                 {
					  fprintf(outf,
					  "@vectorlist {x} width=7 color=%s master={%s dots} master={%s}%s\n", 
					 getColor(i), getAtomName(i), mast[j],extraMstr); 
			     }
                 else
                 {
					 fprintf(outf,
					 "@vectorlist {x} color=%s master={%s dots} master={%s}%s\n",
					 getColor(i), getAtomName(i), mast[j],extraMstr); 
			     }
              }
              else
              {
                 if(LOneDotEach)
                 {
					fprintf(outf,
					 "@vectorlist {x} width=7 color=%s master={%s}%s\n",
					 getColor(i), mast[j],extraMstr); 
	             }
                 else
                 {
					 fprintf(outf,
					 "@vectorlist {x} color=%s master={%s}%s\n",
					 getColor(i), mast[j],extraMstr); 
			     }
              }
           }/* bump w/ spike */
           //else if (j == 4)
           else if (j == 6)/* 04/13/2015 SJ changed to reflect new numbers for different categories*/
           {/* bump w/o spike or H-bond */
             if (AtomMasters)
             {
                if(LOneDotEach)
                {
                   fprintf(outf,
                   "@ringlist {x} radius=0.2 color=%s master={%s dots} master={%s}%s\n",
                   getColor(i), getAtomName(i), mast[j],extraMstr);
                }
                else
                {
                   fprintf(outf,
                    "@dotlist {x} color=%s master={%s dots} master={%s}%s\n",
                   getColor(i), getAtomName(i), mast[j],extraMstr);
                }
             }
             else
             {
                if(LOneDotEach)
                {
                   fprintf(outf,
                   "@ringlist {x} radius=0.2 color=%s master={%s}%s\n",
                   getColor(i), mast[j],extraMstr);
                }
                else
                {
                   fprintf(outf,
                  "@dotlist {x} color=%s master={%s}%s\n",
                   getColor(i), mast[j],extraMstr);
                }
             }
           }/* bump w/o spike or H-bond */
           else if (LweakHbonds)
           {/* H-weak */
             if (AtomMasters)
             {
                if(LOneDotEach)
                {
                   fprintf(outf,
                   "@ringlist {x} radius=0.2 color=%s master={%s dots} master={%s}%s\n",
                   getColor(i), getAtomName(i), mast[j],extraMstr);
                }
                else
                {
                   fprintf(outf,
                    "@dotlist {x} color=%s master={%s dots} master={%s}%s\n",
                   getColor(i), getAtomName(i), mast[j],extraMstr);
                }
             }
             else
             {
                if(LOneDotEach)
                {
                   fprintf(outf,
                   "@ringlist {x} radius=0.2 color=%s master={%s}%s\n",
                   getColor(i), mast[j],extraMstr);
                }
                else
                {
                   fprintf(outf,
                  "@dotlist {x} color=%s master={%s}%s\n",
                   getColor(i), mast[j],extraMstr);
                }
             }
           }/* H-weak */

           lastpointid[0] = '\0'; /* reset */
         }/*node exists, write list headers*/
	 while(node)
         {/*kinemage point or point-line for each node in (ij)th dot-node list*/
            a = node->a; /*owner of dot*/
            t = node->t; /*cause of dot, use for OneDotEach 111013dcr*/
            if(Lgapbins) /*111020dcr*/
            {/*two characters for two pointmasters*/
               ptmast[0] = ' '; /*insure leading space*/
               ptmast[1] = '\''; /*surround char pointmasters with*/
               ptmast[2] = node->ptmaster;
               ptmast[3] = ' '; /*space for later gapbin character 111020dcr*/
               ptmast[4] = '\''; /*single quote marks*/
               ptmast[5] = ' '; /*need trailing space*/
               ptmast[6] = '\0'; /*end string*/
            }
            else if(node->ptmaster == ' ') /*dcr041009*/
            {/*blank means NO pointmaster*/
               ptmast[0] = '\0'; /*which will print as a zero-length str*/
            }
            else
            {
               ptmast[0] = ' '; /*insure leading space*/
               ptmast[1] = '\''; /*surround single char pointmaster with*/
               ptmast[2] = node->ptmaster;
               ptmast[3] = '\''; /*single quote marks*/
               ptmast[4] = ' '; /*need trailing space*/
               ptmast[5] = '\0'; /*end string*/
            }
/*     sprintf(pointid, "%s%c %s%d%c", */
/*       a->atomname, a->altConf, */
/*       a->r->resname, a->r->resid, */
/*       a->r->resInsCode); */
            if(LOneDotEach) /*report cause of dot after owner of dot 111013dcr*/
            {
               sprintf(pointid, "%s%c%s%s%c%s %s%c%s%s%c%s GAP %6.3f",
                  a->atomname, a->altConf,
                  a->r->resname, a->r->Hy36resno,
                  a->r->resInsCode, a->r->chain,
                  t->atomname, t->altConf,
                  t->r->resname, t->r->Hy36resno,
                  t->r->resInsCode, t->r->chain, node->gap);
            }
            else
            {/*usual just report owner of dot*/
               sprintf(pointid, "%s%c%s%s%c%s",
               a->atomname, a->altConf,
               a->r->resname, a->r->Hy36resno,
               a->r->resInsCode, a->r->chain);
            }

            if (strcmp(pointid, lastpointid))
            {
               strcpy(lastpointid, pointid);
               fprintf(outf, "{%s}", pointid);
            }
            else
            {
               fprintf(outf, "{\"}");
            }
            if (ColorGap)
            {
               if (node->t)
               {
                  color = assignGapColorForKin(node->gap, j);
                  fprintf(outf, "%s", color);
               }
               else
               {
                  fprintf(outf, "%s ", OutPointColor);
               }
               /* added "%s" string format to color and */
               /* edited "%s " to "%s" in OutPointColor, wba 110909 */
               /* edited "%s" back to "%s " in OutPointColor, gjk 1200907 */
            }
            if(Lgapbins) /*111020dcr gapbins pointmasters*/
            {
               Lgotgapbin = FALSE; /*until identify which gapbin*/
               for(k=0; k<gaplimit ;k++) /*20111220dcr ngapbins, just std range*/
               {
                  if(node->gap < (((float)k-11.0)/20.0)+0.05 )
                  {/*pt master intervals of 0.05 from -0.5 to +0.5*/
                     ptmast[3] = gapnamestr[k];
                     gapcounts[k] = gapcounts[k] +1;
                     if(gapcounts[k] > maxgapcounts){maxgapcounts = gapcounts[k];}
                     if(k < ngapbins)
                     {
                        Lgotgapbin = TRUE; /*got the gapbin idenity for this dot*/
                     }
                     break; /*out of for loop*/
                  }
               }
               if(!Lgotgapbin)
               {/*assign this node, aka dot, to overflow gapbin*/
                  ptmast[3] = gapnamestr[ngapbins-1]; /*the overflow bin-master*/
                  gapcounts[gaplimit-1] = gapcounts[gaplimit-1] +1; /*last gapbin*/
                  /*if(gapcounts[gaplimit-1] > maxgapcounts)
                          {maxgapcounts = gapcounts[gaplimit-1];}*/
               }
            }
            //if ((j == 2 || j == 3) && spike)
            if ((j == 3 || j == 4 || j == 5) && spike)/* 04/13/2015 SJ changed to reflect new numbers for different categories*/
            {/* bump */
              fprintf(outf,
                     "P %s%.3f,%.3f,%.3f {\"}%s %s%.3f,%.3f,%.3f\n", /*dcr041009*/
                      ptmast,node->loc.x, node->loc.y, node->loc.z, /*dcr041009*/
                      color,      /*** note: second color reference */
                      ptmast,node->spike.x, node->spike.y, node->spike.z); /*dcr041009*/
            }
            else
            {/* contact or H-bond */
              fprintf(outf, "%s%.3f,%.3f,%.3f\n", /*dcr041009*/
                      ptmast,node->loc.x, node->loc.y, node->loc.z); /*dcr041009*/
            }
            node = node->next; /*in this (ij)th dot-node list*/
         }/*kinemage point or point-line for each node in (ij)th dot-node list*/
      }/*OLD j: wide contact, close contact, small overlap, bad overlap, H-bonds*/ /* see new def in NODEWIDTH probe.h*/
   }/*i: Hatom varients + rest of periodic table + base types, see atomprops.h*/
   if(Lgapbins) /*111020dcr gapbins pointmasters, 111215dcr count bin entries*/
   {/*Lgapbins TRUE*/
      fprintf(outf, "@text\n");
      for(k=0; k<gaplimit;k++) /*20111220dcr gaplimit*/
      {
          fprintf(outf, "{%5.2f, %8.0ld }\n",
             ((((float)k-11.0)/20.0)+0.05),gapcounts[k]);
      }

/*kinemage 2*/
      fprintf(outf, "@kinemage 2\n");
      fprintf(outf, "@group {gapbins} dominant\n");
      fprintf(outf, "@vectorlist {gapbins}\n");
      for(k=0; k<(gaplimit-1);k++) /*not plot sum of all over-wide dots*/
      {
          fprintf(outf, "{%5.2f, %8.0ld } %5.2f, %8f, 0.00\n",
             ((((float)k-11.0)/20.0)+0.05),gapcounts[k],
             ((((float)k-11.0)/20.0)+0.05)*(float)maxgapcounts,(float)gapcounts[k]);
      }
      fprintf(outf, "@labellist {gapbins}\n");
      fprintf(outf, "{0} 0.0, -1.0, 0.0\n");

/*kinemage 3*/
     if(LXHvector)
     {
      /*write list headers*/
      fprintf(outf, "@kinemage 3\n");
      /*for(k=1; k<10; k++){fprintf(outf, "@pointmaster '%d' {B<=%d}\n",k,10*k);}*/
      fprintf(outf, "@pointmaster 'A' {B<=10}\n");
      fprintf(outf, "@pointmaster 'B' {B<=15}\n");
      fprintf(outf, "@pointmaster 'C' {B<=20}\n");
      fprintf(outf, "@pointmaster 'D' {B<=25}\n");
      fprintf(outf, "@pointmaster 'E' {B<=30}\n");
      fprintf(outf, "@pointmaster 'F' {B<=40}\n");
      fprintf(outf, "@pointmaster 'G' {B<=50}\n");
      fprintf(outf, "@pointmaster 'H' {B >50}\n");

      fprintf(outf, "@group {AXES} dominant\n");
      fprintf(outf, "@labellist {gapbins}\n");
      for(k=0; k<(gaplimit-1);k++)
      {
          fprintf(outf, "{%5.2f} %5.2f, 0.00, 0.00\n",
             ((((float)k-11.0)/20.0)+0.05), 100*(((((float)k-11.0)/20.0)+0.05)) );
      }
      for(k=90; k<=180; k=k+10)
      {
         fprintf(outf, "{%5.2f} -70.0, %5.2f, 0.00\n",(float)k,(float)(k-90));
      }
      if(sizeof(inputfilename) < nendstr) {nendstr = sizeof(inputfilename);}
      for(k=0; k<nendstr; k++) {strcName[k]=inputfilename[k];}

      /*Note first nemdstr char of inputfilename 20120120dcr*/
      fprintf(outf, "@group {%8s} dominant animate\n",strcName);
      fprintf(outf, "@dotlist {gapangle}\n");

      for (i = 0; i < NUMATOMTYPES; i++)
      {/*i: Hatom varients + rest of periodic table + base types, see atomprops.h*/
       for (j = 0; j < NODEWIDTH; j++)
       {/*Old: j: wide contact, close contact, small overlap, bad overlap, H-bonds*/ /* see NODEWIDTH in probe.h for new def. SJ 04/16/2015*/
              /*optionally: H-weaks if LweakHbonds  111215dcr*/
        node = results[i][j]; /*head of list of dot nodes, i.e. first node*/
        if (node)
        {/*node exists*/

         while(node)
         {/*kinemage point or point-line for each node in (ij)th dot-node list*/
           a = node->a; /*owner of dot*/
           t = node->t; /*cause of dot, use for OneDotEach 111013dcr*/
           for(k=0; k<10; k++){if(a->bval < k*5) {break;}}
           /*2:10A,3:15B,4:20C,5:25D,6:30E,8:40F,10:50G*/
           if(k<=2){ptm[0] = 'A';}
           else if(k<=3){ptm[0] = 'B';}
           else if(k<=4){ptm[0] = 'C';}
           else if(k<=5){ptm[0] = 'D';}
           else if(k<=6){ptm[0] = 'E';}
           else if(k<=2){ptm[0] = 'F';}
           else if(k<=2){ptm[0] = 'G';}
           else{ptm[0] = 'H';}
           /*Note first nemdstr char of inputfilename 20120120dcr*/
           fprintf(outf,"{%6s %s%c%s%s%c%s:O=%4.2f, B=%6.2f; %s%c%s%s%c%s:O=%4.2f, B=%6.2f; GAP %6.3f}'%s' %8.3f, %8.3f , %8.3f\n",
                  strcName,
                  a->atomname, a->altConf,
                  a->r->resname, a->r->Hy36resno,
                  a->r->resInsCode, a->r->chain, a->occ, a->bval,
                  t->atomname, t->altConf,
                  t->r->resname, t->r->Hy36resno,
                  t->r->resInsCode, t->r->chain, t->occ, t->bval, node->gap, ptm,
                  100*(node->gap), (node->angle)-90.0, a->bval);

           node = node->next; /*in this (ij)th dot-node list*/
         }/*kinemage point or point-line for each node in (ij)th dot-node list*/
         fprintf(outf, "\n");
        }/*node exists*/
       }/*NODEWIDTH*/
      }/*NUMATOMTYPES*/
     }/*LXHvector*/
  }/*Lgapbins TRUE*/

}/*writeOutput for kinemage*/
/*}}}writeOutput()___________________________________________________________*/

/*{{{writeAltFmtO()***********************************************************/
void writeAltFmtO(FILE *outf, int showBegin, int showEnd,
      char* groupname, dotNode *results[][NODEWIDTH], int spike)
{
   int i, j, numGroups, gn;
   dotNode *node;
   char *color = "", *prevColor = "";
   //Old: char *mast[NODEWIDTH] = {"WC", "CC", "SO", "BO", "HB"};
   char *mast[NODEWIDTH] = {"WC", "CC", "WH", "SO", "BO", "WO", "HB"}; /* 04/15/2015 SJ changed to match new NODEWIDTH, see probe.h*/
    
   numGroups = (showBegin && showEnd) ? 1 : 2;
   gn = (!showBegin && showEnd) ? 2 : 1;

   for (j = 0; j < NODEWIDTH; j++) {
      fprintf(outf, "begin_object %s%d\n", mast[j], gn);
      fprintf(outf, "mode solid\n");

      for (i = 0; i < NUMATOMTYPES; i++) {
	 prevColor = "--none--";
	 node = results[i][j];
	 if (node) {
	    color = convertKinColorToXV(getColor(i));
	 }
	 while(node) {
	    if (ColorGap) {
	       if (node->t) {
		  color = assignGapColorForO(node->gap, j);
	       }
	       else { color = "124";}
	    }
	    if (color != prevColor) {
	       fprintf(outf, "colour %s\n", color);
	       prevColor = color;
	    }
	    if ((j == 3 || j == 4 || j ==5) && spike) {/* bump */ /* 04/15/2015 SJ changed j numbers to match new ones for so, bo, and wo*/
	       fprintf(outf,
		  "move %.3f %.3f %.3f\nline %.3f %.3f %.3f\n",
		  node->loc.x, node->loc.y, node->loc.z,
		  node->spike.x, node->spike.y, node->spike.z);
	    }
	    else {/* contact or H-bond */
	       fprintf(outf, "sphere_xyz %.3f %.3f %.3f 0.03\n",
		  node->loc.x, node->loc.y, node->loc.z);
	    }
	    node = node->next;
	 }
      }
      fprintf(outf, "end_object\n");
   }

   if (showEnd) {
      if (showBegin) { fprintf(outf, "begin_object %s\n", groupname); }
      else           { fprintf(outf, "begin_object contsurf\n"); }
      for (j = 0; j < NODEWIDTH; j++) {
	 for (i = 0; i < numGroups; i++) {
	    fprintf(outf, "instance %s%d\n", mast[j], i+1);
	 }
      }
      fprintf(outf, "end_object\n");
   }
}
/*}}}writeAltFmtO()__________________________________________________________*/

/*{{{writeAltFmtXV()**********************************************************/
void writeAltFmtXV(FILE *outf, int showBegin, int showEnd,
      char* groupname, dotNode *results[][NODEWIDTH], int spike)
{
   int i, j;
   dotNode *node;
   char *color = "";
   //Old: char *mast[NODEWIDTH] = {"wide contact", "close contact", "small overlap", "bad overlap", "H-bonds"};
    /* 04/15/2015 SJ changed to reflect new NODEWIDTH, see probe.h*/
   char *mast[NODEWIDTH] = {"wide contact", "close contact", "weak H-bonds", "small overlap", "bad overlap", "worse overlap", "H-bonds"};
    
   if (showBegin) {
      if (showEnd) { fprintf(outf, "# begin %s\n", groupname); }
      else         { fprintf(outf, "# begin object\n#  group: %s\n", groupname); }
   }
   else            { fprintf(outf, "#  group: %s\n", groupname); }

   for (j = 0; j < NODEWIDTH; j++) {
      for (i = 0; i < NUMATOMTYPES; i++) {
	 node = results[i][j];
	 if (node) {
	    fprintf(outf, "#  (%s %s)\n", getAtomName(i), mast[j]);
	    color = convertKinColorToXV(getColor(i));
	 }
	 while(node) {
	    if (ColorGap) {
	       if (node->t) {
		  color = assignGapColorForXV(node->gap, j);
	       }
	       else { color = "124";}
	    }
	    if ((j == 3 || j == 4 || j == 5) && spike) {/* bump */ /* 04/15/2015 SJ changed j numbers to match new ones for so, bo, and wo*/
	       fprintf(outf,
		  "%.3f %.3f %.3f %.3f %.3f %.3f %s\n",
		  node->loc.x, node->loc.y, node->loc.z,
		  node->spike.x, node->spike.y, node->spike.z, color);
	    }
	    else {/* contact or H-bond */
	       fprintf(outf, "%.3f %.3f %.3f %.3f %.3f %.3f %s\n",
		  node->loc.x, node->loc.y, node->loc.z,
		  node->loc.x, node->loc.y, node->loc.z, color);
	    }
	    node = node->next;
	 }
      }
   }
   if (showEnd) {
      if (showBegin) { fprintf(outf, "# end %s\n", groupname); }
      else           { fprintf(outf, "#  endgroup: %s\n# end object\n", groupname); }
   }
   else              { fprintf(outf, "#  endgroup: %s\n", groupname); }
}
/*}}}writeAltFmtXV()_________________________________________________________*/

/*{{{dot2bullsEye()***********************************************************/
float dot2bullsEye(point3d *dot, atom *src, atom *targ)
{
        point3d targ2srcVec, targSurfacePoint;

        v3sub(&(src->loc), &(targ->loc), &targ2srcVec);
        v3scale(&targ2srcVec, targ->radius);
        v3add(&targ2srcVec, &(targ->loc), &targSurfacePoint);
 	return v3distance(dot, &targSurfacePoint);
}
/*}}}dot2bullsEye()__________________________________________________________*/

/*{{{dot2srcCenter()**********************************************************/
float dot2srcCenter(point3d *dot, atom *src, atom *targ)
{
        point3d src2targVec, srcSurfacePoint;

        v3sub(&(targ->loc), &(src->loc), &src2targVec);
        v3scale(&src2targVec, src->radius);
        v3add(&src2targVec, &(src->loc), &srcSurfacePoint);
 	return v3distance(dot, &srcSurfacePoint);
}
/*}}}dot2srcCenter()_________________________________________________________*/

/*{{{kissEdge2bullsEye()******************************************************/
float kissEdge2bullsEye(float ra, float rb, float rp)
{
 	return 2.0*ra*sqrt(rb*rp/((ra+rb)*(ra+rp)));
}
/*}}}kissEdge2bullsEye()_____________________________________________________*/

/*{{{writeRaw()***************************************************************/
void writeRaw(FILE *outf, char* groupname, dotNode *results[][NODEWIDTH],
	       float probeRad, char* label, float density, int conFlag) //conFlag added SJ 10/07/2011
{
   int i, j;
   dotNode *node;
   atom *a, *t;
  // char *mast[NODEWIDTH] = {"wc", "cc", "so", "bo", "hb"};
   char *mast[NODEWIDTH] = {"wc","cc","wh","so","bo","wo","hb"}; /* 04/10/2015 SJ to changed to match new categories of weak H bonds and worse overlap (these lists are made only if LweakHbonds and LworseOverlap are true respectively). see NODEWIDTH probe.h*/
   
   float gap, sl, dtgp, ke2be, d2be, d2sc, score;
   double scaledGap;
   for (i = 0; i < NUMATOMTYPES; i++)
   {/*loop over NUMATOMTYPES*/

      for (j = 0; j < NODEWIDTH; j++)
      {/*loop over NODEWIDTH*/

         node = results[i][j];
         if(node) // added 10/05/11 - SJ for condensing rawOutput
         {
             node = Condense(node,conFlag);
             results[i][j]=node;
         }
         while(node)
         {/*while node...*/
             fprintf(outf, "%s:%s:%s:", label, groupname, mast[j]);
             a = node->a;
             t = node->t;
             
             fprintf(outf, "%s%4.4s%c%s %s%c:",
                     a->r->chain, a->r->Hy36resno, a->r->resInsCode,
                     a->r->resname,a->atomname, a->altConf);
             
             if(t)
             {/*t node exists*/
                 fprintf(outf, "%s%4.4s%c%s %s%c:",
                         t->r->chain, t->r->Hy36resno, t->r->resInsCode,
                         t->r->resname,t->atomname, t->altConf);
                 gap = gapSize(&(a->loc), &(t->loc),(a->radius + t->radius));
                 dtgp = node->gap;
                 if (OldStyleU)
                 {
                     ke2be = kissEdge2bullsEye(a->radius, t->radius, probeRad);
                     d2be = dot2bullsEye(&(node->loc), a, t);
                     d2sc = dot2srcCenter(&(node->loc), a, t);
                 }
                 sl = gapSize(&(node->loc), &(node->spike), 0.0);
                 score = 0.0;
                 switch(j) /* SJ 04/10/2015 changed to include new categories*/
                 {
                    case 0:
                    case 1:
                           scaledGap = dtgp/GAPweight;
                           score = exp(-scaledGap*scaledGap);
                           break;
                    case 2: /* TODO: SJ weak H bonds, don't know what to do here, because they can be both wc and cc, so will have to check*/
                    case 3: /* small overlap, doing nothing, as before*/
                    case 4: score = -BUMPweight * sl; break;
                    case 5: score = -BUMPweight * sl; break; /* worse overlap, same as bad overlap*/
                    case 6: score =    HBweight * sl; break;
                 }
                 if(conFlag)
                 {
                      fprintf(outf,"%i:",node->dotCount);
                 }
                 else if(LOneDotEach)  /*dcr20120120*/
                 {
                      fprintf(outf,"%.2f:",node->angle);
                 }
                 if (OldStyleU)
                 {
                     fprintf(outf, "%.3f:%.3f:%.3f:%.3f:%.3f:%.3f:%.4f",
                             gap, dtgp, ke2be, d2be, d2sc, sl, score/density);
                 }
                 else
                 {
                     fprintf(outf, "%.3f:%.3f:%.3f:%.3f:%.3f:%.3f:%.4f",
                             gap, dtgp, node->spike.x, node->spike.y,
                             node->spike.z, sl, score/density);
                 }
             }/*t node exists*/
             else
             {
                 fprintf(outf, ":::::::");
             }

             fprintf(outf,":%s:%s:%.3f:%.3f:%.3f",
                         getAtomName(i),(t?getAtomName(t->atomclass):""),
                         node->loc.x,node->loc.y,node->loc.z);

             if(t)
             {
                 fprintf(outf, ":%.2f:%.2f\n", a->bval, t->bval);
             }
             else
             {
                 fprintf(outf, ":%.2f:\n", a->bval);
             }
             node = node->next;
         }/*while node...*/
      }/*loop over NODEWIDTH*/
   }/*loop over NUMATOMTYPES*/
}
/*}}}writeRaw()______________________________________________________________*/

/*{{{Condense()****(SortdotNodes)***** - SJ 10/05/11 for condensing rawOutput*/
dotNode * Condense(dotNode * head, int conFlag)
{

	// variables of the main for loop
	char sourceCurStart[20],sourceCurEnd[20];
	dotNode *curStart = NULL, *curEnd = NULL, *prevEnd = NULL, *nextStart = NULL, *prev = NULL;
	atom *source = NULL;

	//variables for the sorting and removing duplciates loop
	dotNode *i = NULL, *j = NULL, *next = NULL, *iPrev = NULL, *jPrev = NULL, *minPrev = NULL, *minNode = NULL;
	char minTarget[20],jTarget[20],prevTarget[20],nextTarget[20];
	atom *target = NULL;

	for(curStart=head;curStart!=NULL;curStart=nextStart)
	{ /* main for loop which loops through the source atoms*/
		source=curStart->a;
		sprintf(sourceCurStart,"%s%4.4s%c%s %s%c",source->r->chain,source->r->Hy36resno,
			source->r->resInsCode,source->r->resname,source->atomname,source->altConf);
		//sourceCurStart[17]='\0';
		nextStart = NULL;
		prev = curStart;
		curEnd = NULL;
		for(curEnd=curStart->next;curEnd!=NULL;curEnd=curEnd->next)
		{ /* for loop to determine the nodes for the same source atom*/
			source=curEnd->a;
			sprintf(sourceCurEnd,"%s%4.4s%c%s %s%c",source->r->chain,source->r->Hy36resno,
				source->r->resInsCode,source->r->resname,source->atomname,source->altConf);
			//sourceCurEnd[17]='\0';

			if(strcmp(sourceCurStart,sourceCurEnd)!=0) // if the source atom differes, then break out of the loop
				break;
			prev=curEnd;
		}

		nextStart=curEnd; // the start of the next source atoms
		curEnd=prev; // the last dot for the current source atom
		// isolating all the dots of the current source atom from the prev and next source atom
		if(prevEnd)
			prevEnd->next = NULL;
		curEnd->next = NULL;

		/* sorting code - will have a curState and a curEnd pointed to indicate the start and end of the sorted list*/
		iPrev = NULL;
		jPrev = NULL;
		minPrev = NULL;
		i = NULL;
		j = NULL;
		next = NULL;
		minNode = NULL;
		for(i=curStart;i!=NULL;i=i->next)
		{ /*sorting loop - sorting based on the target atom*/
			target=i->t;
			sprintf(minTarget,"%s%4.4s%c%s %s%c",target->r->chain,target->r->Hy36resno,
				target->r->resInsCode,target->r->resname,target->atomname,target->altConf);
			//minTarget[17]='\0';
			minNode = i;
			minPrev = iPrev;
			//using selection sort
			 jPrev = i;
			 for(j=i->next;j!=NULL; j=j->next)
			 {
			 	 target=j->t;
			 	 sprintf(jTarget,"%s%4.4s%c%s %s%c",target->r->chain,target->r->Hy36resno,
			 	 	target->r->resInsCode,target->r->resname,target->atomname,target->altConf);
			 	 //jTarget[17]='\0';
			 	 if(strcmp(jTarget,minTarget) < 0) // if the target atom of j is less than that of the min then swap it
			 	 {
			 	 	 strcpy(minTarget,jTarget);
			 	 	 minNode = j;
			 	 	 minPrev = jPrev;
			 	 }
			 	 jPrev = j;
			 }

			 if(minNode != i) // if the nodes have to be swapped
			 {
			 	 next = minNode->next;
			 	 if(iPrev)
			 	 	 iPrev->next = minPrev -> next;
			 	 else
			 	 	 curStart=minNode;
			 	 minPrev->next = i;
			 	 minNode->next = i->next;
			 	 i->next = next;
			 	 i = minNode;
			 }
			 iPrev = i;
		}/* sorting loop ends */
		curEnd = iPrev; //updating the end of the dots for the source atom

		if(conFlag)
		{/* code for removing duplicates - will update curEnd */
			i=curStart;
			next = NULL;
			i->dotCount=1; //initialising the dot count

			while(i->next!=NULL)
			{//removing duplicates
				target=i->t;
				sprintf(prevTarget,"%s%4.4s%c%s %s%c",target->r->chain,target->r->Hy36resno,
					target->r->resInsCode,target->r->resname,target->atomname,target->altConf);
				//prevTarget[17]='\0';

				target=i->next->t;
				sprintf(nextTarget,"%s%4.4s%c%s %s%c",target->r->chain,target->r->Hy36resno,
					target->r->resInsCode,target->r->resname,target->atomname,target->altConf);
				//nextTarget[17]='\0';

				if(strcmp(prevTarget,nextTarget) == 0) //free the memory for the duplicate dot
				{
					next = i->next;
					i->next = next->next;
					free(next);
					next=NULL;
					i->dotCount++;
				}
				else
				{
					i=i->next;
					i->dotCount=1;
				}
			}//removing duplicate loop ends
			curEnd = i; //updating the end of the dots of the source atom
		}/*code for removing dots ends*/

		//integrating the condensed list of dots of the current source atom to prev and next source atoms
		curEnd->next = nextStart;
		if(prevEnd)
			prevEnd->next = curStart;
		else
			head=curStart;

		prevEnd=curEnd; //making the current source atom as the prev source atom

	}/* main for loop*/

	return head;
}
/*}}}SortdotNodes()____________________________________________________________*/
/*{{{enumerate()**************************************************************/
void enumerate(FILE *outf, char* groupname, dotNode *results[][NODEWIDTH],
               float probeRad, int method,
	       int nsel, int spike, int outdots, int numSkinDots,
	       float density)
{
   int i, j, doit;
   float gs, hs, bs, hslen, bslen, tgs, ths, tbs, thslen, tbslen, psas, tsas;
   float dtgp, score, tGscore, tHscore, tBscore, tscore, scoreValue, a_radius;
   double scaledGap, slen;
   dotNode *node;
  /* char *mast[NODEWIDTH] = {"wide_contact  ", "close_contact ",
                            "small_overlap ", "bad_overlap   ", "H-bond        "};*/
    /* 04/10/2015 SJ changed to match the new categories, see NODEWIDTH in probe.h*/
    char *mast[NODEWIDTH] = {"wide_contact","close_contact","weak_H-bond","small_overlap","bad_overlap","worse_overlap","H-bond"};

   /* psas and tsas are the partial and total solvent accessible surface */

   fprintf(outf, "        \nsubgroup: %s\n", groupname);
   fprintf(outf, "atoms selected: %d\npotential dots: %d\npotential area: %.1f A^2\n",
				 nsel, numSkinDots, numSkinDots/density);

   if (nsel <= 0 || numSkinDots <= 0) {
      fprintf(outf, "empty selection\n");
      return;
   }
   if (spike) {
      fprintf(outf, "  type                 #      %%       score score/A^2 x 1000\n");
   }
   else {
      fprintf(outf, "  type                 #      %%\n");
   }

   tgs = ths = thslen = tbs = tbslen = tsas = 0.0;
   tGscore = tHscore = tBscore = tscore = 0.0;
   for (i = 0; i < NUMATOMTYPES; i++)
   {
       for (j = 0; j < NODEWIDTH; j++)
       {
           gs = hs = hslen = bs = bslen = score = psas = 0.0;
           node = results[i][j];
           doit = (node != NULL);
           if (doit) {
               fprintf(outf, "%3s %s ", getAtomName(i),
	                   outdots?"external_dots ":mast[j]);
           }
           while(node)
           {
               if (spike)/* 04/10/2015 SJ changed the if statements in this block to reflect change in dots ordering, see NODEWIDTH in probe.h*/
               {
                  // if (j == 0 || j == 1) { /* contact dot */
                   if(j==0 || j ==1 || j ==2) { /* contact dot, and weak H bond*/
                       gs += 1.0;
                       dtgp = node->gap;
                       scaledGap = dtgp/GAPweight;
                       scoreValue = exp(-scaledGap*scaledGap);
                       score   += scoreValue;
                       tGscore += scoreValue;
                   }
                  // else if (j == 2 || j == 3) { /* bump */
                   else if (j ==3 || j ==4 || j ==5){ /* bump*/
                       bs += 1.0;
                       slen = 0.5*FABS(node->gap);
                       bslen += slen;
                       scoreValue = - BUMPweight * slen;
                       score   += scoreValue;
                       tBscore += scoreValue;
                   }
                   else { /* H-bond */
                       hs += 1.0;
                       slen = 0.5*FABS(node->gap);
                       hslen += slen;
                       scoreValue = HBweight * slen;
                       score   += scoreValue;
                       tHscore += scoreValue;
                   }
               }
               else { gs += 1.0; }

               if (method == EXTERNALSURFACE) {
                   a_radius = node->a->radius;
                   psas += (a_radius + probeRad)*(a_radius + probeRad)/(a_radius * a_radius);
               }
               node = node->next;
           }
           if (doit)/* 04/10/2015 SJ changed the if statements in this block to reflect change in dots ordering, see NODEWIDTH in probe.h*/
           {
               if (spike) {
                  // if (j == 0 || j == 1) { /* contact dot */
                   if(j==0 || j ==1 || j ==2) { /* contact dot, and weak H bond*/
                        fprintf(outf, "%7.0f %5.1f%% %9.1f %9.2f\n",
                                gs, 100.0*gs/numSkinDots, score/density,
                                1000.0*score/numSkinDots);
                   }
                  // else if (j == 2 || j == 3) { /* bump */
                   else if (j ==3 || j ==4 || j ==5){ /* bump*/
                        fprintf(outf, "%7.0f %5.1f%% %9.1f %9.2f\n",
                                bs, 100.0*bs/numSkinDots, score/density,
                                1000.0*score/numSkinDots);
                   }
                   else { /* H-bond */
                        fprintf(outf, "%7.0f %5.1f%% %9.1f %9.2f\n",
                                hs, 100.0*hs/numSkinDots, score/density,
                                1000.0*score/numSkinDots);
                   }
               }
               else {
                   fprintf(outf, "%7.0f %5.1f%%\n",
                           gs, 100.0*gs/numSkinDots);
               }
               tgs += gs;
               ths += hs;
               thslen += hslen;
               tbs += bs;
               tbslen += bslen;
               tscore += score;
               if (method == EXTERNALSURFACE) {
                   tsas += psas; /* tally the solvent accessible surface */
               }
           }
       }
   }
   if (spike) {
/*      fprintf(outf, "                                   score/A^2 x 1000\n");*/
      fprintf(outf, "\n     tot contact:  %7.0f %5.1f%% %9.1f %9.2f\n",
		    tgs, 100.0*tgs/numSkinDots, tGscore/density,
		     1000.0*tGscore/numSkinDots);
      fprintf(outf,   "     tot overlap:  %7.0f %5.1f%% %9.1f %9.2f\n",
		    tbs,    100.0*tbs/numSkinDots, tBscore/density,
		     1000.0*tBscore/numSkinDots);
      fprintf(outf,   "     tot  H-bond:  %7.0f %5.1f%% %9.1f %9.2f\n",
		    ths,    100.0*ths/numSkinDots, tHscore/density,
		     1000.0*tHscore/numSkinDots);

      fprintf(outf, "\n       grand tot:  %7.0f %5.1f%% %9.1f %9.2f\n",
		    (tgs+tbs+ths),
		    100.0*(tgs+tbs+ths)/numSkinDots,
		    tscore/density, 1000.0*tscore/numSkinDots);
      fprintf(outf, "\ncontact surface area: %.1f A^2\n",
		     (tgs+tbs+ths)/density);
   }
   else {
      fprintf(outf, "             tot:  %7.0f %5.1f%%\n\n",
		     tgs, 100.0*tgs/numSkinDots);
      fprintf(outf, "   contact surface area: %.1f A^2\n",
		     tgs/density);
      if (method == EXTERNALSURFACE) {
	 fprintf(outf, "accessible surface area: %.1f A^2\n\n",
		     tsas/density);
      }
   }
}
/*}}}enumerate()_____________________________________________________________*/

/*{{{rawEnumerate()***********************************************************/
void rawEnumerate(FILE *outf, char* groupname, dotNode *results[][NODEWIDTH],
               int method, int nsel, int spike, int outdots, int numSkinDots,
	       float density, char *namestring, char *rawname, double scoreBias)
{/*rawEnumerate*/
   int i, j, doit;
   float gs, hs, bs, hslen, bslen, tgs, ths, tbs, thslen, tbslen;
   float dtgp, score, tGscore, tHscore, tBscore, tscore, scoreValue;
   double scaledGap, slen;
   dotNode *node;

   /*autobondrot: countDots==1, rawOutput==1  */
   /*  *rawname holds autobondrot angle values as space_char deliniated string*/

   if (nsel <= 0 || numSkinDots <= 0) { /* empty selection */
      if (method == EXTERNALSURFACE) {
	 fprintf(outf, "%9.3f", 0.0);
      }
      else if (spike) {
	 fprintf(outf, "%9.3f", scoreBias);
      }

      if (*rawname) { fprintf(outf, " %s", rawname); }
      if (*namestring || *groupname) {
	 fprintf(outf, "%s", RAW_HEADER_COMMENT);
	 if (*namestring) { fprintf(outf, " %s", namestring); }
	 if (*groupname)  { fprintf(outf, " %s", groupname);  }
      }
      fprintf(outf, "\n");
      return;
   }

   tgs = ths = thslen = tbs = tbslen = 0.0;
   tGscore = tHscore = tBscore = tscore = 0.0;
   for (i = 0; i < NUMATOMTYPES; i++) {
   for (j = 0; j < NODEWIDTH; j++)
   {
      gs = hs = hslen = bs = bslen = score = 0.0;
      node = results[i][j];
      doit = (node != NULL);
      while(node)
      {
         if (spike)/* 04/10/2015 SJ changed the if statements in this block to reflect change in dots ordering, see NODEWIDTH in probe.h*/
         {
            // if (j == 0 || j == 1) { /* contact dot */
             if(j ==0 || j ==1 || j ==2){ /* contacts, and weak H bonds*/
                 gs += 1.0;
                 dtgp = node->gap;
                 scaledGap = dtgp/GAPweight;
                 scoreValue = exp(-scaledGap*scaledGap);
                 score   += scoreValue;
                 tGscore += scoreValue;
             }
            // else if (j == 2 || j == 3) { /* bump */
             else if(j ==3 || j ==4 || j ==5){ /* bump */
                 bs += 1.0;
                 slen = 0.5*FABS(node->gap);
                 bslen += slen;
                 scoreValue = - BUMPweight * slen;
                 score   += scoreValue;
                 tBscore += scoreValue;
             }
             else { /* H-bond */
                 hs += 1.0;
                 slen = 0.5*FABS(node->gap);
                 hslen += slen;
                 scoreValue = HBweight * slen;
                 score   += scoreValue;
                 tHscore += scoreValue;
             }
         }
         else {
             gs += 1.0;
         }
         node = node->next;
      }
      if (doit) {
         tgs += gs;
         ths += hs;
         thslen += hslen;
         tbs += bs;
         tbslen += bslen;
         tscore += score;
      }
   }
   }

   /*output one line of information */
   if (method == EXTERNALSURFACE) {
      fprintf(outf, "%9.3f", (tgs+tbs+ths)/density);
   }
   else if (spike) /*autobondrot: spike==1 at least in all examples dcr tested*/
   {
      fprintf(outf, "%9.3f", scoreBias + (tscore/density));  /*value # ...*/
   }
   else {
      fprintf(outf, "%9.3f", scoreBias + tgs);
   }
   if (*rawname) { fprintf(outf, " %s", rawname); } /*autobondrot angle values*/
   if (*namestring || *groupname) {
      fprintf(outf, "%s", RAW_HEADER_COMMENT);
      if (*namestring) { fprintf(outf, " %s", namestring); }
      if (*groupname)  { fprintf(outf, " %s", groupname);  }
   }
   fprintf(outf, "\n"); /* NOTE this is the LF for above print statements*/

   /* char* rawname == buf that is filled in autobondrot.c/runnameAndTorsion()*/
   /* contains the autobondrot angle values in the order in which processed*/

}/*rawEnumerate*/
/*}}}rawEnumerate()__________________________________________________________*/

/*{{{countSelected()**********************************************************/
int countSelected(atom *theAtoms, int srcFlag)
{
   atom *a = NULL;
   int ns = 0;
   for(a = theAtoms; a; a = a->next) {
      if (a->flags & srcFlag) { ns++; }
   }
   return ns;
}
/*}}}countSelected()_________________________________________________________*/

/*{{{enumDotSkin()************************************************************/

/* NOTE: allMainAtoms/abins & allMovingAtoms/bbins must be disjoint */
/*       sets of atoms (none in common) or allMovingAtoms/bbins can */
/*       be NULL.                                                   */
/*       allMovingAtoms refers to autobondrot set of atoms          */

int enumDotSkin(atom *allMainAtoms, atomBins *abins,
	    atom *allMovingAtoms, atomBins *bbins,
	    pointSet dots[], int srcFlag)
{
   atom *src = NULL, *atomList = NULL, *atomList2 = NULL;
   int dotTotal = 0, usesMovingAtoms = FALSE;
   int dummy = TRUE;

         /*numSkinDots used to normalize output score*/

   usesMovingAtoms = ((allMovingAtoms != NULL) && (bbins != NULL));

   for(src = allMainAtoms; src; src = src->next) {
      if (src->flags & srcFlag) {

	 atomList  = findTouchingAtoms(src, NULL,     abins, 0.0, 0, &dummy);
	 if (usesMovingAtoms) {
	    atomList2 = findTouchingAtoms(src, atomList, bbins, 0.0, 0, &dummy);
	 }
	 else { atomList2 = atomList; }

	 if (atomList2) {
	    markBonds(src, atomList2, 1, Maxbonded); /*in enumDotSkin()*/
	    if (Maxbonded > 3) { fixupLongBondChains(src, atomList2, 3); }
	 }
	 dotTotal += countSkin(src, atomList2, dots);

	 if(Verbose && ShowTicks) {
	    fprintf(stderr, "%s%s   \r",
	               src->r->resname, src->r->Hy36resno);
	 }
      }
   }

   if (usesMovingAtoms) {
      for(src = allMovingAtoms; src; src = src->next) {
	 if (src->flags & srcFlag) {

	    atomList  = findTouchingAtoms(src, NULL,     abins, 0.0, 0, &dummy);
	    atomList2 = findTouchingAtoms(src, atomList, bbins, 0.0, 0, &dummy);
	    if (atomList2) {
	       markBonds(src, atomList2, 1, Maxbonded); /*in enumDotSkin()*/
	       if (Maxbonded > 3) { fixupLongBondChains(src, atomList2, 3); }
	    }
	    dotTotal += countSkin(src, atomList2, dots);

	    if(Verbose && ShowTicks) {
	       fprintf(stderr, "%s%s   \r",
	               src->r->resname, src->r->Hy36resno);
	    }
	 }
      }
   }
   return dotTotal;
}
/*}}}enumDotSkin()___________________________________________________________*/

/*{{{countSkin()**************************************************************/
int countSkin(atom *src, atom *scratch, pointSet dots[])
{
   atom *targ = NULL;
   int i = 0, within = 0, ok = 0;
   point3d dotloc, dotvect;
   pointSet *srcDots;
   int dotCnt = 0;

   if (src->elem == ignoreAtom
    || src->elem == atomHOd)   { return 0; } /*hb-only-dummy, phantom H atom*/

   srcDots = &(dots[src->elem]);
   if (src->elem == atomC && isCarbonylAtom(src)) {
      srcDots = &COdots;
   }

   for(i=0; i < srcDots->n; i++) {
      dotvect = srcDots->p[i];

      v3add(&(src->loc), &dotvect, &dotloc);

      ok = TRUE;

      /*targ is an atom*, scratch is an atom*,  atom* defined in abin.h */
      /* atom* scratch is a member of atom*: which has member scratch...*/
      /* so for-loop moves through all these atoms while an atom has a scratch*/

      for(targ = scratch; targ; targ = targ->scratch) {
	 if (targ->elem == ignoreAtom) {continue;}

	 /* eliminate if within bonded atom! */
	 if (targ->mark && (targ->elem != atomHOd)) {
#ifdef INLINE_FOR_SPEED
	    within = INRANGEINLINE(dotloc, (targ->loc), targ->radius);
#else
	    within = inRange(&dotloc, &(targ->loc), targ->radius);
#endif

	    if (within) { ok = FALSE; break; }
	 }
      }
      if (ok) { dotCnt++; }
   }
   return dotCnt;
}
/*}}}countSkin()_____________________________________________________________*/

/*{{{updateHydrogenInfo()*****************************************************/
/* Make polor hydrogens HB donors.                                  */
/* Waters without hydrogens are protonated with "phantoms".         */
/* If we are using moving atoms, all the new water H atoms go there */

/* can return a list of new atom "clones" which are just the MainAtoms waters */

/* NOTE: allMainAtoms/abins & allMovingAtoms/bbins must be disjoint */
/*       sets of atoms (none in common) or allMovingAtoms/bbins can */
/*       be NULL.                                                   */
/*       allMovingAtoms refers to autobondrot set of atoms          */

/*usual call: (i.e. NOT autobondrot)
      updateHydrogenInfo(outf,         allMainAtoms,        abins,
                                       NULL,                NULL,
                                       SET1|SET2,           FALSE);
*/
atom* updateHydrogenInfo(FILE *outf, atom *allMainAtoms,   atomBins *abins,
                    atom *allMovingAtoms, atomBins *bbins,
                    int selectedFlag, int mustSaveMainWater)
{/*updateHydrogenInfo*/
   atom *src = NULL, *orig = NULL, *atomList = NULL, *a = NULL;
   atom *newH = NULL, *tempStaticList = NULL, *mainWaterStorage = NULL;
   int type = 0, newHcount = 0, i = 0, whichList = 0;
   int nProx = 0, usesMovingAtoms = FALSE;
   int dummy = FALSE;

   /* Information specific to making sure only one H2O Hydrogen is     */
   /* generated which points into ring atoms on a given aromatic ring. */
   struct {atom  *a; float gap;} proximity[20];

   usesMovingAtoms = ((allMovingAtoms != NULL) && (bbins != NULL));
   if (! usesMovingAtoms) { allMovingAtoms = NULL; bbins = NULL; }

   if (DumpNewHO) { fprintf(outf,"@vectorlist {water H?} color= gray\n"); }

   tempStaticList = NULL; /* new hydrogens for MainAtoms when usesMovingAtoms */
   mainWaterStorage = NULL; /* MainAtoms waters when usesMovingAtoms          */

   for(whichList = 0; whichList <= 1; whichList++)
   { /* loop over both lists of atoms */
      for(src = ((whichList==0)?allMainAtoms:allMovingAtoms); src; src = src->next)
      {/*loop over src atoms*/

         /* Note: we have to examine ALL the atoms at this point.        */
         /* Not testing (src->flags & selectedFlag) for H will slow down */
         /* the code but it is neccessary to get the atom types correct. */

         if (isHatom(src->elem))
         { /* fixup Hs in the input */

            atomList = NULL;
            if (src->flags & selectedFlag)
            {
               atomList = findTouchingAtoms(src, NULL, abins, 0.0, 0, &dummy);
               if (usesMovingAtoms)
               {
                  atomList = findTouchingAtoms(src, atomList, bbins, 0.0, 0, &dummy);
               }
            }
            if (atomList)
            {
               markBonds(src, atomList, 1, 1); /*in updateHydrogenInfo()*/
               type = dotType(src, atomList, ByNABaseColor);

               if ( type == atomN || type == atomO || type == atomS )
               {/* bgn: if (atom types) */
                  if(Verbose && ShowTicks)
                  {
                     fprintf(stderr, "%s%d   \r",src->r->resname, src->r->resid);
                  }
                  src->props |= DONOR_PROP; /* polar hydrogens are donors */

                  if (UsePolarH)
                  {
                     src->elem = atomHpolar; /* very important */
                     src->radius = getRadius(src->elem, 0);
                  }
/* there are a more restrictive set of cutoffs in reduce: 0.7 & 40 */
#define H2OoccCUTTOFF   0.25
#define H2ObvalCUTTOFF 80.0
                  if (  (   FABS(src->occ) < H2OoccCUTTOFF
                         || src->bval >= H2ObvalCUTTOFF   )
                      &&(src->props & WATER_PROP)          )
                  {
                     src->elem = ignoreAtom;
                     src->radius = 0.0;
                  } /* ignore low occupancy High B water hydrogens */

                  for(a = atomList; a; a = a->scratch)
                  {
                     /* since we found a hydrogen, then we know */
                     /* that the heavy atom is not the donor    */
                     /* so we can disable that property         */
                     if(    (a->mark == 1)
                        && !(isHatom(a->elem))
                        &&  (FABS(src->occ) >= H2OoccCUTTOFF)
                        &&  (src->bval < H2ObvalCUTTOFF)     )
                     {

                        a->props &= ~DONOR_PROP; /*set to one's complement of ...*/
                     }
                  }
               }/* end: if (atom types) */
            }
         }/* fixup Hs in the input */
         else if (   (src->flags & selectedFlag)  /* fixup water oxygen */
                  && (src->props & WATER_PROP)
                  && (src->elem == atomO)       )
         {/* else if water */ /* fixup water oxygen */

            /*if connected H already found, the DONOR prop will not be set*/
            /*if it is, we look further ahead in this residue for Hs      */
            if (src->props & DONOR_PROP)
            {
               for(a = src->next; a && (src->r == a->r); a = a->next)
               {
                  if (  isHatom(a->elem)
                      && (src->altConf == a->altConf)
                      && (a->occ > 0.1))
                  {
                     src->props &= ~DONOR_PROP; /* turn off donor property */
                  }
               }
            }
            /* if donor flag still set...  */
            if (src->props & DONOR_PROP)
            {/* water O without occupied Hs */

               /* begin clone procedure */
               /* special case for when we want to remember just these waters */
               /* in this state                                               */
               if (mustSaveMainWater && usesMovingAtoms && (whichList==0))
               {
                  atom* waterClone = NULL;
                  waterClone = (atom *)malloc(sizeof(atom));
                  if (waterClone)
                  {
                     *waterClone = *src;
                     waterClone->next = mainWaterStorage; /* add clone to list */
                     mainWaterStorage = waterClone;
                  }
               }
               /* end clone procedure */

               src->props |= AMBIGWATER_PROP; /* this is a water which needs H? added */
               src->props |=  ACCEPTOR_PROP;
               src->props &= ~DONOR_PROP;    /* acceptor, not donor */
               orig = src; /* keep track of the original Oxygen atom */
               atomList = findTouchingAtoms(src, NULL, abins, 0.25, 0, &dummy);
               if (usesMovingAtoms)
               {
                  atomList = findTouchingAtoms(src, atomList, bbins, 0.25, 0, &dummy);
               }
               /* the 0.25 comes from:
                  OxygenVDW(1.4) + 2*probeRad(0.3) == HpolVDW(1.0) + OHBondLen(1.0)
                  and we want a minimum overlap of 0.1 before we consider possible Hs
               */

              /* for each nearby HB acceptor, keep only 1 per aromatic ring */
              nProx = 0;
              for(a = atomList; a; a = a->scratch)
              {
                 if (a->props & ACCEPTOR_PROP)
                 {/* bgn: if acceptor_prop */
                    float hbgap = gapSize(&(a->loc), &(orig->loc), a->radius + 2.0);
                    int   foundAromRingAtom = FALSE;

                    if (a->props & AROMATIC_PROP)
                    {
                       for (i=0; i < nProx; i++)
                       {
                          if (   a->r == proximity[i].a->r
                              && proximity[i].a->props & AROMATIC_PROP)
                          {

                             if (hbgap < proximity[i].gap)
                             {
                                proximity[i].a = a;
                                proximity[i].gap = hbgap;
                             }
                             foundAromRingAtom = TRUE;
                             break;
                          }
                       }
                    }
                    if (   (!foundAromRingAtom)
                        && nProx < sizeof(proximity)/sizeof(proximity[0]))
                    {
                       proximity[nProx].a = a;
                       proximity[nProx].gap  = hbgap;
                       ++nProx;
                    }
                 }/* end: if acceptor_prop */
              }

              /* build new HOH hydrogens in the direction of each acceptor */
              newHcount = 0;
              for (i=0; i < nProx; i++)
              {/*scan over neighbors*/

                 a = proximity[i].a;

                 if (a->props & ACCEPTOR_PROP)
                 {/*close acceptor! Make a phantom H in the direction of each Acceptor*/
                    newH = (atom *)malloc(sizeof(atom));
                    if (newH)
                    {/*newH allocated, set its parameters*/
                       point3d o2aVec;
                       double waterOHbondLen;

                       *newH = *orig;
                       if (usesMovingAtoms)
                       {
                          newH->next = tempStaticList; /* hook into temp list */
                          tempStaticList = newH;
                       }
                       else
                       {
                          newH->next = src->next; /* hook into list of atoms */
                          src->next  = newH;
                       }

                       v3sub(&(a->loc), &(orig->loc), &o2aVec);
#define BEST_HBOND_OVERLAP 0.6
                       waterOHbondLen = 1.0 +
                       MAX(-1.0, MIN(0.0, proximity[i].gap + BEST_HBOND_OVERLAP));
                       v3scale(&o2aVec, waterOHbondLen);
                       v3add(&(orig->loc), &o2aVec, &(newH->loc));

                       newH->nextInBin= NULL; /* added to bins below */
                       newH->scratch  = NULL;
                       newH->elem     = atomHOd; /*hb-only-dummy, phantom H atom*/
                       newH->radius   = getRadius(newH->elem, 0);
                       newH->covRad   = getCovRad(newH->elem);
                          /*atomprops.h/atomProp AtomTbl covRad*/
                       strcpy(newH->atomname, " H? ");

                       newH->props |=  DONOR_PROP;
                       newH->props &= ~ACCEPTOR_PROP; /* donor, not acceptor */

                       /* *** adding an atom not in input! *** */
                       addNeighbor(newH, (usesMovingAtoms ? bbins : abins));

                       if (DumpNewHO)
                       {
                          if (DebugLevel>3)
                          {
                             fprintf(stderr,
"HETATM%5d  H%-2d%c%3.3s%s%4.4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f      new  H\n",
                                0, ++newHcount, newH->altConf,
                                newH->r->resname, newH->r->chain,
                                newH->r->Hy36resno, newH->r->resInsCode,
                                newH->loc.x, newH->loc.y, newH->loc.z,
                                newH->occ, newH->bval);
                          }
                          fprintf(outf, "{%4.4s%c%3.3s%s%4.4s%c}P %8.3f%8.3f%8.3f\n",
                             orig->atomname, orig->altConf,
                             orig->r->resname, orig->r->chain,
                             orig->r->Hy36resno, orig->r->resInsCode,
                             orig->loc.x, orig->loc.y, orig->loc.z);
                          fprintf(outf, "{%4.4s%c%3.3s%s%4s%c}L %8.3f%8.3f%8.3f\n",
                             newH->atomname, newH->altConf,
                             newH->r->resname, newH->r->chain,
                             newH->r->Hy36resno, newH->r->resInsCode,
                             newH->loc.x, newH->loc.y, newH->loc.z);
                       }

                       if (! usesMovingAtoms)
                       {
                          src = src->next; /* bypass new atom when we iterate */
                                           /* when inserted after src atom    */
                       }
                    }/*newH allocated, set its parameters*/
                 }/*close acceptor! Make a phantom H in the direction of each Acceptor*/
              }/*scan over neighbors*/

            }/* water O without occupied Hs */
         }/* else if water */ /* fixup water oxygen */

         /* non-water polar (N/O/S) atoms in selection list */
         else if (    (src->flags & selectedFlag) && !(src->props & WATER_PROP)
                  &&  (src->elem == atomN || src->elem == atomO || src->elem == atomS))
         {
            /* Remove any DONOR property for non-water non-H            */
            /* polar atoms since the Hydrogen atoms will be the DONORs. */
            /* This was assigned in select.c to facilitate selections.  */

            src->props &= ~DONOR_PROP;
         }
      }/*loop over src atoms */ /* end loop over list of src atoms*/
   }/* loop over both lists of atoms */

   if (tempStaticList)
   {/* new hydrogens for MainAtoms when usesMovingAtoms */
      for(src = allMovingAtoms; src; src = src->next)
      {
         if (src->next == NULL)
         {
            src->next = tempStaticList;
            break; /* make sure we don't continually loop */
         }
      }
   }
   return mainWaterStorage;
}/*updateHydrogenInfo*/
/*}}}updateHydrogenInfo()____________________________________________________*/

/*{{{dump_changes()***********************************************************/
void dump_changes(FILE *outf)
{
fprintf(outf, "Probe change log:\n");
fprintf(outf, "\n");
fprintf(outf, "Note: Not captured prior to Aug 1998.\n");
fprintf(outf, "      JMW = J Michael Word\n");
fprintf(outf, "\n");
fprintf(outf, " 8/ 7/98 -        JMW - added segID parsing for XPLOR compatibility\n");
fprintf(outf, " 8/10/98 -        JMW - added formatting for O and XtalView output\n");
fprintf(outf, "11/12/98 -        JMW - use table to parent H atoms in std. residues,\n");
fprintf(outf, "                        made GAPweight expand with expanded radius,\n");
fprintf(outf, "                        and reorganized O format output in chunks\n");
fprintf(outf, "                        and changed xview format colors and modified\n");
fprintf(outf, "                        hydrogen name parsing for HG## and HE##\n");
fprintf(outf, "11/15/98 -        JMW - added -stdbonds and -noparent flags and built\n");
fprintf(outf, "                        a complete std AA and NA connectivity table\n");
fprintf(outf, "11/23/98 -        JMW - extended std bond table for truncated H names\n");
fprintf(outf, " 3/ 8/99 -        JMW - atom names can include * and ' characters\n");
fprintf(outf, " 3/11/99 -        JMW - fixed HBond problem with implicit hydrogens\n");
fprintf(outf, "                        and altered the skipping of MCMC to skip them\n");
fprintf(outf, "                        only if within a single chain and not a HET\n");
fprintf(outf, " 3/16/99 -        JMW - Made water donor/acceptor insensitive to names\n");
fprintf(outf, "                        and extended atom name parsing for wierd HETS\n");
fprintf(outf, " 4/05/99 -        JMW - Fixed typo in error checking in parse.c for\n");
fprintf(outf, "                        processing within _ of _,_,_\n");
fprintf(outf, " 4/06/99 -        JMW - Updated stdbond table for NT and OXT\n");
fprintf(outf, " 4/07/99 -        JMW - Cleaned up compiler warnings on initialization\n");
fprintf(outf, "                        and unused variables\n");
fprintf(outf, " 5/10/99 -        JMW - Added AIB,ABU and MSE to mc hbond NH & O list\n");
fprintf(outf, "                        and added HN to the list of aliases for H\n");
fprintf(outf, " 7/12/99 -        JMW - Begin changes for version 2\n");
fprintf(outf, " 7/17/99 -        JMW - Added support for autobondrot,  display ref.\n");
fprintf(outf, " 8/ 3/99 -        JMW - Fixed typo in UpdateHydrogenInfo about bins\n");
fprintf(outf, " 8/ 5/99 -        JMW - Sort autobondrot atoms in residue order,\n");
fprintf(outf, "                        add include file processing to autobondrot,\n");
fprintf(outf, "                        add -notick and -kinemage probe options\n");
fprintf(outf, "                        list process Nt/Ct/res for moving atoms\n");
fprintf(outf, " 8/16/99 -        JMW - Flag -kinemage forces kin format (for -auto)\n");
fprintf(outf, "                        added CA-CB for gly in stdbonds for mutants\n");
fprintf(outf, " 9/ 2/99 -        JMW - Added stdbond info to support alternate name\n");
fprintf(outf, "                        C5A for the methyl on a THY nucleic acid base\n");
fprintf(outf, "10/ 5/99 -        JMW - Dropped unused variables flagged by compiler\n");
fprintf(outf, "12/ 6/99 -        JMW - -OUT shows the contact surface area and the\n");
fprintf(outf, "                        solvent accessible surface area\n");
fprintf(outf, "                        -ADDVDW for radius offset(-SCALEVDW for scale)\n");
fprintf(outf, "                        -OUTColor sets the point color for -OUT\n");
fprintf(outf, "12/13/99 -        JMW - added shortcuts: -SCSurface, -CONTACT,\n");
fprintf(outf, "                                         -SASurface, -ACCESS\n");
fprintf(outf, "12/13/99 -        JMW - renamed: -SASurface to -ASurface\n");
fprintf(outf, " 1/ 3/00 -        JMW - added names to surface shortcuts\n");
fprintf(outf, " 1/13/00 -        JMW - renamed: -CONTACT to -EXPOsed\n");
fprintf(outf, " 1/14/00 -        JMW - fixed TRANS keyword in autobondrot\n");
fprintf(outf, " 5/ 4/00 -        JMW - changed -u output to print spike coords and\n");
fprintf(outf, "                        added -oldu flag to give the old output\n");
fprintf(outf, "                        added -ignore flag to drop data from input\n");
fprintf(outf, " 5/31/00 -        JMW - Added modified bases used in tRNAs\n");
fprintf(outf, "                        Included new properties TEST_ACCEPT_ANGLE and\n");
fprintf(outf, "                        CH_DONOR\n");
fprintf(outf, " 6/ 9/00 -        JMW - constructed chain of residues view of atomlist\n");
fprintf(outf, "                        in order to support ring normals\n");
fprintf(outf, " 7/ 5/00 -        JMW - inlined inRange for speed, remembered inosine,\n");
fprintf(outf, "                        dropped dummy atoms, added -segid flag,\n");
fprintf(outf, "                        fixed bug when water was used with -auto\n");
fprintf(outf, " 7/27/00 -        JMW - added code to free extra memory at end of run\n");
fprintf(outf, "10/31/00 -        JMW - added O2* & O2' to NABackbone list (oversite)\n");
fprintf(outf, "11/ 6/00 -        JMW - changed malloc.h and memory.h to stdlib.h\n");
fprintf(outf, "11/14/00 -        JMW - added -basecolor & -colorbase to color DNA/RNA\n");
fprintf(outf, "11/19/00 -        JMW - fixed naBase typo for 7MG and added -nofaceHB\n");
fprintf(outf, "11/27/00 -        JMW - updated properties of odd bases in select.c\n");
fprintf(outf, "11/28/00 -        JMW - added H2U to base class table and moved PSU\n");
fprintf(outf, "12/01/00 -        JMW - fixed H atomclass bug with -colorbase\n");
fprintf(outf, " 2/14/01 -        JMW - changed LowGoodCut from -0.2 to -0.4\n");
fprintf(outf, " 7/11/01 -(v2.4?) JMW - dots no longer broken out by element in .kin,\n");
fprintf(outf, "                        added -element to generate old style output\n");
fprintf(outf, "                        which made kinemage buttons for elements. \n");
fprintf(outf, " 7/25/01 -(v2.5?) JMW - negative residue numbers in patterns   \n");
fprintf(outf, "          (Sometimes need to include extra parentheses or a    \n");
fprintf(outf, "           space to prevent the selection from being treated as\n");
fprintf(outf, "           a command line option.)  \n");
fprintf(outf, "10/ 1/01 -        JMW - allow Nterminal fragnents to Hbond\n");
fprintf(outf, "10/ 9/01 -        JMW - take out Nterm code and make CHX HBs an option\n");
fprintf(outf, " 1/23/03 - 2.09 - JMW - Run default now:-3 -mc -het -self \"altA ogt33\"\n");
fprintf(outf, "                                    was:        -het -self \"altA\"\n");
fprintf(outf, "                        Default for -ignore is now \"olt1\"\n");
fprintf(outf, "                        Recognise insert codes like 123a, 2B-3, & insC\n");
fprintf(outf, "                        Fixed bounding box water-H bug in autobondrot!\n");
fprintf(outf, "10/14/03 - 2.10 - JMW - Fixed bug (found by Patrick Reardon) caused by\n");
fprintf(outf, "                        filtering out non-target contacts too early \n");
fprintf(outf, "                        Also made sure atoms in different models\n");
fprintf(outf, "                        can not interact in any way.\n");
fprintf(outf, "                        Also changed how alternate conformations in\n");
fprintf(outf, "                        in different residues interact. Now altA and\n");
fprintf(outf, "                        altB never see one another. Previously, they\n");
fprintf(outf, "                        could if alts A and B weren't the same res.\n");
fprintf(outf, "                        Added -minocc flag with a default of 0.02\n");
fprintf(outf, "                        as the minimum non-zero occupancy.\n");
fprintf(outf, "                        Added -changes flag to dump change log.\n");
fprintf(outf, "Oct 2004: DCR code annotations, fixups, default modifiations\n");
fprintf(outf, "10/07/04 - 2.10 dcr041007 now presumes atom named _SE_ is Selenium\n");
fprintf(outf, "                       i.e. refmac or cns missadjusted name.\n");
fprintf(outf, "halides flagged as negative ions and H-bond acceptors\n");
fprintf(outf, "default: probe infile > outfile \n");
fprintf(outf, "same as: probe  -3 -mc -het -self \"altA ogt33\" infile > outfile \n");
fprintf(outf, "try    : probe -4H -mc -het -self \"altA ogt33\" infile > outfile \n");
fprintf(outf,"dcr041009 pointmasters M==McMc, S==ScSc, O==Other(McSc,Het--)\n");
fprintf(outf,"dcr041010 -onlybadout option for Bad clashes only\n");
fprintf(outf,"dcr041017 default -4H, pointtmaster P==McSc, report counts\n");
fprintf(outf,"dcr041020 master= {surface}, more work on count reports\n");
fprintf(outf,"dcr041020 annotated (NOT implemented) writeHbonds JACKnANDREA\n");
fprintf(outf,"dcr041023 2D array of count information\n");
fprintf(outf,"dcr041026 NOT report count IntersectBothWays re probe update\n");
fprintf(outf,"dcr041101 -summary colon deliminated multi-count report\n");
fprintf(outf,"dcr041101 -oneline summary all on oneline\n");
fprintf(outf,"dcr041101 -DEFAULTs  same as: <<NO FLAGS>>, \n");
fprintf(outf,"         but allows flags like -summary or -oneline\n");

   /*jmw & dcr agreement of 041110 on version name and maintenance by dcr*/

fprintf(outf,"2.11.041112 by agreement of 041110: maintained now by dcr\n");
fprintf(outf,"041112 -nomodeltest for mage/probe (Lmodeltest) see 050121\n");
fprintf(outf,"041114 more fussing with NMR models\n");
fprintf(outf,"050111,17,18,19 flow, esp. autobondrot, annotations...\n");
fprintf(outf,"050119 very low occ atoms have presence but not show contacts\n");
fprintf(outf,"050121 remove -nomodeltest stuff, mage now sends model # \n");
fprintf(outf,"060129 jEdit type fold-comments on each subroutine \n");
fprintf(outf,"  single vdw contact button replaces both wide and small \n");
fprintf(outf,"  -mastername flags extra master={name} (default: dots)\n");
fprintf(outf,"060212 something about hets also marked prot,rna,dna... \n");
fprintf(outf,"060831 -NOOCC atoms NOT filtered by occ value \n");
fprintf(outf,"061018 not treat HIS as aromatic ACCEPTOR \n");
fprintf(outf,"2.12.0708xx ... version number change (approx date) \n");
fprintf(outf,"070801 rmi Updated for wwPDB remediation and PDB v3.0 \n");
fprintf(outf,"070801 Recognizes both new and old heavy atom names \n");
fprintf(outf,"       Recognizes both new and old residue names for DNA \n");
fprintf(outf,"       Builds atom connectivity and is able to generate contact dots on \n");
fprintf(outf,"   Xplor, PDBv2.3, PDBv3.0, or files with mixed format even within a \n");
fprintf(outf,"   single residue \n");
fprintf(outf,"070821 add strcasecmp to utility.c    rwgk \n");
fprintf(outf,"070905 add alternate names of thymine methyl \n");
fprintf(outf,"070913 add detail on probe unformatted to help \n");
fprintf(outf,"070913 fixed handling of mixed RNA files \n");
fprintf(outf,"071010 added support for Hybrid36 atom and residue numbers \n");
fprintf(outf,"071025 added support of Coot style RNA names \n");
fprintf(outf,"071101 added support for two character chain names \n");
fprintf(outf,"071107 added support for OT1 and  OT2 as C-term carboxylate oxygens \n");
fprintf(outf,"071128 bug fix in parsing of command line chainIds \n");
fprintf(outf,"110413 bug fix onlybadout outputs only bo, not cc     gjk\n");
fprintf(outf,"110830 MacOSX10.6 GCC finicky re. no format string in \n");
fprintf(outf,"       writeOutput function for color \n");
fprintf(outf,"110909 changed the elseif portion of coloroutput to have no extra space \n");
fprintf(outf,"111028 (v2.6?) Swati Jain (SJ) -condense flag, used with -u option \n");
fprintf(outf,"   gives one line per source atom - target atom pair, also  \n");
fprintf(outf,"   gives number of dots in that interaction, as an additional column.\n");
fprintf(outf,"2.14.111204 DCR OneDotEach for atoms across gap on line of sight.\n");
fprintf(outf,"   Logic is quite different so separate examineDots subroutines.\n");
fprintf(outf,"   LOneDotEach option calls examineOneDotEach() subroutine.\n");
fprintf(outf,"   General cleanup of notes and rearrangements for code readability.\n");
fprintf(outf,"111205-120109 DCR tuning ODE, methylene defined in select.c\n");
fprintf(outf,"111028 - 1303xx tweaks and new param to trunk, mainly Jeff \n");
fprintf(outf,"130427 dcr reconciliation of trunk <- probeVector incl ODE\n");
fprintf(outf,"2.15.130427 version number change for merged code! \n");
fprintf(outf,"2.16.130509 jjh added support for segid instead of chaind\n");
fprintf(outf,"2.16.130520 jjh fixed bug in segid handling\n");
fprintf(outf,"04/16/2015 - SJ added the -sepworse flag, if true will seperate the overlaps of >= 0.4 and overlaps of >=0.5. This is default by false. Had to change NODEWIDTH value (see probe.h)\n One can change the 0.5 cutoff for worse overlap by specifying the flag DIVWorse\n");

exit(0);

}/*dump_changes()*/
/*}}}dump_changes()__________________________________________________________*/

/*{{{countsummary()***********************************************************/
/*count arrays indecies: Nth pass, contact width (wide to overlap, Hbond,...)*/
/*N either 1 or 2 depending on direction of relationship                     */
/* Old: contact width                                                  20111215dcr*/
/*    0            1            2           3           4            5
mcmc wide:  mcmc close:   mcmc small : mcmc bad:  mcmc H-bond:  mcmc H-weak
scsc wide:  scsc close:   scsc small:  scsc bad:  scsc H-bond:  scsc H-weak
mcsc wide:  mcsc close:   mcsc small:  mcsc bad:  mcsc H-bond:  mcsc H-weak
other wide: other close:  other small: other bad: other H-bond  other H-weak
wide sum:   close sum:    small sum:   bad sum:   H-bond sum    H-weak sum
*/
/*5 contact width criteria include optional 6th index 5 weakHbond 20111215dcr*/
/* New: contact width 04/15/2015 SJ*/
/*    0            1            2           3           4            5             6
 mcmc wide:  mcmc close:   mcmc H-weak:  mcmc small:  mcmc bad:  mcmc worse:  mcmc H-bond:
 scsc wide:  scsc close:   scsc H-weak:  scsc small:  scsc bad:  scsc worse:  scsc H-bond:
 mcsc wide:  mcsc close:   mcsc H-weak:  mcsc small:  mcsc bad:  mcsc worse:  mcsc H-bond:
 other wide: other close:  other H-weak: other small: other bad: other worse: other H-bond
 wide sum:   close sum:    H-weak sum:   small sum:   bad sum:   worse sum:   H-bond sum:
 */

void countsummary(FILE *outf, char* modestr, int Lines, int Pass) /*dcr041101*/
{/* 04/21/2015 SJ: don't think this needs to check flags before printing, it would just print 0 is LworseOverlap, and LweakHbond is false*/
   char message[200];
   int  N=0; /*index Number*/
   int  k=0;

 if(Pass == 1 || Pass == 0) {N = 0;}
 else                       {N = 1;}

   summcont[N][0] = mcmccont[N][0]+scsccont[N][0]+mcsccont[N][0]+othrcont[N][0];
   summcont[N][1] = mcmccont[N][1]+scsccont[N][1]+mcsccont[N][1]+othrcont[N][1];
   summcont[N][2] = mcmccont[N][2]+scsccont[N][2]+mcsccont[N][2]+othrcont[N][2];
   summcont[N][3] = mcmccont[N][3]+scsccont[N][3]+mcsccont[N][3]+othrcont[N][3];
   summcont[N][4] = mcmccont[N][4]+scsccont[N][4]+mcsccont[N][4]+othrcont[N][4];
   summcont[N][5] = mcmccont[N][5]+scsccont[N][5]+mcsccont[N][5]+othrcont[N][5];
   summcont[N][6] = mcmccont[N][6]+scsccont[N][6]+mcsccont[N][6]+othrcont[N][6]; /* added SJ 04/15/2015*/
   /*20111215dcr weakHbond optional Nodewith criteria*/

 if(Pass == 0)
 {/*expect another pass, so save these values, and clear arrays for next pass*/
      for(k=0 ; k<=6 ; k++) /* changed from k<=5 to k<=6 SJ 04/15/2015*/
      {/*copy values*/
         mcmccont[1][k] = mcmccont[0][k];
         scsccont[1][k] = scsccont[0][k];
         mcsccont[1][k] = mcsccont[0][k];
         othrcont[1][k] = othrcont[0][k];
         summcont[1][k] = summcont[0][k];
      }
      for(k=0 ; k<=6 ; k++)/* changed from k<=5 to k<=6 SJ 04/15/2015*/
      {/*clear arrays*/
         mcmccont[0][k] = 0;
         scsccont[0][k] = 0;
         mcsccont[0][k] = 0;
         othrcont[0][k] = 0;
         summcont[0][k] = 0;
      }
 }
 if(Pass == 2)
 {/*should have been a previous pass with values saved*/
      for(k=0 ; k<=6 ; k++)/* changed from k<=5 to k<=6 SJ 04/15/2015*/
      {
         mcmccont[0][k] = mcmccont[0][k] + mcmccont[1][k];
         scsccont[0][k] = scsccont[0][k] + scsccont[1][k];
         mcsccont[0][k] = mcsccont[0][k] + mcsccont[1][k];
         othrcont[0][k] = othrcont[0][k] + othrcont[1][k];
         summcont[0][k] = summcont[0][k] + summcont[1][k];
      }
 }

 if(Pass == 1 || Pass == 2)
 {/*output*/
  if(Lines == 1) /*for accummulating oneliners for comparisons dcr041101*/
  {
   /*Old: 0 wide contact,1 close contact,2 small overlap,3 bad overlap,4 H-bonds*/
   /*Old: INDEX 5 NOW USED FOR OPTIONAL 6TH CONTACT CRITERIA OF WEAKHBONDS*/
   /* New:0 wide contact,1 close contact,2 weak H-bonds, 3 small overlap,4 bad overlap,5 worse overlap, 6 H-bonds - SJ 04/15/2015*/
   // added the [0][6] to all SJ 04/15/2015
      fprintf(outf,": %s "
                ":%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld "
                ":%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld "
                ":%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld "
                ":%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld "
                ":%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :\n"
     ,inputfilename
     ,mcmccont[0][0],mcmccont[0][1],mcmccont[0][2],mcmccont[0][3],mcmccont[0][4]
     ,mcmccont[0][5],mcmccont[0][6]
     ,scsccont[0][0],scsccont[0][1],scsccont[0][2],scsccont[0][3],scsccont[0][4]
     ,scsccont[0][5],scsccont[0][6]
     ,mcsccont[0][0],mcsccont[0][1],mcsccont[0][2],mcsccont[0][3],mcsccont[0][4]
     ,mcsccont[0][5],mcsccont[0][6]
     ,othrcont[0][0],othrcont[0][1],othrcont[0][2],othrcont[0][3],othrcont[0][4]
     ,othrcont[0][5],othrcont[0][6]
     ,summcont[0][0],summcont[0][1],summcont[0][2],summcont[0][3],summcont[0][4]
     ,summcont[0][5],summcont[0][6]
   );

  }
  else /*for kinemage*/
  {
   fprintf(outf,"@text\n");
   fprintf(outf,"probe: %s\n",modestr);
   fprintf(outf,"%s\n",inputfilename);
   /*Old: 0 wide contact,1 close contact,2 small overlap,3 bad overlap,4 H-bonds*/
   //fprintf(outf,":CONTACT:   WIDE   :  CLOSE   :  SMALL   :   BAD    :  H-BOND  : weakHbond :\n");
   /*New: SJ 04/15/2015 0 wide contact,1 close contact,2 weak H-bonds, 3 small overlap,4 bad overlap,5 worse overlap, 6 H-bonds*/
   /* added [0][6] for all SJ 04/21/2015*/
   fprintf(outf,":CONTACT:   WIDE   :  CLOSE   :  weak H-bonds  : SMALL   :   BAD    :  WORSE  :  H-BOND  :\n");
   fprintf(outf,":MCMC   :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :\n"
     ,mcmccont[0][0],mcmccont[0][1],mcmccont[0][2],mcmccont[0][3],mcmccont[0][4]
     ,mcmccont[0][5],mcmccont[0][6]);
   fprintf(outf,":SCSC   :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :\n"
     ,scsccont[0][0],scsccont[0][1],scsccont[0][2],scsccont[0][3],scsccont[0][4]
     ,scsccont[0][5],scsccont[0][6]);
   fprintf(outf,":MCSC   :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :\n"
     ,mcsccont[0][0],mcsccont[0][1],mcsccont[0][2],mcsccont[0][3],mcsccont[0][4]
     ,mcsccont[0][5],mcsccont[0][6]);
   fprintf(outf,":OTHER  :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :\n"
     ,othrcont[0][0],othrcont[0][1],othrcont[0][2],othrcont[0][3],othrcont[0][4]
     ,othrcont[0][5],othrcont[0][6]);
   fprintf(outf,":SUM    :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :\n"
     ,summcont[0][0],summcont[0][1],summcont[0][2],summcont[0][3],summcont[0][4]
     ,summcont[0][5],summcont[0][6]);
  }

  if (Verbose)
  {
   sprintf(message,"%s",inputfilename);
   note(message);
   /*Old: 0 wide contact,1 close contact,2 small overlap,3 bad overlap,4 H-bonds*/
   //sprintf(message,":CONTACT:   WIDE   :  CLOSE   :  SMALL   :   BAD    :  H-BOND  : weakHbond :");
   /*New: SJ 04/15/2015 0 wide contact,1 close contact,2 weak H-bonds, 3 small overlap,4 bad overlap,5 worse overlap, 6 H-bonds*/
   /* added [0][6] for all SJ 04/21/2015*/
   sprintf(message,":CONTACT:   WIDE   :  CLOSE   :  weak H-bonds  : SMALL   :   BAD    :  WORSE  :  H-BOND  :\n");
   note(message);
   sprintf(message,":MCMC   :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :"
     ,mcmccont[0][0],mcmccont[0][1],mcmccont[0][2],mcmccont[0][3],mcmccont[0][4]
     ,mcmccont[0][5],mcmccont[0][6]);
   note(message);
   sprintf(message,":SCSC   :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :"
     ,scsccont[0][0],scsccont[0][1],scsccont[0][2],scsccont[0][3],scsccont[0][4]
     ,scsccont[0][5],scsccont[0][6]);
   note(message);
   sprintf(message,":MCSC   :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :"
     ,mcsccont[0][0],mcsccont[0][1],mcsccont[0][2],mcsccont[0][3],mcsccont[0][4]
     ,mcsccont[0][5],mcsccont[0][6]);
   note(message);
   sprintf(message,":OTHER  :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :"
     ,othrcont[0][0],othrcont[0][1],othrcont[0][2],othrcont[0][3],othrcont[0][4]
     ,othrcont[0][5],othrcont[0][6]);
   note(message);
   sprintf(message,":SUM    :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :%9ld :"
     ,summcont[0][0],summcont[0][1],summcont[0][2],summcont[0][3],summcont[0][4]
     ,summcont[0][5],summcont[0][6]);
   note(message);
  }
 }/*output*/
}
/*}}}countsummary()__________________________________________________________*/
