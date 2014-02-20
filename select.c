/*{{{select.c general comments                                               */
/* name: select.c                                */
/* author: J. Michael Word                       */
/* date written: 2/20/96                         */
/* purpose: parse atom selections                */

/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 1999 J. Michael Word                            */
/*****************************************************************/
/*}}}select.c general comments  */

/*{{{--includes                                                              */
#include <string.h>  /*060902 needs this for strstr() */
#include "select.h"
#include "atomprops.h"
#include "stdconntable.h"
#include "utility.h"

/*}}}--includes  */

/*{{{getPat() ****************************************************************/
/*select.h includes parse.h, so why is getPat here instead of in parse.c ? */
pattern* getPat(char *line, char *which, int verbose)
{
   pattern *pat;

   if (!line) return NULL;

   pat = parseArg(line);  /* parse.c/parseArg() */

   if (verbose && pat) {
      fprintf(stderr, "%s: ", which);
      printPattern(stderr, pat); /* parse.c/printPattern() */
   }

   return pat;
}
/*}}}getPat() _______________________________________________________________*/

/*{{{declarations ...           global AA and Atom properties                */

/*{{{AromaticAtomsTbl[]      ***/
static ResidueAndAtomPair AromaticAtomsTbl[] = {
":PHE:", ": HD1: HD2: HE1: HE2: HZ : DD1: DD2: DE1: DE2: DZ :", 0,
":HIS:", ": HD1: HD2: HE1: HE2: DD1: DD2: DE1: DE2:",           0,
":TYR:", ": HD1: HD2: HE1: HE2: DD1: DD2: DE1: DE2:",           0,
":TRP:", ": HD1: HE1: HE3: HZ2: HZ3: HH2: DD1: DE1: DE3: DZ2: DZ3: DH2:", 0,
":  U:URA:UTP:UDP:UMP: UR:",    ": H3 : HN3: H5 : H6 : D3 : DN3: D5 : D6 :", 0,
":  T:THY:TTP:TDP:TMP:5MU: DT: TR:",          ": H3 : HN3: H6 : D3 : DN3: D6 :", 0,
":  A:ADE:ATP:ADP:AMP:1MA:RIA:T6A: DA: AR:",            ": H8 : H2 : D8 : D2 :", 0,
":  C:CYT:CTP:CDP:CMP:5MC:OMC: DC: CR:",                ": H5 : H6 : D5 : D6 :", 0,
":  G:GUA:GTP:GDP:GMP:GSP:2MG:M2G:7MG:OMG: DG: GR:",
                                      ": H8 : H1 : HN1: D8 : D1 : DN1:", 0,
": YG:1MG:",                                    ": H8 : D8 :",           0,
":PSU:",   ": H6 : D6 : H1 : HN1: D1 : DN1: H3 : HN3: D3 : DN3:",        0,
":  I: DI:",                    ": H8 : H2 : H1 : HN1: D8 : D2 : D1 : DN1:", 0,
":PHE:", ": CG : CD1: CD2: CE1: CE2: CZ :",                TEST_ACCEPT_ANGLE_PROP,
":HIS:", ": ND1: CD2: CE1: NE2: CG :", 0,
":TYR:", ": CG : CD1: CD2: CE1: CE2: CZ :",                TEST_ACCEPT_ANGLE_PROP,
":TRP:", ": CG : CD1: CD2: NE1: CE2: CE3: CZ2: CZ3: CH2:", TEST_ACCEPT_ANGLE_PROP,
":  U:URA:UTP:UDP:UMP:PSU: UR:",
               ": N1 : C2 : N3 : C4 : C5 : C6 :",               TEST_ACCEPT_ANGLE_PROP,
":  T:THY:TTP:TDP:TMP:5MU: DT: TR:",
               ": N1 : C2 : N3 : C4 : C5 : C6 :",               TEST_ACCEPT_ANGLE_PROP,
":  A:ADE:ATP:ADP:AMP:1MA:RIA:T6A:  I: DA: DI: AR:",
               ": N1 : C2 : N3 : C4 : C5 : C6 : N7 : C8 : N9 :",TEST_ACCEPT_ANGLE_PROP,
":  C:CYT:CTP:CDP:CMP:5MC:OMC: DC: CR:",
               ": N1 : C2 : N3 : C4 : C5 : C6 :",               TEST_ACCEPT_ANGLE_PROP,
":  G:GUA:GTP:GDP:GMP:GSP:1MG:2MG:M2G:7MG:OMG: DG: GR:",
               ": N1 : C2 : N3 : C4 : C5 : C6 : N7 : C8 : N9 :",TEST_ACCEPT_ANGLE_PROP,
": YG:",       ": N1 : C2 : N3 : C4 : C5 : C6 : N7 : C8 : N9 :\
                                              : N2 : C11: C12:",TEST_ACCEPT_ANGLE_PROP,
":HEM:", ": N A: C1A: C2A: C3A: C4A: N B: C1B: C2B: C3B: C4B:", TEST_ACCEPT_ANGLE_PROP,
":HEM:", ": N C: C1C: C2C: C3C: C4C: N D: C1D: C2D: C3D: C4D:", TEST_ACCEPT_ANGLE_PROP,
0, 0, 0};

/*061018 not allow HIS to accept H-bonds as aromatic ring system
":HIS:", ": ND1: CD2: CE1: NE2: CG :", TEST_ACCEPT_ANGLE_PROP,      */

/*}}}AromaticAtomsTbl[] ___________________________________________*/
/*{{{AAList                  ***/
static char *AAList = ":GLY:ALA:VAL:PHE:PRO:MET:ILE:LEU:ASP:GLU:LYS:ARG:\
SER:THR:TYR:HIS:CYS:ASN:GLN:TRP:ASX:GLX:ACE:FOR:NH2:NME:MSE:AIB:ABU:PCA:";
/*}}}AAList  */
/*{{{HphobicAAList           ***/
static char *HphobicAAList = ":ALA:VAL:PHE:MET:ILE:LEU:TYR:CYS:TRP:MSE:AIB:ABU:";
/*}}}HphobicAAList */
/*{{{HphilicAAList           ***/
static char *HphilicAAList = ":SER:THR:ASN:GLN:";
/*}}}HphilicAAList */
/*{{{AromaticAAList          ***/
static char *AromaticAAList = ":PHE:HIS:TYR:TRP:HEM:";
/*}}}AromaticAAList */
/*{{{NAList                  ***/
static char *NAList = ":  C:  G:  A:  T:  U:CYT:GUA:ADE:THY:URA:\
CTP:CDP:CMP:GTP:GDP:GMP:ATP:ADP:AMP:TTP:TDP:TMP:UTP:UDP:UMP:\
GSP:H2U:PSU:1MG:2MG:M2G:7MG:5MC:5MU:T6A:1MA:RIA:\
OMC:OMG: YG:  I: DA: DT: DC: DG: DI: AR: UR: TR: GR: CR:";
/*}}}NAList */
/*{{{NAbackboneList          ***/
static char *NAbackboneList = ": P  : O1P: O2P: OP1: OP2: PA : PB : PG :\
: O1A: O2A: O3A: O1B: O2B: O3B: O1G: O2G: O3G: S1G:\
: O5*: C5*: C4*: O4*: C3*: O3*: C2*: O2*: C1*:\
: O5': C5': C4': O4': C3': O3': C2': O2': C1':\
: H1*: H1': H3*: H3': H4*: H4':1H5*:\
:2H5*:*H51:*H52:H5'': H5':1H2*:2H2*:\
: H2*:H2'':2HO*:*HO2: H2': H3T: H5T:\
:3HO*:HO3':*HO3:5HO*:*HO5:HO5':";
/*}}}NAbackboneList */
/*{{{NAbaseGrouping[]         **/
/* used to treat similar nucleic acid bases similarly */
static ResidueAndAtomPair NAbaseGrouping[] = {
":  U:URA:UTP:UDP:UMP: UR:  T:THY:TTP:TDP:TMP:5MU: DT: TR:",     "", baseTU,
":  A:ADE:ATP:ADP:AMP:1MA:RIA:T6A: DA: AR:",                 "", baseA,
":  C:CYT:CTP:CDP:CMP:5MC:OMC: DC: CR:",                     "", baseC,
":  G:GUA:GTP:GDP:GMP:GSP:1MG:2MG:M2G:7MG:OMG: DG: GR:",     "", baseG,
":  I: YG:H2U:PSU: DI:",                                 "", baseOther,
0, 0, 0};
/*}}}NAbaseGrouping[] */
/*{{{MethylResList           ***/
/* for identifying methyl groups */
static char *MethylResList =
":THR:ALA:MET:LEU:VAL:ILE:  T:THY: DT:AIB:ABU:ACE:MSE:NME:HEM:\
 :OMG:OMC:1MA:1MG:2MG:M2G:5MU:5MC:7MG: YG:T6A:";
 /*}}}MethylResList */
 /*{{{MethylAtomsTbl[]       ***/
static ResidueAndAtomPair MethylAtomsTbl[] = { /* including xplor names */
":THR:", ": CG2:1HG2:2HG2:3HG2:HG21:HG22:HG23:",      0,
":ALA:", ": CB :1HB :2HB :3HB : HB1: HB2: HB3:",      0,
":MET:", ": CE :1HE :2HE :3HE : HE1: HE2: HE3:",      0,
":LEU:", ": CD1:1HD1:2HD1:3HD1:HD11:HD12:HD13:\
          : CD2:1HD2:2HD2:3HD2:HD21:HD22:HD23:",      0,
":VAL:", ": CG1:1HG1:2HG1:3HG1:HG11:HG12:HG13:\
          : CG2:1HG2:2HG2:3HG2:HG21:HG22:HG23:",      0,
":ILE:", ": CG2:1HG2:2HG2:3HG2:HG21:HG22:HG23:\
          : CD1:1HD1:2HD1:3HD1:HD11:HD12:HD13:",      0,
":  T:", ": C5M:1H5M:2H5M:3H5M:H5M1:H5M2:H5M3:",      0,
":THY:", ": C5A:1H5 :2H5 :3H5 : H51: H52: H53:",      0,
": DT:", ": C7 : H71: H72: H73:", 		      0,
":AIB:", ": CB1:1HB1:2HB1:3HB1:HB11:HB12:HB13:\
          : CB2:1HB2:2HB2:3HB2:HB21:HB22:HB23:",      0,
":ABU:", ": CG :1HG :2HG :3HG : HG1: HG2: HG3:",      0,
":ACE:", ": CH3:1HH3:2HH3:3HH3:HH31:HH32:HH33:",      0,
":MSE:", ": CE :1HE :2HE :3HE : HE1: HE2: HE3:",      0,
":NME:", ": CH3:1HH3:2HH3:3HH3:HH31:HH32:HH33:",      0,
":HEM:", ": CMA:1HMA:2HMA:3HMA:HMA1:HMA2:HMA3:\
          : CMB:1HMB:2HMB:3HMB:HMB1:HMB2:HMB3:\
          : CMC:1HMC:2HMC:3HMC:HMC1:HMC2:HMC3:\
          : CMD:1HMD:2HMD:3HMD:HMD1:HMD2:HMD3:",      0,

":OMG:OMC:",
         ": CM2:1HM2:2HM2:3HM2:HM22:HM22:HM23:", 0,
":1MA:1MG:",
         ": CM1:1HM1:2HM1:3HM1:HM12:HM12:HM13: C1A:1H1A:2H1A:3H1A:", 0,
":2MG:M2G:", ": CM1:1HM1:2HM1:3HM1:HM12:HM12:HM13: C1A:1H1A:2H1A:3H1A:\
              : CM2:1HM2:2HM2:3HM2:HM22:HM22:HM23: C2A:1H2A:2H2A:3H2A:", 0,
":5MU:5MC:",
         ": CM5:1HM5:2HM5:3HM5:HM52:HM52:HM53: C5A:1H5A:2H5A:3H5A:", 0,
":7MG:", ": CM7:1HM7:2HM7:3HM7:HM72:HM72:HM73: C7A:1H7A:2H7A:3H7A:", 0,
": YG:", ": C3 :1H3 :2H3 :3H3 : H32: H32: H33:\
          : C10:1H10:2H10:3H10:H102:H102:H103:\
          : C19:1H19:2H19:3H19:H192:H192:H193:\
          : C24:1H24:2H24:3H24:H242:H242:H243:", 0,
":T6A:", ": C15:1H15:2H15:3H15:H152:H152:H153:", 0,
0, 0, 0};
/*}}}MethylAtomsTbl[] */

/*{{{MethyleneResList, protein only ***/
/* for identifying methylene groups on which to find methylene H 20111211dcr*/
static char *MethyleneResList =
":LYS:ILE:ARG:PRO:GLN:GLU:MET:CIS:SER:ASN:ASP:PHE:TYR:TRP:LEU:";
 /*}}}MethyleneResList */

/*{{{MethyleneAtomsTbl[]        ***/
static ResidueAndAtomPair MethyleneAtomsTbl[] = {
":LYS:", ": HE2: HE3: HD2: HD3: HG2: HG3: HB2: HB3:",      0,
":ILE:", ":HG12:HG13:",      0,
":PRO:", ": HD2: HD3: HG2: HG3: HB2: HB3:",      0,
":ARG:", ": HD2: HD3: HG2: HG3: HB2: HB3:",      0,
":GLN:", ": HG2: HG3: HB2: HB3:",      0,
":GLU:", ": HG2: HG3: HB2: HB3:",      0,
":MET:", ": HG2: HG3: HB2: HB3:",      0,
":CIS:", ": HB2: HB3:",      0,
":SER:", ": HB2: HB3:",      0,
":ASN:", ": HB2: HB3:",      0,
":ASP:", ": HB2: HB3:",      0,
":PHE:", ": HB2: HB3:",      0,
":TYR:", ": HB2: HB3:",      0,
":HIS:", ": HB2: HB3:",      0,
":TRP:", ": HB2: HB3:",      0,
":LEU:", ": HB2: HB3:",      0,
0, 0, 0};
/*}}}MethyleneAtomsTbl[] */

/*{{{ChargedAAList           ***/
/* for computing charge state (currently treating HIS as charged) */
static char *ChargedAAList = ":ASP:GLU:LYS:ARG:HIS:HEM:";
/*}}}ChargedAAList */
/*{{{ChargedAAAtomsTbl[]     ***/
static ResidueAndAtomPair ChargedAAAtomsTbl[] = {
":ASP:", ": OD1: OD2:",                              NEGATIVE_PROP,
":GLU:", ": OE1: OE2:",                              NEGATIVE_PROP,
":LYS:", ":1HZ :2HZ :3HZ :1DZ :2DZ :3DZ : NZ :\
          : HZ1: HZ2: HZ3: DZ1: DZ2: DZ3:",	     POSITIVE_PROP,
":ARG:", ": HE :1HH1:2HH1:1HH2:2HH2:\
          : DE :1DH1:2DH1:1DH2:2DH2: NE : NH1: NH2:\
          :HH11:HH12:HH21:HH22:\
          :DH11:DH12:DH21:DH22:",		     POSITIVE_PROP,
":HIS:", ": HD1: HE2: DD1: DE2: ND1: NE2:\
          : CG : CD2: CE1:",                         POSITIVE_PROP,
":HEM:", ": O1A: O2A: O1D: O2D:",                    NEGATIVE_PROP,
0, 0, 0};
/*}}}ChargedAAAtomsTbl[] */
/*{{{AlwaysChargedAtomsList  ***/
static char *AlwaysChargedAtomsList =
":1H  :2H  :3H  :1D  :2D  :3D  : HT1: HT2: HT3: DT1: DT2: DT3:\
 : NT : OXT:1OXT:2OXT:";
 /*}}}AlwaysChargedAtomsList */
 /*{{{AmbigChargedAtomsList  ***/
static char *AmbigChargedAtomsList = ": O  : N  :";
/*}}}AmbigChargedAtomsList */
/*{{{ChargedNucAcidAtomsList ***/
static char *ChargedNucAcidAtomsList = ": P  : O1P: O2P:\
 PA : PB : PG : O1A: O2A: O3A: O1B: O2B: O3B: O1G: O2G: O3G:\
 S1G:";
 /*}}}ChargedNucAcidAtomsList */
 /*{{{WaterList              ***/
static char *WaterList = ":HOH:DOD:H2O:WAT:TIP:SOL:MTO:";
/*}}}WaterList */
/*{{{DonorAcceptorAtomTbl[]  ***/
static ResidueAndAtomPair DonorAcceptorAtomTbl[] = {
":GLY:ALA:VAL:PHE:PRO:MET:ILE:LEU:ASP:GLU:LYS:ARG:\
 :SER:THR:TYR:HIS:CYS:ASN:GLN:TRP:ASX:GLX:NH2:NME:MSE:AIB:ABU:PCA:",
   ": N  : NT : H  : HT1: HT2: HT3:1H  :2H  :3H  :\
         : HN : D  : DT1: DT2: DT3:1D  :2D  :3D  :",
				       DONOR_PROP,

":GLY:ALA:VAL:PHE:PRO:MET:ILE:LEU:ASP:GLU:LYS:ARG:\
 :SER:THR:TYR:HIS:CYS:ASN:GLN:TRP:ASX:GLX:ACE:FOR:MSE:AIB:ABU:PCA:",
         ": O  : OXT:1OXT:2OXT:",      ACCEPTOR_PROP,

":ASN:ASX:", ": AD1: AD2:",            DONOR_PROP|ACCEPTOR_PROP,
":ASN:ASX:", ":1HD1:2HD1:HD11:HD12:\
              :1DD1:2DD1:DD11:DD12:",  DONOR_PROP, /* if ass backwards */
":GLN:GLX:", ": AE1: AE2:",            DONOR_PROP|ACCEPTOR_PROP,
":GLN:GLX:", ":1HE1:2HE1:HE11:HE12:\
              :1DE1:2DE1:DE11:DE12:",  DONOR_PROP, /* if ass backwards */

":ASN:", ": OD1:",                     ACCEPTOR_PROP,
":ASN:", ": ND2:1HD2:2HD2:HD21:HD22:\
               :1DD2:2DD2:DD21:DD22:", DONOR_PROP,
":GLN:", ": OE1:",                     ACCEPTOR_PROP,
":GLN:", ": NE2:1HE2:2HE2:HE21:HE22:\
               :1DE2:2DE2:DE21:DE22:", DONOR_PROP,

":ASP:", ": OD1: OD2:",                ACCEPTOR_PROP,
":GLU:", ": OE1: OE2:",                ACCEPTOR_PROP,
":LYS:", ": NZ :1HZ :2HZ :3HZ : HZ1: HZ2: HZ3:\
               :1DZ :2DZ :3DZ : DZ1: DZ2: DZ3:",
				       DONOR_PROP,
":ARG:", ": NE : NH1: NH2:\
          : HE :1HH1:2HH1:1HH2:2HH2:HH11:HH12:HH21:HH22:\
	  : DE :1DH1:2DH1:1DH2:2DH2:DH11:DH12:DH21:DH22:",
				       DONOR_PROP,

":HIS:", ": ND1: NE2:",                DONOR_PROP|ACCEPTOR_PROP,
":HIS:", ": HD1: HE2: DD1: DE2:",      DONOR_PROP,
":HIS:", ": HD2: HE1: DD2: DE1:",      DONOR_PROP|CH_DONOR_PROP,
":SER:", ": OG :",                     DONOR_PROP|ACCEPTOR_PROP,
":SER:", ": HG : DG :",                DONOR_PROP,
":THR:", ": OG1:",                     DONOR_PROP|ACCEPTOR_PROP,
":THR:", ": HG1: DG1:",                DONOR_PROP,
":TYR:", ": OH :",                     DONOR_PROP|ACCEPTOR_PROP,
":TYR:", ": HH : DH :",                DONOR_PROP,
":CYS:", ": SG :",                     ACCEPTOR_PROP,           /* usually */
":CYS:", ": HG : DG :",                DONOR_PROP,
":TRP:", ": NE1: HE1: DE1:",           DONOR_PROP,
":MET:", ": SD :",                     ACCEPTOR_PROP,
":MSE:", ":SED :",                     ACCEPTOR_PROP,
":NH2:", ": HN1: HN2:1HN :2HN : DN1: DN2:1DN :2DN :", DONOR_PROP,

":  C:  G:  A:  T:  U:CYT:GUA:ADE:THY:URA:\
 :CTP:CDP:CMP:GTP:GDP:GMP:ATP:ADP:AMP:TTP:TDP:TMP:UTP:UDP:UMP:GSP:\
 :H2U:PSU:1MG:2MG:M2G:5MC:5MU:T6A:1MA: YG:  I: DA: DT: DC: DG: DI: AR: UR: TR: CR: GR:",
          ": O2*: O2':",                         DONOR_PROP|ACCEPTOR_PROP,
":  C:  G:  A:  T:  U:CYT:GUA:ADE:THY:URA:\
 :CTP:CDP:CMP:GTP:GDP:GMP:ATP:ADP:AMP:TTP:TDP:TMP:UTP:UDP:UMP:GSP:\
 :H2U:PSU:1MG:2MG:M2G:5MC:5MU:T6A:1MA:RIA:OMC:OMG: YG:  I: DA: DT: DC: DG: DI: AR: UR: TR: CR: GR:",
  ":2HO*:3HO*:5HO*: H3T: H5T:HO3':HO5':\
   :2DO*:3DO*:5DO*: D3T: D5T:DO3':DO5':",
/*  ":2HO*:3HO*:5HO*: H2': H3': H5': H3T: H5T:HO3':HO5':\  RMI removed H2', H3' and H5' 070723 */
/*   :2DO*:3DO*:5DO*: D2': D3': D5': D3T: D5T:DO3':DO5':", */
                                                         DONOR_PROP,
":  C:  G:  A:  T:  U:CYT:GUA:ADE:THY:URA:\
 :CTP:CDP:CMP:GTP:GDP:GMP:ATP:ADP:AMP:TTP:TDP:TMP:UTP:UDP:UMP:GSP:\
 :H2U:PSU:1MG:2MG:M2G:5MC:5MU:T6A:1MA:RIA:OMC:OMG: YG:  I: DA: DT: DC: DG: DI: AR: UR: TR: CR: GR:",
  ": O1P: O2P: OP1: OP2: O1A: O2A: O3A: O1B: O2B: O3B: O1G: O2G: O3G: S1G:\
   : O3*: O5*: O3': O5':",                               ACCEPTOR_PROP,

":  U:URA:UTP:UDP:UMP:H2U: DU: UR:",": N3 : H3 : HN3: D3 : DN3:",DONOR_PROP,
":PSU:",       ": N1 : H1 : D1 : HN1: DN1:\
                : N3 : H3 : D3 : HN3: DN3:",             DONOR_PROP,
":  U:URA:UTP:UDP:UMP:H2U:PSU: DU: UR:",   ": O2 : O4 :",        ACCEPTOR_PROP,
":  T:THY:TTP:TDP:TMP:5MU: DT: TR:",": N3 : H3 : HN3: D3 : DN3:",DONOR_PROP,
":  T:THY:TTP:TDP:TMP:5MU: DT: TR:",  ": O2 : O4 :",             ACCEPTOR_PROP,
":  A:ADE:ATP:ADP:AMP:1MA:RIA: DA: AR:",
         ": N6 :1H6 :2H6 : H61: H62:1HN6:2HN6:\
               :1D6 :2D6 : D61: D62:1DN6:2DN6:",         DONOR_PROP,
":  A:ADE:ATP:ADP:AMP:RIA: DA: AR:",     ": N1 : N3 : N7 :",     ACCEPTOR_PROP,
                ":1MA:  I: DI:",          ": N3 : N7 :",     ACCEPTOR_PROP,

":T6A:", ": N6 : HN6: DN6: N11: HN1: DN1: H14: D14:\
                              : H11: D11: HO4: DO4:",    DONOR_PROP,
":T6A:", ": N1 : N3 : N7 : O10:AO13:BO13: O14:",         ACCEPTOR_PROP,

":  C:CYT:CTP:CDP:CMP:OMC:5MC: DC: CR:",
         ": N4 :1H4 :2H4 : H41: H42:1HN4:2HN4:\
               :1D4 :2D4 : D41: D42:1DN4:2DN4:",         DONOR_PROP,
":  C:CYT:CTP:CDP:CMP:OMC:5MC: DC: CR:",  ": O2 : N3 :",         ACCEPTOR_PROP,
":  G:GUA:GTP:GDP:GMP:GSP:OMG:7MG: DG: GR:",
         ": N1 : N2 : H1 :1H2 :2H2 : H21: H22: HN1:1HN2:2HN2:\
                    : D1 :1D2 :2D2 : D21: D22: DN1:1DN2:2DN2:",DONOR_PROP,
":2MG:", ": N1 : N2 : H1 : H2 : HN1: HN2:\
                    : D1 : D2 : DN1: DN2:",              DONOR_PROP,
":  G:GUA:GTP:GDP:GMP:GSP:OMG:1MG:2MG:M2G: DG: GR:",
                                 ": O6 : N3 : N7 :",     ACCEPTOR_PROP,
": YG:",     ": O6 : N2 : N7 : O17: O17: O22: O23:",     ACCEPTOR_PROP,

":1MG:", ": N2 :1H2 :2H2 : H21: H22:1HN2:2HN2:\
               :1D2 :2D2 : D21: D22:1DN2:2DN2:",         DONOR_PROP,
":M2G:", ": N1 : H1 : HN1: D1 : DN1:",                   DONOR_PROP,
":7MG:", ": O6 : N3 :",                                  ACCEPTOR_PROP,

":  A:ADE:ATP:ADP:AMP:1MA:RIA:T6A:  I: DA: DI: AR:",        ": H2 : H8 :", DONOR_PROP|CH_DONOR_PROP,
":  G:GUA:GTP:GDP:GMP:GSP:OMG:1MG:2MG:M2G:7MG: YG: DG: GR:", ": H8 :", DONOR_PROP|CH_DONOR_PROP,

#ifdef EXPLICIT_WATER_ATOM_NAMES
":HOH:DOD:H2O:WAT:TIP:SOL:MTO:", ": O  : OH2: OD2: OW :", DONOR_PROP|ACCEPTOR_PROP,
":HOH:DOD:H2O:WAT:TIP:SOL:MTO:",
  ": H  : H1 : H2 :1H  :2H  : D  : D1 : D2 :1D  :2D  : H? :", DONOR_PROP,
#endif

/* also, the aromatic heavy atoms are considered H bond acceptors (see setAromaticProp) */

0, 0, 0};
/*}}}DonorAcceptorAtomTbl[] */
/*}}}declarations */

/*{{{naBaseCategory() ********************************************************/
/* used to assign a category to atoms to allow coloring */
/* dots by base rather than atom type */
/* naBaseCategory() must be called after setProperties() */

int naBaseCategory(atom *a) {
   ResidueAndAtomPair *pair;

   if ((a->props & DNA_PROP) && (a->props & SC_PROP)) {
      if(strstr(NAList, a->r->resname)) {
         for (pair = NAbaseGrouping; pair->rlist; pair++){
	    if(strstr(pair->rlist, a->r->resname)) {
	       return pair->bits; /* bits contain the atom type */
	    }
         }
         return baseOther; /* not in table */
      }
     else return baseOther; /* not in list */
   }
   else return nonBase; /* backbone or not recognized as a nucleic acid */
}
/*}}}naBaseCategory() _______________________________________________________*/

/*{{{setAromaticProp() *******************************************************/
/* hunt for aromatic atoms and set AROMATIC property & atomHarom */
/* also find aromatic carbon atoms and set atomCarom */

void setAromaticProp(atom *a) {
   ResidueAndAtomPair *pair;

   if((strstr(AromaticAAList, a->r->resname))
   || (strstr(NAList,         a->r->resname))) {
      for (pair = AromaticAtomsTbl; pair->rlist && pair->alist; pair++){
	 if((strstr(pair->rlist, a->r->resname))
	 && (strstr(pair->alist, a->atomname))) {
	    a->props |= AROMATIC_PROP;
	    a->props |= pair->bits; /* selectively set TEST_ACCEPT_ANGLE prop */
	    if(isHatom(a->elem)) { a->elem = atomHarom; }
	    if(isCatom(a->elem)) { a->elem = atomCarom; }
	    break;
	 }
      }
   }
}
/*}}}setAromaticProp() ______________________________________________________*/

/*{{{setMethylProp() *********************************************************/
/* hunt for Methyl atoms and set METHYL property */

void setMethylProp(atom *a) {
   ResidueAndAtomPair *pair;

   if(strstr(MethylResList, a->r->resname)) {
      for (pair = MethylAtomsTbl; pair->rlist && pair->alist; pair++){
	 if((strstr(pair->rlist, a->r->resname))
	 && (strstr(pair->alist, a->atomname))) {
	    a->props |= METHYL_PROP;
	    break;
	 }
      }
   }
}
/*}}}setMethylProp()_________________________________________________________*/

/*{{{setMethyleneProp() ******************************************************/
/* hunt for Methyene atoms and set METHYLENE property */

void setMethyleneProp(atom *a) {
   ResidueAndAtomPair *pair;

   if(strstr(MethyleneResList, a->r->resname)) {
      for(pair = MethyleneAtomsTbl; pair->rlist && pair->alist; pair++) {
         if(   (strstr(pair->rlist, a->r->resname))
            && (strstr(pair->alist, a->atomname))   ) {
            a->props |= METHYLENE_PROP;
            break;
         }
      }
   }
}
/*}}}setMethyleneProp()______________________________________________________*/

/*{{{isCarbonylAtom() ********************************************************/
int isCarbonylAtom(atom *a) { /* limitation: this will not handle het groups properly */
   return (strstr(AAList, a->r->resname) && (!strcmp(a->atomname,  " C  ")))
   ||((strstr(":ASP:ASN:ASX:",a->r->resname))&&(strstr(": CG :",a->atomname)))
   ||((strstr(":GLU:GLN:GLX:",a->r->resname))&&(strstr(": CD :",a->atomname)));
}
/*}}}isCarbonylAtom() _______________________________________________________*/

/*{{{setProperties() *********************************************************/
/* set property flags for this atom and return T/F as an indicator     */
/*        that this atom is important in figuring out if N and O atoms */
/*        are in the N or C terminus                                   */
/* later functions can look through a list of these atoms to see if    */
/* atoms marked MAYBECHG_PROP are in a terminal residue                */

/*060212 Properties as specified in current pdb-format files, */
/*are not logically complete.  HETATM records and residue names that are*/
/*known valid protein or nucleic acid names can occur both in regions that */
/*are protein/nucleic chain and regions that are het-chains */
/*   The only way to beat this is to accummulate evidence of the type of */
/*region as atoms are read in, then set a new and different flag for region*/

void setProperties(atom *a, int hetflag, int hb2aromaticFace, int chohb) {
   char *s;
   ResidueAndAtomPair *pair;
   int addCharge = FALSE, checkTerminal = FALSE;

   a->props = 0;

   setAromaticProp(a);
   setMethylProp(a);
   setMethyleneProp(a); /*20111211dcr*/

   if (strstr(WaterList, a->r->resname)) {
      a->props |= WATER_PROP;
   }
   else if(hetflag) { /* non water hets */
      if (! strstr(":MSE:", a->r->resname)) {
	 /* SelenoMet is not a HET for our purposes here */

	 a->props |= HET_PROP;
      }
   }

   if(strstr(AAList, a->r->resname)) { /* proteins and peptides */
      a->props |= PROT_PROP;
           if(strstr(HphobicAAList, a->r->resname)){a->props |= RHPHOBIC_PROP;}
      else if(strstr(HphilicAAList, a->r->resname)){a->props |= RHPHILIC_PROP;}
      else if(strstr(ChargedAAList, a->r->resname)){a->props |= RCHARGED_PROP;}
      s = a->atomname;
      if (a->elem == atomC || isHatom(a->elem)) {
             if (s[2] == 'A') {
	       if (!strcmp(s, "2HA ")){ a->props |= ALPHA_PROP; }
	       else                   { a->props |= ALPHA_PROP|MC_PROP; }
	}
        else if (s[2] == 'B')         { a->props |=  BETA_PROP; }
        else if (!strcmp(s,  " C  ")) { a->props |=    MC_PROP; }
        else if (!strcmp(s+1, "H  ")) { a->props |=    MC_PROP; }
        else if (!strcmp(s+1, "D  ")) { a->props |=    MC_PROP; }
        else if (strstr(": HXT: DXT: HN :", s)){ a->props|= MC_PROP; }
      }
      else if (!strcmp(s, " N  ")) { a->props |= MC_PROP; }
      else if (!strcmp(s, " O  ")) { a->props |= MC_PROP; }

      else if (strstr(AlwaysChargedAtomsList, s)) { a->props |= MC_PROP; }

      if (!(a->props & MC_PROP)) { a->props |= SC_PROP;}
   }
   else if(strstr(NAList, a->r->resname)) { /* DNA and RNA */
      a->props |= DNA_PROP;
      if (strstr(NAbackboneList, a->atomname)) { a->props |= MC_PROP; }
      else { a->props |= SC_PROP;}
   }
   switch(a->elem) {
   case atomN: a->props |= N_PROP; break;
   case atomC: a->props |= C_PROP; break;
   case atomO: a->props |= O_PROP; break;
   case atomS: a->props |= S_PROP; break;
   case atomP: a->props |= P_PROP; break;
   }

   if (isHatom(a->elem)) { a->props |= H_PROP; }

   if (atomHasProp(a->elem, METALIC_ATOM_FLAG)) { a->props |= METAL_PROP; }
   if (atomHasProp(a->elem, IONIC_ATOM_FLAG)) { a->props |= ION_PROP; }/*dcr041007*/

/* HBond Donor/Acceptor status  */
   if (hb2aromaticFace
    && (a->props & TEST_ACCEPT_ANGLE_PROP)) { a->props |= ACCEPTOR_PROP; }

   if (a->props & WATER_PROP) { /* water O or H */
      if (a->elem == atomO) { a->props |= DONOR_PROP|ACCEPTOR_PROP; }
      else if (isHatom(a->elem)) { a->props |= DONOR_PROP; }
   }
   else if ((a->elem == atomN)||(a->elem == atomO)||(a->elem == atomS)||isHatom(a->elem)){
      if((strstr(AAList, a->r->resname))
      || (strstr(NAList, a->r->resname))
      || (strstr(WaterList, a->r->resname))) {
	 for (pair = DonorAcceptorAtomTbl; pair->rlist && pair->alist; pair++){
	    if((strstr(pair->rlist, a->r->resname))
	    && (strstr(pair->alist, a->atomname))) {

               if ((! chohb) && (pair->bits & CH_DONOR_PROP)) {
	          /* do nothing --- not doing CH..O type hbonds */
               }
	       else {        /* we generally fall in this category */
	          a->props |= pair->bits; /* set donor/acceptor properties */
	       }
	       break;
	    }
	 }
	 if ( !( isHatom(a->elem) ||
	         (a->props & (ACCEPTOR_PROP|DONOR_PROP)) ) ) {
	    a->props |= ACCEPTOR_PROP; /* default for N,O,S */
	 }
      }
      else if (!isHatom(a->elem)){/* mostly handles het groups we are unsure about */
	 a->props |= DONOR_PROP|ACCEPTOR_PROP;
      }
   }
/* note: final H atom DONOR assignment done later in updateHydrogenInfo()     */
/*       but the above code will set (for selections) most HB donor hydrogens */

   /* determining charge state of atoms */
   addCharge = FALSE;
   checkTerminal = FALSE;
   if (a->props & PROT_PROP) {
      if(strstr(AmbigChargedAtomsList, a->atomname)) {
           if (!strcmp(a->atomname, " N  ")) { /* *** N+ must be first residue *** */
	       a->props |= MAYBECHG_PROP;
	   }
           else if (!strcmp(a->atomname, " O  ")) { /* *** O- must be last residue  *** */
	       a->props |= MAYBECHG_PROP;
	       checkTerminal = TRUE; /* O used to indicate end of chain */
	   }
      }
      else if(strstr(AlwaysChargedAtomsList, a->atomname)) {
	 addCharge = TRUE;
	 checkTerminal = TRUE; /* 1H,2H,3H or OXT  indicate end of chain */
      }
      else if (strstr(ChargedAAList, a->r->resname)) {
	 for (pair = ChargedAAAtomsTbl; pair->rlist && pair->alist; pair++){
	    if((strstr(pair->rlist, a->r->resname))
	    && (strstr(pair->alist, a->atomname))) {
	       /* addCharge = TRUE; */

	       a->props |= pair->bits;
	       break;
	    }
	 }
      }
   }
   else if ((a->props & DNA_PROP)
        &&  (strstr(ChargedNucAcidAtomsList, a->atomname))) {
      addCharge = TRUE;
   }
   else if (a->props & METAL_PROP) { addCharge = TRUE; }
   else if (a->props & ION_PROP) { addCharge = TRUE; } /*dcr041007*/

   if(addCharge) {
      if (a->props &  METAL_PROP) { a->props |= POSITIVE_PROP; }
      else if (a->props & H_PROP) { a->props |= POSITIVE_PROP; }
      else {
	 switch(a->elem) {
	 case atomN:  a->props |= POSITIVE_PROP; break;
	 case atomO:  a->props |= NEGATIVE_PROP; break;
	 case atomS:  a->props |= NEGATIVE_PROP; break;
	 case atomP:  a->props |= NEGATIVE_PROP; break;
         /*dcr041007 should O,S,Se be allowed to be acceptors ???? */
         /*dcr041007 halides get (NEGATIVE_PROP|ACCEPTOR_PROP) */
	 case atomF:   a->props |= (NEGATIVE_PROP|ACCEPTOR_PROP); break;
	 case atomCl:  a->props |= (NEGATIVE_PROP|ACCEPTOR_PROP); break;
	 case atomBr:  a->props |= (NEGATIVE_PROP|ACCEPTOR_PROP); break;
	 case atomI:   a->props |= (NEGATIVE_PROP|ACCEPTOR_PROP); break;
         /*dcr041007 halides only need ACCEPTOR_PROP to not clash with polar H*/
	 }
      }
   }

   if (checkTerminal) {
      /* now we have to fixup some of the MAYBECHG_PROP atoms elsewhere */
      a->props |= CHECK_ENDS_PROP; /* these atoms will help us do that  */
   }

}/*setProperties*/
/*}}}setProperties() ________________________________________________________*/

/*{{{matchPat() ************* only called from probe.c/selectSource() ********/
/* recursive at OR,AND,NOT nodes; TRUE,FALSE nodes allow those logic choices*/
int matchPat(atom *a, pattern *pat)
{
  int lp, rc = FALSE;

  if (pat)
  {
   switch(pat->type)
   {
   case     OR_NODE: rc = matchPat(a, pat->lhs) ?
			   TRUE : matchPat(a, pat->rhs); break;
   case    AND_NODE: rc = matchPat(a, pat->lhs) ?
			   matchPat(a, pat->rhs) : FALSE; break;
   case    NOT_NODE: rc = !matchPat(a, pat->lhs); break;
   case   FILE_NODE: rc = (a->r->file  == pat->val); break;
   case  MODEL_NODE: rc = (a->r->model == pat->val); break;
   case  CHAIN_NODE: lp = strlen(lexString(pat->val));
                     rc = (strncmp(lexString(pat->val), a->r->chain, lp)
                                                            == 0); break;
   case    ALT_NODE: rc = (a->altConf == ' ' || a->altConf == pat->val); break;
   case    RES_NODE: rc = (a->r->resid == pat->val); break;
   case  RANGE_NODE: rc = (pat->lhs->val <= a->r->resid
                             && a->r->resid <= pat->rhs->val); break;
   case  RTYPE_NODE: lp = strlen(lexString(pat->val));
		     rc = (strncmp(lexString(pat->val), a->r->resname, lp)
							    == 0); break;
   case   DIST_NODE: rc = atomWithinDistance(a, pat->fvec); break;
   case   PROP_NODE: rc = (a->props & pat->val); break;
   case  ANAME_NODE: lp = strlen(lexString(pat->val));
		     rc = (strncmp(lexString(pat->val), a->atomname, lp)
							    == 0); break;
   case  SEGID_NODE: lp = strlen(lexString(pat->val));
		     rc = (strncmp(lexString(pat->val), a->r->segid, lp)
							    == 0); break;
   case   TRUE_NODE: rc = TRUE;  break;
   case  FALSE_NODE: rc = FALSE; break;
   case OCC_LT_NODE: rc = (a->occ  < pat->val/100.0); break;
   case OCC_GT_NODE: rc = (a->occ  > pat->val/100.0); break;
   case   B_LT_NODE: rc = (a->bval < pat->val); break;
   case   B_GT_NODE: rc = (a->bval > pat->val); break;
   case    INS_NODE: rc = (a->r->resInsCode == pat->val); break;

 /* ignores the insertion code */
   case INS_RANGE_NODE: rc = (((pat->lhs->val == AND_NODE)
                                ? pat->lhs->lhs->val
				: pat->lhs->val)
				<= a->r->resid
                             && a->r->resid <=
			      ((pat->rhs->val == AND_NODE)
                                ? pat->rhs->lhs->val
				: pat->rhs->val));
			break;
   default: { char msg[100];
         sprintf(msg, "unknown pattern type: %d", pat->type);
         halt(msg);
      }
   }/*switch*/
  }
  return rc;
}
/*}}}matchPat() _____________________________________________________________*/

/*{{{atomWithinDistance() ****************************************************/
/* is the atom within the specified distance from the given point? */
int atomWithinDistance(atom *a, float *fvec) {
   double dx = a->loc.x - fvec[1];
   double dy = a->loc.y - fvec[2];
   double dz = a->loc.z - fvec[3];
   return ((dx*dx)+(dy*dy)+(dz*dz)) < (fvec[0]*fvec[0]);
}
/*}}}atomWithinDistance() ___________________________________________________*/

/*{{{setHydrogenParentName() *************************************************/
/* Set the parent atom name for a hydrogen     */
/* if can be determined from the hydrogen name.*/
/* These parent names are used when a hydrogen */
/* has more than one parent close enough to    */
/* bond. No such bonding table is currently    */
/* setup for non-hydrogens.                    */
char * setHydrogenParentName(char *rname, char *aname) {

   return searchForStdBondingPartner(rname, aname, TRUE);
}
/*}}}setHydrogenParentName() ________________________________________________*/

/*{{{setMainchainBonding() ***************************************************/
/* For heavy atoms, set the name of atoms which would bond with */
/* in standard residues.                                        */
/* This is used only when the -STDBOND flag                     */
/* is on the command line (UseStdBond == TRUE).                 */

char * setMainchainBonding(char *rname, char *aname) {

   return searchForStdBondingPartner(rname, aname, FALSE);
}
/*}}}setMainchainBonding() __________________________________________________*/
