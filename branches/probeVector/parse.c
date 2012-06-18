/* name: parse.c                                         */
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

#include <stdlib.h> /*060902 needs this for malloc() */
#include <string.h>
#include <ctype.h>
#include "utility.h"
#include "parse.h"

static symbolEntry SymTable[LEXSYMMAX] = {0};
static char *inputStr = 0;
static int tokenval = 0;
static float realtokval = 0.0;
static int lookahead = 0;
static char lexbuf[LXBUFSZ+1] = {0};
static char Lexemes[LEXSTRINGMAX] = {0};
static int LastLexChar = -1;
static int LastSymEntry = 0;
static int ParseError = FALSE;

pattern* parseArg(char *s) { /*called by select.c/getPat() */
   pattern *pat;
   inputStr = s;
   lookahead = lexan();

   ParseError = FALSE;

   pat = exprItem();
   matchTok(DONE_TOK);

   if (ParseError == TRUE) {
      pat = freePattern(pat);
   }

   return pat;
}

pattern* freePattern(pattern *pat) {
    if (pat) {
      pat->lhs = freePattern(pat->lhs);
      pat->rhs = freePattern(pat->rhs);
      free(pat);
    }
    return NULL;
}

pattern* exprItem() {
   pattern *lhs, *rhs;
   
   lhs = setItem();
   for (;;) {
      if (lookahead == '|') {
	 matchTok('|');
	 rhs = setItem();
	 lhs = makeNode(lhs, rhs, OR_NODE);
      }
      else {
	 break;
      }
   }
   return lhs;
}

pattern* setItem() {
   pattern *lhs, *rhs;
   
   lhs = limitItem();
   for (;;) {
      if((lookahead == '(') /* a limitItem must begin with a feature! */
      || (lookahead == NOT_TOK)
      || (lookahead == NUM_TOK)
      || (lookahead ==  ID_TOK)
      || (lookahead == WITHIN_TOK)
      || (lookahead == '-')) {
	 rhs = limitItem();
	 lhs = makeNode(lhs, rhs, AND_NODE);
      }
      else {
	 break;
      }
   }
   return lhs;
}

pattern* limitItem() {
   pattern *lhs, *rhs;
   
   lhs = featureItem();
   for (;;) {
      if (lookahead == ',') {
	 matchTok(',');
	 rhs = featureItem();
	 lhs = makeNode(lhs, rhs, OR_NODE);
      }
      else {
	 break;
      }
   }
   return lhs;
}

pattern* featureItem() {
   pattern *f = NULL;

   if (lookahead == '(') {
      matchTok('(');
      f = exprItem();
      matchTok(')');
   }
   else if (lookahead == NOT_TOK) {
      matchTok(NOT_TOK);
      f = featureItem();
      f = makeNode(f, NULL, NOT_NODE);
   }
   else if (lookahead == NUM_TOK) {
      f = numericItem(+1);
   }
   else if (lookahead == ID_TOK) {
      f = idItem();
   }
   else if (lookahead == WITHIN_TOK) {
      f = distItem();
   }
   /* if an item starts with a minus  */
   /* it is a negative residue number */
   else if (lookahead == '-') {
      matchTok('-');
      if (lookahead == NUM_TOK) {
         f = numericItem(-1);
      }
      else {
         char msg[100];
         describeLookahead(msg,"syntax: minus sign requires a number ");
         errmsg(msg);
         ParseError = TRUE;
      }
   }
   else {
     char msg[100];
     describeLookahead(msg,"syntax: improper feature ");
     errmsg(msg);
     ParseError = TRUE;
   }
   return f;
}

pattern* numericItem(int signFactor) {
   pattern *lhs, *rhs, *rrs;
   int inscode = 0;

   lhs = makeTerminal(RES_NODE, signFactor*tokenval);
   inscode = lexbuf[0];
   matchTok(NUM_TOK);

   if (isalpha(inscode) || inscode == ' ') { /* combine res# with insert code */
      rhs = makeTerminal(INS_NODE, inscode);
      lhs = makeNode(lhs, rhs, AND_NODE);
   }

   if (lookahead == '-') { /* first dash is the range separator */
      matchTok('-');

      signFactor = +1;
      if (lookahead == '-') { /* this dash specifies a negative num */
         matchTok('-');
         signFactor = -1;
      }
      rhs = makeTerminal(RES_NODE, signFactor*tokenval);
      inscode = lexbuf[0];
      matchTok(NUM_TOK);

      if (isalpha(inscode) || inscode == ' ') {  /* combine res# with insert code */
         rrs = makeTerminal(INS_NODE, inscode);
         rhs = makeNode(rhs, rrs, AND_NODE);
      }

      if ((lhs->type == AND_NODE)|| (rhs->type == AND_NODE)) {
        /* more than just a span of numbers */
	/* includes insert codes            */
        lhs = makeNode(lhs, rhs, INS_RANGE_NODE);
      }
      else {
        lhs = makeNode(lhs, rhs, RANGE_NODE);
      }
   }
   return lhs;
}

pattern* idItem() {
   pattern *lhs = NULL;
   int type = 0, val = 0;
   char *s, *tmpstr;

   s = SymTable[tokenval].lexptr;

   if ((strncmp(s, "FILE", 4) == 0) && isdigit(s[4])) {
      type = FILE_NODE; val = parseInteger(s,4,10);
   }
   else if ((strncmp(s, "CHAIN", 5) == 0) && (isalnum(s[5])||(s[5]=='_'))) {
      type = CHAIN_NODE;
      for(tmpstr = s+5;  *tmpstr; tmpstr++) { /* substitute blank for _ */
         if (*tmpstr == '_') {*tmpstr = ' ';}
      }
      if (s[6] == '\0') { 
         s[4] = ' ';
         val = lookup(s+4); if(val == 0) val = insert(s+4, ID_TOK);
      }
      else {
         val = lookup(s+5); if(val == 0) val = insert(s+5, ID_TOK);
      }
     /* type = CHAIN_NODE; val = (unsigned char)s[5]; */
   }
   else if ((strncmp(s, "ALT", 3) == 0) && (isalnum(s[3])||(s[3]=='_'))) {
      if (s[3] == '_') {s[3] = ' ';}
      type = ALT_NODE; val = (unsigned char)s[3];
   }
   else if ((strncmp(s, "MODEL", 5) == 0) && isdigit(s[5])) {
      type = MODEL_NODE; val = parseInteger(s,5,10);
   }
   else if ((strncmp(s, "ATOM", 4) == 0) && (isalnum(s[4])||(s[4]=='_'))) {
      type = ANAME_NODE;
      for(tmpstr = s+4; *tmpstr; tmpstr++) { /* substitute blank for _ */
	 if (*tmpstr == '_') {*tmpstr = ' ';}
      }
      val = lookup(s+4); if(val == 0) val = insert(s+4, ID_TOK);
   }
   else if ((strncmp(s, "SEG", 3) == 0) && (isalnum(s[3])||(s[3]=='_'))) {
      type = SEGID_NODE;
      for(tmpstr = s+3; *tmpstr; tmpstr++) { /* substitute blank for _ */
	 if (*tmpstr == '_') {*tmpstr = ' ';}
      }
      val = lookup(s+3); if(val == 0) val = insert(s+3, ID_TOK);
   }
   else if (strcmp(s, "ALL") == 0) {
      type = TRUE_NODE; val = 0;
   }
   else if (strcmp(s, "NONE") == 0) {
      type = FALSE_NODE; val = 0;
   }
   else if (strcmp(s, "ALPHA") == 0) {
      type = PROP_NODE; val = ALPHA_PROP;
   }
   else if (strcmp(s, "BETA") == 0) {
      type = PROP_NODE; val = BETA_PROP;
   }
   else if ((strcmp(s, "MAINCHAIN") == 0) || (strcmp(s, "MC") == 0)
        ||  (strcmp(s, "BACKBONE") == 0)) {
      type = PROP_NODE; val = MC_PROP;
   }
   else if ((strcmp(s, "SIDECHAIN") == 0) || (strcmp(s, "SC") == 0)
        ||  (strcmp(s, "BASE") == 0)) {
      type = PROP_NODE; val = SC_PROP;
   }
   else if (strcmp(s, "NITROGEN") == 0) {
      type = PROP_NODE; val = N_PROP;
   }
   else if (strcmp(s, "CARBON") == 0) {
      type = PROP_NODE; val = C_PROP;
   }
   else if (strcmp(s, "OXYGEN") == 0) {
      type = PROP_NODE; val = O_PROP;
   }
   else if (strcmp(s, "SULFUR") == 0) {
      type = PROP_NODE; val = S_PROP;
   }
   else if (strcmp(s, "PHOSPHORUS") == 0) {
      type = PROP_NODE; val = P_PROP;
   }
   else if (strcmp(s, "HYDROGEN") == 0) {
      type = PROP_NODE; val = H_PROP;
   }
   else if ((strcmp(s, "H2O") == 0) || (strcmp(s, "WATER") == 0)) {
      type = PROP_NODE; val = WATER_PROP;
   }
   else if ((strcmp(s, "HET") == 0) || (strcmp(s, "HETATM") == 0)) {
      type = PROP_NODE; val = HET_PROP;
   }
   else if ((strcmp(s, "PROTEIN") == 0) || (strcmp(s, "PEPTIDE") == 0)) {
      type = PROP_NODE; val = PROT_PROP;
   }
   else if ((strcmp(s, "DNA") == 0) || (strcmp(s, "RNA") == 0)) {
      type = PROP_NODE; val = DNA_PROP;
   }
   else if (strcmp(s, "METAL") == 0) {
      type = PROP_NODE; val = METAL_PROP;
   }
   else if (strcmp(s, "METHYLENE") == 0) {
      type = PROP_NODE; val = METHYLENE_PROP;
   }/*longer name first*/
      else if (strcmp(s, "METHYL") == 0) {
      type = PROP_NODE; val = METHYL_PROP;
   }
   else if (strcmp(s, "DONOR") == 0) {
      type = PROP_NODE; val = DONOR_PROP;
   }
   else if (strcmp(s, "ACCEPTOR") == 0) {
      type = PROP_NODE; val = ACCEPTOR_PROP;
   }
   else if ((strcmp(s, "NONPOLAR") == 0) || (strcmp(s, "HYDROPHOBIC") == 0)
        ||  (strcmp(s, "HPHOBIC")== 0)|| (strcmp(s, "PHOBIC") == 0)) {
      type = PROP_NODE; val = RHPHOBIC_PROP;
   }
   else if ((strcmp(s, "POLAR") == 0) || (strcmp(s, "HYDROPHILIC") == 0)
        ||  (strcmp(s, "HPHILIC")== 0)|| (strcmp(s, "PHILIC") == 0)) {
      type = PROP_NODE; val = RHPHILIC_PROP;
   }
   else if ((strcmp(s, "CHARGED") == 0) || (strcmp(s, "CHARGE") == 0)) {
      type = PROP_NODE; val = RCHARGED_PROP;
   }
   else if (strcmp(s, "AROMATIC") == 0) {
      type = PROP_NODE; val = AROMATIC_PROP;
   }
   else if ((strcmp(s, "GLY") == 0) || (strcmp(s, "G") == 0)) {
      type = RTYPE_NODE;
      val = lookup("GLY"); if(val == 0) val = insert("GLY", ID_TOK);
   }
   else if ((strcmp(s, "ALA") == 0) || (strcmp(s, "A") == 0)) {
      type = RTYPE_NODE;
      val = lookup("ALA"); if(val == 0) val = insert("ALA", ID_TOK);
   }
   else if ((strcmp(s, "VAL") == 0) || (strcmp(s, "V") == 0)) {
      type = RTYPE_NODE;
      val = lookup("VAL"); if(val == 0) val = insert("VAL", ID_TOK);
   }
   else if ((strcmp(s, "PHE") == 0) || (strcmp(s, "F") == 0)) {
      type = RTYPE_NODE;
      val = lookup("PHE"); if(val == 0) val = insert("PHE", ID_TOK);
   }
   else if ((strcmp(s, "PRO") == 0) || (strcmp(s, "P") == 0)) {
      type = RTYPE_NODE;
      val = lookup("PRO"); if(val == 0) val = insert("PRO", ID_TOK);
   }
   else if ((strcmp(s, "MET") == 0) || (strcmp(s, "M") == 0)) {
      type = RTYPE_NODE;
      val = lookup("MET"); if(val == 0) val = insert("MET", ID_TOK);
   }
   else if ((strcmp(s, "ILE") == 0) || (strcmp(s, "I") == 0)) {
      type = RTYPE_NODE;
      val = lookup("ILE"); if(val == 0) val = insert("ILE", ID_TOK);
   }
   else if ((strcmp(s, "LEU") == 0) || (strcmp(s, "L") == 0)) {
      type = RTYPE_NODE;
      val = lookup("LEU"); if(val == 0) val = insert("LEU", ID_TOK);
   }
   else if ((strcmp(s, "ASP") == 0) || (strcmp(s, "D") == 0)) {
      type = RTYPE_NODE;
      val = lookup("ASP"); if(val == 0) val = insert("ASP", ID_TOK);
   }
   else if ((strcmp(s, "GLU") == 0) || (strcmp(s, "E") == 0)) {
      type = RTYPE_NODE;
      val = lookup("GLU"); if(val == 0) val = insert("GLU", ID_TOK);
   }
   else if ((strcmp(s, "LYS") == 0) || (strcmp(s, "K") == 0)) {
      type = RTYPE_NODE;
      val = lookup("LYS"); if(val == 0) val = insert("LYS", ID_TOK);
   }
   else if ((strcmp(s, "ARG") == 0) || (strcmp(s, "R") == 0)) {
      type = RTYPE_NODE;
      val = lookup("ARG"); if(val == 0) val = insert("ARG", ID_TOK);
   }
   else if ((strcmp(s, "SER") == 0) || (strcmp(s, "S") == 0)) {
      type = RTYPE_NODE;
      val = lookup("SER"); if(val == 0) val = insert("SER", ID_TOK);
   }
   else if ((strcmp(s, "THR") == 0) || (strcmp(s, "T") == 0)) {
      type = RTYPE_NODE;
      val = lookup("THR"); if(val == 0) val = insert("THR", ID_TOK);
   }
   else if ((strcmp(s, "TYR") == 0) || (strcmp(s, "Y") == 0)) {
      type = RTYPE_NODE;
      val = lookup("TYR"); if(val == 0) val = insert("TYR", ID_TOK);
   }
   else if ((strcmp(s, "HIS") == 0) || (strcmp(s, "H") == 0)) {
      type = RTYPE_NODE;
      val = lookup("HIS"); if(val == 0) val = insert("HIS", ID_TOK);
   }
   else if ((strcmp(s, "CYS") == 0) || (strcmp(s, "C") == 0)) {
      type = RTYPE_NODE;
      val = lookup("CYS"); if(val == 0) val = insert("CYS", ID_TOK);
   }
   else if ((strcmp(s, "ASN") == 0) || (strcmp(s, "N") == 0)) {
      type = RTYPE_NODE;
      val = lookup("ASN"); if(val == 0) val = insert("ASN", ID_TOK);
   }
   else if ((strcmp(s, "GLN") == 0) || (strcmp(s, "Q") == 0)) {
      type = RTYPE_NODE;
      val = lookup("GLN"); if(val == 0) val = insert("GLN", ID_TOK);
   }
   else if ((strcmp(s, "TRP") == 0) || (strcmp(s, "W") == 0)) {
      type = RTYPE_NODE;
      val = lookup("TRP"); if(val == 0) val = insert("TRP", ID_TOK);
   }
   else if (strcmp(s, "ACE") == 0) {
      type = RTYPE_NODE;
      val = lookup("ACE"); if(val == 0) val = insert("ACE", ID_TOK);
   }
   else if (strcmp(s, "ASX") == 0) {
      type = RTYPE_NODE;
      val = lookup("ASX"); if(val == 0) val = insert("ASX", ID_TOK);
   }
   else if (strcmp(s, "GLX") == 0) {
      type = RTYPE_NODE;
      val = lookup("GLX"); if(val == 0) val = insert("GLX", ID_TOK);
   }
   else if (strcmp(s, "MSE") == 0) {
      type = RTYPE_NODE;
      val = lookup("MSE"); if(val == 0) val = insert("MSE", ID_TOK);
   }
   else if (strcmp(s, "PCA") == 0) {
      type = RTYPE_NODE;
      val = lookup("PCA"); if(val == 0) val = insert("PCA", ID_TOK);
   }
   else if (strcmp(s, "NH2") == 0) {
      type = RTYPE_NODE;
      val = lookup("NH2"); if(val == 0) val = insert("NH2", ID_TOK);
   }
   else if (strcmp(s, "NME") == 0) {
      type = RTYPE_NODE;
      val = lookup("NME"); if(val == 0) val = insert("NME", ID_TOK);
   }
   else if (strcmp(s, "FOR") == 0) {
      type = RTYPE_NODE;
      val = lookup("FOR"); if(val == 0) val = insert("FOR", ID_TOK);
   }
   else if ((strncmp(s, "RES", 3) == 0) && (isalnum(s[3])||(s[3]=='_'))) {
      type = RTYPE_NODE;
      for(tmpstr = s+3; *tmpstr; tmpstr++) { /* substitute blank for _ */
	 if (*tmpstr == '_') {*tmpstr = ' ';}
      }
      val = lookup(s+3); if(val == 0) val = insert(s+3, ID_TOK);
   }
   else if ((strncmp(s, "OLT", 3) == 0) && isdigit(s[3])) {
      type = OCC_LT_NODE; val = parseInteger(s,3,10);
   }
   else if ((strncmp(s, "OGT", 3) == 0) && isdigit(s[3])) {
      type = OCC_GT_NODE; val = parseInteger(s,3,10);
   }
   else if ((strncmp(s, "BLT", 3) == 0) && isdigit(s[3])) {
      type = B_LT_NODE; val = parseInteger(s,3,10);
   }
   else if ((strncmp(s, "BGT", 3) == 0) && isdigit(s[3])) {
      type = B_GT_NODE; val = parseInteger(s,3,10);
   }
   else if ((strncmp(s, "INS", 3) == 0) && (isalnum(s[3])||(s[3]=='_'))) {
      if (s[3] == '_') {s[3] = ' ';}
      type = INS_NODE; val = (unsigned char)s[3];
   }
   else {
        char msg[200];
	sprintf(msg,"syntax: unknown identifier: %s", s);
	errmsg(msg);
	ParseError = TRUE;
        type = FALSE_NODE;
        val = 0;
   }

   lhs = makeTerminal(type, val);
   matchTok(ID_TOK);
   return lhs;
}

pattern* distItem() {
   pattern *f;
   float dst = 0.0, px = 0.0, py = 0.0, pz = 0.0, signFactor = 0.0;
   int ok = 1;

/* within */
   f = makeTerminal(DIST_NODE, 0);
   matchTok(WITHIN_TOK);
/* distance */
   if (ok && lookahead == REALNUM_TOK) {
      dst = realtokval;
      matchTok(REALNUM_TOK);
   }
   else if (ok && lookahead == NUM_TOK) {
      dst = tokenval;
      matchTok(NUM_TOK);
   }
   else { ok = 0; }
/* of */
   if (ok && lookahead == ID_TOK
          && strcmp(SymTable[tokenval].lexptr, "OF") == 0) {
      matchTok(ID_TOK);
   }
   else { ok = 0; }
/* x */
   signFactor = 1.0;
   if (ok && lookahead == '-') {
      matchTok('-');
      signFactor = -1.0;
   }

   if (ok && lookahead == REALNUM_TOK) {
      px = signFactor * realtokval;
      matchTok(REALNUM_TOK);
   }
   else if (ok && lookahead == NUM_TOK) {
      px = signFactor * tokenval;
      matchTok(NUM_TOK);
   }
   else { ok = 0; }
/* , */
   if (ok && lookahead == ',') {
      matchTok(',');
   }
   else { ok = 0; }
/* y */
   signFactor = 1.0;
   if (ok && lookahead == '-') {
      matchTok('-');
      signFactor = -1.0;
   }

   if (ok && lookahead == REALNUM_TOK) {
      py = signFactor * realtokval;
      matchTok(REALNUM_TOK);
   }
   else if (ok && lookahead == NUM_TOK) {
      py = signFactor * tokenval;
      matchTok(NUM_TOK);
   }
   else { ok = 0; }
/* , */
   if (ok && lookahead == ',') {
      matchTok(',');
   }
   else { ok = 0; }
/* z */
   signFactor = 1.0;
   if (ok && lookahead == '-') {
      matchTok('-');
      signFactor = -1.0;
   }

   if (ok && lookahead == REALNUM_TOK) {
      pz = signFactor * realtokval;
      matchTok(REALNUM_TOK);
   }
   else if (ok && lookahead == NUM_TOK) {
      pz = signFactor * tokenval;
      matchTok(NUM_TOK);
   }
   else { ok = 0; }

/* info is packed into the floating-point vector */
   if (ok) {
      f->fvec[0] = dst;
      f->fvec[1] = px;
      f->fvec[2] = py;
      f->fvec[3] = pz;
   }
   else {
      char msg[100];
      describeLookahead(msg,"syntax: parsing WITHIN: ");
      errmsg(msg);
      ParseError = TRUE;
   }
   return f;
}

pattern* makeNode(pattern *lhs, pattern *rhs, int type) {
   pattern *node;
   node = (pattern *)malloc(sizeof(pattern));
   if (node) {
      node->lhs = lhs;
      node->rhs = rhs;
      node->type = type;
      node->val  = 0;
      node->fvec[0] = 0.0;
      node->fvec[1] = 0.0;
      node->fvec[2] = 0.0;
      node->fvec[3] = 0.0;
   }
   else {
      halt("out of memory in makeNode()");
   }
   return node;
}

pattern* makeTerminal(int type, int val) {
   pattern *term;
   term = (pattern *)malloc(sizeof(pattern));
   if (term) {
      term->lhs  = NULL;
      term->rhs  = NULL;
      term->type = type;
      term->val  = val;
      term->fvec[0] = 0.0;
      term->fvec[1] = 0.0;
      term->fvec[2] = 0.0;
      term->fvec[3] = 0.0;
   }
   else {
      halt("out of memory in makeTerminal()");
   }
   return term;
}

void matchTok(int t) {
   if (lookahead == t) {
      lookahead = lexan();
   }
   else {
        char msg[100];
        describeLookahead(msg,"improper syntax ");
	errmsg(msg);
        sprintf(msg,"expect match code: %d", t);
        errmsg(msg);
	ParseError = TRUE;
   }
}

char* lexString(int p) {

   return SymTable[p].lexptr;
}

int lexan() {
   int p, t, rc, cnt;
   
   for (;;) {
      t = toupper(*inputStr++);
      if (t == ' ' || t == '\t' || t == '\n' || t == '\r') {
	 /* skip over whitespace  */
      }
      else if (isdigit(t) || t == '.') {
	 int insCode = 0;
	 int fractional = (t == '.');
	 lexbuf[cnt = 0] = t;
	 while ( isdigit(*inputStr)
	    || (  (*inputStr == '.' ||
	            isalpha(*inputStr) || *inputStr == '_')
	      && ! (insCode || fractional)) ) {

	         if (*inputStr == '.')   { fractional = 1;               }
	    else if (isalpha(*inputStr)) { insCode = toupper(*inputStr); }
	    else if (*inputStr == '_')   { insCode = ' ';                }
	    if (cnt >= LXBUFSZ-1) { halt("number too long"); }
	    lexbuf[++cnt] = toupper(*inputStr++);
	 }
	 lexbuf[++cnt] = '\0';

	 if (fractional) {
	    realtokval = parseReal(lexbuf, 0, cnt);
	    tokenval = 0;
	    rc = REALNUM_TOK;
	 }
	 else { /* need to save insCode */
	    tokenval = parseInteger(lexbuf, 0, cnt);
            cnt = 0;
	    if (insCode) { lexbuf[cnt++] = insCode; }
	    lexbuf[cnt] = '\0';
	    rc = NUM_TOK;
	 }
	 break;
      }
      else if (t == 'N' && toupper(inputStr[0]) == 'O'
           && toupper(inputStr[1]) == 'T' && !isalnum(inputStr[2]) ) {
	 inputStr++;
	 inputStr++;
	 tokenval = 0;
	 rc = NOT_TOK;
	 break;
      }
      else if (t == 'W' && toupper(inputStr[0]) == 'I'
           && toupper(inputStr[1]) == 'T' && toupper(inputStr[2]) == 'H'
           && toupper(inputStr[3]) == 'I' && toupper(inputStr[4]) == 'N'
	   && !isalnum(inputStr[5]) ) {
	 inputStr += 5;
	 tokenval = 0;
	 rc = WITHIN_TOK;
	 break;
      }
      else if (isalpha(t)) {
	 lexbuf[cnt = 0] = t;
	 while (isalnum(*inputStr) || (*inputStr == '_')
	     || (*inputStr == '*') || (*inputStr == '\'')
	     || (*inputStr == '"')) {
	    if (cnt >= LXBUFSZ-1) { halt("symbol too long"); }
	    lexbuf[++cnt] = toupper(*inputStr++);
	 }
	 lexbuf[++cnt] = '\0';

	 p = lookup(lexbuf);
	 if (p == 0) {
	    p = insert(lexbuf, ID_TOK);
	 }
	 tokenval = p;
	 rc = ID_TOK;
	 break;
      }
      else if (t == '\0') {
	 tokenval = 0;
	 rc = DONE_TOK;
	 break;
      }
      else {
	 tokenval = 0;
	 rc = t;
	 break;
      }
   }
   return rc;
}

int lookup(char *s) {
   int p;
   for (p = LastSymEntry; p > 0; p--) {
      if (strcmp(SymTable[p].lexptr, s) == 0) {
	 return p;
      }
   }
   return 0;
}

int insert(char *s, int tok) {
   int len;
   len = strlen(s);

   if (LastSymEntry + 1 >= LEXSYMMAX) {
      halt("symbol table full");
   }
   if (LastLexChar + len + 1 >= LEXSTRINGMAX) {
      halt("symbol name array full");
   }
   LastSymEntry++;

   SymTable[LastSymEntry].token = tok;
   SymTable[LastSymEntry].lexptr = &Lexemes[LastLexChar + 1];
   LastLexChar += len + 1;
   strcpy(SymTable[LastSymEntry].lexptr, s);

   return LastSymEntry;
}

int  modelInPattern(pattern *pat) 
{
 static  int model=0;
 static  int ireturn=0;

   if(pat == NULL) {model=0;}
   else
   {
      switch(pat->type) 
      {
         case OR_NODE:
         case AND_NODE:
            ireturn = modelInPattern(pat->rhs);
         case NOT_NODE:
            ireturn = modelInPattern(pat->lhs);
         break;
         case MODEL_NODE:
            model = pat->val;
         break;
      }/*switch*/
   }
   if(model>0){ireturn = model;}
   return(ireturn);
}/*modelInPattern()*/

void printPattern(FILE *outf, pattern *pat) {
   if (pat) {
      recPrint(outf, pat);
      fprintf(outf, "\n");
   }
}

void recPrint(FILE *outf, pattern *pat) {

   if (pat == NULL) return;

   switch(pat->type) {
   case OR_NODE:
      fprintf(outf, "(");
      recPrint(outf, pat->lhs);
      fprintf(outf, " or ");
      recPrint(outf, pat->rhs);
      fprintf(outf, ")");
      break;
   case AND_NODE:
      fprintf(outf, "(");
      recPrint(outf, pat->lhs);
      fprintf(outf, " and ");
      recPrint(outf, pat->rhs);
      fprintf(outf, ")");
      break;
   case NOT_NODE:
      fprintf(outf,"not(");
      recPrint(outf, pat->lhs);
      fprintf(outf,")");
      break;
   case RANGE_NODE:
      fprintf(outf,"%d <= resno <= %d", pat->lhs->val, pat->rhs->val);
      break;
   case FILE_NODE:
      fprintf(outf,"file=%d", pat->val);
      break;
   case MODEL_NODE:
      fprintf(outf,"model=%d", pat->val);
      break;
   case CHAIN_NODE:
      fprintf(outf,"chain=\"%s%s\"", lexString(pat->val), 
                              ((strlen(lexString(pat->val))<3)?"*":""));
      break;
   case ALT_NODE:
      fprintf(outf,"alt='%c'", pat->val);
      break;
   case ANAME_NODE:
      fprintf(outf,"atom=\"%s%s\"", lexString(pat->val),
			      ((strlen(lexString(pat->val))<4)?"*":""));
      break;
   case SEGID_NODE:
      fprintf(outf,"segID=\"%s%s\"", lexString(pat->val),
			      ((strlen(lexString(pat->val))<4)?"*":""));
      break;
   case RES_NODE:
      fprintf(outf,"resno=%d", pat->val);
      break;
   case RTYPE_NODE:
      fprintf(outf,"res=\"%s%s\"", lexString(pat->val),
			      ((strlen(lexString(pat->val))<3)?"*":""));
      break;
   case DIST_NODE:
      fprintf(outf,"atoms within %.1f of %.3f, %.3f, %.3f", pat->fvec[0],
	 pat->fvec[1], pat->fvec[2], pat->fvec[3]);
      break;
   case PROP_NODE:
      switch(pat->val) {
      case ALPHA_PROP: fprintf(outf,"alpha");           break;
      case  BETA_PROP: fprintf(outf,"beta");            break;
      case    MC_PROP: fprintf(outf,"mainchain");       break;
      case    SC_PROP: fprintf(outf,"sidechain");       break;
      case     N_PROP: fprintf(outf,"Nitrogen");        break;
      case     C_PROP: fprintf(outf,"Carbon");          break;
      case     O_PROP: fprintf(outf,"Oxygen");          break;
      case     S_PROP: fprintf(outf,"Sulfur");          break;
      case     P_PROP: fprintf(outf,"Phosphorus");      break;
      case     H_PROP: fprintf(outf,"Hydrogen");        break;
      case WATER_PROP: fprintf(outf,"H2O");             break;
      case   HET_PROP: fprintf(outf,"het-group");       break;
      case  PROT_PROP: fprintf(outf,"protein");         break;
      case   DNA_PROP: fprintf(outf,"dna/rna");         break;
      case METAL_PROP: fprintf(outf,"metal");	        break;
      case METHYL_PROP: fprintf(outf,"methyl");	        break;
      case METHYLENE_PROP: fprintf(outf,"methylene");   break; /*20111210dcr*/
      case DONOR_PROP: fprintf(outf,"donor");	        break;
      case ACCEPTOR_PROP: fprintf(outf,"acceptor");     break;
      case RHPHOBIC_PROP:
      fprintf(outf,"A,V,L,I,M,C,F,W,Y (but not G,P!)"); break;
      case RHPHILIC_PROP: fprintf(outf,"S,T,N,Q");      break;
      case RCHARGED_PROP: fprintf(outf,"D,E,K,R(+H)");  break;
      case AROMATIC_PROP:fprintf(outf,"aromatic-atoms");break;
      default: fprintf(outf,"property=%d",pat->val);    break;
      }
      break;
   case TRUE_NODE:
      fprintf(outf,"all");
      break;
   case FALSE_NODE:
      fprintf(outf,"none");
      break;
   case OCC_LT_NODE:
      fprintf(outf,"occupancy < %.2f", pat->val/100.0);
      break;
   case OCC_GT_NODE:
      fprintf(outf,"occupancy > %.2f", pat->val/100.0);
      break;
   case B_LT_NODE:
      fprintf(outf,"Bval < %d", pat->val);
      break;
   case B_GT_NODE:
      fprintf(outf,"Bval > %d", pat->val);
      break;
   case INS_NODE:
      fprintf(outf,"insert='%c'", pat->val);
      break;
   case INS_RANGE_NODE:
      if (pat->lhs->type == RES_NODE) {
        fprintf(outf,"%d", pat->rhs->val);
      }
      else if (pat->lhs->type == AND_NODE) {
        fprintf(outf,"%d%c", pat->lhs->lhs->val,
	                     pat->lhs->rhs->val);
      }
      else {
        errmsg("***INTERNAL DEFECT - EXPECT AND***");
        fprintf(outf, "(");
        recPrint(outf, pat->lhs);
        fprintf(outf, ")");
      }
      fprintf(outf," <= resno <= ");
      if (pat->rhs->type == RES_NODE) {
        fprintf(outf,"%d", pat->rhs->val);
      }
      else if (pat->rhs->type == AND_NODE) {
        fprintf(outf,"%d%c", pat->rhs->lhs->val,
	                     pat->rhs->rhs->val);
      }
      else {
        errmsg("***INTERNAL DEFECT - EXPECT AND***");
        fprintf(outf, "(");
        recPrint(outf, pat->rhs);
        fprintf(outf, ")");
      }
      fprintf(outf," /* warning: insert codes ignored in ranges */");
      break;
   default:
      fprintf(outf,"<<bad type: %d>>", pat->type);
   }
}/*recPrint()*/

char* describeLookahead(char formatstring[], char *s) {

   if (lookahead == ID_TOK) {
     sprintf(formatstring,"%s\"%s\"", s, lexString(tokenval));
   }
   else if (lookahead == NUM_TOK) {
     sprintf(formatstring,"%s%d", s, tokenval);
   }
   else if (lookahead == REALNUM_TOK) {
     sprintf(formatstring,"%s%g", s, realtokval);
   }
   else if (lookahead == NOT_TOK) {
     sprintf(formatstring,"%sNOT", s);
   }
   else if (lookahead == WITHIN_TOK) {
     sprintf(formatstring,"%sWITHIN", s);
   }
   else if (lookahead == DONE_TOK) {
     sprintf(formatstring,"%s(end of pattern)", s);
   }
   else {
     sprintf(formatstring,"%s\'%c\'", s, lookahead);
   }
   return formatstring;
}
