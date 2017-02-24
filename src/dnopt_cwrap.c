#include "dnopt_cwrap.h"

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnBegin(dnProblem *prob, char* name, char* prtfile, int summOn) {
  int leniw, lenrw, len;

  dninit2zero(prob);

  leniw = 500;
  lenrw = 500;

  dnallocI(prob, leniw);
  dnallocR(prob, lenrw);

  prob->name = name;

  len = strlen(prtfile);

  if (summOn != 0) {
    printf(" ==============================\n");
    printf("   DNOPT  C interface  1.0.0   ");
    fflush(stdout);
    //------123456789|123456789|123456789|
  }

  f_dnbegin(prtfile, len, summOn,
	    prob->iw, prob->leniw,
	    prob->rw, prob->lenrw);
  prob->initCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dninit2zero(dnProblem *prob) {
  prob->name       = NULL;

  prob->memCalled  = 0;
  prob->initCalled = 0;

  prob->dnLog      = NULL;
  prob->dnLogQP    = NULL;
  prob->dnSTOP     = NULL;

  prob->leniw      = 0;
  prob->lenrw      = 0;
  prob->iw         = NULL;
  prob->rw         = NULL;

  prob->leniu      = 0;
  prob->lenru      = 0;
  prob->iu         = NULL;
  prob->ru         = NULL;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnallocI(dnProblem *prob, int len) {
  prob->leniw = len;
  prob->iw    = calloc(len, sizeof(int));
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnallocR(dnProblem *prob, int len) {
  prob->lenrw = len;
  prob->rw    = calloc(len, sizeof(double));
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnreallocI(dnProblem *prob, int len) {
  prob->leniw = len;
  prob->iw    = realloc(prob->iw, sizeof(int)*prob->leniw);

  dnSetIntParameter(prob, (char*)"Total int workspace", prob->leniw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnreallocR(dnProblem *prob, int len) {
  prob->lenrw = len;
  prob->rw = realloc(prob->rw, sizeof(double)*prob->lenrw);

  dnSetIntParameter(prob, (char*)"Total real workspace", prob->lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnPrintfile(dnProblem *prob, char *prtname) {
  int len = strlen(prtname);

  assert(prob->initCalled == 1);
  f_dnsetprint(prtname, len,
	       prob->iw, prob->leniw, prob->rw, prob->lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnSpecsfile(dnProblem *prob, char *spcname) {
  int inform;
  int len = strlen(spcname);

  assert(prob->initCalled == 1);
  f_dnspec(spcname, len, &inform,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnSetParameter(dnProblem *prob, char stropt[]) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dnset(stropt, len, &errors,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnGetParameter(dnProblem *prob, char stropt[], char strout[]) {
  int errors;
  int inlen  = strlen(stropt);
  int outlen = strlen(strout);

  assert(prob->initCalled == 1);
  f_dngetc(stropt, inlen, strout, outlen, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnSetIntParameter(dnProblem *prob, char stropt[], int opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dnseti(stropt, len, opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnGetIntParameter(dnProblem *prob, char stropt[], int opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dngeti(stropt, len, &opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnSetRealParameter(dnProblem *prob, char stropt[], double opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dnsetr(stropt, len, opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnGetRealParameter(dnProblem *prob, char stropt[], double opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dngetr(stropt, len, &opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnSetUserI(dnProblem *prob, int *iu, int leniu) {
  prob->iu    = iu;
  prob->leniu = leniu;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnSetUserR(dnProblem *prob, double *ru, int lenru) {
  prob->ru    = ru;
  prob->lenru = lenru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnSetUserspace(dnProblem *prob, int *iu, int leniu,
		    double *ru, int lenru) {
  prob->iu    = iu;
  prob->leniu = leniu;

  prob->ru    = ru;
  prob->lenru = lenru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnSetWorkspace(dnProblem *prob, int mLCon, int mNCon, int n,
		    int nnJac, int nnObj, int iObj) {
  int miniw, minrw, iiObj, inform;
  int memGuess = 0;

  assert(prob->initCalled == 1);

  iiObj = iObj+1;

  f_dnmem (&inform, mLCon, mNCon, n, nnJac,
	   nnObj, iiObj, &miniw, &minrw,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  if (miniw > prob->leniw) { dnreallocI(prob, miniw); }
  if (minrw > prob->lenrw) { dnreallocR(prob, minrw); }

  prob->memCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnSetLog(dnProblem *prob, idnLog dnLog, idnLogQP dnLogQP) {
  prob->dnLog   = dnLog;
  prob->dnLogQP = dnLogQP;
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dnSetSTOP(dnProblem *prob, idnSTOP dnSTOP) {
  prob->dnSTOP = dnSTOP;
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnopt(dnProblem *prob, int start,
	  int n, int mLCon, int mNCon, int nnJac, int nnObj,
	  int iObj, double ObjAdd,
	  dnCon funcon, dnObj funobj,
	  int *state, double *A, int ldA, double *bl, double *bu,
	  double *J, int ldJ, double *H, int ldH,
	  double *objective, int *nS, int* nInf, double* sInf,
	  double *x, double *y) {

  int    inform, dniObj, miniw, minrw, nb;
  dnHes  funhes;

  assert(prob->initCalled == 1);

  dniObj = iObj+1;
  nb     = n + mLCon + mNCon;

  if (prob->memCalled == 0) {
    dnSetWorkspace(prob, mLCon, mNCon, n, nnJac, nnObj, dniObj);
  }

  funhes = NULL;

  f_dnker(start, n, nb, mLCon, mNCon, nnJac, nnObj,
	  prob->name, dniObj, ObjAdd,
	  funcon, funobj, funhes,
	  prob->dnLog, prob->dnLogQP, prob->dnSTOP,
	  state, A, ldA, bl, bu,
	  J, ldJ, H, ldH,
	  objective, nInf, sInf,
	  x, y,
	  &inform, &miniw, &minrw,
	  prob->iu, prob->leniu, prob->ru, prob->lenru,
	  prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnopth(dnProblem *prob, int start,
	   int n, int mLCon, int mNCon, int nnJac, int nnObj,
	   int iObj, double ObjAdd,
	   dnCon funcon, dnObj funobj, dnHes funhes,
	   int *state, double *A, int ldA, double *bl, double *bu,
	   double *J, int ldJ, double *H, int ldH,
	   double *objective, int *nS, int* nInf, double* sInf,
	   double *x, double *y) {

  int    inform, dniObj, miniw, minrw, nb;

  assert(prob->initCalled == 1);

  dniObj = iObj+1;
  nb     = n + mLCon + mNCon;

  if (prob->memCalled == 0) {
    dnSetWorkspace(prob, mLCon, mNCon, n, nnJac, nnObj, dniObj);
  }

  f_dnker(start, n, nb, mLCon, mNCon, nnJac, nnObj,
	  prob->name, dniObj, ObjAdd,
	  funcon, funobj, funhes,
	  prob->dnLog, prob->dnLogQP, prob->dnSTOP,
	  state, A, ldA, bl, bu,
	  J, ldJ, H, ldH,
	  objective, nInf, sInf,
	  x, y,
	  &inform, &miniw, &minrw,
	  prob->iu, prob->leniu, prob->ru, prob->lenru,
	  prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void deleteDNOPT(dnProblem *prob) {
  f_dnend(prob->iw, prob->leniw, prob->rw, prob->lenrw);

  free(prob->iw);
  free(prob->rw);

  dninit2zero(prob);
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
