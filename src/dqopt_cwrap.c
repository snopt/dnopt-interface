#include "dqopt_cwrap.h"

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dqBegin(dqProblem *prob, char* name, char* prtfile, int summOn) {
  int leniw, lenrw, len;

  init2zeroQ(prob);

  leniw = 500;
  lenrw = 500;

  allocIQ(prob, leniw);
  allocRQ(prob, lenrw);

  prob->name = name;

  len = strlen(prtfile);

  if (summOn != 0) {
    printf(" ==============================\n");
    printf("   DQOPT  C interface  1.0.0   ");
    fflush(stdout);
    //------123456789|123456789|123456789|
  }

  f_dqbegin(prtfile, len, summOn,
	    prob->iw, prob->leniw,
	    prob->rw, prob->lenrw);
  prob->initCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void init2zeroQ(dqProblem *prob) {
  prob->name       = NULL;

  prob->memCalled  = 0;
  prob->initCalled = 0;

  prob->dqLog      = NULL;

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

void allocIQ(dqProblem *prob, int len) {
  prob->leniw = len;
  prob->iw    = calloc(len, sizeof(int));
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void allocRQ(dqProblem *prob, int len) {
  prob->lenrw = len;
  prob->rw    = calloc(len, sizeof(double));
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void reallocIQ(dqProblem *prob, int len) {
  prob->leniw = len;
  prob->iw    = realloc(prob->iw, sizeof(int)*prob->leniw);

  dqSetIntParameter(prob, (char*)"Total int workspace", prob->leniw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void reallocRQ(dqProblem *prob, int len) {
  prob->lenrw = len;
  prob->rw = realloc(prob->rw, sizeof(double)*prob->lenrw);

  dqSetIntParameter(prob, (char*)"Total real workspace", prob->lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dqPrintfile(dqProblem *prob, char *prtname) {
  int len = strlen(prtname);

  assert(prob->initCalled == 1);
  f_dqsetprint(prtname, len,
	       prob->iw, prob->leniw, prob->rw, prob->lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqSpecsfile(dqProblem *prob, char *spcname) {
  int inform;
  int len = strlen(spcname);

  assert(prob->initCalled == 1);
  f_dqspec(spcname, len, &inform,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqSetParameter(dqProblem *prob, char stropt[]) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dqset(stropt, len, &errors,
	  prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqGetParameter(dqProblem *prob, char stropt[], char strout[]) {
  int errors;
  int inlen  = strlen(stropt);
  int outlen = strlen(strout);

  assert(prob->initCalled == 1);
  f_dqgetc(stropt, inlen, strout, outlen, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqSetIntParameter(dqProblem *prob, char stropt[], int opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dqseti(stropt, len, opt, &errors,
	   prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqGetIntParameter(dqProblem *prob, char stropt[], int opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dqgeti(stropt, len, &opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqSetRealParameter(dqProblem *prob, char stropt[], double opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dqsetr(stropt, len, opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqGetRealParameter(dqProblem *prob, char stropt[], double opt) {
  int errors, len = strlen(stropt);

  assert(prob->initCalled == 1);
  f_dqgetr(stropt, len, &opt, &errors,
	    prob->iw, prob->leniw, prob->rw, prob->lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dqSetUserI(dqProblem *prob, int *iu, int leniu) {
  prob->iu    = iu;
  prob->leniu = leniu;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dqSetUserR(dqProblem *prob, double *ru, int lenru) {
  prob->ru    = ru;
  prob->lenru = lenru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dqSetUserspace(dqProblem *prob, int *iu, int leniu,
		    double *ru, int lenru) {
  prob->iu    = iu;
  prob->leniu = leniu;

  prob->ru    = ru;
  prob->lenru = lenru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void dqSetLog(dqProblem *prob, idqLog dqLog) {
  prob->dqLog  = dqLog;
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqopt(dqProblem *prob, int start, int n,
	  int mCon, int nnH, int ncObj,
	  int iObj, double ObjAdd,
	  double *A, int ldA, double *bl, double *bu,
	  double *cObj, double *H, int ldH, dqHx userHx,
	  int *eType, int *state, double *x, double *y,
	  double *objQP, int *nInf, double *sInf) {

  int    inform, dniObj, miniw, minrw, nb;

  assert(prob->initCalled == 1);

  dniObj = iObj+1;

  /*
  if (prob->memCalled == 0) {
    dnSetWorkspace(prob, mLCon, mNCon, n, nnJac, nnObj, dniObj);
  }
  */

  nb = n + mCon;
  f_dqker(start, n, nb, mCon, nnH, ncObj,
	  dniObj, ObjAdd, prob->name,
	  A, ldA, bl, bu, cObj, H, ldH,
	  userHx, prob->dqLog,
	  eType, state, x, y,
	  &inform, &miniw, &minrw,
	  objQP, nInf, sInf,
	  prob->iu, prob->leniu, prob->ru, prob->lenru,
	  prob->iw, prob->leniw, prob->rw, prob->lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

void deleteDQOPT(dqProblem *prob) {
  f_dqend(prob->iw, prob->leniw, prob->rw, prob->lenrw);

  free(prob->iw);
  free(prob->rw);

  init2zeroQ(prob);
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
