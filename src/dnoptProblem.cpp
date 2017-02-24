#include <assert.h>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include "dnopt.h"
#include "dnoptProblem.hpp"

using namespace std;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblem::dnoptProblem() {
  init2zero();
  sprintf(Prob, "%8s", "        ");

  allocI(500);
  allocR(500);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblem::dnoptProblem(const char *name) {
  init2zero();

  sprintf(Prob, "%8s", name);

  allocI(500);
  allocR(500);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblem::~dnoptProblem() {
  f_dnend(iw, leniw, rw, lenrw);

  delete []rw;  delete []iw;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::init2zero() {
  initCalled = 0; memCalled = 0;

  leniw = 0; lenrw = 0;
  iw    = 0; rw    = 0;

  leniu = 0; lenru = 0;
  iu    = 0; ru    = 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::allocI(int aleniw) {
  // Reset work array lengths.
  // Allocate new memory for work arrays.
  leniw = aleniw;
  iw    = new int[leniw];
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::allocR(int alenrw) {
  // Reset work array lengths.
  // Allocate new memory for work arrays.
  lenrw = alenrw;
  rw    = new double[lenrw];
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::reallocI(int aleniw) {
  int  tleniw = leniw;
  int    *tiw = iw;

  // Allocate new memory
  allocI(aleniw);

  // Copy old workspace into new.
  int mleniw = leniw < tleniw ? leniw : tleniw;
  memcpy(iw, tiw, mleniw*sizeof(int));

  // Delete temporary work arrays
  delete []tiw;

  setIntParameter((char*)"Total int workspace", leniw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::reallocR(int alenrw) {
  int  tlenrw = lenrw;
  double *trw = rw;

  // Allocate new memory
  allocR(alenrw);

  // Copy old workspace into new.
  int mlenrw = lenrw < tlenrw ? lenrw : tlenrw;
  memcpy(rw, trw, mlenrw*sizeof(double));

  // Delete temporary work arrays
  delete []trw;

  setIntParameter((char*)"Total real workspace   ", lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::setProbName(const char *name) {
  sprintf(Prob, "%8s", name);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::setPrintFile(const char *prtname) {
  assert(initCalled == 1);

  int len = strlen(prtname);
  f_dnsetprint(prtname, len, iw, leniw, rw, lenrw);
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dnoptProblem::setParameter(const char *stropt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);
  f_dnset(stropt, stropt_len, &errors, iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dnoptProblem::getParameter(const char *stroptin, char *stroptout) {
  assert(initCalled == 1);

  int errors;
  int stroptin_len  = strlen(stroptin);
  int stroptout_len = strlen(stroptout);

  f_dngetc(stroptin, stroptin_len, stroptout, stroptout_len,
	    &errors, iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dnoptProblem::setIntParameter(const char *stropt, int opt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);

  f_dnseti(stropt, stropt_len, opt, &errors,
	    iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dnoptProblem::getIntParameter(const char *stropt, int &opt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);
  f_dngeti(stropt, stropt_len, &opt, &errors,
	    iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dnoptProblem::setRealParameter(const char *stropt, double opt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);
  f_dnsetr(stropt, stropt_len, opt, &errors,
	    iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dnoptProblem::getRealParameter(const char *stropt, double &opt) {
  assert(initCalled == 1);

  int errors, stropt_len = strlen(stropt);
  f_dngetr(stropt, stropt_len, &opt, &errors,
	    iw, leniw, rw, lenrw);
  return errors;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::setUserI(int *aiu, int aleniu) {
  leniu = aleniu;
  iu    = aiu;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::setUserR(double *aru, int alenru) {
  lenru = alenru;
  ru    = aru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblem::setUserspace(int*aiu,     int aleniu,
				double *aru, int alenru) {
  leniu = aleniu;
  iu    = aiu;

  lenru = alenru;
  ru    = aru;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblemNLP::dnoptProblemNLP() {
  init2zero();
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblemNLP::dnoptProblemNLP(const char *name) : dnoptProblem(name) {
  init2zero();
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblemNLP::~dnoptProblemNLP() {
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblemNLP::init2zero() {
  dnLog = 0;  dnLogQP = 0;  dnSTOP = 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblemNLP::initialize(const char *prtfile, int summOn) {
  int len = strlen(prtfile);

  if (summOn != 0) {
    printf(" ==============================\n");
    printf("  DNOPT  C++ interface  1.0.0  ");
    fflush(stdout);
    //------123456789|123456789|123456789|
  }

  f_dnbegin(prtfile, len, summOn, iw, leniw, rw, lenrw);
  initCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dnoptProblemNLP::setSpecsFile(const char *specname) {
  assert(initCalled == 1);

  int inform, len = strlen(specname);

  f_dnspec(specname, len, &inform, iw, leniw, rw, lenrw);
  if (inform != 101) {
    printf("Warning: unable to find specs file %s \n", specname);
  }
  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblemNLP::setLog(idnLog adnLog, idnLogQP adnLogQP) {
  dnLog   = adnLog;  dnLogQP = adnLogQP;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblemNLP::setSTOP(idnSTOP adnSTOP) {
  dnSTOP = adnSTOP;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblemB::dnoptProblemB() {
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblemB::dnoptProblemB(const char *name) : dnoptProblemNLP(name) {
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dnoptProblemB::~dnoptProblemB() {
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dnoptProblemB::setWorkspace(int mLCon, int mNCon, int n,
				 int nnJac, int nnObj, int iObj) {
  assert(initCalled == 1);

  int dniObj, inform, miniw, minrw;

  dniObj = iObj + 1;

  f_dnmem(&inform, mLCon, mNCon, n, nnJac, nnObj, dniObj,
	  &miniw, &minrw, iw, leniw, rw, lenrw);

  if (miniw > leniw) { reallocI (miniw); }
  if (minrw > lenrw) { reallocR (minrw); }

  memCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dnoptProblemB::solve(int starttype, int n, int mLCon, int mNCon,
			 int nnJac, int nnObj, int iObj, double ObjAdd,
			 dnCon funcon, dnObj funobj,
			 int *state, double *A, int ldA,
			 double *bl, double *bu,
			 double *J, int ldJ, double *H, int ldH,
			 double &objective, int &nInf, double &sInf,
			 double *x, double *y) {
  assert(initCalled == 1);

  int  dniObj, inform, miniw, minrw, nb;
  dnHes funhes = 0;

  if (memCalled == 0) {
    setWorkspace(mLCon, mNCon, n, nnJac, nnObj, iObj);
  }

  nb     = n + mLCon + mNCon;
  dniObj = iObj+1;

  f_dnker(starttype, n, nb, mLCon, mNCon, nnJac, nnObj,
	  Prob, dniObj, ObjAdd,
	  funcon, funobj, funhes, dnLog, dnLogQP, dnSTOP,
	  state, A, ldA, bl, bu, J, ldJ, H, ldH,
	  &objective, &nInf, &sInf, x, y,
	  &inform, &miniw, &minrw, iu, leniu, ru, lenru,
	  iw, leniw, rw, lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dnoptProblemB::solveh(int starttype, int n, int mLCon, int mNCon,
			  int nnJac, int nnObj, int iObj, double ObjAdd,
			  dnCon funcon, dnObj funobj, dnHes funhes,
			  int *state, double *A, int ldA,
			  double *bl, double *bu,
			  double *J, int ldJ, double *H, int ldH,
			  double &objective, int &nInf, double &sInf,
			  double *x, double *y) {
  assert(initCalled == 1);

  int  dniObj, inform, miniw, minrw, nb;

  if (memCalled == 0) {
    setWorkspace(mLCon, mNCon, n, nnJac, nnObj, iObj);
  }

  nb     = n + mLCon + mNCon;
  dniObj = iObj+1;

  f_dnker(starttype, n, nb, mLCon, mNCon, nnJac, nnObj,
	  Prob, dniObj, ObjAdd,
	  funcon, funobj, funhes, dnLog, dnLogQP, dnSTOP,
	  state, A, ldA, bl, bu, J, ldJ, H, ldH,
	  &objective, &nInf, &sInf, x, y,
	  &inform, &miniw, &minrw, iu, leniu, ru, lenru,
	  iw, leniw, rw, lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dqoptProblem::dqoptProblem() {
  init2zero();
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dqoptProblem::dqoptProblem(const char *name) : dnoptProblem(name) {
  init2zero();
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
dqoptProblem::~dqoptProblem() {
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dqoptProblem::init2zero() {
  dqLog = 0;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dqoptProblem::initialize(const char*prtfile, int summOn) {
  int len = strlen(prtfile);

  if (summOn != 0) {
    printf(" ==============================\n");
    printf("  DQOPT  C++ interface  1.0.0  ");
    fflush(stdout);
    //------123456789|123456789|123456789|
  }

  f_dqbegin(prtfile, len, summOn, iw, leniw, rw, lenrw);
  initCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int dqoptProblem::setSpecsFile(const char *specname) {
  assert(initCalled == 1);

  int inform, len = strlen(specname);
  f_dqspec(specname, len, &inform, iw, leniw, rw, lenrw);
  if (inform != 101) {
    printf("Warning: unable to find specs file %s \n", specname);
  }

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dqoptProblem::setLog(idqLog adqLog) {
  dqLog  = adqLog;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void dqoptProblem::setWorkspace() {
  assert(initCalled == 1);


  memCalled = 1;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int dqoptProblem::solve(int starttype, int mCon, int n, int nnH, int ncObj,
			int iObj, double ObjAdd,
			double *A, int ldA, double *bl, double *bu,
			double *cObj, double *H, int ldH, dqHx qpHx,
			int *eType, int *state, double *x, double *y,
			double &objective, int &nInf, double &sInf) {

  assert(initCalled == 1);

  int dniObj, inform, miniw, minrw, nb;

  if (memCalled == 0) { setWorkspace(); }

  dniObj = iObj+1;
  nb     = n + mCon;

  f_dqker(starttype, n, nb, mCon, nnH,
	  ncObj, dniObj, ObjAdd, Prob,
	  A, ldA, bl, bu, cObj, H, ldH,
	  qpHx, dqLog,
	  eType, state, x, y,
	  &inform, &miniw, &minrw,
	  &objective, &nInf, &sInf,
	  iu, leniu, ru, lenru,
	  iw, leniw, rw, lenrw);

  return inform;
}

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
