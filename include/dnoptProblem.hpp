#ifndef DNOPTPROBLEM_H
#define DNOPTPROBLEM_H

#include "dnopt.h"

/* File dnoptProblem.hpp
 *   C++ interface for DNOPT and DQOPT
 */

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class dnoptProblem {
protected:
  dnoptProblem();
  dnoptProblem(const char*name);
  ~dnoptProblem();

  void init2zero();

  char    Prob[30];

  int     initCalled, memCalled;

  int     leniw, lenrw;
  double *rw;
  int    *iw;

  int     lenru, leniu;
  double *ru;
  int    *iu;

  void allocI    (int leniw);
  void allocR    (int lenrw);
  void reallocI  (int leniw);
  void reallocR  (int lenrw);
  int  errMsgExit(const char *var);

public:
  void setProbName    (const char *Prob);
  void setPrintFile   (const char *prtname);

  int getParameter    (const char *stroptin, char *stroptout);
  int getIntParameter (const char *stropt,   int    &opt);
  int getRealParameter(const char *stropt,   double &opt);
  int setParameter    (const char *stroptin);
  int setIntParameter (const char *stropt,   int     opt);
  int setRealParameter(const char *stropt,   double  opt);

  void setUserI       (int    *iu, int leniu);
  void setUserR       (double *ru, int lenru);
  void setUserspace   (int    *iu, int leniu,
		       double *ru, int lenru);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class dnoptProblemNLP : public dnoptProblem {
protected:
  dnoptProblemNLP();
  dnoptProblemNLP(const char*name);
  ~dnoptProblemNLP();

  void init2zero();

  idnLog   dnLog;
  idnLogQP dnLogQP;
  idnSTOP  dnSTOP;

public:
  void initialize  (const char *prtfile, int summOn);
  int  setSpecsFile(const char *specname);

  void setLog      (idnLog dnLog, idnLogQP dnLogQP);
  void setSTOP     (idnSTOP dnSTOP);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class dnoptProblemB : public dnoptProblemNLP {
public:
  dnoptProblemB();
  dnoptProblemB(const char*name);
  ~dnoptProblemB();

  void setWorkspace(int mLCon, int mNCon, int n,
		    int nnJac, int nnObj, int iObj);

  int solve(int starttype, int n, int mLCon, int nMCon,
	    int nnJac, int nnObj, int iObj, double ObjAdd,
	    dnCon funcon, dnObj funobj,
	    int *state, double *A, int ldA,
	    double *bl, double *bu,
	    double *J, int ldJ, double *H, int ldH,
	    double &objective, int &nInf, double &sInf,
	    double *x, double *y);

  int solveh(int starttype, int n, int mLCon, int nMCon,
	     int nnJac, int nnObj, int iObj, double ObjAdd,
	     dnCon funcon, dnObj funobj, dnHes funhes,
	     int *state, double *A, int ldA,
	     double *bl, double *bu,
	     double *J, int ldJ, double *H, int ldH,
	     double &objective, int &nInf, double &sInf,
	     double *x, double *y);
};

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

class dqoptProblem : public dnoptProblem {
private:
  idqLog dqLog;
  void init2zero();

public:
  dqoptProblem();
  dqoptProblem(const char* name);
  ~dqoptProblem();

  void initialize  (const char *prtfile, int summOn);
  int  setSpecsFile(const char *specname);

  void setLog      (idqLog dqLog);
  void setSTOP     (idnSTOP dnSTOP);

  void setWorkspace();
  int  solve       (int starttype, int mCon, int n, int nnH, int ncObj,
		    int iObj, double ObjAdd,
		    double *A, int ldA, double *bl, double *bu,
		    double *cObj, double *H, int ldH, dqHx qpHx,
		    int *eType, int *state, double *x, double *y,
		    double &objective, int &nInf, double &sInf);
};


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#endif /* DNOPTPROBLEM_H */
