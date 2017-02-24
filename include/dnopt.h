#ifndef DNOPT_H
#define DNOPT_H

/*
 * File:  dnopt.h
 *   Header file for the DNOPT functions.
 *
 *   10 Jul 2014: First version.
 */

#ifdef __cplusplus
extern "C" {
#endif

  typedef void (*idnLog)
  (int *iAbort, int Tags[], int KTcond[], int *minimz,
   int *n, int *nb, int *mnCon0, int *nZ,
   int *itn, int *nMajor, int *nMinor,
   double *condHz, double *ObjAdd, double *fMrt, double PenParm[], double *step,
   double *prInf, double *duInf, double *maxVi, double *maxViRel, int state[],
   double scales[], double bl[], double bu[], double fCon[], double fmul[],
   double *JQP, int *ldJQP, double x[],
   char cu[], int *lencu, int iu[], int *leniu, double ru[], int *lenru,
   char cw[], int *lencw, int iw[], int *leniw, double rw[], int *lenrw);

  typedef void (*idnLogQP)
  (int *qpProbType, int *phase, char *runInfo, char *qpProbTag,
   int *Elastic, int *gotR, int *Named, int *Optimal, int *nNames, char Names[],
   int *stateDel, int *jAdd, int *jDel, int *printDel, int *minimize,
   int *n, int *nactiv, int *nb, int *mLCon, int *nfree, int *nnH,
   int *nonOpt, int *nZ, int *nZr,
   int *nlineSP, int *nlinesS, int *itn, int *itQP,
   int *jOj, double *scaleObj, double *ObjAdd,
   double *objc, int *nInf, double *sInf, int *nInfE, double *sInfE, double *wtInf,
   double *condHz, double *yDel, double *normgZr, double *step,
   int state[], int kx[], double *R, int *ldR, double *T, int ldT, double x[],
   double wrk1[], int iw[], int *leniw);

  typedef void (*idqLog)
  (int *qpProbType, int *phase, char *runInfo, char *qpProbTag,
   int *Elastic, int *gotR, int *Named, int *Optimal, int *nNames, char Names[],
   int *stateDel, int *jAdd, int *jDel, int *printDel, int *minimize,
   int *n, int *nactiv, int *nb, int *mLCon, int *nfree, int *nnH,
   int *nonOpt, int *nZ, int *nZr,
   int *nlineSP, int *nlinesS, int *itn, int *itQP,
   int *jOj, double *scaleObj, double *ObjAdd,
   double *objc, int *nInf, double *sInf, int *nInfE, double *sInfE, double *wtInf,
   double *condHz, double *yDel, double *normgZr, double *step,
   int state[], int kx[], double *R, int *ldR, double *T, int ldT, double x[],
   double wrk1[], int iw[], int *leniw);

  typedef void (*idnSTOP)
  (int *iAbort, int KTcond[], int *minimz,
   int *mCon0, int *n, int *nb, int *mnCon0, int *mnCon, int *nnObj0, int *nnObj, int *nZ,
   int *itn, int *nMajor, int *nMinor,
   double *condHz, double *ObjAdd,
   double *fMrt, double PenParm[], double *step,
   double *prInf, double *duInf, double *maxVi, double *maxViRel, int state[],
   double scale[], double bl[], double bu[],
   double *fObj, double gObj[], double fCon[], double fmul[],
   double *JQP, int *ldJQP, double gQ[], double x[], double yCon[],
   char cu[], int *lencu, int iu[], int *leniu, double ru[], int *lenru,
   char cw[], int *lencw, int iw[], int *leniw, double rw[], int *lenrw);

  typedef void (*dnCon)
  (int *modeC, int *mnCon, int *nnJac, double x[],
   double fCon[], double *JCon, int *nState,
   char     cu[], int   *lencu,
   int      iu[], int   *leniu,
   double   ru[], int   *lenru);

  typedef void (*dnObj)
  (int *modeF, int *nnObj, double x[],
   double *fObj, double gObj[], int *nState,
   char    cu[], int   *lencu,
   int     iu[], int   *leniu,
   double  ru[], int   *lenru);

  typedef void (*dnHes)
  (int *modeH, int *nnH, int *mnCon,
   double x[], double *fMul,
   double *H, int *ldH, int *status,
   char    cu[], int   *lencu,
   int     iu[], int   *leniu,
   double  ru[], int   *lenru);

  typedef void (*dqHx)
  (int ncolH, double *H, int ldH,
   double x[], double Hx[], int qpState,
   char    cu[], int   *lencu,
   int     iu[], int   *leniu,
   double  ru[], int   *lenru);


  /* DNOPT */
  void f_dnbegin(const char *name, int len, int summOn,
		 int iw[], int leniw, double rw[], int lenrw);
  void f_dnsetprint(const char *name, int len,
		    int iw[], int leniw, double rw[], int lenrw);

  void f_dnspec(const char *specfile, int len, int *inform,
		int iw[], int leniw, double rw[], int lenrw);

  void f_dngetc(const char *buffer, int lenb, char *ivalue,
		int lenc, int *errors,
		int iw[], int leniw, double rw[], int lenrw);
  void f_dngeti(const char *buffer, int len, int  *ivalue, int *errors,
		int iw[], int leniw, double rw[], int lenrw);
  void f_dngetr(const char *buffer, int len, double *rvalue, int *errors,
		int iw[], int leniw, double rw[], int lenrw);

  void f_dnset(const char *buffer, int len, int *errors,
	       int iw[], int leniw, double rw[], int lenrw);
  void f_dnseti(const char *buffer, int len, int iopt, int *errors,
		int iw[], int leniw, double rw[], int lenrw);
  void f_dnsetr(const char *buffer, int len, double rvalue, int *errors,
		int iw[], int leniw, double rw[], int lenrw);
  void f_dnend(int iw[], int leniw, double rw[], int lenrw);

  void f_dnmem(int *iExit, int mLCon, int mNCon, int n,
	       int nnJac, int nnObj, int iObj, int *miniw, int *minrw,
	       int iw[], int leniw, double rw[], int lenrw);

  void f_dnopt(int start, int n, int nb,
	       int mLCon, int mNcon, int nnJac, int nnObj,
	       const char *prob,
	       int iobj, double objadd,
	       dnCon funcon, dnObj funobj,
	       int *state, double *A, int ldA, double *bl, double *bu,
	       double *fObj, double gObj[], double fCon[],
	       double *Jcon, int ldJ, double *H, int ldH,
	       double *objNP, int *nInf, double *sInf, double *x, double *y,
	       int *inform, int *miniw, int *minrw,
	       int iu[], int leniu, double ru[], int lenru,
	       int iw[], int leniw, double rw[], int lenrw);

  void f_dnopth(int start, int n, int nb,
		int mLCon, int mNcon, int nnJac, int nnObj,
		const char *prob,
		int iobj, double objadd,
		dnCon funcon, dnObj funobj, dnHes funhes,
		int *state, double *A, int ldA, double *bl, double *bu,
		double *fObj, double gObj[], double fCon[],
		double *Jcon, int ldJ, double *H, int ldH,
		double *objNP, int *nInf, double *sInf, double *x, double *y,
		int *inform, int *miniw, int *minrw,
		int iu[], int leniu, double ru[], int lenru,
		int iw[], int leniw, double rw[], int lenrw);

  void f_dnker(int start, int n, int nb,
	       int mLCon, int mNcon, int nnJac, int nnObj,
	       const char *prob,
	       int iobj, double objadd,
	       dnCon funcon, dnObj funobj, dnHes funhes,
	       idnLog dnLog, idnLogQP dnLogQP, idnSTOP dnSTOP,
	       int *state, double *A, int ldA,
	       double *bl, double *bu,
	       double *Jcon, int ldJ, double *H, int ldH,
	       double *objNP, int *nInf, double *sInf, double *x, double *y,
	       int *inform, int *miniw, int *minrw,
	       int iu[], int leniu, double ru[], int lenru,
	       int iw[], int leniw, double rw[], int lenrw);

  /* DQOPT */
  void f_dqbegin(const char *name, int len, int summOn,
		 int iw[], int leniw, double rw[], int lenrw);
  void f_dqsetprint(const char *name, int len,
		    int iw[], int leniw, double rw[], int lenrw);

  void f_dqspec(const char *specfile, int len, int *inform,
		int iw[], int leniw, double rw[], int lenrw);

  void f_dqgetc(const char *buffer, int lenb, char *ivalue,
		int lenc, int *errors,
		int iw[], int leniw, double rw[], int lenrw);
  void f_dqgeti(const char *buffer, int len, int  *ivalue, int *errors,
		int iw[], int leniw, double rw[], int lenrw);
  void f_dqgetr(const char *buffer, int len, double *rvalue, int *errors,
		int iw[], int leniw, double rw[], int lenrw);

  void f_dqset(const char *buffer, int len, int *errors,
	       int iw[], int leniw, double rw[], int lenrw);
  void f_dqseti(const char *buffer, int len, int iopt, int *errors,
		int iw[], int leniw, double rw[], int lenrw);
  void f_dqsetr(const char *buffer, int len, double rvalue, int *errors,
		int iw[], int leniw, double rw[], int lenrw);
  void f_dqend(int iw[], int leniw, double rw[], int lenrw);

  void f_dqopt(int start, int n, int nb, int mCon, int nnH, int ncObj,
	       int iobj, double objadd, const char *prob,
	       double *A, int ldA, double *bl, double *bu,
	       double *cObj, double *H, int ldH, dqHx userHx,
	       int *eType, int *state, double *x, double *y,
	       int *inform, int *miniw, int *minrw,
	       double *objQP, int *nInf, double *sInf,
	       int iu[], int leniu, double ru[], int lenru,
	       int iw[], int leniw, double rw[], int lenrw);

  void f_dqker(int start, int n, int nb, int mCon, int nnH, int ncObj,
	       int iobj, double objadd, const char *prob,
	       double *A, int ldA, double *bl, double *bu,
	       double *cObj, double *H, int ldH,
	       dqHx userHx, idqLog dqLog,
	       int *eType, int *state, double *x, double *y,
	       int *inform, int *miniw, int *minrw,
	       double *objQP, int *nInf, double *sInf,
	       int iu[], int leniu, double ru[], int lenru,
	       int iw[], int leniw, double rw[], int lenrw);


#ifdef __cplusplus
}
#endif

#endif /* DNOPT_H */
