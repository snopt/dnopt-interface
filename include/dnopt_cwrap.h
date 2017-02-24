#ifndef DNOPT_C_H
#define DNOPT_C_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "dnopt.h"

/* File: dnopt_cwrap.h
 *   C interface for DNOPT.
 *
 * 10 Jul 2014: First version (based on cwrap).
 */
typedef struct {
  char   *name;

  int     memCalled;
  int     initCalled;

  idnSTOP  dnSTOP;
  idnLog   dnLog;
  idnLogQP dnLogQP;

  int     lenrw, leniw;
  int     *iw;
  double  *rw;

  int     lenru, leniu;
  int    *iu;
  double  *ru;

} dnProblem;

void dnBegin        (dnProblem *prob, char* name, char* prtfile, int summaryOn);
void dninit2zero    (dnProblem *prob);

void dnallocI       (dnProblem *prob, int len);
void dnallocR       (dnProblem *prob, int len);
void dnreallocI     (dnProblem *prob, int len);
void dnreallocR     (dnProblem *prob, int len);

void dnPrintfile    (dnProblem *prob, char* prtname);
int  dnSpecsfile    (dnProblem *prob, char* spcname);

int dnSetParameter    (dnProblem *prob, char stropt[]);
int dnGetParameter    (dnProblem *prob, char stropt[], char strout[]);
int dnSetIntParameter (dnProblem *prob, char stropt[], int opt);
int dnGetIntParameter (dnProblem *prob, char stropt[], int opt);
int dnSetRealParameter(dnProblem *prob, char stropt[], double opt);
int dnGetRealParameter(dnProblem *prob, char stropt[], double opt);

void dnSetUserI       (dnProblem *prob, int *iu, int leniu);
void dnSetUserR       (dnProblem *prob, double *ru, int lenru);
void dnSetUserspace   (dnProblem *prob, int *iu, int leniu,
		      double *ru, int lenru);

void dnSetLog         (dnProblem *prob, idnLog dnLog, idnLogQP dnLogQP);
void dnSetSTOP        (dnProblem *prob, idnSTOP dnSTOP);

void dnSetWorkspace   (dnProblem *prob, int mLCon, int mNCon, int n,
		       int nnJac, int nnObj, int iObj);

int dnopt(dnProblem *prob, int start, int n, int mLCon, int mNCon,
	  int nnJac, int nnObj, int iObj, double ObjAdd,
	  dnCon funcon, dnObj funobj,
	  int *state, double *A, int ldA, double *bl, double *bu,
	  double *J, int ldJ, double *H, int ldH,
	  double *objective, int *nS, int* nInf, double* sInf,
	  double *x, double *y);

int dnopth(dnProblem *prob, int start, int n, int mLCon, int mNCon,
	   int nnJac, int nnObj, int iObj, double ObjAdd,
	   dnCon funcon, dnObj funobj, dnHes funhes,
	   int *state, double *A, int ldA, double *bl, double *bu,
	   double *J, int ldJ, double *H, int ldH,
	   double *objective, int *nS, int* nInf, double* sInf,
	   double *x, double *y);

void deleteDNOPT    (dnProblem *prob);

#endif /* DNOPT_C_H */
