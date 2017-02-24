#ifndef DQOPT_C_H
#define DQOPT_C_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "dnopt.h"

/* File: dqopt_cwrap.h
 *   C interface for DQOPT.
 *
 */
typedef struct {
  char   *name;

  int     memCalled;
  int     initCalled;

  idqLog  dqLog;

  int     lenrw, leniw;
  int     *iw;
  double  *rw;

  int     lenru, leniu;
  int    *iu;
  double  *ru;

} dqProblem;

void dqBegin        (dqProblem *prob, char* name, char* prtfile, int summaryOn);
void init2zeroQ     (dqProblem *prob);

void allocIQ        (dqProblem *prob, int len);
void allocRQ        (dqProblem *prob, int len);
void reallocIQ      (dqProblem *prob, int len);
void reallocRQ      (dqProblem *prob, int len);

void dqPrintfile    (dqProblem *prob, char* prtname);
int  dqSpecsfile    (dqProblem *prob, char* spcname);

int dqSetParameter    (dqProblem *prob, char stropt[]);
int dqGetParameter    (dqProblem *prob, char stropt[], char strout[]);
int dqSetIntParameter (dqProblem *prob, char stropt[], int opt);
int dqGetIntParameter (dqProblem *prob, char stropt[], int opt);
int dqSetRealParameter(dqProblem *prob, char stropt[], double opt);
int dqGetRealParameter(dqProblem *prob, char stropt[], double opt);

void dqSetUserI       (dqProblem *prob, int *iu, int leniu);
void dqSetUserR       (dqProblem *prob, double *ru, int lenru);
void dqSetUserspace   (dqProblem *prob, int *iu, int leniu,
		       double *ru, int lenru);

void dqSetLog         (dqProblem *prob, idqLog dqLog);

void dqSetWorkspace   (dqProblem *prob);

int dqopt(dqProblem *prob, int start, int n,
	  int mCon, int nnH, int ncObj,
	  int iobj, double objadd,
	  double *A, int ldA, double *bl, double *bu,
	  double *cObj, double *H, int ldH, dqHx userHx,
	  int *eType, int *state, double *x, double *y,
	  double *objQP, int *nInf, double *sInf);

void deleteDQOPT(dqProblem *prob);

#endif /* DQOPT_C_H */
