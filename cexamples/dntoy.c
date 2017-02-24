/*------------------------------------------------------------------
  File dntoy.f: implements the sample nonlinear problem defined
  in the dnOpt User's Guide.
  ------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dnopt_cwrap.h"

void toycon(int *modeC, int *mnCon, int *nnJac, double x[],
	     double fCon[], double JCon[*nnJac][*mnCon], int *nState,
	     char     cu[], int   lencu,
	     int      iu[], int   leniu,
	     double   ru[], int   lenru) {

  if (*modeC == 0  ||  *modeC == 2) {
    fCon[0]  = pow(x[0],2) +  pow(x[1],2);
    fCon[1]  =                pow(x[1],4);
  }

  if (*modeC >= 1) {
    JCon[0][0] = 2.0*x[0];
    JCon[1][0] = 2.0*x[1];

    JCon[0][1] = 0.0;
    JCon[1][1] = 4.0*pow(x[1],3);
  }
}

void toyobj(int *modeF, int nnObj, double x[],
	     double *fObj, double gObj[], int nState,
	     char    cu[], int   lencu,
	     int     iu[], int   leniu,
	     double  ru[], int   lenru) {
  double sum;

  sum = x[0] + x[1] + x[2];

  if (*modeF == 0  || *modeF == 2) {
    *fObj    = sum*sum;
  }

  if (*modeF == 1  ||  *modeF ==  2) {
    sum     = 2.0*sum;
    gObj[0] = sum;
    gObj[1] = sum;
    gObj[2] = sum;
  }

}

int main(int argc, char *argv[]) {
  int    cold   =  0;
  int    n      =  4;
  int    mLCon  =  2;
  int    mNCon  =  2;
  int    nb     =  8; // n + mLCon + mNCon

  int    nnJac  =  2;
  int    nnObj  =  3;
  int    iObj   =  1;

  double ObjAdd = 0.0;
  double infBnd = 1e20;

  int    state[8];
  double bl[8], bu[8], x[8], y[8], A[2*4], J[2*4], H[4*4];

  int i, info, nInf, nS;
  double objective, sInf;

  dnProblem toy;

  dnBegin(&toy, "dntoy", "dntoy.out", 1);
  dnSpecsfile(&toy, "dntoy.spc");
  dnSetIntParameter(&toy, "Verify level", 3);

  for (i = 0; i < mNCon*n; i++) {
    J[i] = 0.0;
  }

  J[2*mNCon+0]  = 1.0;
  J[3*mNCon+1]  = 1.0;

  for (i = 0; i < mLCon*n; i++) {
    A[i] = 0.0;
  }

  A[0*mLCon+0   ] = 2.0;  // A(0,0)
  A[1*mLCon+0   ] = 4.0;  // A(0,1)
  A[2*mLCon+iObj] = 3.0;  // A(iObj,2)
  A[3*mLCon+iObj] = 5.0;  // A(iObj,3)

  bl[n] =  2.0;
  bu[n] =  2.0;

  bl[n+1] =  4.0;
  bu[n+1] =  4.0;

  bl[n+2] =  0.0;
  bu[n+2] =  infBnd;

  bl[n+3] = -infBnd;
  bu[n+3] =  infBnd;


  for (i = 0; i < n; i++) {
    bl[i] = -infBnd;
    bu[i] =  infBnd;
  }

  bl[2] =  0.0;
  bl[3] =  0.0;

  for (i = 0; i < nb; i++) {
    state[i] = 0;
    x[i]     = 0.0;
    y[i]     = 0.0;
  }

  info = dnopt(&toy, cold, n, mLCon, mNCon,
	       nnJac, nnObj, iObj, ObjAdd,
	       (dnCon) toycon, (dnObj) toyobj,
	       state, A, mLCon, bl, bu, J, mNCon, H, n,
	       &objective, &nS, &nInf, &sInf, x, y);

  deleteDNOPT(&toy);

}
