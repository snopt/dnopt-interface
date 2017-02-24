/*------------------------------------------------------------------
  File dnmain.cpp
  Illustrates the use of subroutine dnopt.
  ------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <cmath>

#include "dnoptProblem.hpp"

void hexcon(int *modeC, int *mnCon, int *nnJac, double x[],
	     double fCon[], double *JCon, int *nState,
	     char     cu[], int   *lencu,
	     int      iu[], int   *leniu,
	     double   ru[], int   *lenru) {

  fCon[0]    =   pow(x[0],2)  +  pow(x[5],2);
  fCon[1]    =   pow(x[1] - x[0],2)  +  pow(x[6] - x[5],2);
  fCon[2]    =   pow(x[2] - x[0],2)  +  pow(x[5],2);
  fCon[3]    =   pow(x[0] - x[3],2)  +  pow(x[5] - x[7],2);
  fCon[4]    =   pow(x[0] - x[4],2)  +  pow(x[5] - x[8],2);
  fCon[5]    =   pow(x[1],2)  +  pow(x[6],2);
  fCon[6]    =   pow(x[2] - x[1],2)  +  pow(x[6],2);
  fCon[7]    =   pow(x[3] - x[1],2)  +  pow(x[7] - x[6],2);
  fCon[8]    =   pow(x[1] - x[4],2)  +  pow(x[6] - x[8],2);
  fCon[9]    =   pow(x[3] - x[2],2)  +  pow(x[7],2);
  fCon[10]   =   pow(x[4] - x[2],2)  +  pow(x[8],2);
  fCon[11]   =   pow(x[3],2)  +  pow(x[7],2);
  fCon[12]   =   pow(x[3] - x[4],2)  +  pow(x[8] - x[7],2);
  fCon[13]   =   pow(x[4],2)  +  pow(x[8],2);

  JCon[0*(*mnCon)+0]  =   2.*x[0];
  JCon[5*(*mnCon)+0]  =   2.*x[5];

  JCon[0*(*mnCon)+1]  = - 2.*(x[1] - x[0]);
  JCon[1*(*mnCon)+1]  =   2.*(x[1] - x[0]);
  JCon[5*(*mnCon)+1]  = - 2.*(x[6] - x[5]);
  JCon[6*(*mnCon)+1]  =   2.*(x[6] - x[5]);

  JCon[0*(*mnCon)+2]  = - 2.*(x[2] - x[0]);
  JCon[2*(*mnCon)+2]  =   2.*(x[2] - x[0]);
  JCon[5*(*mnCon)+2]  =   2.*x[5];

  JCon[0*(*mnCon)+3]  =   2.*(x[0] - x[3]);
  JCon[3*(*mnCon)+3]  = - 2.*(x[0] - x[3]);
  JCon[5*(*mnCon)+3]  =   2.*(x[5] - x[7]);
  JCon[7*(*mnCon)+3]  = - 2.*(x[5] - x[7]);

  JCon[0*(*mnCon)+4]  =   2.*(x[0] - x[4]);
  JCon[4*(*mnCon)+4]  = - 2.*(x[0] - x[4]);
  JCon[5*(*mnCon)+4]  =   2.*(x[5] - x[8]);
  JCon[8*(*mnCon)+4]  = - 2.*(x[5] - x[8]);

  JCon[1*(*mnCon)+5]  =   2.*x[1];
  JCon[6*(*mnCon)+5]  =   2.*x[6];

  JCon[1*(*mnCon)+6]  = - 2.*(x[2] - x[1]);
  JCon[2*(*mnCon)+6]  =   2.*(x[2] - x[1]);
  JCon[6*(*mnCon)+6]  =   2.*x[6];

  JCon[1*(*mnCon)+7]  = - 2.*(x[3] - x[1]);
  JCon[3*(*mnCon)+7]  =   2.*(x[3] - x[1]);
  JCon[6*(*mnCon)+7]  = - 2.*(x[7] - x[6]);
  JCon[7*(*mnCon)+7]  =   2.*(x[7] - x[6]);

  JCon[1*(*mnCon)+8]  =   2.*(x[1] - x[4]);
  JCon[4*(*mnCon)+8]  = - 2.*(x[1] - x[4]);
  JCon[6*(*mnCon)+8]  =   2.*(x[6] - x[8]);
  JCon[8*(*mnCon)+8]  = - 2.*(x[6] - x[8]);

  JCon[2*(*mnCon)+9] = - 2.*(x[3] - x[2]);
  JCon[3*(*mnCon)+9] =   2.*(x[3] - x[2]);
  JCon[7*(*mnCon)+9] =   2.*x[7];

  JCon[2*(*mnCon)+10] = - 2.*(x[4] - x[2]);
  JCon[4*(*mnCon)+10] =   2.*(x[4] - x[2]);
  JCon[8*(*mnCon)+10] =   2.*x[8];

  JCon[3*(*mnCon)+11] =   2.*x[3];
  JCon[7*(*mnCon)+11] =   2.*x[7];

  JCon[3*(*mnCon)+12] =   2.*(x[3] - x[4]);
  JCon[4*(*mnCon)+12] = - 2.*(x[3] - x[4]);
  JCon[7*(*mnCon)+12] = - 2.*(x[8] - x[7]);
  JCon[8*(*mnCon)+12] =   2.*(x[8] - x[7]);

  JCon[4*(*mnCon)+13] =   2.*x[4];
  JCon[8*(*mnCon)+13] =   2.*x[8];
}

void hexobj(int *modeF, int *nnObj, double x[],
	     double *fObj, double gObj[], int *nState,
	     char    cu[], int   *lencu,
	     int     iu[], int   *leniu,
	     double  ru[], int   *lenru) {

  *fObj =  - x[1]*x[5] + x[0]*x[6] - x[2]*x[6] - x[4]*x[7]
    + x[3]*x[8] + x[2]*x[7];

  gObj[0] =   x[6];
  gObj[1] = - x[5];
  gObj[2] = - x[6] + x[7];
  gObj[3] =   x[8];
  gObj[4] = - x[7];
  gObj[5] = - x[1];
  gObj[6] = - x[2] + x[0];
  gObj[7] = - x[4] + x[2];
  gObj[8] =   x[3];
}


int main(int argc, char *argv[]) {
  int    info, j, nInf;
  double objective, sInf;

  int    cold   =  0;
  int    n      =  9;
  int    mLCon  =  4;
  int    mNCon  =  14;
  int    nnObj  =  n;
  int    nnJac  =  n;

  int    iObj   = -1;
  double ObjAdd =  0;

  int     nb = n + mLCon + mNCon;
  double inf = 1e20;

  int state[nb];
  double x[nb], y[nb], bl[nb], bu[nb], A[mLCon*n], J[mNCon*n], H[n*n];

  dnoptProblemB myProb("dnmain");

  myProb.initialize("dnmain.out", 1);


  // Initialize states, x and multipliers
  for (j = 0; j < nb; j++) {
    state[j] = 0;
    x[j]     = 0;
    y[j]     = 0;
  }

  x[0] = 0.1;
  x[1] = 0.125;
  x[2] =  .666666;
  x[3] =  .142857;
  x[4] =  .111111;
  x[5] =  .2;
  x[6] =  .25;
  x[7] = -.2;
  x[8] = -.25;


  // Set up the bounds of the problem.
  for (j = 0; j < nb; j++) {
    bl[j] = -inf;
    bu[j] =  inf;
  }
  for (j = n; j < n+mNCon; j++) {
    bl[j] = -inf;
    bu[j] =  1.0;
  }
  for (j = n+mNCon; j < nb; j++) {
    bl[j] =  0.0;
    bu[j] =  inf;
  }

  bl[0] =  0.0;
  bl[2] = -1.0;
  bl[4] =  0.0;
  bl[5] =  0.0;
  bl[6] =  0.0;

  bu[2] =  1.0;
  bu[7] =  0.0;
  bu[8] =  0.0;


  // Set up the Jacobian matrix
  for (j = 0; j < mLCon*n; j++) {
    A[j] = 0.0;
  }

  A[0*mLCon+0] = -1.0;
  A[1*mLCon+0] =  1.0;
  A[1*mLCon+1] = -1.0;
  A[2*mLCon+1] =  1.0;
  A[2*mLCon+2] =  1.0;
  A[3*mLCon+2] = -1.0;
  A[3*mLCon+3] =  1.0;
  A[4*mLCon+3] = -1.0;

  for (j = 0; j < mNCon*n; j++) {
    J[j] = 0.0;
  }

  for (j = 0; j < n*n; j++) {
    H[j] = 0.0;
  }

  myProb.setSpecsFile("dnmain.spc");
  myProb.setIntParameter("Derivative level", 3);
  myProb.setParameter("Sticky parameters  yes");

  info = myProb.solve(cold, n, mLCon, mNCon, nnJac, nnObj,
		      iObj, ObjAdd,
		      (dnCon) hexcon, (dnObj) hexobj,
		      state, A, mLCon, bl, bu, J, mNCon, H, n,
		      objective, nInf, sInf, x, y);

  printf("Final objective value: %f\n", objective);
  printf("Final x:\n");
  for (j = 0; j < n; j++) {
    printf("x[%d]:  %f\n", j, x[j]);
  }

}
