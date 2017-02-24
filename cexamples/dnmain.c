/*------------------------------------------------------------------
  File: dnmain.f
  ------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dnopt_cwrap.h"


void hexcon(int *modeC, int *mnCon, int *nnJac, double x[],
	     double fCon[], double JCon[*nnJac][*mnCon], int *nState,
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

  JCon[0][0]  =   2.*x[0];
  JCon[5][0]  =   2.*x[5];

  JCon[0][1]  = - 2.*(x[1] - x[0]);
  JCon[1][1]  =   2.*(x[1] - x[0]);
  JCon[5][1]  = - 2.*(x[6] - x[5]);
  JCon[6][1]  =   2.*(x[6] - x[5]);

  JCon[0][2]  = - 2.*(x[2] - x[0]);
  JCon[2][2]  =   2.*(x[2] - x[0]);
  JCon[5][2]  =   2.*x[5];

  JCon[0][3]  =   2.*(x[0] - x[3]);
  JCon[3][3]  = - 2.*(x[0] - x[3]);
  JCon[5][3]  =   2.*(x[5] - x[7]);
  JCon[7][3]  = - 2.*(x[5] - x[7]);

  JCon[0][4]  =   2.*(x[0] - x[4]);
  JCon[4][4]  = - 2.*(x[0] - x[4]);
  JCon[5][4]  =   2.*(x[5] - x[8]);
  JCon[8][4]  = - 2.*(x[5] - x[8]);

  JCon[1][5]  =   2.*x[1];
  JCon[6][5]  =   2.*x[6];

  JCon[1][6]  = - 2.*(x[2] - x[1]);
  JCon[2][6]  =   2.*(x[2] - x[1]);
  JCon[6][6]  =   2.*x[6];

  JCon[1][7]  = - 2.*(x[3] - x[1]);
  JCon[3][7]  =   2.*(x[3] - x[1]);
  JCon[6][7]  = - 2.*(x[7] - x[6]);
  JCon[7][7]  =   2.*(x[7] - x[6]);

  JCon[1][8]  =   2.*(x[1] - x[4]);
  JCon[4][8]  = - 2.*(x[1] - x[4]);
  JCon[6][8]  =   2.*(x[6] - x[8]);
  JCon[8][8]  = - 2.*(x[6] - x[8]);

  JCon[2][9] = - 2.*(x[3] - x[2]);
  JCon[3][9] =   2.*(x[3] - x[2]);
  JCon[7][9] =   2.*x[7];

  JCon[2][10] = - 2.*(x[4] - x[2]);
  JCon[4][10] =   2.*(x[4] - x[2]);
  JCon[8][10] =   2.*x[8];

  JCon[3][11] =   2.*x[3];
  JCon[7][11] =   2.*x[7];

  JCon[3][12] =   2.*(x[3] - x[4]);
  JCon[4][12] = - 2.*(x[3] - x[4]);
  JCon[7][12] = - 2.*(x[8] - x[7]);
  JCon[8][12] =   2.*(x[8] - x[7]);

  JCon[4][13] =   2.*x[4];
  JCon[8][13] =   2.*x[8];
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
  dnProblem toy;

  int    i, info, nS, nInf;
  double objective, sInf;

  int    cold   =  0;
  int    n      =  9;
  int    mLCon  =  4;
  int    mNCon  =  14;
  int    nb     =  27; // n + mLCon + mNCon
  int    nnObj  =  n;
  int    nnJac  =  n;
  int    iObj   = -1;
  double ObjAdd =  0;

  double inf = 1e20;

  int    state[27];
  double bl[27], bu[27], x[27], y[27], A[4*9], J[14*9], H[9*9];

  // Print file hexagon.out
  // Summary output turned off
  dnBegin(&toy, "hexagon", "hexagon.out", 1);

  // Set some options
  dnSpecsfile(&toy, "dnmain.spc");
  dnSetIntParameter(&toy, "Derivative level", 3);
  dnSetParameter   (&toy, "Sticky parameters  yes");


  // Set up the bounds of the problem.
  for (i = 0; i < nb; i++) {
    bl[i] = -inf;
    bu[i] =  inf;
  }
  for (i = n; i < n+mNCon; i++) {
    bl[i] = -inf;
    bu[i] =  1.0;
  }
  for (i = n+mNCon; i < nb; i++) {
    bl[i] =  0.0;
    bu[i] =  inf;
  }

  bl[0] =  0.0;
  bl[2] = -1.0;
  bl[4] =  0.0;
  bl[5] =  0.0;
  bl[6] =  0.0;

  bu[2] =  1.0;
  bu[7] =  0.0;
  bu[8] =  0.0;


  // Initialize states, x and multipliers
  for (i = 0; i < nb; i++) {
    state[i] = 0;
    x[i]     = 0;
    y[i]     = 0;
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


  // Set up the Jacobian matrix
  for (i = 0; i < mLCon*n; i++) {
    A[i]= 0.0;
  }

  A[0*mLCon+0] = -1.0;
  A[1*mLCon+0] =  1.0;
  A[1*mLCon+1] = -1.0;
  A[2*mLCon+1] =  1.0;
  A[2*mLCon+2] =  1.0;
  A[3*mLCon+2] = -1.0;
  A[3*mLCon+3] =  1.0;
  A[4*mLCon+3] = -1.0;

  for (i = 0; i < mNCon*n; i++) {
    J[i] = 0.0;
  }

  info = dnopt(&toy, cold, n, mLCon, mNCon,
	       nnJac, nnObj, iObj, ObjAdd,
	       (dnCon) hexcon, (dnObj) hexobj,
	       state, A, mLCon, bl, bu, J, mNCon, H, n,
	       &objective, &nS, &nInf, &sInf, x, y);

  printf("Final objective value: %f\n", objective);
  printf("Final x:\n");
  for (i = 0; i < n; i++) {
    printf("x[%d]:  %f\n", i, x[i]);
  }

  deleteDNOPT(&toy);
}
