/*------------------------------------------------------------------
  File dntoy.cpp: implements the sample nonlinear problem defined
  in the dnOpt User's Guide.
  ------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <cmath>

#include "dnoptProblem.hpp"

void toycon(int *modeC, int *mnCon, int *nnJac, double x[],
	     double fCon[], double *JCon, int *nState,
	     char     cu[], int   lencu,
	     int      iu[], int   leniu,
	     double   ru[], int   lenru) {

  if (*modeC == 0  ||  *modeC == 2) {
    fCon[0]  = x[0]*x[0] +  x[1]*x[1];
    fCon[1]  =            pow(x[1],4);
  }

  if (*modeC >= 1) {
    JCon[0+0*(*mnCon)] = 2.0*x[0];
    JCon[1+0*(*mnCon)] = 2.0*x[1];

    JCon[0+1*(*mnCon)] = 0.0;
    JCon[1+1*(*mnCon)] = 4.0*pow(x[1],3);
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

int main(int argc, char **argv) {

  int    cold   =  0;
  int    n      =  4;
  int    mLCon  =  2;
  int    mNCon  =  2;
  int    nb     = n + mLCon + mNCon;

  int    nnJac  =  2;
  int    nnObj  =  3;
  int    iObj   =  1;

  double ObjAdd = 0.0;
  double infBnd = 1e20;

  int state[nb];
  double x[nb], y[nb], bl[nb], bu[nb], A[mLCon*n], J[mNCon*n], H[n*n];

  int    j, info, nInf;
  double  objective, sInf;

  dnoptProblemB ToyProb("ToyB");

  ToyProb.initialize("toyB.out", 1);


  // Initialize states, x and multipliers
  for (j = 0; j < nb; j++) {
    state[j] = 0;
    x[j]     = 0;
    y[j]     = 0;
  }

  // Initialize Jacobian matrices
  for(j = 0; j < mNCon*n; j++) {
    J[j] = 0.0;
  }

  J[2*mNCon+0]  = 1.0;
  J[3*mNCon+1]  = 1.0;

  for(j = 0; j < mLCon*n; j++) {
    A[j] = 0.0;
  }

  A[0*mLCon+0   ] = 2.0;  // A(0,0)
  A[1*mLCon+0   ] = 4.0;  // A(0,1)
  A[2*mLCon+iObj] = 3.0;  // A(iObj,2)
  A[3*mLCon+iObj] = 5.0;  // A(iObj,3)

  // Set bound sfor the problem
  for(j = 0; j < n; j++) {
    bl[j] = -infBnd;
    bu[j] =  infBnd;
  }

  bl[2] =  0.0;
  bl[3] =  0.0;

  bl[n] =  2.0;
  bu[n] =  2.0;

  bl[n+1] =  4.0;
  bu[n+1] =  4.0;

  bl[n+2] =  0.0;
  bu[n+2] =  infBnd;

  bl[n+3] = -infBnd;
  bu[n+3] =  infBnd;

  ToyProb.setSpecsFile("dntoy.spc");
  ToyProb.setIntParameter("Verify level", 3);

  info = ToyProb.solve(cold, n, mLCon, mNCon, nnJac, nnObj,
		       iObj, ObjAdd,
		       (dnCon) toycon, (dnObj) toyobj,
		       state, A, mLCon, bl, bu, J, mNCon, H, n,
		       objective, nInf, sInf, x, y);
}
