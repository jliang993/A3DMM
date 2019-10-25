/* -*- mode: c++; c-basic-offset: 3 -*-
 *
 * TV-minimization 
 * this can be linked with C code (tested only on linux) provided the
 * flag -lstdc++ is given to the linker
 * TV4 = with nearest neighbours interaction
 * TV8 = with also next nearest
 * TV16 = with 16 neighbours
 * 
 * Based on the code by 
 *	A. Chambolle and J. Darbon: On total variation
 * 	minimization and surface evolution using parametric maximum flows,
 *	preprint (2008).
 * Their code implements Dorit Hochbaum's algorithm:
 *    D. S. Hochbaum: An efficient algorithm for image segmentation,
 *    Markov random fields and related problems. J. ACM, 48(4):686--701,
 *    2001.     
 * 
 *
 * 
 * 
 * 
 */

/**
 * @authors
 *  Jalal Fadili
 * @date 2008-04-22
 *
 */

/*
 * @file tvmin_mex.cpp
 * @brief Fast TV-minimization
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <iostream>
#include "mex.h"
#include "matrix.h"

extern "C" void TV1D(double *,int,int,double,float);
extern "C" void TV4(double *,int,int,double,float);
extern "C" void TV8(double *,int,int,double,float);
extern "C" void TV16(double *,int,int,double,float);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  double* tmp1, *tmp2;
  int m, n;
  double lambda, lmin, lmax;
  int numdeep;
  
  if(nrhs!=2)
	 mexErrMsgTxt("2 inputs required");
  if(nlhs!=1)
	 mexErrMsgTxt("1 output required");
  tmp1	= mxGetPr(prhs[0]);
  m 	= mxGetM(prhs[0]); 
  n 	= mxGetN(prhs[0]);
  lambda  = *(mxGetPr(prhs[1]));

  tmp2 = (double*)malloc(m*n*sizeof(double));
  memcpy(tmp2,tmp1,m*n*sizeof(double));
  if(m==1 || n==1)	
  	TV1D(tmp2, m, n, lambda, 0);
  else			
  	TV8(tmp2, m, n, lambda, 0);
     // TV4(tmp2, m, n, lambda, 0);
  
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  tmp1 = mxGetPr(plhs[0]);
  memcpy(tmp1,tmp2,m*n*sizeof(double));
  
  free(tmp2);
  
  return;
}
