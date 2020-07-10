#ifdef STANDALONE
# include <stdio.h>
# include <stdlib.h>
# include "area.h"
# define CALLOC calloc
#else
# include "mex.h"
# define CALLOC mxCalloc
#endif

#define TIMING

#ifdef TIMING
#include <sys/time.h>
#endif
#include <math.h>
#include <string.h>

typedef unsigned char PIX;
double time_diff(struct timeval *tv1, struct timeval *tv2);
void compute_means(float *L, float *R, int width, int height, int wx, int wy, float *meanL, float *meanR);

/* matrix reference macros Matlab element order (column wise) */
#define REF(p,x,y)  p[(y)*width+x]
#define REF_ML(p,x,y)  p[(x)*height+(y)]
#define REF3(p,x,y,z)  p[((z)*width+(x))*height+(y)]


// function implementing ZNCC between patches
inline double match_patch(float *L, float *meanL, int xL, int yL, float *R, float *meanR, int xR, int yR, int wx, int wy,
        int width, int height) {
    int i, j;
    double sumL, sumR, sumLR, dL, dR, num, den, muL, muR;

    sumL = 0;
    sumR = 0;
    sumLR = 0;
    muL = REF(meanL, xL, yL);
    muR = REF(meanR, xR, yR);
    
    // iterate over window and calculate sums
    for (i = -wx; i <= wx; ++i) {
        for (j = -wy; j <= wy; ++j) {
            dL = ((double)REF(L, xL+i, yL+j) - muL);
            dR = ((double)REF(R, xR+i, yR+j) - muR);
            sumL += dL * dL;
            sumR += dR * dR;
            sumLR += dL * dR;
        }
    }
    
    num = sumLR;
    den = sqrt(sumL * sumR);
    
    // if denominator is zero, zncc is undefined
    if (den < 1e-6) {
        return mxGetInf();
    }
    // ro = 1 - (num / den)
    return (den - num) / den;
}


// function implementing classic stereo matching strategy using ZNCC
void stereo_matching_zncc(float *L, float *R, double *dsi, int width, int height, int wx, int wy, int dmin, int dmax) {
    int x, y, disp, disp_lo, disp_hi;
    float *meanL, *meanR;
    double zncc;
    
    // allocate mean matrices and calculate means of all patches
    meanL = (float*) mxCalloc(height * width, sizeof(float));
    meanR = (float*) mxCalloc(height * width, sizeof(float));
    compute_means(L, R, width, height, wx, wy, meanL, meanR);

    for (y = wy; y < height - wy; ++y) {
        for (x = wx; x < width - wx; ++x) {
            disp_lo = dmin;
            disp_hi = (dmax <= x - wx) ? dmax : (x - wx);

            for (disp = disp_lo; disp <= disp_hi; ++disp) {
                zncc = match_patch(L, meanL, x, y, R, meanR, x-disp, y, wx, wy, width, height);
                REF3(dsi, x, y, disp - dmin) = zncc;
            } 
        }
    }

    // free the mean matrices
    mxFree(meanL);
	mxFree(meanR);
}


// compute means for all possible patches
void compute_means(float *L, float *R, int width, int height, int wx, int wy, float *meanL, float *meanR) {
    int x, y, i, j;
    float sumL, sumR, area;
    area = (2*wx+1) * (2*wy+1);

    // simple brute force
    for (y = wy; y < height - wy; ++y) {
        for (x = wx; x < width - wx; ++x) {
            sumL = 0;
            sumR = 0;
            for (i = -wx; i <= wx; ++i) {
                for (j = -wy; j <= wy; ++j) {
                    sumL += REF(L, x+i, y+j);
                    sumR += REF(R, x+i, y+j);
                }
            }
            REF(meanL, x, y) = sumL / area;
            REF(meanR, x, y) = sumR / area;
        }
    }
}

/*
 * disp = stereo(L, R, w, disprange, metric)
 *
 *  L and R are images (uint8 or double) of size NxM
 *  w is a 1-vector wx=wy=w or a 2-vector [wx wy] window size
 *  disprange is a 2-vector, D = abs(disprange(2)-disprange(1))
 *
 *  disp is a NxMxD matrix of int16
 */
void 
mexFunction(
		 int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  unsigned int width_l, height_l, height_r, width_r, width, height;
  int wx, wy, i, j, dispmin, dispmax;
  float *leftI, *rightI, *leftI_ptr, *rightI_ptr; 
  double *p, *scoresD;
  mwSize    dims[3];
  double Z_NAN = mxGetNaN();
  double Z_INFINITY = mxGetInf();
#ifdef  TIMING
  struct timeval	t0, t1, t2, t3, t4;
#endif
  char  metric[14];


#ifdef  TIMING
  gettimeofday (&t0, NULL);
#endif

  /* Check for proper number of arguments */

  switch (nrhs) {
  case 4:
    strcpy(metric, "Classic");
    break;
  case 5:
    if (!mxIsChar(prhs[4]))
        mexErrMsgTxt("approach must be specified by a string");
    mxGetString(prhs[4], metric, 8);
    break;
  default:
    mexErrMsgTxt("expecting 4 or 5 arguments");
    return;
  }

  /* Check that input images are same size */
  height_l = mxGetM(prhs[0]);
  width_l = mxGetN(prhs[0]);
  height_r = mxGetM(prhs[1]);
  width_r = mxGetN(prhs[1]);

  if ((height_l != height_r) || (width_l != width_r))
    mexErrMsgTxt("Left and right images must be the same size");

  height = height_l;
  width = width_l;

  /* get window size */
  switch (mxGetNumberOfElements(prhs[2])) {
  case 1:
    wx = wy = mxGetScalar(prhs[2]);
    break;
  case 2: {
    double  *t = mxGetPr(prhs[2]);

    wx = t[0];
    wy = t[1];
    break;
  }
  default:
    mexErrMsgTxt("Window size must be 1- or 2-vector");
    return;
  }

  /* get disparity range */
  switch (mxGetNumberOfElements(prhs[3])) {
  case 1:
    dispmin = 0;
    dispmax = (int) mxGetScalar(prhs[3]);
    break;
  case 2: {
    double  *d = mxGetPr(prhs[3]);

	dispmin = (int)d[0];
	dispmax = (int)d[1];
    break;
  }
  default:
    mexErrMsgTxt("Disparities must be 1- or 2-vector");
    return;
  }

#ifdef  TIMING
  gettimeofday (&t1, NULL);
#endif

  /* allocate float images to hold copies of the images */
  leftI = (float*) mxCalloc (width * height, sizeof(float));
  rightI = (float*) mxCalloc (width * height, sizeof(float));

  leftI_ptr = leftI;
  rightI_ptr = rightI;

  if (mxGetClassID(prhs[0]) != mxGetClassID(prhs[1])) 
    mexErrMsgTxt("Images must be of the same class");

  switch (mxGetClassID(prhs[0])) {
  case mxDOUBLE_CLASS: {
      double   *l, *r;

      l = (double *) mxGetData(prhs[0]);
      r = (double *) mxGetData(prhs[1]);

      for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
          
          *(leftI_ptr++) = (float) REF_ML(l,j,i);
          *(rightI_ptr++) = (float) REF_ML(r,j,i);
        }
      }
    }
    break;

  case mxUINT8_CLASS: {
      PIX   *l, *r;

      l = (PIX *) mxGetData(prhs[0]);
      r = (PIX *) mxGetData(prhs[1]);

      for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
          
          *(leftI_ptr++) = (float) REF_ML(l,j,i);
          *(rightI_ptr++) = (float) REF_ML(r,j,i);
        }
      }
    }
    break;

  default:
    mexErrMsgTxt("Images must be double or uint8 class");
  }

  
#ifdef  TIMING
  gettimeofday (&t2, NULL);
#endif

  /* Create an NxMxD matrix for the return arguments */
  dims[0] = height_l;
  dims[1] = width_l;
  dims[2] = (int) (dispmax - dispmin);
  if (dims[2] < 0)
    dims[2] = -dims[2];
  dims[2] += 1;


  /* create 3D array to hold similarity scores */
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  scoresD = (double *) mxGetData(plhs[0]);

  /* set all values to NaN now */
  for (i=0; i<dims[0]*dims[1]*dims[2]; i++)
    scoresD[i] = Z_NAN;

#ifdef  TIMING
  gettimeofday (&t3, NULL);
#endif

  if (strcmp(metric, "Classic") == 0) {
      stereo_matching_zncc(leftI, rightI, scoresD, width, height, wx, wy, dispmin, dispmax);
  }
  else {
      mexErrMsgTxt("Unknown approach");
  }

#ifdef  TIMING
  gettimeofday (&t4, NULL);

  mexPrintf("arg checking     %.1fms\n", time_diff(&t1, &t0)*1000);
  mexPrintf("image conversion %.1fms\n", time_diff(&t2, &t1)*1000);
  mexPrintf("score array      %.1fms\n", time_diff(&t3, &t2)*1000);
  mexPrintf("stereo           %.1fms\n", time_diff(&t4, &t3)*1000);
#endif

  /* free the temporary float images */
  mxFree(leftI);
  mxFree(rightI);
}

#ifdef  TIMING
double
time_diff(struct timeval *tv1, struct timeval *tv2)
{
    double  dt;
    int     du;

    dt = tv1->tv_sec - tv2->tv_sec;
    du = tv1->tv_usec - tv2->tv_usec;

    if (du < 0) {
        dt--;
        du += 1000000;
    }

    dt += du / 1e6;

    return dt;
}
#endif