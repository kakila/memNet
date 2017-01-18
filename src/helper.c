/*
 * Copyright (C) 2013 - Juan Pablo Carbajal
 *
 * This progrm is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "helper.h"
#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

#include <stdlib.h>
#include <time.h>


int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1);
  }

  return(0);
}

void PrintInfo(void *mem, realtype t)
{
  int retval, kused;
  long int nst;
  realtype hused;

  retval = IDAGetLastOrder(mem, &kused);
  check_flag(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_flag(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.4Le %3ld  %1d %12.4Le\n",
         t, nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%10.4le %3ld  %1d %12.4le\n",
         t, nst, kused, hused);
#else
  printf("%10.4e %3ld  %1d %12.4e\n",
         t, nst, kused, hused);
#endif
}

void Print2File(FILE* fout, int N, realtype t, N_Vector y)
{
  realtype *yval;
  int i;

  yval  = NV_DATA_S(y);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fout, "%12.8Le ",t);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(fout, "%12.8le ",t);
#else
  fprintf(fout,"%12.8e ",t);
#endif

  for (i=0; i<N; ++i){
#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(fout,"%12.8Le ",yval[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(fout,"%12.8le ",yval[i]);
#else
  fprintf(fout,"%12.8e ",yval[i]);
#endif
  }

  fprintf(fout,"\n");

  fflush (fout);
  
}

void PrintFinalStats(void *mem)
{
  int retval;
  long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;

  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_flag(&retval, "IDAGetNumResEvals", 1);
  retval = IDADlsGetNumJacEvals(mem, &nje);
  check_flag(&retval, "IDADlsGetNumJacEvals", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDADlsGetNumResEvals(mem, &nreLS);
  check_flag(&retval, "IDADlsGetNumResEvals", 1);
  retval = IDAGetNumGEvals(mem, &nge);
  check_flag(&retval, "IDAGetNumGEvals", 1);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);
  printf("Number of root fn. evaluations     = %ld\n", nge);
}

// Input signals for learning tasks

#define TAU 6.283185307179586231995926937088370323181152343750 //two pi
#define PI  3.141592653589793115997963468544185161590576171875 // pi
realtype freqSequence (realtype t, realtype* ts, realtype* freq, size_t nT)
{
  unsigned int i, idx;
  realtype x, dur, phas[nT];

  // Time is before first time stamp
  if (t < ts[0]) return 0.0;

  // Calculate durations and phases
  phas[0] = 0.0;
  for (i=0; i<nT-1; i++)
  {
    dur  = ts[i+1] - ts[i];
    phas[i+1] = phas[i] + TAU*freq[i]*dur;
  }

  // Find where is t in the table ts
  i = 0;
  while (t >= ts[i] && i < nT) i++;
  idx = i-1;

  // Evaluate the corresponding phase
  x = TAU * freq[idx] * (t-ts[idx]) + phas[idx];

  return sin (x);
}

realtype __ts__[]= {0,0.3,0.7,1};
realtype __f__[] = {5,16,5};
size_t __nT__    = sizeof __f__ / sizeof __f__[0];

realtype fsweep (realtype t)
{
  return freqSequence (t, __ts__, __f__, __nT__);
}

realtype trapesoidPulse (realtype t, realtype t_ini, realtype t_end)
{
  realtype t_r, T, beta=1e-3;

  // Map time to [0,T] interval
  T   = t_end - t_ini;
  t_r = t - t_ini;

  if (t_r <= 0)
    t_r = 0.0;
  else if (t_r < beta)
    t_r = t_r/beta;
  else if (t_r >= beta && t_r <= T - beta)
    t_r = 1.0;
  else if (t_r > T - beta)
    t_r = (T - t_r)/beta;
  else if (t_r >= T)
    t_r = 0.0;

  return t_r;
}

realtype trapesoidTrain (realtype t, realtype t_cycle, realtype duty, realtype t_max)
{

  realtype t_end;

  if (t >= t_max) return 0.0;

  t_end = t_cycle;
  while (t_end < t)
    t_end += t_cycle;

  return trapesoidPulse (t, t_end - duty*t_cycle, t_end);

}

/*********************************************/
/*************** Random numbers **************/
realtype gaussrand()
{
  // Knuth polar method on p 122 of the 2nd volume of TAOCP
	static realtype V1, V2, S;
	static int phase = 0;
	realtype X;

	if(phase == 0) {
		do {
			realtype U1 = (realtype)rand() / RAND_MAX;
			realtype U2 = (realtype)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}

int rand_lim(int limit) {
/* return a random number between 0 and limit inclusive.
 */

    int divisor = RAND_MAX/(limit+1);
    int retval;

    do {
        retval = rand() / divisor;
    } while (retval > limit);

    return retval;
}
/*********************************************/


enum states {CTE,RUP,RDW,RND,OFF};
#define MAX_IN 10 // Maximum number of independent inputs

realtype trapesoidTrain_N (realtype t,unsigned int N, realtype t_cycle, realtype duty, realtype t_max)
{
  FILE *fout;
  
  static realtype t_last=-1;
  static unsigned char init = FALSE;
  static char fname[80];
  static unsigned int nsamples, A;
  static int fpos;


  if (!init)
  {
    sprintf (fname,"trapesoidTrain%d.dat",N);
    t_last = t;
    A = rand_lim(N-1) + 1;
  
    fout = fopen (fname,"w");

    fprintf (fout,"# Duration of symbol:\n");
    fprintf (fout,"# name: dT\n# type: scalar\n %.8g",t_cycle);

    fprintf (fout,"# Random symbols for each signal (time of swap, symbol)\n");
    fprintf (fout,"# name: symbols\n# type: matrix\n# rows: ");
    fpos = ftell(fout);
    fprintf (fout,"00000000\n# columns: 2\n");
    // Length of the matrix is not know a priori in this function
    // We will update by the length using nsamples and fpos
    fprintf(fout," %g %d\n",t_last,A);

    fclose (fout);

    nsamples = 1;
    init = TRUE;

  } 
  
  if (t >= t_max) return 0.0;
  if ((t - t_last) >= duty*t_cycle && (t - t_last) < t_cycle) return 0.0;
  
  if ((t - t_last) >= t_cycle)
  {

    t_last = t;
    A = rand_lim(N-1) + 1;
    
    nsamples += 1;

    fout = fopen(fname,"r+");
    // Update length
    fseek (fout, fpos, SEEK_SET);
    fprintf (fout, "%-8d\n",nsamples);

    // Append to file
    fseek (fout, 0, SEEK_END);
    fprintf(fout," %g %d\n",t_last,A);

    fclose(fout);
  }

  return A*trapesoidPulse (t-t_last, 0.0, duty*t_cycle);

}
realtype symbol_osczerom (realtype t, unsigned int N, realtype t_cycle, realtype t_max, unsigned int i, const char* fname)
{
  // Encode N symbols as oscillatory signals with zero mean.
  static realtype t_last[MAX_IN]={0.0};
  static unsigned int index[MAX_IN]={0};
  static unsigned char init=FALSE, init2[MAX_IN] ={FALSE};
  static realtype** A;
  static realtype* wS;
  static unsigned int nsamples;
  static int fpos;


  FILE *fout;

  unsigned int j,k;
  int d;
  
  realtype retval=0.0;
  realtype scA = 1.5;

  if (!init)
  {
    // Initialize weight matrix
    A = (realtype **)malloc(N * sizeof(realtype*));
    for(j = 0; j < N; j++)
      A[j] = (realtype *)malloc(N * sizeof(realtype));
    if (NULL == A)
    {
      printf ("Error in allocation A\n");
      exit(1);
    }

// EXPONENTIAL MATRIX
    for (j=0; j < N; j++)
    {
      retval = 0.0;
      for (k=0; k < N; k++)
        {
          d = k-j;
          //A[j][k] = exp (- (abs (k - j) + 1.0) / scA) / exp (-1.0 / scA);
          A[j][k] = sin (TAU/N * (scA*d + N/4.0) ) * exp (- (abs (d) + 1.0) / scA);
          retval += fabs(A[j][k]);
        }
      for (k=0; k < N; k++)
          A[j][k] = A[j][k] / retval;
    }


    // Initialize frequencies
    wS = (realtype *)malloc(N * sizeof(realtype));
    if (NULL == wS) {
      printf ("Error in allocation wS\n");
      exit(1);
    }
    for (j=0; j < N; j++)
      wS[j] = (TAU/t_cycle) * (j+1);

    fout = fopen (fname,"w");

    fprintf (fout,"# Duration of symbol:\n");
    fprintf (fout,"# name: dT\n# type: scalar\n %.8g",t_cycle);

    fprintf (fout,"# Frequencies in radians:\n");
    fprintf (fout,"# name: wS\n# type: matrix\n# rows: 1\n# columns: %d\n",N);
    for (j=0; j < N; j++)
      fprintf (fout," %.8g\n",wS[j]);
    fprintf (fout,"\n");

    fprintf (fout,"# Mixing coeffs:\n");
    fprintf (fout,"# name: A\n# type: matrix\n# rows: %d\n# columns: %d\n",N,N);
    for (j=0; j < N; j++)
    {
      for (k=0; k < N; k++)
        fprintf (fout," %.8g ",A[j][k]);
      fprintf (fout,"\n");
    }
    fprintf (fout,"\n");

    fprintf (fout,"# Random symbols for each signal (time of swap, symbol, signal)\n");
    fprintf (fout,"# name: symbols\n# type: matrix\n# rows: ");
    fpos = ftell(fout);
    fprintf (fout,"00000000\n# columns: 3\n");
    // Length of the matrix is not know a priori in this function
    // We will update by the length using nsamples and fpos
    fclose (fout);

    nsamples = 0;
    init = TRUE;
  }

  if (!init2[i])
  {
    index[i] = rand_lim (N-1);
    t_last[i] = t;

    nsamples += 1;

    fout = fopen(fname,"r+");
    // Update length
    fseek (fout, fpos, SEEK_SET);
    fprintf (fout, "%-8d\n",nsamples);

    // Append to file
    fseek (fout, 0, SEEK_END);
    fprintf(fout," %g %d %d\n",t_last[i],index[i],i+1);

    fclose(fout);

    init2[i] = TRUE;
  }

  if (t >= t_max)
    return 0.0;

  if ((t - t_last[i]) >= t_cycle)
  {
    index[i] = rand_lim (N-1);
    t_last[i] = t;

    nsamples += 1;

    fout = fopen(fname,"r+");
    // Update length
    fseek (fout, fpos, SEEK_SET);
    fprintf (fout, "%-8d\n",nsamples);

    // Append to file
    fseek (fout, 0, SEEK_END);
    fprintf(fout," %g %d %d\n",t_last[i],index[i],i+1);

    fclose(fout);

  }

  retval = 0.0;
  for (j=0; j<N; j++)
    retval += A[index[i]][j] * sin (wS[j] * (t - t_last[i]));

  return retval;
}

realtype symbol_osczeromS (realtype t, unsigned int N, realtype t_cycle, realtype t_max, unsigned int i)
{
  // Encode N symbols as 2 oscillatory signals with zero mean.
  static realtype t_last[2]={0.0};
  static unsigned int index[2]={0};
  static unsigned char init=FALSE, init2[2] ={FALSE};
  static realtype** A;
  static realtype* wS;
  static unsigned int symbol[2];

  FILE *fout;
  char fname[80];
  unsigned int j,k;
  realtype retval=0.0;

  realtype scA = 0.1;
  sprintf (fname,"symbol_osczerom%d.dat",N);

  if (!init)
  {
    // Initialize weight matrix
    A = (realtype **)malloc(N * sizeof(realtype*));
    for(j = 0; j < N; j++)
      A[j] = (realtype *)malloc(N * sizeof(realtype));
    if (NULL == A)
    {
      printf ("Error in allocation A\n");
      exit(1);
    }

   // EXPONENTIAL MATRIX
    for (j=0; j < N; j++)
    {
      for (k=0; k < N; k++)
          A[j][k] = exp (- (abs (k - j) + 1.0) / scA) / exp (-1.0 / scA);
    }

    // Initialize frequencies
    wS = (realtype *)malloc(N * sizeof(realtype));
    if (NULL == wS) {
      printf ("Error in allocation wS\n");
      exit(1);
    }
    for (j=0; j < N; j++)
      wS[j] = (TAU/t_cycle) * (j+1);

    fout = fopen (fname,"w");

    fprintf (fout,"# Duration of symbol:\n");
    fprintf (fout,"# name: dT\n# type: scalar\n %.8g",t_cycle);

    fprintf (fout,"# Frequencies in radians:\n");
    fprintf (fout,"# name: wS\n# type: matrix\n# rows: 1\n# columns: %d\n",N);
    for (j=0; j < N; j++)
      fprintf (fout," %.8g\n",wS[j]);
    fprintf (fout,"\n");

    fprintf (fout,"# Mixing coeffs:\n");
    fprintf (fout,"# name: A\n# type: matrix\n# rows: %d\n# columns: %d\n",N,N);
    for (j=0; j < N; j++)
    {
      for (k=0; k < N; k++)
        fprintf (fout," %.8g ",A[j][k]);
      fprintf (fout,"\n");
    }
    fprintf (fout,"\n");

    fprintf (fout,"# Sequence of symbols for each signal (time of swap, symbol, signal)\n");
//    fprintf (fout,"# name: symbols\n# type: matrix\n# rows: %d\n# columns: 3\n",2*(int)pow(N,2.0));
    fprintf (fout,"# name: symbols\n# type: matrix\n# rows: %d\n# columns: 3\n",2*(int)(t_max/t_cycle));
    fclose (fout);

    init = TRUE;
  }

  if (!init2[i])
  {

    t_last[i] = t;
    symbol[i] = index[i]/N;
    symbol[i] = (i == 0) ? symbol[i] : index[i]-N*symbol[i];

    fout = fopen(fname,"a");
    // Append to file
    fprintf(fout," %g %d %d\n",t_last[i],symbol[i],i+1);

    fclose(fout);

    index[i] = index[i]+1;
    init2[i] = TRUE;
  }

  if (t >= t_max)
    return 0.0;

  if ((t - t_last[i]) >= t_cycle)
  {
    t_last[i] = t;
    symbol[i] = index[i]/N;
    symbol[i] = (i == 0) ? symbol[i] : index[i]-N*symbol[i];

    fout = fopen(fname,"a");
    // Append to file
    fprintf(fout," %g %d %d\n",t_last[i],symbol[i],i+1);

    fclose(fout);

    index[i] = index[i]+1;

    if (index[i] >= (int)pow(N,2.0))
      index[i] = 0;
    
  }

  retval = 0.0;
  for (j=0; j<N; j++)
    retval += A[symbol[i]][j] * sin (wS[j] * (t - t_last[i]));

  return retval;
}

realtype symbol_osczerom_Servo (realtype t, unsigned int N, realtype t_cycle, realtype t_max, unsigned int i)
{
  // Read symbols from file.
  static realtype t_last[4]={0.0};
  static unsigned char init=FALSE, init2[4] ={FALSE};
  static realtype** A;
  static realtype* wS;
  static unsigned int symbol[4];
  static int fpos;

  FILE *fout;
  FILE *fin;
  char fname[80];
  char fname_in[]="servo_train.dat";
  unsigned int j,k;
  int d;
  realtype retval=0.0;

  realtype scA = 0.1;
  sprintf (fname,"symbol_osczerom%d.dat",N);

  if (!init)
  {
    // Initialize weight matrix
    A = (realtype **)malloc(N * sizeof(realtype*));
    for(j = 0; j < N; j++)
      A[j] = (realtype *)malloc(N * sizeof(realtype));
    if (NULL == A)
    {
      printf ("Error in allocation A\n");
      exit(1);
    }

// EXPONENTIAL MATRIX
    for (j=0; j < N; j++)
    {
      retval = 0.0;
      for (k=0; k < N; k++)
        {
          d = k-j;
          //A[j][k] = exp (- (abs (k - j) + 1.0) / scA) / exp (-1.0 / scA);
          A[j][k] = sin (TAU/N * (scA*d + N/4.0) ) * exp (- (abs (d) + 1.0) / scA);
          retval += fabs(A[j][k]);
        }
      for (k=0; k < N; k++)
          A[j][k] = A[j][k] / retval;
    }


    // Initialize frequencies
    wS = (realtype *)malloc(N * sizeof(realtype));
    if (NULL == wS) {
      printf ("Error in allocation wS\n");
      exit(1);
    }
    for (j=0; j < N; j++)
      wS[j] = (TAU/t_cycle) * (j+1);

    fout = fopen (fname,"w");

    fprintf (fout,"# Duration of symbol:\n");
    fprintf (fout,"# name: dT\n# type: scalar\n %.8g",t_cycle);

    fprintf (fout,"# Frequencies in radians:\n");
    fprintf (fout,"# name: wS\n# type: matrix\n# rows: 1\n# columns: %d\n",N);
    for (j=0; j < N; j++)
      fprintf (fout," %.8g\n",wS[j]);
    fprintf (fout,"\n");

    fprintf (fout,"# Mixing coeffs:\n");
    fprintf (fout,"# name: A\n# type: matrix\n# rows: %d\n# columns: %d\n",N,N);
    for (j=0; j < N; j++)
    {
      for (k=0; k < N; k++)
        fprintf (fout," %.8g ",A[j][k]);
      fprintf (fout,"\n");
    }
    fprintf (fout,"\n");

    fprintf (fout,"# Sequence of symbols for each signal (time of swap, symbol, signal)\n");

    // Read number of data points from file
    fin = fopen (fname_in,"r");
    fscanf (fin,"%d",&k);
    fpos = ftell(fin);
    fclose (fin);

    fprintf (fout,"# name: symbols\n# type: matrix\n# rows: %d\n# columns: 3\n",4*k);
    fclose (fout);

    init = TRUE;

  }

  if (!init2[i])
  {
    fin = fopen (fname_in,"r");
    fseek (fin, fpos, SEEK_SET);
    fscanf (fin,"%d",&symbol[i]);
    fpos = ftell (fin);
    fclose (fin);

    t_last[i] = t;

    fout = fopen(fname,"a");
    // Append to file
    fprintf(fout," %g %d %d\n",t_last[i],symbol[i],i+1);
    fclose(fout);

    init2[i] = TRUE;
  }

  if (t >= t_max)
    return 0.0;

  if ((t - t_last[i]) >= t_cycle)
  {
    fin = fopen (fname_in,"r");
    fseek (fin, fpos, SEEK_SET);
    fscanf (fin,"%d",&symbol[i]);
    fpos = ftell (fin);
    fclose (fin);

    t_last[i] = t;

    fout = fopen(fname,"a");
    // Append to file
    fprintf(fout," %g %d %d\n",t_last[i],symbol[i],i+1);

    fclose(fout);

  }

  retval = 0.0;
  for (j=0; j<N; j++)
    retval += A[symbol[i]][j] * sin (wS[j] * (t - t_last[i]));

  return retval;
}

realtype symbol_cascade (realtype t, unsigned int N, realtype t_cycle, realtype t_max, unsigned int i)
{
  // Generates a sequence of N values in [0,1]. Each value remains constant for a random period of time.
  unsigned int j;
  realtype retval=0.0;

  for (j=0; j<N-1; j++)
    retval += trapesoidTrain_rndu (t, t_cycle, 0.49, t_max, j+N*i);
  retval = retval / (N-1);

  return retval;
}

realtype trapesoidTrain_rndu (realtype t, realtype t_cycle, realtype duty_var, realtype t_max, unsigned int i)
{

  static realtype t_end[MAX_IN]={0.0}, duty[MAX_IN]={0.0};
  static enum states state[MAX_IN]={[0 ... MAX_IN-1] = RND};

  if (t >= t_max) return 0.0;

  if (t >= t_end[i]) state[i] = RND;

  if (state[i] == RND)
  {
    while (t_end[i] < t)
      t_end[i] += t_cycle;

    duty[i] = 0.5 + duty_var * (2.0 * rand () / (RAND_MAX + 1.0) - 1.0);

    state[i] = CTE;
  }

  return trapesoidPulse (t, t_end[i] - duty[i]*t_cycle, t_end[i]);

}

realtype trapesoidTrain_rnde (realtype t, realtype d_min, realtype beta, realtype t_max)
{
  static enum states state=CTE;
  static realtype t_trig=0.0, value=0.0, t0=0.0, dt = 1e-3;

  if (state == CTE)
    if (t-t0 >= d_min)
      state = RND;

  if (state == RND)
  {
    if (rand()/(RAND_MAX + 1.0) > exp( -(t-t0 - d_min)/beta ))
    {
      state = (value > 0.0)? RDW : RUP;
      t_trig = t;
    }
  }

  if (state == RUP)
  {
    value = (t - t_trig) / dt;
    if (value > 1.0)
    {
      value = 1.0;
      state = CTE;
      t0 = t;
    }
  }

  if (state == RDW)
  {
    value = 1.0 - (t - t_trig)/ dt;
    if (value < 0.0)
    {
      value = 0.0;
      state = CTE;
      t0 = t;
    }
  }

  if (t > t_max)
  {
    if (value > 0.0)
      value = 1.0 - (t - t_max)/ dt;
    if (value < 0.0)
    {
      value = 0.0;
      state = OFF;
    }
  }

  return value;

}

realtype stcproc (realtype t, realtype mean, realtype d_min, unsigned int i)
{
  static realtype t0[MAX_IN]={0.0},
                  y0[MAX_IN]={0.0},
                  r1[MAX_IN]={0.0},
                  r2[MAX_IN]={0.0},
                  t_prev[MAX_IN]={0.0};
  realtype duration, dt;

  duration = t-t0[i];

  if (duration > d_min)
  {
    t0[i]  = t;
    r1[i]  = 2.0 * rand()/(RAND_MAX+1.0) - 1.0;
    r1[i]  = r1[i]*fabs(r1[i]);
    r2[i]  = 2.0 * rand()/(RAND_MAX+1.0) -1.0;
  }

  dt     = t-t_prev[i];
  t_prev[i] = t;

  y0[i] = y0[i] + (r2[i] - (mean + r1[i]) * y0[i] * (1 - y0[i]*y0[i])) * dt;

  return y0[i];
}

/****** Detrended Normalized Finite Wiener process in t [0,1]*******/ 
realtype defiWiener (realtype t, int Nmax, int n0, realtype s, int id)
{
  static char init[MAX_IN] = {FALSE};
  static realtype* A[MAX_IN];

  int n;
  realtype retval;

  // Sample random vector of weights only once
  if (!init[id])
  {
    if (id > MAX_IN)
    {
      printf ("Too many inputs, increase MAX_IN\n");
      exit(1);
    }

    // Initialize weight vector
    A[id] = (realtype *)malloc(Nmax * sizeof(realtype));

    if (NULL == A[id])
    {
      printf ("Error in allocation A\n");
      exit(1);
    }

    for(n = 0; n < Nmax; n++)
      A[id][n] = gaussrand()/sqrt(2);

    // Balance spectrum
    if (s > 0)
    {
      for(n = 0; n < Nmax; n++)
      {
        retval = -log (fabs (n-n0)/(2.0*s) + 1.0);
        A[id][n] = exp ( retval * retval) * A[id][n];
      }
    }
    
    init[id] = TRUE;
  }
  
  retval = 0;
  for (n = 1; n <= Nmax; n++)
  {
    if (n%2 == 0)
      retval += A[id][n-1] * sin (PI*n*t) / ( PI * n);
    else
      retval += A[id][n-1] * (sin (PI*n*t) - 2.0 / (PI * n)) / ( PI * n);
  }

  return retval;

}
/*********************************************/


/*********************************************/
/******** General Frequency Encoding *********/

// Mean value of encoding in interval [0,T]
realtype mean_FE(realtype w, realtype T, realtype b)
{
  return (1.0/b/w/w) * (
         sin(w*b) * (1.0 - cos(w*T)) 
       - sin(T*w) * (1.0 - cos(b*w))
       );     
}

// Read data from file.
void readData_FE(const char* fname, realtype** DATA, unsigned int len, unsigned int dim)
{
    FILE *fin;
    unsigned int i,j;
    float tmp;

    fin = fopen (fname,"r");
    if (NULL == fin)
    {
      printf ("Error opening %s\n",fname);
      exit(1);
    }
    // Read first line
    fscanf (fin,"%f %f",&tmp, &tmp);
    
    //Read data
    for (i=0; i < len; i++)
    {
      for (j=0; j < dim; j++)
      {
        fscanf (fin,"%f",&tmp);
        DATA[i][j] = tmp;
      }
    }
    
    fclose (fin);
}

// Zero-ends trapezoid envelope in [0,1]
realtype Trpz (realtype t, realtype beta)
{
  if (t <= 0 || t >= 1.0)
    return 0.0;
  else if (t < beta)
    return t / beta;
  else if (t >= beta && t <= (1.0 - beta))
    return 1.0;
  else if (t > (1.0 - beta))
    return (1.0 - t) / beta;
  
  return 0.0;
}

// Encoding
realtype FE (realtype t, realtype t_cycle, realtype b, unsigned int idx, const char* fname)
{
  static unsigned char init=FALSE;
  static realtype* t_last;
  static realtype* w;
  static realtype** DATA;

  static unsigned int *ndim;
  static unsigned int nsample;
  static int len, dim;
  static char fname_out[80];
  
  FILE *fout, *fin;
  unsigned int i, k;
  
  if (!init)
  {
    sprintf (fname_out,"%s_FE.dat",fname);

    // Read length of data, inputs size 
    fin = fopen (fname,"r");    
    if (NULL == fin)
    {
      printf ("Error opening %s\n",fname);
      exit(1);
    }
    fscanf (fin,"%d %d",&len, &dim);
    fclose (fin);

    // Initialize DATA matrix
    DATA = (realtype **)malloc(len * sizeof(realtype*));
    for(i = 0; i < len; i++)
      DATA[i] = (realtype *)malloc(dim * sizeof(realtype));
    if (NULL == DATA)
    {
      printf ("Error in allocation A\n");
      exit(1);
    }

    // Read all data from file
    readData_FE(fname, DATA, len, dim);

    if (idx > dim)
    {
      printf ("init Error: idx(%d) bigger than dim(%d)\n",idx,dim);
      exit(1);
    }

    w = (realtype *)malloc(dim * sizeof(realtype));
    if (NULL == w)
    {
      printf ("init Error in allocation w\n");
      exit(1);
    }
    for (i=0; i<dim; i++)
      w[i] = TAU * DATA[nsample][i] / t_cycle;

    t_last = (realtype *)malloc(dim * sizeof(realtype));
    if (NULL == t_last)
    {
      printf ("init Error in allocation t_last\n");
      exit(1);
    }
    for (i=0; i<dim; i++)
      t_last[i] = t;

    ndim = (unsigned int*)malloc(dim * sizeof(unsigned int));
    if (NULL == ndim)
    {
      printf ("init Error in allocation ndim\n");
      exit(1);
    }
    for (i=0; i<dim; i++)
      ndim[i] = 0;

    // Write file header & first data
    fout = fopen (fname_out,"w+");
    if (NULL == fout)
    {
      printf ("Error opening %s\n",fname_out);
      exit(1);
    }
    
    fprintf (fout,"# Duration of symbol:\n");
    fprintf (fout,"# name: dT\n# type: scalar\n %.8g\n",t_cycle);
    fprintf (fout,"# Sequence of symbols for each signal (time of swap, freqs)\n");
    fprintf (fout,"# name: encode_event\n# type: matrix\n# rows: %d\n# columns: %d\n",len,dim+1);

    fprintf(fout," %g ",t_last[idx]);
    for (i=0; i<dim; i++)
      fprintf(fout,"%g ",w[i]);
    fprintf(fout,"\n");

    fclose (fout);

    printf ("\r%d/%d",nsample+1,len);
    fflush(stdout);

    nsample = 1;
 
    init = TRUE;
    
  }

  if (idx >= dim)
  {
    printf ("Error: idx(%d) bigger than dim(%d)\n",idx,dim);
    exit(1);
  }

  if ((t - t_last[idx]) >= t_cycle)
  {
    if(nsample >= len)
      return 0.0;

    t_last[idx] = t;
      
    // Encode value in freq
    w[idx] = TAU * DATA[nsample][idx] / t_cycle;

    ndim[idx] = 1;
    k = 0;
    for (i=0;i<dim;i++)
      k +=ndim[i];
       
    if (k == dim)
    {
      // Append data to file
      fout = fopen (fname_out,"a");
      if (NULL == fout)
      {
        printf ("Error opening %s\n",fname_out);
        exit(1);
      }
      
      fprintf(fout," %g ",t_last[idx]);
      for (i=0; i<dim; i++)
        fprintf(fout,"%g ",w[i]);
      fprintf(fout,"\n");
      fclose (fout);

      for (i=0;i<dim;i++)
        ndim[i] = 0;
      nsample +=1;

      printf ("\r%d/%d",nsample,len);
      fflush(stdout);

    }
  }
 
  return ( sin( w[idx]*t ) - (1.0 - cos(w[idx]*t_cycle))/w[idx] )* Trpz( (t-t_last[idx])/t_cycle , b );
}

/*********************************************/
/********* General Phase Encoding ************/
#define N_FREQ 2
void readData(FILE* fin, realtype*** A, realtype nlvl, unsigned int len, unsigned int dim)
{
    // FIXME Works for N_FREQ == 2
    unsigned int i,j;
    float tmp;

    //Read data
    for (i=0; i < len; i++)
    {
      for (j=0; j < dim; j++)
      {
        fscanf (fin,"%f",&tmp);
        // Make sure A is in [0,1]. Clip.
        tmp = tmp + nlvl*gaussrand();
        if (tmp < 0)
          tmp = 0;
        else if (tmp>1)
          tmp = 1; 
        
        A[0][i][j] = tmp;
        A[1][i][j] = 1 - A[0][i][j];
      }
    }
}


// Encoding
// Encodes real valued input data in [0,1] as mixing amplitudes of 
// N_FREQ given frequencies
realtype PE (realtype t, realtype t_cycle, realtype nlvl, unsigned int idx, const char* fname)
{
  static unsigned char init=FALSE;
  static realtype* t_last;
  static realtype w[N_FREQ];

  static realtype*** A;

  static unsigned int *ndim;
  static unsigned int nsample;
  static int len, dim;
  static char fname_out[80];
  
  FILE *fout, *fin;
  unsigned int i, k;
  realtype retval;
  
  
  if (!init)
  {
    sprintf (fname_out,"%s_PE.oct.dat",fname);

    // Read length of data, inputs size and number of symbols(ignored)
    fin = fopen (fname,"r");    
    if (NULL == fin)
    {
      printf ("Error opening %s\n",fname);
      exit(1);
    }
    fscanf (fin,"%d %d %d",&len, &dim, &i);
    printf ("samples: %d, inputs: %d, values: %d\n",len,dim,i);
    
    // Initialize Amplitude matrix
    A = (realtype ***)malloc(N_FREQ * sizeof(realtype**));
    for(i = 0; i < N_FREQ; i++)
    {
      A[i] = (realtype **)malloc(len * sizeof(realtype*));
      for(k = 0; k < len; k++)
        A[i][k] = (realtype *)malloc(dim * sizeof(realtype));
    }
    if (NULL == A)
    {
      printf ("Error in allocation A\n");
      exit(1);
    }
    // Read all data from file as amplitudes
    // Add noise if given
    readData(fin, A, nlvl, len, dim);
    fclose (fin);

    nsample = 0;
        
        
    // Define frequencies
    w[0] = TAU * 2.0 / t_cycle;
    w[1] = TAU * 3.0 / t_cycle;
    
    t_last = (realtype *)malloc(dim * sizeof(realtype));
    if (NULL == t_last)
    {
      printf ("init Error in allocation t_last\n");
      exit(1);
    }
    for (i=0; i<dim; i++)
      t_last[i] = t;

    ndim = (unsigned int*)malloc(dim * sizeof(unsigned int));
    if (NULL == ndim)
    {
      printf ("init Error in allocation ndim\n");
      exit(1);
    }
    for (i=0; i<dim; i++)
      ndim[i] = 0;

    // Write file header & first data
    fout = fopen (fname_out,"w+");
    if (NULL == fout)
    {
      printf ("Error opening %s\n",fname_out);
      exit(1);
    }
    
    fprintf (fout,"# Duration of symbol:\n");
    fprintf (fout,"# name: dT\n# type: scalar\n %.8g\n",t_cycle);

    fprintf (fout,"# Encoding frequencies:\n");
    fprintf (fout,"# name: wS\n# type: matrix\n# rows: 1\n# columns: %d\n",N_FREQ);
    for (i=0; i<N_FREQ; i++)
      fprintf(fout," %g",w[i]);
    fprintf(fout,"\n");
    
    
    fprintf (fout,"# Sequence of symbols for each signal (time of swap, freqs)\n");
    fprintf (fout,"# name: encode_event\n# type: matrix\n# rows: %d\n# columns: %d\n",len,N_FREQ*dim+1);

    fprintf(fout," %g ",t_last[idx]);
    for (i=0; i<dim; i++)
      fprintf(fout,"%g %g ",A[0][nsample][i],A[1][nsample][i]);
    fprintf(fout,"\n");
    fclose (fout);

    printf ("\r%d/%d",nsample+1,len);
    fflush(stdout);

    init = TRUE;
    
  }

  if (idx >= dim)
  {
    printf ("Error: idx(%d) bigger than dim(%d)\n",idx,dim);
    exit(1);
  }

  if ((t - t_last[idx]) >= t_cycle)
  {
    t_last[idx] = t;
    
    ndim[idx] = 1;
    k = 0;
    for (i=0;i<dim;i++)
      k +=ndim[i];

    if (k == dim)
    {
      // Reset counter
      for (i=0;i<dim;i++)
        ndim[i] = 0;
      nsample +=1;

      if(nsample >= len)
        nsample = len-1;

      // Append data to file
      fout = fopen (fname_out,"a");
      if (NULL == fout)
      {
        printf ("Error opening %s\n",fname_out);
        exit(1);
      }
      
      fprintf(fout," %g ",t_last[idx]);
      for (i=0; i<dim; i++)
        fprintf(fout,"%g %g ",A[0][nsample][i],A[1][nsample][i]);
      fprintf(fout,"\n");
      fclose (fout);

      printf ("\r%d/%d",nsample+1,len);
      fflush(stdout);

    }

  }
  
  retval = 0.0;
  for (i=0; i<N_FREQ; i++)
    retval += A[i][nsample][idx]*sin(w[i]*(t - t_last[idx]));
 
  return retval;
}
