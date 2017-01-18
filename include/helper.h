#ifndef HELPER_H
#define HELPER_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <math.h>

// MACROS
#define positive(a) (fmax (a,0.0))
     
// Function prototypes
int check_flag(void *flagvalue, char *funcname, int opt);
void PrintInfo(void *mem, realtype t);
void Print2File(FILE* fout, int N, realtype t, N_Vector y);
void PrintFinalStats(void *mem);
int rand_lim(int limit);

// Signal prototypes
realtype freqSequence (realtype t, realtype* ts, realtype* freq, size_t n);
realtype fsweep (realtype t);

realtype trapesoidPulse (realtype t, realtype t_ini, realtype t_end);
realtype trapesoidTrain (realtype t, realtype t_cycle, realtype duty, realtype t_max);
realtype trapesoidTrain_N (realtype t,unsigned int N, realtype t_cycle, realtype duty, realtype t_max);
realtype trapesoidTrain_rndu (realtype t, realtype t_cycle, realtype duty_var, realtype t_max, unsigned int i);
realtype trapesoidTrain_rnde (realtype t, realtype d_min, realtype beta, realtype t_max);
realtype stcproc (realtype t, realtype mean, realtype d_min, unsigned int i);
realtype symbol_cascade (realtype t, unsigned int N, realtype t_cycle, realtype t_max, unsigned int i);
realtype symbol_osczerom (realtype t, unsigned int N, realtype t_cycle, realtype t_max, unsigned int i, const char* fname);
realtype symbol_osczeromS (realtype t, unsigned int N, realtype t_cycle, realtype t_max, unsigned int i);
realtype symbol_osczerom_Servo (realtype t, unsigned int N, realtype t_cycle, realtype t_max, unsigned int i);

// Detrended band-limited balanced Wiener process
realtype defiWiener (realtype t, int Nmax, int n0, realtype s, int id);

// Frequency encoding
realtype FE (realtype t, realtype t_cycle, realtype b, unsigned int idx, const char* fname);

// Phase encoding
realtype PE (realtype t, realtype t_cycle, realtype nlvl, unsigned int idx, const char* fname);

#endif /* HELPER_H */
