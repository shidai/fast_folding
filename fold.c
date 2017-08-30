#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
#include "omp.h"
#include "read_fits.h"

void fold_main (search_mode *s, char *pred_name)
{
    int i,j,k;
    int jj;

    int ret;
    T2Predictor pred;

    long double phase0, phase;
    long double mjd;
    long double freq;

    if (ret=T2Predictor_Read(&pred, pred_name))
    //if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
    {
	    printf("Error: unable to read predictor\n");
	    exit(1);
    }

    s->nbin = 128;
    s->prof = (double *)malloc(sizeof(double)*s->nbin);
    for (i=0; i<s->nbin; i++)
    {
	    s->prof[i] = 0.0;
    }

    float tc = 0.0; 
    float t0 = 0.0; 
    int ncyc = 0;
    long double *freq_phase;
    freq_phase = (long double *)malloc(sizeof(long double)*s->nchan);
    long double *freq_period;
    freq_period = (long double *)malloc(sizeof(long double)*s->nchan);

    for (i=1; i<=s->nsub; i++)
    {
			for (j = 0; j<s->nsblk; j++)
			{
				tc = (s->nsblk*(i-1) + j)*s->tsample;
				mjd = s->mjd0 + (s->nsblk*(i-1) + j)*s->tsample/86400.0L;
				
				if (tc >= ncyc)
				{
					printf ("test %f %d %d %f\n", tc, ncyc, s->nsblk, s->tsample);
					t0 = (s->nsblk*(i-1) + j)*s->tsample;
					for (jj = 0; jj<s->nchan; jj++)
					{
						freq = (long double)s->freqs[jj];
						freq_phase[jj] = T2Predictor_GetPhase(&pred,mjd,freq);
						freq_period[jj] = 1.0/T2Predictor_GetFrequency(&pred,mjd,freq);   // second
						//freq_phase[jj] = 0.1;
						//freq_period[jj] = 1.0/218.81184385;   // second
					}
					ncyc += 2;
				}
				
				for (k = 0; k<s->nchan; k++)
				{
					phase0 = freq_phase[k] + (tc-t0)/freq_period[k];
					phase = (phase0 - floor(phase0));
					temp = phase*s->nbin;
					index = (int)(temp+0.5)>(int)temp?(int)temp+1:(int)temp;
					s->prof[index] += s->fdata[i-1][s->nchan*j+k];
				}
			}
		}

		for (i=0; i<s->nbin; i++)
		{
			printf ("%f\n", s->prof[i]);
		}

    free(freq_phase);
    free(freq_period);
}

