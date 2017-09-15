#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
#include "read_fits.h"

void fold_main (search_mode *s, char *pred_name)
{
    int i,j,k;

    int ret;
    T2Predictor pred;

    long double phase0, phase;
    long double mjd;
    long double freq;

    double temp;
    int index;

    int ncyc = 0;
    long double **freq_phase;
    long double **freq_period;
    //long double *freq_phase;
    //long double *freq_period;
    float t0;
    float tc;
    float t0_temp;

    if (ret=T2Predictor_Read(&pred, pred_name))
    //if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
    {
	    printf("Error: unable to read predictor\n");
	    exit(1);
    }

    s->nbin = 128;
    s->prof = (double *)malloc(sizeof(double)*s->nbin);
    memset(s->prof, 0.0, sizeof(double)*s->nbin);

		double prof[s->nsub][s->nbin];
		memset(prof, 0.0, sizeof(prof));

		ncyc = (int)(s->nsub*s->nsblk*s->tsample/s->tcyc) + 1;
    freq_phase = (long double **)malloc(sizeof(long double)*ncyc);
    freq_period = (long double **)malloc(sizeof(long double)*ncyc);
		for (i=0; i<ncyc; i++)
		{
			freq_phase[i] = (long double *)malloc(sizeof(long double)*s->nchan);
			freq_period[i] = (long double *)malloc(sizeof(long double)*s->nchan);
		}

    ///////////////////////////////////////////////
		omp_set_num_threads(10);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(s,ncyc,pred,freq_phase,freq_period)
#endif
		for (i=0; i<ncyc; i++)
		{
				int k;
				long double mjd, freq;
		    mjd = s->mjd0 + s->tcyc*i/86400.0L;
			  for (k = 0; k<s->nchan; k++)
			  {
				  freq = (long double)s->freqs[k];
				  freq_phase[i][k] = T2Predictor_GetPhase(&pred,mjd,freq);
				  freq_period[i][k] = 1.0/T2Predictor_GetFrequency(&pred,mjd,freq);   // second
			  }
		}
		printf ("number of cycles: %d\n", ncyc);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(s,freq_phase,freq_period,prof)
#endif
    for (i = 0; i<s->nsub; i++)
    {
			int j,k,ncyc,index;
			float tc,t0;
			long double phase0,phase;
			double temp;

	    for (j = 0; j<s->nsblk; j++)
	    {
		    tc = (s->nsblk*i + j)*s->tsample;
		    ncyc = (int)(tc/s->tcyc);
				t0 = ncyc*s->tcyc;
		    for (k = 0; k<s->nchan; k++)
		    {
					phase0 = freq_phase[ncyc][k] + (tc-t0)/freq_period[ncyc][k];
	    		phase = (phase0 - floor(phase0));
	    		temp = phase*s->nbin;
	    		index = (int)(temp+0.5)>(int)temp?(int)temp+1:(int)temp;
	    		//index = (int)(temp);
	    		//s->prof[index] += s->fdata[i][s->nchan*j + k];
	    		prof[i][index] += s->fdata[i][s->nchan*j + k];
		    }
			}
		}

    for (i = 0; i<s->nsub; i++)
    {
			for (j = 0; j<s->nbin; j++)
			{
				s->prof[j] += prof[i][j];
			}
		}

    //for (i=0; i<s->nbin; i++)
    //{
	  //  printf ("%f\n", s->prof[i]);
    //}

    for (i=0; i<ncyc; i++)
    {
			free(freq_phase[i]);
			free(freq_period[i]);
		}

    free(freq_phase);
    free(freq_period);
}

