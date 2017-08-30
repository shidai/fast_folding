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

    //float tc = 0.0; 
    //float t0 = 0.0; 
    int ncyc = 0;
    long double **freq_phase;
    long double **freq_period;
    long double *freq_phase_temp;
    long double *freq_period_temp;
    float **t0;
    float t0_temp;
    float **tc;

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

    freq_phase = (long double **)malloc(sizeof(long double *)*s->nsub);
    freq_period = (long double **)malloc(sizeof(long double *)*s->nsub);
    t0 = (float **)malloc(sizeof(float *)*s->nsub);
    tc = (float **)malloc(sizeof(float *)*s->nsub);
    for (i=0; i<s->nsub; i++)
    {
    	freq_phase[i] = (long double *)malloc(sizeof(long double)*s->nchan*s->nsblk);
    	freq_period[i] = (long double *)malloc(sizeof(long double)*s->nchan*s->nsblk);

	t0[i] = (float *)malloc(sizeof(float)*s->nsblk);
	tc[i] = (float *)malloc(sizeof(float)*s->nsblk);
    }

    freq_phase_temp = (long double *)malloc(sizeof(long double)*s->nchan);
    freq_period_temp = (long double *)malloc(sizeof(long double)*s->nchan);

    ///////////////////////////////////////////////

    for (i=0; i<100; i++)
    //for (i=0; i<s->nsub; i++)
    {
	    for (j = 0; j<s->nsblk; j++)
	    {
		    mjd = s->mjd0 + (s->nsblk*i + j)*s->tsample/86400.0L;
		    tc[i][j] = (s->nsblk*i + j)*s->tsample;
		    
		    if (tc[i][j] >= ncyc)
		    {
			    //t0[i][j] = (s->nsblk*i + j)*s->tsample;
			    t0_temp = (s->nsblk*i + j)*s->tsample;
			    for (k = 0; k<s->nchan; k++)
			    {
				    freq = (long double)s->freqs[k];
				    freq_phase_temp[k] = T2Predictor_GetPhase(&pred,mjd,freq);
				    freq_period_temp[k] = 1.0/T2Predictor_GetFrequency(&pred,mjd,freq);   // second
				    printf ("test %f %d\n", tc[i][j], ncyc);
			    }
			    ncyc += 2;
		    }
		    
		    t0[i][j] = t0_temp;
		    for (k = 0; k<s->nchan; k++)
		    {
			    freq_phase[i][j*s->nchan + k] = freq_phase_temp[k];
			    freq_period[i][j*s->nchan + k] = freq_period_temp[k];   // second
		    }
			    //printf ("test %f %d %f %Lf\n", tc[i][j*s->nchan + k], ncyc, t0[i][j*s->nchan + k], freq_period[i][j*s->nchan + k]);
	    }
    }
    //printf ("I am here\n");

    for (i=0; i<100; i++)
    //for (i=0; i<s->nsub; i++)
    {
	    printf ("%d\n", i);
	    for (j = 0; j<s->nsblk; j++)
	    {
		    for (k = 0; k<s->nchan; k++)
    		    {
			    printf ("I am here %d %d %d %Lf %Lf\n", i, j, k, freq_phase[i][j*s->nchan + k], freq_period[i][j*s->nchan + k]);
			    phase0 = freq_phase[i][j*s->nchan + k] + (tc[i][j]-t0[i][j])/freq_period[i][j*s->nchan + k];
			    //printf ("I am here here %d %Lf\n", k, phase0);
	    		    phase = (phase0 - floor(phase0));
			    //printf ("I am here here here %d %Lf\n", k, phase);
	    		    temp = phase*s->nbin;
	    		    index = (int)(temp+0.5)>(int)temp?(int)temp+1:(int)temp;
			    //printf ("I am here here here here %d %Lf %d\n", k, phase, index);
	    		    s->prof[index] += s->fdata[i][s->nchan*j + k];
		    }
	    }
    }
    /*
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
		*/

    for (i=0; i<s->nbin; i++)
    {
	    printf ("%f\n", s->prof[i]);
    }

    for (i=0; i<s->nsub; i++)
    {
    	free(freq_phase[i]);
    	free(freq_period[i]);
    	free(t0[i]);
    	free(tc[i]);
    }
    free(freq_phase);
    free(freq_period);
    free(freq_phase_temp);
    free(freq_period_temp);
    free(t0);
    free(tc);
}

