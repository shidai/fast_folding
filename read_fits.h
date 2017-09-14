#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

typedef struct search_mode{
	long double imjd, smjd, mjd0;
	int nsub, nchan, nsblk, nbit;
	int nbin;
	double tsample;
	float tcyc;

	double *freqs;
	double *prof;
	float **fdata;
}search_mode;

void get_PSRFITS_subint(float *fdata, fitsfile *fp, int isub, int nbit, int nchan, int nsblk);

void read_PSRFITS_files(search_mode *s, char *fname);

void demalloc (search_mode *s);
