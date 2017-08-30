#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "T2toolkit.h"
#include "tempo2pred.h"

typedef struct search_mode{
	long double imjd, smjd, mjd0;
	int nsub, nchan, nsblk, nbit;
	int nbin;
	double tsample;

	double *freqs;
	double *prof;
	float *fdata;
}search_mode;

void get_PSRFITS_subint(float *fdata, fitsfile *fp, int isub, int nbit, int nchan, int nsblk);

void read_PSRFITS_files(search_mode *s, char *fname, char *pred_name);

void demalloc (search_mode *s);
