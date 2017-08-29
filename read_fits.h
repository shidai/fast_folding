#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "T2toolkit.h"
#include "tempo2pred.h"

void get_PSRFITS_subint(float *fdata, fitsfile *fp, int isub, int nbit, int nchan, int nsblk);

void read_PSRFITS_files(char *fname, char *pred_name);

