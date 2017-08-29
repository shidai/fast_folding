#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "T2toolkit.h"
#include "tempo2pred.h"

void get_PSRFITS_subint(float *fdata, fitsfile *fp, int isub, int nbit, int nchan, int nsblk)
{
    int ii, status = 0, anynull;
    int colnum;

    int FITS_typecode = 11;
    //long cur_subint = 1L;
    long cur_subint = isub;

    int bits_per_sample = nbit;
    int num_channels = nchan;
    int num_pols = 1;
    int samples_per_spectra = nchan;
    int spectra_per_subint = nsblk;
    int samples_per_subint = nsblk*nchan;
    int numtoread = samples_per_subint;

    unsigned char *ctmp, *cdata;

    // The following allows us to read byte-packed data
    if (bits_per_sample < 8) {
        numtoread = (samples_per_subint * bits_per_sample) / 8;
        //ctmp = gen_bvect(numtoread);
    }
    // or 16-bit data that is listed as being bytes
    if (bits_per_sample == 16 && FITS_typecode == 11)
        numtoread = samples_per_subint * 2;

    ctmp = (unsigned char *)malloc(numtoread*sizeof(unsigned char));
    cdata = (unsigned char *)malloc((8/bits_per_sample)*numtoread*sizeof(unsigned char));

    // Read the weights, offsets, and scales if required
    //if (s->apply_weight)
    //    fits_read_col(s->fitsfiles[cur_file], TFLOAT, s->dat_wts_col, cur_subint, 1L,
    //                  s->num_channels, 0, weights, &anynull, &status);
    //if (s->apply_offset)
    //    fits_read_col(s->fitsfiles[cur_file], TFLOAT, s->dat_offs_col, cur_subint,
    //                  1L, s->num_channels * s->num_polns, 0, offsets, &anynull,
    //                  &status);
    //if (s->apply_scale)
    //    fits_read_col(s->fitsfiles[cur_file], TFLOAT, s->dat_scl_col, cur_subint, 1L,
    //                  s->num_channels * s->num_polns, 0, scales, &anynull, &status);

    // Now actually read the subint into the temporary buffer
    fits_get_colnum(fp, 0, "DATA", &colnum, &status);
    fits_read_col(fp, FITS_typecode, colnum, cur_subint, 1L, numtoread, 0, ctmp, &anynull, &status);
    //printf ("First element: %u\n", ctmp[512]);

    if (status) {
        fprintf(stderr, "Error!:  Problem reading record from PSRFITS data file\n"
               "subint = %ld.  FITS status = %d.  Exiting.\n",
                cur_subint, status);
        exit(1);
    }
    // The following converts that byte-packed data into bytes
    if (bits_per_sample == 4) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numtoread,cdata,ctmp)
#endif
        for (ii = 0; ii < numtoread; ii++) {
            const unsigned char uctmp = ctmp[ii];
            const int jj = 2 * ii;
            cdata[jj] = uctmp >> 4;
            cdata[jj + 1] = uctmp & 0x0F;
    	    //printf ("First element: %u %u\n", cdata[jj], cdata[jj+1]);
        }
    } else if (bits_per_sample == 2) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numtoread,cdata,ctmp)
#endif
        for (ii = 0; ii < numtoread; ii++) {
            const unsigned char uctmp = ctmp[ii];
            const int jj = 4 * ii;
            cdata[jj] = ((uctmp >> 0x06) & 0x03);
            cdata[jj + 1] = ((uctmp >> 0x04) & 0x03);
            cdata[jj + 2] = ((uctmp >> 0x02) & 0x03);
            cdata[jj + 3] = (uctmp & 0x03);
        }
    } else if (bits_per_sample == 1) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(numtoread,cdata,ctmp)
#endif
        for (ii = 0; ii < numtoread; ii++) {
            const unsigned char uctmp = ctmp[ii];
            const int jj = 8 * ii;
            cdata[jj] = ((uctmp >> 0x07) & 0x01);
            cdata[jj + 1] = ((uctmp >> 0x06) & 0x01);
            cdata[jj + 2] = ((uctmp >> 0x05) & 0x01);
            cdata[jj + 3] = ((uctmp >> 0x04) & 0x01);
            cdata[jj + 4] = ((uctmp >> 0x03) & 0x01);
            cdata[jj + 5] = ((uctmp >> 0x02) & 0x01);
            cdata[jj + 6] = ((uctmp >> 0x01) & 0x01);
            cdata[jj + 7] = (uctmp & 0x01);
        }
    }

    //printf ("First element: %u %u\n", cdata[1025*2], cdata[1025*2+1]);
    //printf ("First element: %u\n", cdata[512*2]);

    if (bits_per_sample < 8)
        free(ctmp);

    //if (s->bits_per_sample == 1 && s->flip_bytes) {
    //    // Hack to flip each byte of data if needed
    //    for (ii = 0; ii < s->bytes_per_subint / 8; ii++) {
    //        int jj;
    //        const int offset = ii * 8;
    //        for (jj = 0; jj < 4; jj++) {
    //            unsigned char uctmp = cdata[offset + jj];
    //            cdata[offset + jj] = cdata[offset + 8 - 1 - jj];
    //            cdata[offset + 8 - 1 - jj] = uctmp;
    //        }
    //    }
    //}

    // Now convert all of the data into floats

    // The following allows us to work with single polns out of many
    // or to sum polarizations if required
    {
        if (bits_per_sample == 16) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(s,cdata,fdata,scales,offsets,weights)
#endif
            for (ii = 0; ii < spectra_per_subint; ii++) {
                int jj;
                float *fptr = fdata + ii * num_channels;
                const short *sptr = (short *) cdata + ii * samples_per_spectra;
                for (jj = 0; jj < num_channels; jj++)
                    //fptr[jj] = (((float) sptr[jj] - s->zero_offset) * scales[jj] + offsets[jj]) * weights[jj];
                    fptr[jj] = (float) sptr[jj];
            }
        } else if (bits_per_sample == 32) {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(s,cdata,fdata,scales,offsets,weights)
#endif
            for (ii = 0; ii < spectra_per_subint; ii++) {
                int jj;
                float *fptr = fdata + ii * num_channels;
                const float *ftptr = (float *) cdata + ii * samples_per_spectra;
                for (jj = 0; jj < num_channels; jj++)
                    //fptr[jj] = (((float) ftptr[jj] - s->zero_offset) * scales[jj] + offsets[jj]) * weights[jj];
                    fptr[jj] = (float) ftptr[jj];
            }
        } else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(s,cdata,fdata,scales,offsets,weights)
#endif
            for (ii = 0; ii < spectra_per_subint; ii++) {
                int jj;
                float *fptr = fdata + ii * num_channels;
                const unsigned char *cptr = cdata + ii * samples_per_spectra;
                for (jj = 0; jj < num_channels; jj++)
                    //fptr[jj] = (((float) cptr[jj] - s->zero_offset) * scales[jj] + offsets[jj]) * weights[jj];
                    fptr[jj] = (float) cptr[jj];
            }
        }
    }
}

void read_PSRFITS_files(char *fname, char pred_name)
{
    int i, nsub;
    int IMJD, SMJD, itmp, ii, status = 0;
    double OFFS, dtmp;
    long double MJDf;
    char ctmp[80], comment[120];

    int ret;
    T2Predictor pred;
    long double phase0, phase;

    long double imjd, smjd;
    long double mjd0;
    long double mjd;
    long double freq;

    double tsample;   // second
    //long double tsample = 0.000128;   // second

    int nchan;
    int nsblk;
    int nbit;

    float ftmp;
    long repeat, width;
    int colnum, anynull;

    if (ret=T2Predictor_Read(&pred, &pred_name))
    //if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
    {
	    printf("Error: unable to read predictor\n");
	    exit(1);
    }

    /*
    i = 0;
    while (i<512*2048)
    {
    	phase0 = T2Predictor_GetPhase(&pred,mjd0,freq);
	phase = (phase0 - floor(phase0));
    	printf ("Predicted phase: %lf\n", phase);
	i++;
    }
    */

    fitsfile *fp;	
    fits_open_file(&fp, fname, READONLY, &status);

    // Is the data in search mode?
    //fits_read_key(s->fitsfiles[ii], TSTRING, "OBS_MODE", ctmp, comment, &status);
    fits_read_key(fp, TSTRING, "OBS_MODE", ctmp, comment, &status);
    printf ("Mode: %s\n", ctmp);

    // Quick fix for Parkes DFB data (SRCH?  why????)...
    if (strcmp("SRCH", ctmp) == 0) {
        strncpy(ctmp, "SEARCH", 40);
    }
    if (strcmp(ctmp, "SEARCH")) {
        fprintf(stderr,
                "\nError!  File '%s' does not contain SEARCH-mode data!\n",
                fname);
        exit(1);
    }

    // Now get the stuff we need from the primary HDU header
    //fits_read_key(s->fitsfiles[ii], TSTRING, "TELESCOP", ctmp, comment, &status);
    fits_read_key(fp, TSTRING, "TELESCOP", ctmp, comment, &status);
    printf ("Telescope: %s\n", ctmp);

    fits_read_key(fp, TSTRING, "STT_IMJD", ctmp, comment, &status);
    imjd = atof(ctmp);
    fits_read_key(fp, TSTRING, "STT_SMJD", ctmp, comment, &status);
    smjd = atof(ctmp);
    mjd0 = imjd + smjd/86400.0L;
    printf ("Start MJD: %Lf %Lf %Lf\n", imjd, smjd, mjd0);

    // Now switch to the SUBINT HDU header
    fits_movnam_hdu(fp, BINARY_TBL, "SUBINT", 0, &status);
    fits_read_key(fp, TSTRING, "NAXIS2", ctmp, comment, &status);
    nsub = atoi(ctmp);

    // Read in nchan, nbit, tsamp
    fits_read_key(fp, TSTRING, "NCHAN", ctmp, comment, &status);
    nchan = atoi(ctmp);

    fits_read_key(fp, TSTRING, "NSBLK", ctmp, comment, &status);
    nsblk = atoi(ctmp);

    fits_read_key(fp, TSTRING, "NBITS", ctmp, comment, &status);
    nbit = atoi(ctmp);

    fits_read_key(fp, TSTRING, "TBIN", ctmp, comment, &status);
    tsample = atof(ctmp);
    //printf("test here %s\n", ctmp);

    printf ("NSUB: %d; NCHAN: %d; NSBLK: %d; NBITS: %d; TSAMP: %lf\n", nsub, nchan, nsblk, nbit, tsample);

    // Observing frequencies
    double *freqs;
    freqs = (double *)malloc(sizeof(double)*nchan);

    fits_get_colnum(fp, 0, "DAT_FREQ", &colnum, &status);
    if (status == COL_NOT_FOUND) {
        printf("Warning!:  Can't find the channel freq column!\n");
        status = 0;     // Reset status
    } else {
        fits_read_col(fp, TDOUBLE, colnum, 1L, 1L, nchan, 0, freqs, &anynull, &status);
        //printf("Frequency: %lf\n", freqs[511]);
    }

    int nbin = 128;
    float *prof;
    prof = (float *)malloc(sizeof(float)*nbin);
    for (i=0; i<nbin; i++)
    {
	    prof[i] = 0.0;
    }

    // Data
    float *fdata;
    int j,k;
    int index;
    float temp;
    fdata = (float *)malloc(sizeof(float)*nsblk*nchan);

    // Now pull stuff from the other columns
    float tc = 0.0; 
    float t0 = 0.0; 
    int ncyc = 0;
    long double *freq_phase;
    freq_phase = (long double *)malloc(sizeof(long double)*nchan);
    long double *freq_period;
    freq_period = (long double *)malloc(sizeof(long double)*nchan);
    int jj;

    for (i=1; i<=nsub; i++)
    //for (i=1; i<=100; i++)
    {
	get_PSRFITS_subint(fdata, fp, i, nbit, nchan, nsblk);
	for (j = 0; j<nsblk; j++)
	//for (j = 0; j<1; j++)
	{
		tc = (nsblk*(i-1) + j)*tsample;
		mjd = mjd0 + (nsblk*(i-1) + j)*tsample/86400.0L;

		if (tc >= ncyc)
		{
			printf ("test %f %d %d %f\n", tc, ncyc, nsblk, tsample);
			t0 = (nsblk*(i-1) + j)*tsample;
			for (jj = 0; jj<nchan; jj++)
			{
				freq = (long double)freqs[jj];
    				freq_phase[jj] = T2Predictor_GetPhase(&pred,mjd,freq);
    				freq_period[jj] = 1.0/T2Predictor_GetFrequency(&pred,mjd,freq);   // second
    				//freq_period[jj] = 1.0/218.81184385;   // second
			}
			//ncyc++;
			ncyc += 2;
		}

		for (k = 0; k<nchan; k++)
		//for (k = 128; k<136; k++)
		//k = 128;
		{
			phase0 = freq_phase[k] + (tc-t0)/freq_period[k];
			phase = (phase0 - floor(phase0));
    			//printf ("Predicted phase: %.10Lf\n", phase);
    			//printf ("Predicted phase: %.10Lf %Lf %.10Lf\n", phase, freq, mjd);

			temp = phase*nbin;
			index = (int)(temp+0.5)>(int)temp?(int)temp+1:(int)temp;
    			//printf ("Predicted index: %d\n", index);

			prof[index] += fdata[nchan*j+k];
		}
	}
	//printf ("%d\n", i);

    }

    for (i=0; i<nbin; i++)
    {
	    printf ("%f\n", prof[i]);
    }

    free(freqs);
    free(fdata);
    free(freq_phase);
    free(freq_period);
    free(prof);
}

