//#include <sys/types.h>
//#include <pwd.h>
//#include "presto.h"
//#include "mask.h"
//#include "psrfits.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include "fitsio.h"
#include "read_fits.h"
#include "fold.h"

int main(int argc, char *argv[])
{
	char in_name[128]; // input file name
	char pred_name[128]; // predictor name
	int i;

	search_mode *s;
	s = (search_mode *)malloc(sizeof(search_mode));

	for (i=0; i<argc; i++)
	{
		if (strcmp(argv[i], "-f") == 0)
		{
			strcpy(in_name, argv[++i]);
			printf ("File name: %s\n", in_name);
		}
		else if (strcmp(argv[i], "-pred") == 0)
		{
			strcpy(pred_name, argv[++i]);
			printf ("Predictor name: %s\n", pred_name);
		}
	}

	read_PSRFITS_files(s, in_name);

	fold_main(s, pred_name);

	demalloc(s);

	return 0;
}


