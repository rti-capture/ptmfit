#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <string.h>

#include <vips/vips8>

#include "computepoly.h"
#include "LinearSystem.h"

#define INTERNAL_FEATURES_ENABLE // for internal release only

#define VERSION_NUMBER 1.01

char lpfile[STRSIZE];
char fname[STRSIZE];

Basis_e base = QUADRATIC_BIVARIATE;

bool outputfilegiven = false;
bool cache = false;

int crop_left = 0;
int crop_top = 0;
int crop_width = 1000;
int crop_height = 1000;

void Usage(char *argv0)
{
	printf("%s usage:\n", argv0);

	printf("  -i filename\n");
	printf("    Full filename for lp file specifing input files and light positions. \n\n");

	printf("  -PTM <path>/<file.ptm>\n");
	printf("  -o <path>/<file.ptm>\n");
	printf("     Output file name.\n\n");

	printf("  -BIVARIATE | -UNIVARIATE\n");
	printf("  -b basis  (basis either 0 for biquadratic or 1 for univariate\n");
	printf("     Calculate a least squares fit of one independent variable:  -UNIVARIATE\n");
	printf("                                   or two independent variables: -BIVARIATE\n");
	printf("     (Default: BIVARIATE)\n\n");

	printf("  -crop LEFT TOP WIDTH HEIGHT\n");
	printf("    Only process part of input frames, crops are 10 x percent\n\n");
	printf("  -cache\n");
	printf("    Cache calculated coefficients (needs lots of mem)\n\n");

	printf("  -version\n");
	printf("    Prints software version\n\n");

	printf("  -h\n");
	printf("     List command line options\n");
}


void Process_parameters(int argc, char **argv)
{
	int i;

	/*
	 * Defaults:
	 */

	strcpy(lpfile,"");

	base = QUADRATIC_BIVARIATE;

	for( i = 1; i < argc; i++ )
	{
		if( strcmp( argv[i], "-i" ) == 0)
		{
			if (strlen(argv[i+1]) > STRSIZE)
			{
				printf("Sorry, input filename too long, please try again\n");
				exit(-1);
			}
			strcpy(lpfile, argv[++i]);
		} else

		if( strcmp( argv[i], "-PTM" ) == 0 || strcmp( argv[i], "-o" ) == 0)
		{
			strcpy(fname, argv[++i]);
			outputfilegiven = true;
		} else

		if( strcmp( argv[i], "-BIVARIATE") == 0)
		{
			base = QUADRATIC_BIVARIATE;
		} else

		if( strcmp( argv[i], "-UNIVARIATE") == 0)
		{
			base = QUADRATIC_UNIVARIATE;
		} else

		if( strcmp( argv[i], "-b") == 0)
		{
			if (atoi(argv[++i]) == 0)
				base = QUADRATIC_BIVARIATE;
			else
				base = QUADRATIC_UNIVARIATE;
		} else

		if( strcmp( argv[i], "-rgb") == 0)
		{
			// ignore this one, some other fitter
		} else

		if( strcmp( argv[i], "-lrgb") == 0)
		{
			// ignore this one, some other fitter
		} else

		if( strcmp( argv[i], "-cache") == 0)
		{
			cache = true;
		} else

		if( strcmp( argv[i], "-crop" ) == 0)
		{
			if( argc - i < 5 ) {
				printf("too few arguments for crop\n");
				exit(-1);
			}
			crop_left = atoi( argv[i + 1] );
			crop_top = atoi( argv[i + 2] );
			crop_width = atoi( argv[i + 3] );
			crop_height = atoi( argv[i + 4] );
			i += 4;
		} else

		if( strcmp( argv[i], "-version") == 0)
		{
			printf("PTM Fitter version %3.2f\n",VERSION_NUMBER);
		} else

		if( strcmp( argv[i], "-h") == 0 || strcmp( argv[i], "/?") == 0)
		{
			Usage(argv[0]);
			exit(0);
		} else

		{
			printf("Error in parameter list: %s\n", argv[i]);
			Usage(argv[0]);
			exit(0);
		}
	}

	if (strlen(lpfile) == 0)
	{

		printf("Error: -i option must be specified.\n");
		Usage(argv[0]);
		exit(0);
	}

	return;
}


int
main(int argc, char* argv[])
{
	printf("Fast Polynomial Texture Map (PTM) Fitter\n");
	printf("based on Hewlett-Packard Code from 2001. All rights reserved.\n");
	
	//Gather user input

	int stat;

	if (VIPS_INIT(argv[0]))
		return(-1);

	compute_polys_get_type();

	if(argc > 1)
		// we have command line arguments, do not need to poll user
	{
		Process_parameters(argc, argv);
	}
	else
	{							 //prompt user

		int polybase;

		std::cout << "Enter filename for lp file (light positions file)" << "\n";
		std::cin >> lpfile;

		std::cout << "Enter basis(polynomials)" <<  "\n";
		std::cout << "0 = quadratic in two variables, 1 = quadratic in 1 variable" << "\n";
		std::cin >> polybase;
		if(polybase != 0 && polybase !=1)
		{
			std::cout << " Wrong input; assume the default (quadratic)" <<"\n";
			polybase =0;
		}

		if(polybase ==0)
			base = QUADRATIC_BIVARIATE;
		else
			base = QUADRATIC_UNIVARIATE;
	}

	LinearSystem lin(base, cache, crop_left, crop_top, crop_width, crop_height);

	stat = lin.FitPTM(lpfile);

	if(stat == -1)
	{
		std::cout << "Error in fitting PTM" << "\n";
		std::cout.flush();
		return (-1 );
	}

	if(outputfilegiven == false)
	{
		std::cout << "Enter filename for output PTM file: " << "\n";
		std::cin >> fname;
	}

	std::cout << "Writing file to : " << fname << "\n";
	std::cout.flush();

	lin.WriteFileVersion1_2(fname);

	return( 0 );
}
