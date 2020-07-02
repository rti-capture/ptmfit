// LinearSystem.cpp: implementation of the LinearSystem class.
//
//////////////////////////////////////////////////////////////////////

/*
#define DEBUG 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <cstdlib>
#include <iostream>
#include <math.h>
#include <string.h>

#include <vips/vips8>

#include "LinearSystem.h"
#include "svd.h"
#include "nrutil.h"
#include "computepoly.h"
#include "writeptm.h"

using namespace vips;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

LinearSystem::LinearSystem(Basis_e b, bool c, int crop_left, int crop_top, 
	int crop_width, int crop_height)
{
  basis_m = b;
  cache = c;
  crop_left_m = crop_left;
  crop_top_m = crop_top;
  crop_width_m = crop_width;
  crop_height_m = crop_height;
  colors = 0;
  Samples_m = 0;
}

void
readLine (FILE * fp, char *buf, int limit)
{
  int c, i;

  for (i = 0; i < limit - 1 && (c = getc (fp)) != EOF && c != '\n'; i++)
    buf[i] = c;

  if (c == '\n')
    {
      buf[i] = c;
      i++;
    }

  if (c != '\n' && c != EOF)
    {
      printf ("File reading warning:\n");
      printf ("attempted to read line and stopped before newline\n");
      printf ("Please use a shorter path or check for newlines\n");
    }

  buf[i] = '\0';
}

// load the files from the filenames that were found in the .lp file
int
LinearSystem::LoadFiles ()
{
  int stat = 1;
  int i;

  /* Output images.  high images has r,g,b = A,B,C repectively */
  /* low images have r,g,b = D,E,F */

  printf ("Reading images:");

  for (i = 0; i < Images_m; i++)
    {
      char buf[5];

      sprintf (buf, "%3i", i + 1);
      if (i)
	printf ("\b\b\b");
      printf ("%s", buf);
      fflush (stdout);

      char *basename = g_path_get_basename(Samples_m[i].filename); 

      VImage im = VImage::new_from_file(basename, 
	VImage::option()->set("access", VIPS_ACCESS_SEQUENTIAL));

      g_free( basename ); 

      int left = crop_left_m / 1000.0 * im.width();
      int top = crop_top_m / 1000.0 * im.height();
      int width = crop_width_m / 1000.0 * im.width();
      int height = crop_height_m / 1000.0 * im.height();

      im = im.extract_area (left, top, width, height);

      Samples_m[i].im = im;
      Samples_m[i].xsize = im.width();
      Samples_m[i].ysize = im.height();

      if (Samples_m[0].xsize != im.width()
	  || Samples_m[0].ysize != im.height())
	{
	  fprintf (stderr, "images differ in size\n");
	  return -1;
	}
    }

  printf ("\n");

  return stat;
}

int
LinearSystem::InitFiles (char *lpfile)
{
  char filedesc[250];
  float vecmag;
  int readargs;

  FILE *infofp;

  infofp = fopen (lpfile, "r");
  if (infofp == NULL)
    {
      fprintf (stderr, "Light Position file: %s, not found.\n", lpfile);
      return -1;
    }

  if (fscanf (infofp, "%i\n", &Images_m) != 1)
    {
      fprintf (stderr, "Light Position file: first line wrong\n");
      return -1;
    }

  Samples_m = new RGB_Image[Images_m];

  for (int i = 0; i < Images_m; i++)
    {
      char inputLine[STRSIZE];

      readLine (infofp, inputLine, STRSIZE);
      readargs = sscanf (inputLine, "%s %f %f %f", (char *) &filedesc,
			 &(Samples_m[i].x), &(Samples_m[i].y),
			 &(Samples_m[i].z));

      if (readargs != 4)
	{
	  printf
	    ("Caution, couldn't find required values: filename  x  y  z in lp file.\n");
	  printf ("For line: %s\n", inputLine);
	}

      vecmag = (float) sqrt (Samples_m[i].x * Samples_m[i].x +
			     Samples_m[i].y * Samples_m[i].y +
			     Samples_m[i].z * Samples_m[i].z);

      if (vecmag > 0)
	{
	  // normalize
	  vecmag = 1.0f / vecmag;
	  Samples_m[i].x *= vecmag;
	  Samples_m[i].y *= vecmag;
	  Samples_m[i].z *= vecmag;
	}
      else
	{
	  printf ("Caution, length zero vector found in lp file: %f %f %f\n",
		  Samples_m[i].x, Samples_m[i].y, Samples_m[i].z);
	}

      Samples_m[i].filename = new char[strlen (filedesc) + 1];
      strcpy (Samples_m[i].filename, filedesc);
    }

  fclose (infofp);

  return 0;
}

int
LinearSystem::FitPTM (char *lpfile)
{
  int stat = 1;

  if (InitFiles (lpfile) == -1)
    return -1;

  if (LoadFiles () == -1)
    return -1;

  // here we do the set up for solving the PTM;
  // here we first set up a matrix and the right hand side
  // then we call the solver to obtain the result
  // Afterwards scale and bias have to be computed and the PTM file has to 
  // be written
  double **M;

  stat = BuildMatrix (M);
  if (stat == -1)
    return stat;

  stat = ComputePolynomials (M);
  if (stat == -1)
    return stat;

  std::vector<VImage> in;
  for (int i = 0; i < Images_m; i++)
    in.push_back(Samples_m[i].im);

  VImage::call("compute_polys", VImage::option()->
	set( "in", in )->
	set( "out", &coeffs )->
	set( "M", vipsM ) );

  /* With cache enabled, write to a huge memory buffer.
   */
  if( cache ) 
	  coeffs = coeffs.write(VImage::new_memory());

  ComputeScaleAndBias ();

  printf ("computation done!\n");

  // if we're not caching, we'll need to scan again for write
  if( !cache ) {
	  // we open our files in streaming mode, so after the pass where we 
	  // calculate scale/bias, we need to reopen for the ptm write and 
	  // regen coeffs

	  LoadFiles ();

	  std::vector<VImage> in2;
	  for (int i = 0; i < Images_m; i++)
	    in2.push_back(Samples_m[i].im);

	  VImage::call("compute_polys", VImage::option()->
		set( "in", in2 )->
		set( "out", &coeffs )->
		set( "M", vipsM ) );
  }

  int coldim;

  if (basis_m == QUADRATIC_BIVARIATE)
    coldim = 6;
  else if (basis_m == QUADRATIC_UNIVARIATE)
    coldim = 3;
  else
    coldim = 6;			// use QUADRATIC_BIVARIATE as default

  free_dmatrix (M, 1, Images_m, 1, coldim);

  return stat;
}

int
LinearSystem::BuildMatrix (double **&M)
{
  int coldim;
  int stat = 1;

  if (basis_m == QUADRATIC_BIVARIATE)
    coldim = 6;
  else if (basis_m == QUADRATIC_UNIVARIATE)
    coldim = 3;
  else
    coldim = 6;			// use QUADRATIC_BIVARIATE as default

  if (Images_m < coldim)
    {
      printf ("Error: Not enough samples for fitting a PTM\n");
      if (coldim == 6)
	printf ("A Bivariate PTM requires 6 or more images\n");
      else if (coldim == 3)
	printf ("A Univariate PTM requires 3 or more images\n");
      stat = -1;
    }

  M = dmatrix (1, Images_m, 1, coldim);

  if (basis_m == QUADRATIC_BIVARIATE)
    {
      for (int k = 1; k <= Images_m; k++)
	{
	  M[k][1] = 1.0;
	  M[k][2] = Samples_m[k - 1].y;
	  M[k][3] = Samples_m[k - 1].x;
	  M[k][4] = M[k][2] * M[k][3];
	  M[k][5] = M[k][2] * M[k][2];
	  M[k][6] = M[k][3] * M[k][3];
	}
    }
  else if (basis_m == QUADRATIC_UNIVARIATE)
    {
      for (int k = 1; k <= Images_m; k++)
	{
	  M[k][1] = 1.0;
	  M[k][2] = Samples_m[k - 1].x;
	  M[k][3] = M[k][2] * M[k][2];
	}
    }
  else
    {
      printf ("Basis not implemented\n");
      stat = -1;
    }

  return stat;
}

int
LinearSystem::ComputePolynomials (double **M)
{
  // R is the right hand side and hence is crucial

  // FILE *fp = fopen("Matrix","w");
  // fclose(fp);

  int i, j, k, l;
  int coldim;

  if (basis_m == QUADRATIC_BIVARIATE)
    coldim = 6;
  else if (basis_m == QUADRATIC_UNIVARIATE)
    coldim = 3;
  else
    coldim = 6;			// use QUADRATIC_BIVARIATE as default

  double *Diag = dvector (1, Images_m);
  double *R = dvector (1, Images_m);
  double **V = dmatrix (1, Images_m, 1, Images_m);

  svdcmp (M, Images_m, coldim, Diag, V);

  for (k = 1; k <= coldim; k++)
    if (fabs (Diag[k]) <= 1.0e-10)
      {
	printf ("System can not be solved, not enough info to "
		"compute coefficients!\n");
	printf ("Most likely cause: sample locations are redundant; "
		"e.g. are colinear\n");
	free_dvector (Diag, 1, Images_m);
	free_dvector (R, 1, Images_m);
	free_dmatrix (V, 1, Images_m, 1, Images_m);
	free_dmatrix (M, 1, Images_m, 1, coldim);

	return -1;
      }

  double **UT = MatrixStyleTranspose (M, Images_m, coldim);

  for (k = 1; k <= coldim; k++)
    for (l = 1; l <= Images_m; l++)
      UT[k][l] = UT[k][l] / Diag[k];

  double **InverseMatrix =
    MatrixStyleMult (V, coldim, coldim, UT, coldim, Images_m);

  // make the matrix image
  vipsM = VImage::new_matrix(Images_m, coldim);
  for (j = 0; j < coldim; j++)
    for (i = 0; i < Images_m; i++)
	*VIPS_MATRIX( vipsM.get_image(), i, j ) = InverseMatrix[j + 1][i + 1];

  free_dvector (Diag, 1, Images_m);
  free_dvector (R, 1, Images_m);
  free_dmatrix (V, 1, Images_m, 1, Images_m);
  free_dmatrix (UT, 1, coldim, 1, Images_m);
  free_dmatrix (InverseMatrix, 1, coldim, 1, Images_m);

  return 1;
}

void
LinearSystem::ComputeScaleAndBias ()
{
  // First compute the minimum and maximum of each coefficient for each channel
  int i;
  int basedim;

  printf ("computing scale and bias ... \n");

  if (basis_m == QUADRATIC_BIVARIATE)
    basedim = 6;
  else if (basis_m == QUADRATIC_UNIVARIATE)
    basedim = 3;
  else
    basedim = 6;

  VImage stats = coeffs.stats ();

#ifdef DEBUG
  std::cout << stats;
#endif /*DEBUG*/

  // the first three channels are RGB, the subsequent 3 or 6 are the poly
  // coeffs

  double lummin[6], lummax[6];

  for (i = 0; i < basedim; i++)
    {
      lummin[i] = *VIPS_MATRIX( stats.get_image(), 0, i + 4);
      lummax[i] = *VIPS_MATRIX( stats.get_image(), 1, i + 4);
    }

#ifdef DEBUG
  printf( "min: " );
  for( i = 0; i < basedim; i++ )
  	printf( "%g ", lummin[i] );
  printf( "\n" );

  printf( "max: " );
  for( i = 0; i < basedim; i++ )
  	printf( "%g ", lummax[i] );
  printf( "\n" );
#endif /*DEBUG*/

  int lumscale[6];

  for (i = 0; i < basedim; i++)
    {
      // check what the next highest 2 power is
      frexp (lummax[i] - lummin[i], &(lumscale[i]));
      //256 or 8 is what we can deal with, higher power requires scaling
      lumscale[i] = lumscale[i] - 8;
    }

  for (i = 0; i < basedim; i++)
    scale[i] = (float) pow (2, lumscale[i]);

  for (i = 0; i < basedim; i++)
    {
      if (lummin[i] < 0)
	bias[i] = (int) (0.0 - lummin[i] / scale[i] + 1);

      while (bias[i] > 255)
	{
	  scale[i] = scale[i] * 2;
	  bias[i] = (int) (0.0 - lummin[i] / scale[i] + 1);
	}
      if (lummin[i] > 0)
	{
	  while (lummax[i] / scale[i] > 255)
	    scale[i] = scale[i] * 2;

	  bias[i] = 0;
	}
    }

#ifdef DEBUG
  for (i = 0; i < basedim; i++)
    printf ("SCALE AND BIAS: %lf  %d\n", scale[i], bias[i]);
#endif /*DEBUG*/
}

void
LinearSystem::WriteFileVersion1_2 (char *fname)
{
	if( writeptm( coeffs.get_image (), fname, scale, bias ) )
	{
		std::cerr << "Error writing file\n"; 
	}
}

LinearSystem::~LinearSystem ()
{
  // clean up the samples as well
  if (Samples_m)
    {
      for (int i = 0; i < Images_m; i++)
	{
	  delete Samples_m[i].filename;
	}
    }
  delete[]Samples_m;
}
