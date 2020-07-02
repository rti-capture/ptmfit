// LinearSystem.h: interface for the LinearSystem class.
//
//////////////////////////////////////////////////////////////////////

#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "RGBImage.h"

enum Basis_e {QUADRATIC_BIVARIATE, QUADRATIC_UNIVARIATE};

#define STRSIZE 256

class LinearSystem
{
	public:
		LinearSystem(Basis_e b, bool c = false, 
				int crop_left = 0, int crop_top = 0, 
				int crop_width = 1000, int crop_height = 1000);
		int FitPTM(char *lpfile);
		virtual ~LinearSystem();
		void WriteFileVersion1_2(char *filename);

	private:
		int InitFiles(char *lpfile);
		int LoadFiles();
		int BuildMatrix(double **  &M);
		int ComputePolynomials(double **M);
		void ComputeScaleAndBias();
		void ComputeQuantizedRGBPolynomials();
		void ComputeQuantizedLumPolynomials();

		int Images_m;			 // number of images read in

		RGB_Image *Samples_m;	 // pointer to actual input images also contains the position
		// should be in 3D though.
		Basis_e basis_m;

		// keep a mem cache
		bool cache;

		// only work on this part of the input images
		// crop expressed as 10 x percent
		int crop_left_m;
		int crop_top_m;
		int crop_width_m;
		int crop_height_m;

		// inverse matrix, as passed to compute_polys()
		vips::VImage vipsM;

		// huge array of computed coefficients
		vips::VImage coeffs;

		double scale[6];
		int bias[6];

		RGB_Image *qredhigh ;
		RGB_Image *qredlow;
		RGB_Image *qgreenhigh;
		RGB_Image *qgreenlow;
		RGB_Image *qbluehigh;
		RGB_Image *qbluelow;

		RGB_Image *qlumlow;
		RGB_Image *qlumhigh;
		RGB_Image *colors;
};

#endif /*LINEARSYSTEM_H*/
