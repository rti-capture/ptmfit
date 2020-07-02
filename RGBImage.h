// RGBImage.h: interface for the RGBImage class.

#ifndef RGBIMAGE_H
#define RGBIMAGE_H

/* Confusingly, this is actually used as a C++ Class, see LinearSystem.cpp.
 */

struct pixelRGB
{
  unsigned char r, g, b;
};

struct RGB_Image
{
  int xsize, ysize;		/* size of image */

  vips::VImage im;

  float x;
  float y;
  float z;			//position of light vector or view vector or H vector or whatever vector

  char *filename;
};

#endif /*RGBIMAGE_H*/
