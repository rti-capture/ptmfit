# ptmfit program with libvips backend
DEVELOPENT OF THE 2014 version by JC
needs: 
new svd.c
some path smarts
easy cross compile for win10 and mac
----------

image are uncompressed to disc temps in /tmp, or in $TMPDIR, if set

memory use is 9 floats for every pixel in the output image



----------------------------


Instructions for new fitting code
This fitting update should make it more user friendly:

By executing main.exe the user gets prompted for several inputs:

1)Input directory where samples are located eg. D:\sample06
2)Input prefix of the actual images files(ppm). If the files are called
ex-001.pppm, ex-002.ppm and so forth one should input “ex”. See remarks for
more info.
3)Select format RGB or LRGB. Code works for both, RGB fits 3 polynomials and
LRGB fits 1 polynomial for luminance and average rgb value.
4)Weighting scheme: right now only the choice 0 (NONE) is supported.
5)UNIVARIATE or BIVARIATE: Typically one selects  BIVARIATE (0). Some of the
animations you might have seen (the Golden gate sequence) uses univariate.
Here still 3 coordinates are specified but only x is used. 
6)Number of images to be used: -1 uses all, a number  N>0 has the effect that
only the first N images are used for fitting. The last M-N images are
discarded (M = total number of images)
7)Output file: Specify name including the entire path. The program appends the
exetnsion .ppm


REMARKS:
The fitter expects a file ex.lp in the directory where the files ex-001.ppm
etc are located. This file contains th 3d positions. I have included an
example directory sample06 which illustrates this (please look at 06.lp)
