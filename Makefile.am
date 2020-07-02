bin_PROGRAMS = ptmfit

ptmfit_SOURCES = \
	computepoly.c \
	computepoly.h \
	LinearSystem.cpp \
	LinearSystem.h \
	main.cpp \
	nrutil.c \
	nrutil.h \
	RGBImage.h \
	svd.c \
	svd.h \
	writeptm.c \
	writeptm.h 

AM_CPPFLAGS = @VIPS_CFLAGS@ @VIPS_INCLUDES@
AM_LDFLAGS = @LDFLAGS@ 
LDADD = @VIPS_CFLAGS@ @VIPS_LIBS@


