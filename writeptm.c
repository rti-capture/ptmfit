/* write ptm file format
 */

/*
#define DEBUG 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdlib.h>
#include <stdio.h>

#include <vips/vips.h>

#include "writeptm.h"

/* What we track during a PTM write.
 */
typedef struct {
	VipsImage *in;
	char *name;
	double *scale;
	int *bias;

	VipsPel *line;
	FILE *fp;

	/* We have to write the file backwards, since PTM files have the origin
	 * at byte 0 and (almost) all other file formats have the top left
	 * corner at byte 0.
	 */

	long coeff_start;
	long rgb_start;
} Write;

static void
write_destroy( Write *write )
{
	VIPS_FREEF( fclose, write->fp );
	VIPS_FREE( write->name );
	VIPS_FREE( write->line );

	vips_free( write );
}

static Write *
write_new( VipsImage *in, const char *name, double *scale, int *bias )
{
	Write *write;

	if( !(write = VIPS_NEW( NULL, Write )) )
		return( NULL );

	write->in = in;
	write->name = vips_strdup( NULL, name );
	write->scale = scale;
	write->bias = bias;
	write->line = VIPS_ARRAY( NULL, in->Xsize * 6, VipsPel );
	write->fp = NULL;

        if( !(write->fp = fopen( name, "wb" )) ) {
		write_destroy( write );
		vips_error( "writeptm", 
			"unable to open \"%s\" for writing", name);
		return( NULL );
	}
	if( !write->name || !write->line ) {
		write_destroy( write );
		return( NULL );
	}
	
        return( write );
}

static int
write_coeff_block( Write *write, VipsRegion *region, VipsRect *area )
{
	double * restrict scale = write->scale;
	int * restrict bias = write->bias;

	int x, y, i;

	for( y = 0; y < area->height; y++ ) {
		float * restrict p;
		VipsPel * restrict q;
		long offset;

		p = (float * restrict) 
			VIPS_REGION_ADDR( region, 0, area->top + y );
		q = write->line;

		for( x = 0; x < area->width; x++ ) {
			for( i = 5; i >= 0; i-- ) {
				float v = p[i + 3];
				VipsPel ch;

				ch = (VipsPel) (v / scale[i] + bias[i] + 0.5);

				q[5 - i] = ch;
			}

			p += 9;
			q += 6;
		}

		offset = write->coeff_start;
		offset += (write->in->Ysize - 1 - (area->top + y)) * 
			area->width * 6;
		fseek( write->fp, offset, SEEK_SET ); 

		if( !fwrite( write->line, area->width * 6, 1, write->fp ) ) {
			vips_error( "writeptm", 
				"%s", _( "write error ... disc full?" ) );
			return( -1 );
		}

#ifdef DEBUG
		printf( "written coeff scanline %d to %zd ",
			area->top + y, offset );
		for( i = 0; i < 5; i++ )
			printf( "%d ", write->line[i] );
		printf( "\n" ); 
#endif /*DEBUG*/
	}

	return( 0 );
}

static int
write_rgb_block( Write *write, VipsRegion *region, VipsRect *area )
{
	int x, y, i;

	for( y = 0; y < area->height; y++ ) {
		float * restrict p;
		VipsPel * restrict q;
		long offset;

		p = (float * restrict) 
			VIPS_REGION_ADDR( region, 0, area->top + y );
		q = write->line;

		for( x = 0; x < area->width; x++ ) {
			for( i = 0; i < 3; i++ ) {
				float v = p[i];
				VipsPel ch;

				ch = (VipsPel) (v + 0.5);

				q[i] = ch;
			}

			p += 9;
			q += 3;
		}

		offset = write->rgb_start;
		offset += (write->in->Ysize - 1 - (area->top + y)) * 
			area->width * 3;
		fseek( write->fp, offset, SEEK_SET ); 

		if( !fwrite( write->line, area->width * 3, 1, write->fp ) ) {
			vips_error( "writeptm", 
				"%s", _( "write error ... disc full?" ) );
			return( -1 );
		}

		/*
		printf( "written rgb scanline %d to %zd ",
			area->top + y, offset );
		for( i = 0; i < 5; i++ )
			printf( "%d ", write->line[i] );
		printf( "\n" ); 
		 */
	}

	return( 0 );
}

static int
write_block( VipsRegion *region, VipsRect *area, void *a )
{
	Write *write = (Write *) a;

	if( write_coeff_block( write, region, area ) ||
		write_rgb_block( write, region, area ) )
		return( -1 ); 

	return( 0 );
}

static int
write_ptm( Write *write )
{
	VipsImage *in = write->in;

	int i;

	fprintf( write->fp, "PTM_1.2\n" );
	fprintf( write->fp, "%s\n", "PTM_FORMAT_LRGB" );

	fprintf( write->fp, "%i\n", in->Xsize );
	fprintf( write->fp, "%i\n", in->Ysize );

	for( i = 0; i <= 5; i++ )
		fprintf( write->fp, "%f ", write->scale[5 - i] );
	fprintf( write->fp, "\n" );

	for( i = 0; i <= 5; i++ )
		fprintf( write->fp, "%i ", write->bias[5 - i] );
	fprintf( write->fp, "\n" );

	/* Set up file layout. The coeff area is 6 bytes per pixel.
	 */
	write->coeff_start = ftell( write->fp ); 
	write->rgb_start = write->coeff_start + 
		((size_t) write->in->Xsize) * write->in->Ysize * 6;

	if( vips_sink_disc( write->in, write_block, write ) )
		return( -1 );

	return( 0 );
}

int
writeptm( VipsImage *in, const char *filename, double *scale, int *bias )
{
	Write *write;

	if( vips_check_format( "writeptm", in, VIPS_FORMAT_FLOAT ) || 
		vips_check_bands( "writeptm", in, 9 ) || 
		vips_check_uncoded( "writeptm", in ) )
		return( -1 );

	if( !(write = write_new( in, filename, scale, bias )) )
		return( -1 );

	if( write_ptm( write ) ) {
		write_destroy( write );
		return( -1 );
	}
	write_destroy( write );

	return( 0 );
}
