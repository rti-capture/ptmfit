/* compute polynomials
 *
 * 2/2/17
 * 	- rewrite as a vips8 class
 */

/*
#define DEBUG 
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdlib.h>

#include <vips/vips.h>
#include <vips/debug.h>

#include "computepoly.h"

typedef struct _ComputePolys {
	VipsOperation parent_instance;

	/* Args.
	 */
	VipsArrayImage *in;
	VipsImage *out;
	VipsImage *M;

	VipsImage **arr;
	int n;

} ComputePolys;

typedef VipsOperationClass ComputePolysClass;

G_DEFINE_TYPE( ComputePolys, compute_polys, VIPS_TYPE_OPERATION );

/* Per-thread state.
 */
typedef struct {
	VipsRegion **ir;
	VipsPel * restrict * restrict p;
	double * restrict R;	
} ComputePolysSeq;

/* Free a sequence value.
 */
static int
compute_polys_stop( void *vseq, void *a, void *b )
{
	ComputePolysSeq *seq = (ComputePolysSeq *) vseq;
	ComputePolys *polys = (ComputePolys *) b;

	int i;

	if( seq->ir )
		for( i = 0; i < polys->n; i++ ) 
			VIPS_UNREF( seq->ir[i] );

	return( 0 );
}

/* Make a sequence value.
 */
static void *
compute_polys_start( VipsImage *out, void *a, void *b )
{
	ComputePolys *polys = (ComputePolys *) b;
	ComputePolysSeq *seq;
	int i;

	if( !(seq = VIPS_NEW( out, ComputePolysSeq )) )
		return( NULL );
	seq->ir = NULL;
	seq->p = NULL;
	seq->R = NULL;

	/* Attach regions and arrays.
	 */
	seq->ir = VIPS_ARRAY( out, polys->n + 1, VipsRegion * );
	seq->p = VIPS_ARRAY( out, polys->n + 1, VipsPel * );
	seq->R = VIPS_ARRAY( out, polys->n, double );
	if( !seq->ir || 
		!seq->p || 
		!seq->R ) {
		compute_polys_stop( seq, a, b );
		return( NULL );
	}

	for( i = 0; i < polys->n; i++ )
		if( !(seq->ir[i] = vips_region_new( polys->arr[i] )) ) {
			compute_polys_stop( seq, a, b );
			return( NULL );
		}
	seq->ir[i] = NULL;

	return( (void *) seq );
}

static int
compute_polys_gen( VipsRegion *or, void *vseq, void *a, void *b, gboolean *stop )
{
	ComputePolysSeq *seq = (ComputePolysSeq *) vseq;
	ComputePolys *polys = (ComputePolys *) b;
	VipsRect *r = &or->valid;

	int x, y, i, j;

	/* 8.5 added this faster thing.
	 */
#if VIPS_MAJOR_VERSION > 8 || \
	(VIPS_MAJOR_VERSION == 8 && VIPS_MINOR_VERSION >= 5)
	if( vips_reorder_prepare_many( or->im, seq->ir, r ) )
#else
	if( vips_region_prepare_many( seq->ir, r ) )
#endif
		return( -1 );

	VIPS_GATE_START( "compute_polys_gen: work" ); 

	for( y = 0; y < r->height; y++ ) {
		float * restrict q = (float * restrict) 
			VIPS_REGION_ADDR( or, r->left, r->top + y );

		for( i = 0; i < polys->n; i++ )
			seq->p[i] = VIPS_REGION_ADDR( seq->ir[i], 
				r->left, r->top + y );

		for( x = 0; x < r->width; x++ ) {
			double maxy;
			double adoty[3];
			double ydoty;

#ifdef DEBUG
			printf( "pix %d, %d: ", x + r->left, y + r->top );
			for( i = 0; i < polys->n; i++ ) {
				printf( "[" );
				for( j = 0; j < 3; j++ )
					printf( "%3d ", seq->p[i][j] );
				printf( "] " );
			}
			printf( "\n" );
#endif /*DEBUG*/

			/* Find normalised luminance for the input images. 
			 */

			maxy = 0;
			for( i = 0; i < polys->n; i++ ) {
				seq->R[i] = 0.2125 * seq->p[i][0] + 
					0.7154 * seq->p[i][1] +
					0.0721 * seq->p[i][2];
				seq->R[i] /= 255.0;

				if( seq->R[i] > maxy )
					maxy = seq->R[i];
			}

			if( maxy != 0 )
				for( i = 0; i < polys->n; i++ ) 
					seq->R[i] /= maxy;

			/* Find average colour component of inputs.
			 */

			for( j = 0; j < 3; j++ )
				adoty[j] = 0.0;
			ydoty = 0.0;
			for( i = 0; i < polys->n; i++ ) {
				for( j = 0; j < 3; j++ )
					adoty[j] += 
						seq->p[i][j] / 255.0 * seq->R[i];
				ydoty += seq->R[i] * seq->R[i];
			}

			for( i = 0; i < polys->n; i++ ) 
				seq->p[i] += 3;

			/* Write normalised RGB to the first three bands.
			 */

			for( j = 0; j < 3; j++ ) {
				double av;

				if( ydoty != 0 ) {
					av = adoty[j] / ydoty;

					if( av < 0.0 )
						av = 0.0;
					else if( av > 1.0 )
						av = 1.0;
				}
				else
					av = 0.0;

				q[j] = av * 255.0;
			}

#ifdef DEBUG
			printf( "RGB: " );
			for( j = 0; j < 3; j++ )
				printf( "%g ", q[j] );
			printf( "\n" );
#endif /*DEBUG*/

			q += 3;
			
			/* Put the luminance signal through the matrix to get
			 * the poly coefficients.
			 */

			g_assert( polys->M->Xsize == polys->n );
			g_assert( polys->M->Ysize == polys->out->Bands - 3 );

			for( j = 0; j < polys->M->Ysize; j++ ) {
				double sum;
				double * restrict coeff = 
					VIPS_MATRIX( polys->M, 0, 0 );
				double * restrict M = 
					coeff + j * polys->M->Xsize;

				sum = 0.0; 
				for( i = 0; i < polys->M->Xsize; i++ )
					sum += seq->R[i] * M[i];

				q[j] = 255.0 * sum;
			}

#ifdef DEBUG
			printf( "Lum poly: " );
			for( j = 0; j < polys->M->Ysize; j++ )
				printf( "%g ", q[j] );
			printf( "\n" );
#endif /*DEBUG*/

			q += polys->M->Ysize; 
		}
	}

	VIPS_GATE_STOP( "vips_bandary_gen: work" ); 

	return( 0 );
}

static int
compute_polys_build( VipsObject *object )
{
	ComputePolys *polys = (ComputePolys *) object;

	int i;

	if( VIPS_OBJECT_CLASS( compute_polys_parent_class )->build( object ) )
		return( -1 );

	polys->arr = vips_array_image_get( polys->in, &polys->n );

	if( polys->n < 1 ) {
		vips_error( "compute_polys", "%s", _( "zero input images!" ) );
		return( -1 );
	}
	if( polys->n != polys->M->Xsize ) {
		vips_error( "compute_polys", "%s", _( "M width != n images" ) );
		return( -1 );
	}

	for( i = 0; i < polys->n; i++ ) 
		if( vips_check_uncoded( "compute_polys", polys->arr[0] ) ||
			vips_check_size_same( "compute_polys", 
				polys->arr[0], polys->arr[i] ) ||
			vips_check_bands( "compute_polys", polys->arr[i], 3 ) ||
			vips_check_format( "compute_polys", 
				polys->arr[i], VIPS_FORMAT_UCHAR ) )
			return( -1 );

	g_object_set( object, "out", vips_image_new(), NULL ); 

	if( vips_image_pipeline_array( polys->out, 
		VIPS_DEMAND_STYLE_FATSTRIP, polys->arr ) )
		return( -1 );

	polys->out->BandFmt = VIPS_FORMAT_FLOAT;
	polys->out->Bands = 3 + polys->M->Ysize;

	if( vips_image_generate( polys->out,
		compute_polys_start, compute_polys_gen, compute_polys_stop, 
		polys->arr, polys ) )
		return( -1 );

	return( 0 );
}

static void
compute_polys_class_init( ComputePolysClass *class )
{
	GObjectClass *gobject_class = G_OBJECT_CLASS( class );
	VipsObjectClass *vobject_class = VIPS_OBJECT_CLASS( class );

	VIPS_DEBUG_MSG( "compute_polys_class_init\n" );

	gobject_class->set_property = vips_object_set_property;
	gobject_class->get_property = vips_object_get_property;

	vobject_class->nickname = "compute_polys";
	vobject_class->description = _( "compute a set of polynomials" );
	vobject_class->build = compute_polys_build;

	VIPS_ARG_BOXED( class, "in", 0, 
		_( "Input" ), 
		_( "Array of input images" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( ComputePolys, in ),
		VIPS_TYPE_ARRAY_IMAGE );

	VIPS_ARG_IMAGE( class, "out", 1, 
		_( "Output" ), 
		_( "Output image" ),
		VIPS_ARGUMENT_REQUIRED_OUTPUT, 
		G_STRUCT_OFFSET( ComputePolys, out ) );

	VIPS_ARG_IMAGE( class, "M", 2, 
		_( "M" ), 
		_( "Coefficient array" ),
		VIPS_ARGUMENT_REQUIRED_INPUT,
		G_STRUCT_OFFSET( ComputePolys, M ) );


}

static void
compute_polys_init( ComputePolys *polys )
{
}
