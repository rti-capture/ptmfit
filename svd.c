#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <math.h>

#include "nrutil.h"

#define SIGN(a,b) ( (b> 0.0) ? fabs(a) : -fabs(a))

static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a), maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a), iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

static double dsqrarg;
#define SQR(a) ((dsqrarg = (a) ) == 0.0 ? 0.0 : dsqrarg * dsqrarg)

static double
pythag (double a, double b)
{
  double absa, absb;

  absa = fabs (a);
  absb = fabs (b);

  if (absa > absb)
    return absa * sqrt (1.0 + SQR (absb / absa));
  else
    return (absb == 0.0 ? 0.0 : absb * sqrt (1.0 + SQR (absa / absb)));
}

void
svdcmp (double **a, int m, int n, double w[], double **v)
{
  int flag, i, its, j, jj, k, l, nm;
  double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

  rv1 = dvector (1, n);
  g = scale = anorm = 0.0;

  for (i = 1; i <= n; i++)
    {
      l = i + 1;
      rv1[i] = scale * g;
      g = s = scale = 0.0;

      if (i <= m)
	{
	  for (k = i; k <= m; k++)
	    scale += fabs (a[k][i]);

	  if (scale)
	    {
	      for (k = i; k <= m; k++)
		{
		  a[k][i] /= scale;
		  s += a[k][i] * a[k][i];
		}

	      f = a[i][i];
	      g = -SIGN (sqrt (s), f);
	      h = f * g - s;
	      a[i][i] = f - g;

	      for (j = l; j <= n; j++)
		{
		  for (s = 0.0, k = i; k <= m; k++)
		    s += a[k][i] * a[k][j];
		  f = s / h;
		  for (k = i; k <= m; k++)
		    a[k][j] += f * a[k][i];
		}

	      for (k = i; k <= m; k++)
		a[k][i] *= scale;
	    }
	}

      w[i] = scale * g;
      g = s = scale = 0.0;

      if (i <= m && i != n)
	{
	  for (k = l; k <= n; k++)
	    scale += fabs (a[i][k]);

	  if (scale)
	    {
	      for (k = l; k <= n; k++)
		{
		  a[i][k] /= scale;
		  s += a[i][k] * a[i][k];
		}

	      f = a[i][l];
	      g = -SIGN (sqrt (s), f);
	      h = f * g - s;
	      a[i][l] = f - g;

	      for (k = l; k <= n; k++)
		rv1[k] = a[i][k] / h;

	      for (j = l; j <= m; j++)
		{
		  for (s = 0.0, k = l; k <= n; k++)
		    s += a[j][k] * a[i][k];
		  for (k = l; k <= n; k++)
		    a[j][k] += s * rv1[k];
		}

	      for (k = l; k <= n; k++)
		a[i][k] *= scale;
	    }
	}

      anorm = FMAX (anorm, (fabs (w[i]) + fabs (rv1[i])));
    }

  for (i = n; i >= 1; i--)
    {
      if (i < n)
	{
	  if (fabs (g) > 0.0)
	    {
	      for (j = l; j <= n; j++)
		v[j][i] = (a[i][j] / a[i][l]) / g;

	      for (j = l; j <= n; j++)
		{
		  for (s = 0.0, k = l; k <= n; k++)
		    s += a[i][k] * v[k][j];

		  for (k = l; k <= n; k++)
		    v[k][j] += s * v[k][i];
		}
	    }

	  for (j = l; j <= n; j++)
	    v[i][j] = v[j][i] = 0.0;
	}
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
    }

  for (i = IMIN (m, n); i >= 1; i--)
    {
      l = i + 1;
      g = w[i];

      for (j = l; j <= n; j++)
	a[i][j] = 0.0;

      if (fabs (g) > 0.0)
	{
	  g = 1.0 / g;

	  for (j = l; j <= n; j++)
	    {
	      for (s = 0.0, k = l; k <= m; k++)
		s += a[k][i] * a[k][j];
	      f = (s / a[i][i]) * g;
	      for (k = i; k <= m; k++)
		a[k][j] += f * a[k][i];
	    }

	  for (j = i; j <= m; j++)
	    a[j][i] *= g;
	}
      else
	for (j = i; j <= m; j++)
	  a[j][i] = 0.0;

      ++a[i][i];
    }

  for (k = n; k >= 1; k--)
    {
      for (its = 1; its <= 30; its++)
	{
	  flag = 1;

	  for (l = k; l >= 1; l--)
	    {
	      nm = l - 1;

	      if ((double) (fabs (rv1[l]) + anorm) == anorm)
		{
		  flag = 0;
		  break;
		}

	      if ((double) (fabs (w[nm]) + anorm) == anorm)
		break;
	    }

	  if (flag)
	    {
	      c = 0.0;
	      s = 1.0;

	      for (i = l; i <= k; i++)
		{
		  f = s * rv1[i];
		  rv1[i] = c * rv1[i];

		  if ((double) (fabs (f) + anorm) == anorm)
		    break;

		  g = w[i];
		  h = pythag (f, g);
		  w[i] = h;
		  h = 1.0 / h;
		  c = g * h;
		  s = -f * h;

		  for (j = 1; j <= m; j++)
		    {
		      y = a[j][nm];
		      z = a[j][i];
		      a[j][nm] = y * c + z * s;
		      a[j][i] = z * c - y * s;
		    }
		}
	    }

	  z = w[k];

	  if (l == k)
	    {
	      if (z < 0.0)
		{
		  w[k] = -z;

		  for (j = 1; j <= n; j++)
		    v[j][k] = -v[j][k];
		}
	      break;
	    }

	  if (its == 30)
	    nrerror ("no convergence in 30 svdcmp iterations");

	  x = w[l];
	  nm = k - 1;
	  y = w[nm];
	  g = rv1[nm];
	  h = rv1[k];
	  f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
	  g = pythag (f, 1.0);
	  f = ((x - z) * (x + z) + h * ((y / (f + SIGN (g, f))) - h)) / x;
	  c = s = 1.0;

	  for (j = l; j <= nm; j++)
	    {
	      i = j + 1;
	      g = rv1[i];
	      y = w[i];
	      h = s * g;
	      g = c * g;
	      z = pythag (f, h);
	      rv1[j] = z;
	      c = f / z;
	      s = h / z;
	      f = x * c + g * s;
	      g = g * c - x * s;
	      h = y * s;
	      y *= c;

	      for (jj = 1; jj <= n; jj++)
		{
		  x = v[jj][j];
		  z = v[jj][i];
		  v[jj][j] = x * c + z * s;
		  v[jj][i] = z * c - x * s;
		}

	      z = pythag (f, h);
	      w[j] = z;

	      if (z)
		{
		  z = 1.0 / z;
		  c = f * z;
		  s = h * z;
		}

	      f = c * g + s * y;
	      x = c * y - s * g;

	      for (jj = 1; jj <= m; jj++)
		{
		  y = a[jj][j];
		  z = a[jj][i];
		  a[jj][j] = y * c + z * s;
		  a[jj][i] = z * c - y * s;
		}
	    }

	  rv1[l] = 0.0;
	  rv1[k] = f;
	  w[k] = x;
	}
    }

  free_dvector (rv1, 1, n);
}

void
PrintMatrix (FILE * file, double **M, int m, int n)
{
  int i, j;

  for (i = 1; i <= m; i++)
    {
      for (j = 1; j <= n; j++)
	fprintf (file, "%lf ", M[i][j]);
      fprintf (file, "\n");
    }
}

double **
MatrixStyleTranspose (double **A, int m, int n)
{
  int i, j;
  double **AT = dmatrix (1, n, 1, m);

  for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++)
      {
	AT[j][i] = A[i][j];
      }

  return AT;
}

double **
MatrixStyleMult (double **A, int m, int n, double **B, int o, int p)
{
  int i, j, k;
  double **M;
  double sum;

  if (n != o)
    {
      fprintf (stderr, "Illegal matrix multiply\n");
      return (0);
    }

  M = dmatrix (1, m, 1, p);

  for (i = 1; i <= m; i++)
    for (j = 1; j <= p; j++)
      {
	sum = 0.0;

	for (k = 1; k <= n; k++)
	  {
	    sum += A[i][k] * B[k][j];
	  }
	M[i][j] = sum;
      }

  return M;
}

double *
MatrixStyleVectorMult (double **A, int m, int n, double *V, int o)
{
  int i, j;
  double sum;
  double *Vec;

  if (n != o)
    {
      fprintf (stderr, "Illegal matrix multiply\n");
      return (0);
    }

  Vec = (double *) malloc (m * sizeof (double));

  for (i = 0; i < m; i++)
    {
      sum = 0.0;

      for (j = 0; j < n; j++)
	sum += A[i + 1][j + 1] * V[j + 1];
      Vec[i] = sum;
    }

  return Vec;
}


/* (C) Copr. 1986-92 Numerical Recipes Software 42,. */
