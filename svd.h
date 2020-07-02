#ifndef SVD_H
#define SVD_H

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

#include <stdlib.h>
#include <stdio.h>

#include "nrutil.h"

void svdcmp(double **a, int m, int n, double w[], double **v);
double **MatrixStyleMult(double ** A, int m, int n, double **B, int o, int p);
double *MatrixStyleVectorMult(double **A,int m, int n, double *V, int o);
void PrintMatrix(FILE *file, double **M, int m, int n);
double **MatrixStyleTranspose(double ** A, int m, int n);

#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif /*SVD_H*/
