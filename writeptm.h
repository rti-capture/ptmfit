#ifndef WRITEPTM_H
#define WRITEPTM_H

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

#include <vips/vips.h>

/* Keep i18n stuff happy.
 */
#define _(S) (S)

int writeptm( VipsImage *in, const char *filename, double *scale, int *bias );

#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif /*WRITEPTM_H*/
