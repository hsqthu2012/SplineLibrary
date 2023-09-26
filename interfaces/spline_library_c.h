#ifndef SPLINE_LIBRARY_C_H
#define SPLINE_LIBRARY_C_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  double x;
  double y;
  double z;
} SL_Vector3;

typedef void* SL_pSpline;

SL_pSpline SL_createSpline(int num_pts, const SL_Vector3 *pts);

void SL_destroySpline(SL_pSpline spline);

double SL_getMaxT(const SL_pSpline spline);

SL_Vector3 SL_getPosition(const SL_pSpline spline, double knot);

SL_Vector3 SL_getPositionDerivative(const SL_pSpline spline, double knot, SL_Vector3 *der1, SL_Vector3 *der2);

double SL_getTotalArclength(const SL_pSpline spline);

void SL_getAccumulateArclengths(const SL_pSpline spline, int num_pts, double *knots, double *arclengths);

double SL_arclengthToKnot(const SL_pSpline spline, double arclength);

#ifdef __cplusplus
}
#endif

#endif // SPLINE_LIBRARY_C_H