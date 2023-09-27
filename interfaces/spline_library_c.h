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

/**
 * @brief Create Spline Library object
 * 
 * @param num_pts Number of data points
 * @param pts Pointer of data points in SL_Vector3 type
 * @return SL_pSpline Pointer of spline object
 */
SL_pSpline SL_createSpline(int num_pts, const SL_Vector3 *pts);

/**
 * @brief Destroy Spline Library object
 * 
 * @param spline Pointer of spline object
 */
void SL_destroySpline(SL_pSpline spline);

/**
 * @brief Spline ranges in [0, Tmax], get Tmax 
 * 
 * @param spline Pointer of spline object
 * @return double Tmax
 */
double SL_getMaxT(const SL_pSpline spline);

/**
 * @brief Calculate points on the spline with given knot in [0, Tmax]
 * 
 * @param spline Pointer of spline object
 * @param knot Given knot in [0, Tmax]
 * @return SL_Vector3 Sampled points
 */
SL_Vector3 SL_getPosition(const SL_pSpline spline, double knot);

/**
 * @brief Calculate points with derivatives on the spline with given knot in [0, Tmax]
 * 
 * @param spline Pointer of spline object
 * @param knot Given knot in [0, Tmax]
 * @param der1 First-order derivatives in XYZ
 * @param der2 Second-order derivatives in XYZ
 * @return SL_Vector3 
 */
SL_Vector3 SL_getPositionDerivative(const SL_pSpline spline, double knot, SL_Vector3 *der1, SL_Vector3 *der2);

/**
 * @brief Calculate total arclength of the spline object
 * 
 * @param spline Pointer of spline object
 * @return double Arclength value
 */
double SL_getTotalArclength(const SL_pSpline spline);

/**
 * @brief Calculate accumulative arclengths at sampled knots
 * 
 * @param spline Pointer of spline object
 * @param num_pts Number of sampled knots
 * @param knots Pointer of sampled knots
 * @param arclengths Pointer of calculated arclengths
 */
void SL_getAccumulateArclengths(const SL_pSpline spline, int num_pts, double *knots, double *arclengths);

/**
 * @brief Calculate knot on spline with given arclength
 * 
 * @param spline Pointer of spline object
 * @param arclength Given arclength value
 * @return double Knot value on the spline
 */
double SL_arclengthToKnot(const SL_pSpline spline, double arclength);

/**
 * @brief Calculate knot_end so that arclength(knot_start, knot_end) ~= delta_arclength
 * 
 * @param spline Pointer of spline object
 * @param knot_start Value of knot_start
 * @param delta_arclength Value of delta_arclength
 * @return double Value of knot_end
 */
double SL_deltaArclengthToKnot(const SL_pSpline spline, double knot_start, double delta_arclength);

#ifdef __cplusplus
}
#endif

#endif // SPLINE_LIBRARY_C_H