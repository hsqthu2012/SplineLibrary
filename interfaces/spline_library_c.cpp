#include "spline_library_c.h"
#include "spline.h"
#include "vector.h"
#include "splines/cubic_hermite_spline.h"
#include "utils/arclength.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef Vector<3, double> Vector3b;

SL_pSpline SL_createSpline(int num_pts, const SL_Vector3 *pts)
{
  std::vector<Vector3b> splinePoints(num_pts);
  for (size_t i = 0; i < num_pts; i++)
  {
    splinePoints[i] = Vector3b({pts[i].x, pts[i].y, pts[i].z});
  }

  CubicHermiteSpline<Vector3b, double> *spline = nullptr;
  try {
    spline = new CubicHermiteSpline<Vector3b, double>(splinePoints, 0.5);
  } catch (std::bad_alloc& ba) {
    throw std::runtime_error("Failed to allocate memory for spline");
  }

  return (SL_pSpline) spline;
}

void SL_destroySpline(SL_pSpline spline)
{
  delete static_cast<CubicHermiteSpline<Vector3b, double> *>(spline);
}

double SL_getMaxT(const SL_pSpline spline)
{
  auto real_spline = static_cast<CubicHermiteSpline<Vector3b, double> *>(spline);

  return real_spline->getMaxT();
}

SL_Vector3 SL_getPosition(const SL_pSpline spline, double knot)
{
  auto real_spline = static_cast<CubicHermiteSpline<Vector3b, double> *>(spline);

  Vector3b pos = real_spline->getPosition(knot);

  SL_Vector3 vec = {pos[0], pos[1], pos[2]};

  return vec;
}

SL_Vector3 SL_getPositionDerivative(const SL_pSpline spline, double knot, SL_Vector3 *der1, SL_Vector3 *der2)
{
  auto real_spline = static_cast<CubicHermiteSpline<Vector3b, double> *>(spline);

  auto ptc = real_spline->getCurvature(knot);

  SL_Vector3 pos = {ptc.position[0], ptc.position[1], ptc.position[2]};

  der1->x = ptc.tangent[0];
  der1->y = ptc.tangent[1];
  der1->z = ptc.tangent[2];

  der2->x = ptc.curvature[0];
  der2->y = ptc.curvature[1];
  der2->z = ptc.curvature[2];

  return pos;
}

double SL_getTotalArclength(const SL_pSpline spline)
{
  auto real_spline = static_cast<CubicHermiteSpline<Vector3b, double> *>(spline);

  return real_spline->totalLength();
}

void SL_getAccumulateArclengths(const SL_pSpline spline, double *arclengths)
{
  auto real_spline = static_cast<CubicHermiteSpline<Vector3b, double> *>(spline);

  size_t segs = real_spline->segmentCount();

  for (size_t i = 0; i < segs; i++)
  {
    double knot = real_spline->segmentForT(i);
    arclengths[i] = real_spline->arcLength(0.0, knot);
  }
}

double SL_arclengthToKnot(const SL_pSpline spline, double arclength)
{
  auto real_spline = static_cast<CubicHermiteSpline<Vector3b, double> *>(spline);

  return ArcLength::solveLength(*real_spline, 0.0, arclength);
}

#ifdef __cplusplus
}
#endif