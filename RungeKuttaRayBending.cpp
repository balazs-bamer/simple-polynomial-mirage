#include "RungeKuttaRayBending.h"
#include <cmath>


RungeKuttaRayBending::Result RungeKuttaRayBending::solve4xFlat(Vertex const &aStart, Vector const &aDir, double const aX) {
  typename Eikonal::Variables start;
  start[0u] = aStart(0u);
  start[1u] = aStart(1u);
  start[2u] = aStart(2u);
  auto slowness = mDiffEq.getSlowness(aStart(1u));  // from height
  start[3u] = aDir(0u) * slowness;
  start[4u] = aDir(1u) * slowness;
  start[5u] = aDir(2u) * slowness;
  auto solution = mSolver.solve(start,
      [aX](double const, typename Eikonal::Variables const& aY){ return aY[0] >= aX; },
    [this](typename Eikonal::Variables const& aYprev, typename Eikonal::Variables const& aYnow) { return decide2resetBigStep(aYprev, aYnow); });
  Result result;
  result.mValid = solution.mValid;
  result.mValue(0u) = solution.mValue[0u];
  result.mValue(1u) = solution.mValue[1u];
  result.mValue(2u) = solution.mValue[2u];
  return result;
}

RungeKuttaRayBending::Result RungeKuttaRayBending::solve4xRound(Vertex const &aStart, Vector const &aDir, double const aX) {
  typename Eikonal::Variables start;
  start[0u] = aStart(0u);
  start[1u] = aStart(1u) + mDiffEq.getEarthRadius();
  start[2u] = aStart(2u);
  auto slowness = mDiffEq.getSlowness(aStart(1u));  // from height
  start[3u] = aDir(0u) * slowness;
  start[4u] = aDir(1u) * slowness;
  start[5u] = aDir(2u) * slowness;
  auto solution = mSolver.solve(start,
      [aX](double const, typename Eikonal::Variables const& aY){ return aY[0] >= aX; },     // We now neglect the variation in perpendicular along the travelled distance.
    [this](typename Eikonal::Variables const& aYprev, typename Eikonal::Variables const& aYnow) { return decide2resetBigStep(aYprev, aYnow); });
  Result result;
  result.mValid = solution.mValid;
  result.mValue(0u) = solution.mValue[0u];
  result.mValue(1u) = solution.mValue[1u] - mDiffEq.getEarthRadius();
  result.mValue(2u) = solution.mValue[2u];
  return result;
}
