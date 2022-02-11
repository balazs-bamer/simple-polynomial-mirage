#include "Angle2apparentMirrorDepth.h"

extern "C" {
#include <quadmath.h>
}
#include <functional>
#include <cmath>

#include <iostream> // TODO remove
#include <vector>

Angle2apparentMirrorDepth::Static::Static() {
  mHeightLimit = std::move(std::make_unique<PolynomApprox>(csTplate, csHeightLimit, csTempProfilePointCount, csTempProfileDegree));
  mB           = std::move(std::make_unique<PolynomApprox>(csTplate, csB,           csTempProfilePointCount, csTempProfileDegree));
  mDelta       = std::move(std::make_unique<PolynomApprox>(csTplate, csDelta,       csTempProfilePointCount, csTempProfileDegree));
}

Angle2apparentMirrorDepth::Angle2apparentMirrorDepth(double const aTempDiffSurface)
  : mTempDiffSurface(std::max(0.0, aTempDiffSurface)) {
  static Static tempCoeffProfiles;

  mHeightLimit = tempCoeffProfiles.mHeightLimit->eval(csTempAmbient + mTempDiffSurface);
  mB           = tempCoeffProfiles.mB->eval(csTempAmbient + mTempDiffSurface);
  mDelta       = tempCoeffProfiles.mDelta->eval(csTempAmbient + mTempDiffSurface);
  initReflection();

  std::cout << "Tplate = " << csTempAmbient + mTempDiffSurface;
  std::cout << "\nheightLimit = " << mHeightLimit;
  std::cout << "\nB = " << mB;
  std::cout << "\ndelta = " << mDelta << '\n';
}

double Angle2apparentMirrorDepth::getTempAtHeight(double const aHeight) const {
  auto height = 100.0 * aHeight;
//  auto delta = (height > mHeightLimit ? mDelta : csDeltaFallback);
  return csTempAmbient * ::exp(mB * ::pow(height, 1.0 - mDelta));
}

double Angle2apparentMirrorDepth::getHeightAtTemp(double const aTemp) const {
  return ::exp(::log(::log(aTemp / csTempAmbient) / mB) / (1.0 - mDelta)) / 100.0;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

std::vector<uint32_t> gIters;
bool                  gLog = false;

void Angle2apparentMirrorDepth::initReflection() {
  auto const minInclination = binarySearch(csAlmostVertical, csAlmostHorizontal, csEpsilon, [this](double const aInclination) {
    return getReflectionDepth(aInclination) ? -1.0f : 1.0f;
  }) + csEpsilon;
gIters.clear();
  std::vector<double> inclinations;
  std::vector<double> depths;
  auto increment = (csAlmostHorizontal - minInclination) / (csInclinationProfilePointCount - 1u);
  auto inclination = minInclination;
gLog = true;
  for(uint32_t i = 0u; i < csInclinationProfilePointCount; ++i) {
    inclinations.push_back(inclination);
    depths.push_back(getReflectionDepth(inclination).value());
    inclination += increment;
  }
std::cout << "====================\n==================\n==================\nincl = [";
for(auto d : inclinations) {
std::cout << d * 180.0 / cgPi << ", ";
}
std::cout << "\ndepth = [";
for(auto d : depths) {
std::cout << d << ", ";
}
std::cout << "\niter = [";
for(auto d : gIters) {
std::cout << d << ", ";
}
std::cout << "\n";
  mInclinationProfile = std::move(std::make_unique<PolynomApprox>(inclinations, depths, csInclinationProfileDegree));
std::cout << "incl error: " << mInclinationProfile->getRrmsError() << "\n";
}

std::optional<double> Angle2apparentMirrorDepth::getReflectionDepth(double const aInclination0) const {
  double currentHeight = csInitialHeight;
  double currentTemp = getTempAtHeight(currentHeight);
  __float128 currentRefractionIndex = getRefractionAtTemp(currentTemp);
  double sinInclination0 = ::sin(aInclination0);
  __float128 sinInclination = sinInclination0;
  __float128 horizDisp = 0.0f;

  // TODO check reflection calculation
  // TODO unit test for temp2height and inverse

  std::optional<double> result;
uint32_t iter = 0u;
std::vector<double> x;
std::vector<double> y;
x.push_back(horizDisp);
y.push_back(currentHeight);
  while(currentHeight > csMinimalHeight) {
    double nextTemp = currentTemp + static_cast<long double>(csLayerDeltaTemp);
    double nextHeight = getHeightAtTemp(nextTemp);
    horizDisp += static_cast<__float128>(currentHeight - nextHeight) * sinInclination / sqrtq(static_cast<__float128>(1.0) - sinInclination * sinInclination);
    __float128 nextRefractionIndex = getRefractionAtHeight(nextHeight);
    __float128 factor = currentRefractionIndex / nextRefractionIndex;
    sinInclination *= factor;
x.push_back(horizDisp);
y.push_back(nextHeight);
    if(sinInclination >= static_cast<__float128>(1.0)) { // The problem here is how much sinInclanation is greater than 1. Optimal would be to have it just touch 1, so the last temperature step canát be uniform.
      auto r = csInitialHeight - horizDisp * ::sqrt(1.0 - sinInclination0 * sinInclination0) / sinInclination0;
      result = r;
      break;
    }
    else {
      currentTemp = nextTemp;
      currentHeight = nextHeight;
      currentRefractionIndex = nextRefractionIndex;
    }
++iter;
//std::cout << iter << "   " << currentHeight << "   " << currentTemp << " - " << sinInclination << "   " << horizDisp << "\n";
  }
//std::cout << "\n";
if(result && gLog) {
gIters.push_back(iter);
/*int const degrees1m = aInclination0 * 180000000.0 / cgPi;
double const tangentLength = 1000.0;
std::cout << "x" << degrees1m << " = [";
for(auto d : x) {
  std::cout << d << ", ";
}
std::cout << "\ny" << degrees1m << " = [";
for(auto d : y) {
  std::cout << d << ", ";
}
  std::cout << "\nxt" << degrees1m << " = [" << -::sin(aInclination0) * tangentLength << ", " << ::sin(aInclination0) * tangentLength << "]\n";
  std::cout << "yt" << degrees1m << " = [" << csInitialHeight + ::cos(aInclination0) * tangentLength << ", " << csInitialHeight - ::cos(aInclination0) * tangentLength << "]\n";
*/}
  return result;
}
