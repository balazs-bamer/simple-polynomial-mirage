#include "Angle2apparentMirrorDepth.h"

#include <functional>
#include <cmath>

#include <iostream> // TODO remove
#include <vector>

Angle2apparentMirrorDepth::Static::Static() {
  mHeightLimit = std::move(std::make_unique<PolynomApprox>(csTplate, csHeightLimit, csTempProfilePointCount, csTempProfileDegree));
  mB           = std::move(std::make_unique<PolynomApprox>(csTplate, csB,           csTempProfilePointCount, csTempProfileDegree));
  mDelta       = std::move(std::make_unique<PolynomApprox>(csTplate, csDelta,       csTempProfilePointCount, csTempProfileDegree));
/*std::cout << i << '\n';
for(int j = 0; j < 3; ++j) {
  std::cout << mStuff->at(i)->eval(csReferenceTempProfiles[0u][j]) << "   ";
}
std::cout << '\n';*/
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
  auto delta = (height > mHeightLimit ? mDelta : csDeltaFallback);
  return csTempAmbient * ::exp(mB * ::pow(height, 1.0 - delta));
}

double Angle2apparentMirrorDepth::getHeightAtTemp(double const aTemp) const {
  return ::exp(::log(::log(aTemp / csTempAmbient) / mB) / (1.0 - mDelta));  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

std::vector<uint32_t> gIters;

void Angle2apparentMirrorDepth::initReflection() {
  auto const minInclination = binarySearch(csAlmostVertical, csAlmostHorizontal, csEpsilon, [this](double const aInclination) {
    return getReflectionDepth(aInclination) ? -1.0f : 1.0f;
  }) + csEpsilon;
gIters.clear();
  std::vector<double> inclinations;
  std::vector<double> depths;
  auto increment = (csAlmostHorizontal - minInclination) / (csInclinationProfilePointCount - 1u);
  auto inclination = minInclination;
std::cout << "====================\n==================\n==================\nincl: ";
  for(uint32_t i = 0u; i < csInclinationProfilePointCount; ++i) {
    inclinations.push_back(inclination);
    depths.push_back(getReflectionDepth(inclination).value());
std::cout << inclination * 180.0 / cgPi << ", ";
    inclination += increment;
  }
std::cout << "\ndepth: ";
for(auto d : depths) {
std::cout << d << ", ";
}
std::cout << "\niter: ";
for(auto d : gIters) {
std::cout << d << ", ";
}
std::cout << "\n";
  mInclinationProfile = std::move(std::make_unique<PolynomApprox>(inclinations, depths, csInclinationProfileDegree));
std::cout << "incl error: " << mInclinationProfile->getRrmsError() << "\nincl coeff: ";
for(int i = 0; i < mInclinationProfile->size(); ++i)
  std::cout << mInclinationProfile->operator[](i) << ' ';
std::cout << '\n';
}

std::optional<double> Angle2apparentMirrorDepth::getReflectionDepth(double const aInclination0) const {
  double currentHeight = csInitialHeight;
  double currentTemp = getTempAtHeight(currentHeight);
  double currentRefractionIndex = getRefractionAtTemp(currentTemp);
  double sinInclination0 = ::sin(aInclination0);
  double sinInclination = sinInclination0;
  double horizDisp = 0.0f;

  // TODO check reflection calculation
  // TODO unit test for temp2height and inverse

uint32_t iter = 0u;
//std::cout << "--------------- " << aInclination0 * 180.0/cgPi << " ---------------- " << iter << "   " << currentHeight << "   " << currentTemp << " - " << sinInclination0 << '\n';
  std::optional<double> result;
  while(currentHeight > csMinimalHeight) {
    auto nextTemp = currentTemp + csLayerDeltaTemp;
    auto nextHeight = getHeightAtTemp(nextTemp);
    horizDisp += (currentHeight - nextHeight) * sinInclination / ::sqrt(1 - sinInclination * sinInclination);
    double nextRefractionIndex = getRefractionAtHeight(nextHeight);
    sinInclination *= currentRefractionIndex / nextRefractionIndex;
    if(sinInclination >= 1.0f) {
      result = horizDisp * ::sqrt(1 - sinInclination0 * sinInclination0) / sinInclination0;
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
gIters.push_back(iter);
  return result;
}
