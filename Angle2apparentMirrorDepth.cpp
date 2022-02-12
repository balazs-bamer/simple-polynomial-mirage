#include "Angle2apparentMirrorDepth.h"

#include <functional>
#include <vector>
#include <cmath>

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
}

double Angle2apparentMirrorDepth::getTempRiseAtHeight(double const aHeight) const {
  auto height = 100.0 * aHeight;
  return csTempAmbient * ::exp(mB * ::pow(height, 1.0 - mDelta)) - csTempAmbient;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double Angle2apparentMirrorDepth::getHeightAtTempRise(double const aTempRise) const {
  return ::exp(::log(::log((aTempRise + csTempAmbient) / csTempAmbient) / mB) / (1.0 - mDelta)) / 100.0;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double Angle2apparentMirrorDepth::getRefractionAtTempRise(double const aTempRise) const {
  auto tempKelvin = csTempAmbient + aTempRise;
  auto tempCelsius = tempKelvin - csCelsius2kelvin;
  auto r = 1.0 + 7.86e-4 * csAtmosphericPressureKpa / tempKelvin - 1.5e-11 * csRelativeHumidityPercent * (tempCelsius * tempCelsius + 160.0);
  return r;
}

void Angle2apparentMirrorDepth::initReflection() {
  mCriticalInclination = binarySearch(csAlmostVertical, csAlmostHorizontal, csEpsilon, [this](double const aInclination) {
    return getReflectionDepth(aInclination) ? -1.0 : 1.0;
  }) + csEpsilon;
  std::vector<double> inclinations;
  std::vector<double> depths;
  auto increment = (csAlmostHorizontal - mCriticalInclination) / (csInclDepthProfilePointCount - 1u);
  auto inclination = mCriticalInclination;
  for(uint32_t i = 0u; i < csInclDepthProfilePointCount; ++i) {
    inclinations.push_back(inclination);
    depths.push_back(getReflectionDepth(inclination).value());
    inclination += increment;
  }
  mInclinationProfile = std::move(std::make_unique<PolynomApprox>(inclinations, depths, csInclinationProfileDegree));
}

std::optional<double> Angle2apparentMirrorDepth::getReflectionDepth(double const aInclination0) const {
  double currentHeight = csInitialHeight;
  double currentTemp = getTempRiseAtHeight(currentHeight);
  double currentRefractionIndex = getRefractionAtTempRise(currentTemp);
  double sinInclination0 = ::sin(aInclination0);
  double currentSinInclination = sinInclination0;
  double horizDisp = 0.0;

  std::optional<double> result;
  while(currentHeight > csMinimalHeight) {
    double nextTemp = currentTemp + static_cast<long double>(csLayerDeltaTemp);
    double nextHeight = getHeightAtTempRise(nextTemp);
    double nextRefractionIndex = getRefractionAtHeight(nextHeight);
    double factor = currentRefractionIndex / nextRefractionIndex;
    double nextSinInclination = currentSinInclination * factor;
    if(nextSinInclination >= static_cast<double>(1.0)) { // The problem here is how much sinInclanation is greater than 1. Optimal would be to have it just touch 1, so the last temperature step can't be uniform.
      auto const criticalTemp = binarySearch(currentTemp, nextTemp, csEpsilon, [this, currentRefractionIndex, currentSinInclination](double const aTemp) {
        double height = getHeightAtTempRise(aTemp);
        double refractionIndex = getRefractionAtHeight(height);
        double factor = currentRefractionIndex / refractionIndex;
        double sinInclination = currentSinInclination * factor;
        return (sinInclination > static_cast<double>(1.0) ? 1.0 : -1.0);
      });
      double criticalHeight = getHeightAtTempRise(criticalTemp);
      horizDisp += static_cast<double>(currentHeight - criticalHeight) * currentSinInclination / ::sqrt(static_cast<double>(1.0) - currentSinInclination * currentSinInclination);
      auto r = csInitialHeight - horizDisp * ::sqrt(1.0 - sinInclination0 * sinInclination0) / sinInclination0;
      result = r;
      break;
    }
    else {
      horizDisp += static_cast<double>(currentHeight - nextHeight) * currentSinInclination / ::sqrt(static_cast<double>(1.0) - currentSinInclination * currentSinInclination);
      currentTemp = nextTemp;
      currentHeight = nextHeight;
      currentRefractionIndex = nextRefractionIndex;
      currentSinInclination = nextSinInclination;
    }
  }
  return result;
}
