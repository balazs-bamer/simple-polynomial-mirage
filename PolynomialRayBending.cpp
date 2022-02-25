#include "PolynomialRayBending.h"

#include <functional>
#include <vector>
#include <cmath>


PolynomialRayBending::Static::Static() {
  mHeightLimit = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csHeightLimit, std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
  mB           = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csB,           std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
  mDelta       = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csDelta,       std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
}

PolynomialRayBending::PolynomialRayBending(double const aTempDiffSurface)
  : mTempDiffSurface(std::max(0.0, aTempDiffSurface)) {
  static Static tempCoeffProfiles;

  mHeightLimit = tempCoeffProfiles.mHeightLimit->eval(csTempAmbient + mTempDiffSurface);
  mB           = tempCoeffProfiles.mB->eval(csTempAmbient + mTempDiffSurface);
  mDelta       = tempCoeffProfiles.mDelta->eval(csTempAmbient + mTempDiffSurface);
  initReflection();
}

double PolynomialRayBending::getTempRiseAtHeight(double const aHeight) const {
  auto height = 100.0 * aHeight;
  return csTempAmbient * ::exp(mB * ::pow(height, 1.0 - mDelta)) - csTempAmbient;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double PolynomialRayBending::getHeightAtTempRise(double const aTempRise) const {
  return ::exp(::log(::log((aTempRise + csTempAmbient) / csTempAmbient) / mB) / (1.0 - mDelta)) / 100.0;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double PolynomialRayBending::getRefractionAtTempRise(double const aTempRise) const {
  auto tempKelvin = csTempAmbient + aTempRise;
  auto tempCelsius = tempKelvin - csCelsius2kelvin;
  auto r = 1.0 + 7.86e-4 * csAtmosphericPressureKpa / tempKelvin - 1.5e-11 * csRelativeHumidityPercent * (tempCelsius * tempCelsius + 160.0);
  return r;
}

void PolynomialRayBending::initReflection() {
  mCriticalInclination = binarySearch(csAlmostVertical, csAlmostHorizontal, csEpsilon, [this](double const aInclination) {
    return getReflectionDispDepth(aInclination).mAsphalt ? 1.0 : -1.0;
  }) + csEpsilon;
  auto increment = (csAlmostHorizontal - mCriticalInclination) / (csInclDepthProfilePointCount - 1u);
  auto inclination = mCriticalInclination;
  std::deque<Sample> samplesBend;
  for(uint32_t i = 0u; i < csInclDepthProfilePointCount; ++i) {
    auto processed = process(getReflectionDispDepth(inclination));
    if(!processed.mAsphalt) {
      add(samplesBend, processed.mCollection);
    }
    else {
      throw "Should not hit asphalt in bending region.";
    }
    inclination += increment;
  }
  mCriticalInclination -= 2.0 * csEpsilon;
  std::deque<Sample> samplesAsphaltDown;
  std::deque<Sample> samplesAsphaltUp;
  // TODO sample asphalt hit region and gather these
  // TODO process samples*
}

PolynomialRayBending::Gather PolynomialRayBending::getReflectionDispDepth(double const aInclination0) const {
  double currentHeight = csInitialHeight;
  double currentTemp = getTempRiseAtHeight(currentHeight);
  double currentRefractionIndex = getRefractionAtTempRise(currentTemp);
  double sinInclination0 = ::sin(aInclination0);
  double currentSinInclination = sinInclination0;

  Gather result;
  uint32_t iterations = 0u;
  result.mAsphalt = true;
  while(currentHeight > csMinimalHeight) {
    ++iterations;
    double nextTemp = currentTemp + static_cast<long double>(csLayerDeltaTemp);
    double nextHeight = getHeightAtTempRise(nextTemp);
    double nextRefractionIndex = getRefractionAtHeight(nextHeight);
    double factor = currentRefractionIndex / nextRefractionIndex;
    double nextSinInclination = currentSinInclination * factor;
    double disp;
    if(nextSinInclination >= static_cast<double>(1.0)) { // The problem here is how much sinInclanation is greater than 1. Optimal would be to have it just touch 1, so the last temperature step can't be uniform.
      auto const criticalTemp = binarySearch(currentTemp, nextTemp, csEpsilon, [this, currentRefractionIndex, currentSinInclination](double const aTemp) {
        double height = getHeightAtTempRise(aTemp);
        double refractionIndex = getRefractionAtHeight(height);
        double factor = currentRefractionIndex / refractionIndex;
        double sinInclination = currentSinInclination * factor;
        return (sinInclination > static_cast<double>(1.0) ? 1.0 : -1.0);
      });
      double criticalHeight = getHeightAtTempRise(criticalTemp);
      disp = static_cast<double>(currentHeight - criticalHeight) * currentSinInclination / ::sqrt(static_cast<double>(1.0) - currentSinInclination * currentSinInclination);
      result.mCollection.emplace_back(currentHeight, 0.0, std::nan("1"));
      result.mAsphalt = false;
      break;
    }
    else {
      disp = static_cast<double>(currentHeight - nextHeight) * currentSinInclination / ::sqrt(static_cast<double>(1.0) - currentSinInclination * currentSinInclination);
      result.mCollection.emplace_back(currentHeight, disp, currentSinInclination);
      currentTemp = nextTemp;
      currentHeight = nextHeight;
      currentRefractionIndex = nextRefractionIndex;
      currentSinInclination = nextSinInclination;
    }
  }
  return result;
}

PolynomialRayBending::Intermediate PolynomialRayBending::process(Gather const aRaws) {
  Intermediate result;
  result.mAsphalt = aRaws.mAsphalt;
  if(!aRaws.mAsphalt) {
    double previousSinInclination;
    double sumDisp = 0.0;
    result.mCollection.reserve(aRaws.mCollection.size() * 2u - 1u);
    for(auto const &rawSample : aRaws.mCollection) {
      if(!std::isnan(rawSample.mSinInclination)) {
        previousSinInclination = rawSample.mSinInclination;
        result.mCollection.emplace_back(rawSample.mHeight, sumDisp, std::asin(rawSample.mSinInclination) - cgPi / 2.0);
      }
      else {
        result.mCollection.emplace_back(rawSample.mHeight, sumDisp, cgPi / 2.0 - std::asin(previousSinInclination));
      }
      sumDisp += rawSample.mHorizDisp;
    }
    auto i = aRaws.mCollection.crbegin();
    ++i;
    sumDisp += i->mHorizDisp;
    auto previousHeight = i->mHeight;
    previousSinInclination = i->mSinInclination;
    ++i;
    while(i != aRaws.mCollection.crend()) {
      result.mCollection.emplace_back(previousHeight, sumDisp, cgPi / 2.0 - std::asin(i->mSinInclination));
      sumDisp += i->mHorizDisp;
      previousHeight = i->mHeight;
      previousSinInclination = i->mSinInclination;
      ++i;
    }
    result.mCollection.emplace_back(previousHeight, sumDisp, cgPi / 2.0 - std::asin(previousSinInclination));
  }
  else {
    double sumDisp = 0.0;
    result.mCollection.reserve(aRaws.mCollection.size());
    for(auto const &rawSample : aRaws.mCollection) {
      result.mCollection.emplace_back(rawSample.mHeight, sumDisp, std::asin(rawSample.mSinInclination) - cgPi / 2.0);
      sumDisp += rawSample.mHorizDisp;
    }
  }
  return result;
}
