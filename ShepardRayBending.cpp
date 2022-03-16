#include "ShepardRayBending.h"

#include <functional>
#include <stdexcept>
#include <vector>
#include <cmath>

#include <iostream> // TODO remove

ShepardRayBending::Static::Static() {
  mHeightLimit = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csHeightLimit, std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
  mB           = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csB,           std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
  mDelta       = std::move(std::make_unique<PolynomApprox>(csTempProfilePointCount, csDelta,       std::initializer_list<PolynomApprox::Var>{PolynomApprox::Var{csTplate, csTempProfileDegree}}));
}

ShepardRayBending::ShepardRayBending(double const aTempDiffSurface)
  : mTempDiffSurface(std::max(0.0, aTempDiffSurface))
  , mRandomDistribution(-1.0, 1.0) {
  static Static tempCoeffProfiles;

  mHeightLimit = tempCoeffProfiles.mHeightLimit->eval(csTempAmbient + mTempDiffSurface);
  mB           = tempCoeffProfiles.mB->eval(csTempAmbient + mTempDiffSurface);
  mDelta       = tempCoeffProfiles.mDelta->eval(csTempAmbient + mTempDiffSurface);
  initReflection();
}

double ShepardRayBending::getTempRiseAtHeight(double const aHeight) const {
  auto height = 100.0 * aHeight;
  return csTempAmbient * ::exp(mB * ::pow(height, 1.0 - mDelta)) - csTempAmbient;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double ShepardRayBending::getHeightAtTempRise(double const aTempRise) const {
  return ::exp(::log(::log((aTempRise + csTempAmbient) / csTempAmbient) / mB) / (1.0 - mDelta)) / 100.0;  // Here I neglect the difference between mDelta and csDeltaFallback, because this result will only be used indirectly.
}

double ShepardRayBending::getRefractionAtTempRise(double const aTempRise) const {
  auto tempKelvin = csTempAmbient + aTempRise;
  auto tempCelsius = tempKelvin - csCelsius2kelvin;
  auto r = 1.0 + 7.86e-4 * csAtmosphericPressureKpa / tempKelvin - 1.5e-11 * csRelativeHumidityPercent * (tempCelsius * tempCelsius + 160.0);
  return r;
}


void ShepardRayBending::initReflection() {
  mCriticalInclination = binarySearch(csAlmostVertical, csAlmostHorizontal, csEpsilon, [this](double const aInclination) {
    return traceHalf(aInclination).mAsphalt ? 1.0 : -1.0;
  }) + csEpsilon;
std::cout << "ready 1: critical inclination\n";
  auto increment = (csAlmostHorizontal - mCriticalInclination) / (csRayTraceCountBending - 1u);
  auto inclination = mCriticalInclination;
  typename ActualShepard::DataTransfer samplesForward;
  for(uint32_t i = 0u; i < csRayTraceCountBending; ++i) {
    auto rayPath = toRayPath(traceHalf(inclination));
    if(!rayPath.mAsphalt) {
      addForward(samplesForward, rayPath.mCollection, csDispSampleFactorBending);
    }
    else {
      throw std::runtime_error("Should not hit asphalt in bending region.");
    }
    inclination += increment;
  }
  mShepardBending = std::make_unique<ActualShepard>(samplesForward, csSamplesToConsider, csAverageRelativeSize, csShepardExponent);
  samplesForward.clear();
std::cout << "ready 2: mShepardBending\n";

  typename ActualShepard::DataTransfer samplesBackward;
  mCriticalInclination -= 2.0 * csEpsilon;
  increment = (mCriticalInclination - csAsphaltRayAngleLimit) / (csRayTraceCountAsphalt - 1u);
  inclination = csAsphaltRayAngleLimit;
  for(uint32_t i = 0u; i < csRayTraceCountAsphalt; ++i) {
    auto rayPath = toRayPath(traceHalf(inclination));
    if(rayPath.mAsphalt) {
      addForward(samplesForward, rayPath.mCollection, csDispSampleFactorAsphalt);
      addReverse(samplesBackward, rayPath.mCollection, csDispSampleFactorAsphalt);
    }
    else {
      throw std::runtime_error("Should hit asphalt outside bending region.");
    }
    inclination += increment;
  }
  mShepardAsphaltDown = std::make_unique<ActualShepard>(samplesForward, csSamplesToConsider, csAverageRelativeSize, csShepardExponent);
  mShepardAsphaltUp   = std::make_unique<ActualShepard>(samplesBackward, csSamplesToConsider, csAverageRelativeSize, csShepardExponent);
  samplesForward.clear();
  samplesBackward.clear();
std::cout << "ready 3: mShepardAsphaltDown\n";
std::cout << "ready 4: mShepardAsphaltUp\n";
}

ShepardRayBending::Gather ShepardRayBending::traceHalf(double const aInclination0) const {
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

ShepardRayBending::Intermediate ShepardRayBending::toRayPath(Gather const aRaws) {
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

#include<iostream>

void ShepardRayBending::addForward(typename ActualShepard::DataTransfer &aCollector, std::vector<Sample> const &aLot, double const aDispSampleFactor) const {
  auto valid = std::find_if(aLot.begin() + 1u, aLot.end(), [](auto &x){return x.mHorizDisp == 0.0; }) - aLot.begin();
  if(valid > 3u) {
    auto indices = getRandomIndices(valid, static_cast<uint32_t>(aLot[valid - 1u].mHorizDisp * aDispSampleFactor));
    for(uint32_t i = 1u; i < indices.size(); ++i) {
      for(uint32_t j = 0u; j < i; ++j) {          // j -> i
        auto from = indices[j];
        auto to = indices[i];
        typename ActualShepard::Data item;
        item.mLocation[csIndexLocationStartHeight] = aLot[from].mHeight;
        item.mLocation[csIndexLocationStartDir] = aLot[from].mAngleFromHoriz;
        item.mLocation[csIndexLocationHorizDisp] = aLot[to].mHorizDisp - aLot[from].mHorizDisp;
        item.mPayload[csIndexPayloadHeight]  = aLot[to].mHeight;
        item.mPayload[csIndexPayloadDir]  = aLot[to].mAngleFromHoriz;
        aCollector.push_back(item);
      }
    }
  }
  else {} // nothing to do
}

void ShepardRayBending::addReverse(typename ActualShepard::DataTransfer &aCollector, std::vector<Sample> const &aLot, double const aDispSampleFactor) const {
  auto valid = std::find_if(aLot.begin() + 1u, aLot.end(), [](auto &x){return x.mHorizDisp == 0.0; }) - aLot.begin();
  if(valid > 3u) {
    auto indices = getRandomIndices(valid, static_cast<uint32_t>(aLot[valid - 1u].mHorizDisp * aDispSampleFactor));
    for(uint32_t i = 1u; i < indices.size(); ++i) {
      for(uint32_t j = 0u; j < i; ++j) {          // j -> i
        auto to = indices[j];
        auto from = indices[i];
        typename ActualShepard::Data item;
        item.mLocation[csIndexLocationStartHeight] = aLot[from].mHeight;
        item.mLocation[csIndexLocationStartDir] = aLot[from].mAngleFromHoriz;
        item.mLocation[csIndexLocationHorizDisp] = aLot[from].mHorizDisp - aLot[to].mHorizDisp;
        item.mPayload[csIndexPayloadHeight]  = aLot[to].mHeight;
        item.mPayload[csIndexPayloadDir]  = aLot[to].mAngleFromHoriz;
        aCollector.push_back(item);
      }
    }
  }
  else {} // nothing to do
}

std::vector<uint32_t> ShepardRayBending::getRandomIndices(uint32_t const aFromCount, uint32_t const aChosenCount) const {
  std::vector<uint32_t> result;
  auto count = std::min(aFromCount, aChosenCount);
  result.reserve(count);
  if(aChosenCount >= aFromCount) {
     result.resize(count);
     std::iota(result.begin(), result.end(), 0u);
  }
  else {
    double step = static_cast<double>(aFromCount - 1u) / (aChosenCount - 1u);
    auto randomRadius = (step * csRelativeRandomRadius > 1.0 ? csRelativeRandomRadius : 0.0 );
    for(int32_t i = 0; i < aChosenCount; ++i) {
      auto index = static_cast<int32_t>(::round(step * (i + randomRadius * mRandomDistribution(mRandomEngine))));
      index = std::max(0, index);
      index = std::min(index, static_cast<int32_t>(aFromCount - 1));
      result.push_back(index);
    }
  }
  return result;
}
