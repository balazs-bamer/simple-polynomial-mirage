#include "Angle2apparentMirrorDepth.h"

#include <functional>
#include <cmath>

#include <iostream> // TODO remove


Angle2apparentMirrorDepth::Static::Static() {
  mStuff = std::move(std::make_unique<TempTempProfiles>());
  for(uint32_t i = 0u; i < csTempProfilePointCount; ++i) {
    mStuff->at(i) = std::move(std::make_unique<PolynomApprox>(csReferenceTempProfiles[0u].data(), csReferenceTempProfiles[i].data(), csTempProfileCount, 0u, csReferenceTempProfileDegree));
  }
}

Angle2apparentMirrorDepth::Angle2apparentMirrorDepth(double const aTempAmbient, double const aTempDiffSurface)
  : mTempAmbient(aTempAmbient)
  , mTempDiffSurface(std::max(0.0, aTempDiffSurface)) {
  static Static tempTempProfiles;
  std::array<double, csTempProfilePointCount> temperatures;

  for(uint32_t i = 0u; i < csTempProfilePointCount; ++i) {
    temperatures[i] = tempTempProfiles.eval(i, mTempDiffSurface + csReferenceTempAmbient);
  }
  mTempProfile = std::move(std::make_unique<PolynomApprox>(csTempProfileHeights.data(), temperatures.data(), csTempProfilePointCount, csTempProfileDegreeMinus, csTempProfileDegreePlus));

std::cout << "temp error: " << mTempProfile->getRrmsError() << "\ntemp coeff: ";
for(int i = 0; i < mTempProfile->size(); ++i)
  std::cout << mTempProfile->operator[](i) << ' ';
std::cout << '\n';

  initReflection();
}

void Angle2apparentMirrorDepth::initReflection() {
  auto const minInclination = binarySearch(csAlmostVertical, csAlmostHorizontal, csEpsilon, [this](double const aInclination) {
    return getReflectionDepth(aInclination) ? -1.0f : 1.0f;
  }) + csEpsilon;
  std::vector<double> inclinations;
  std::vector<double> depths;
  auto increment = (csAlmostHorizontal - minInclination) / (csInclinationProfilePointCount - 1u);
  auto inclination = minInclination;
  for(uint32_t i = 0u; i < csInclinationProfilePointCount; ++i) {
    inclinations.push_back(inclination);
    depths.push_back(getReflectionDepth(inclination).value());
std::cout << "incl: " << inclination * 180.0 / cgPi << "  depth: " << depths.back() << '\n';
    inclination += increment;
  }
  mInclinationProfile = std::move(std::make_unique<PolynomApprox>(inclinations, depths, csInclinationProfileDegreeMinus, csInclinationProfileDegreePlus));
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

//std::cout << currentHeight << "   " << currentTemp << " - " << sinInclination0 << '\n';
  std::optional<double> result;
  while(currentHeight > csMinimalHeight * 2.0f) {
    auto nextTemp = currentTemp + csLayerDeltaTemp;
    auto nextHeight = binarySearch(csMinimalHeight, currentHeight, csEpsilon, [this, nextTemp](double const aHeight) {
      return getTempAtHeight(aHeight) - nextTemp;
    });
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
//std::cout << currentHeight << "   " << currentTemp << " - " << sinInclination << "   " << horizDisp << '\n';
  }
  return result;
}
