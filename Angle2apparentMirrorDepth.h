#ifndef ANGLE2APPARENTMIRRORDEPTH
#define ANGLE2APPARENTMIRRORDEPTH

#include "3dGeomUtil.h"
#include "mathUtil.h"
#include <memory>
#include <optional>

// These calculations do not take relative humidity in account, since it has less, than 0.5% the effect on air refractive index as temperature and pressure.
class Angle2apparentMirrorDepth final {
private:
  static constexpr uint32_t csReferenceTempProfileDegree    =   2u;
  static constexpr uint32_t csTempProfileDegreeMinus        =   2u;
  static constexpr uint32_t csTempProfileDegreePlus         =   0u;  // 1 would be better, but then the height-temp function is not monotonous.
  static constexpr uint32_t csTempProfilePointCount         =   6u;
  static constexpr uint32_t csTempProfileCount              =   3u;
  static constexpr uint32_t csInclinationProfilePointCount  =   6u;
  static constexpr uint32_t csInclinationProfileDegreeMinus =   0u;
  static constexpr uint32_t csInclinationProfileDegreePlus  =   3u;

  static constexpr double    csRelativeHumidityPercent      =  50.0;
  static constexpr double    csAtmosphericPressureKpa       = 101.0;
  static constexpr double    csReferenceTempAmbient         =  20.0;
  static constexpr double    csLayerDeltaTemp               =   0.05;
  static constexpr double    csInitialHeight                =   0.9;  // Just a little bit over the regression range.
  static constexpr double    csMinimalHeight                =   0.0005;
  static constexpr double    csAlmostVertical               =   0.01;
  static constexpr double    csAlmostHorizontal             =   (cgPi / 2.0) * 0.999;
  static constexpr double    csEpsilon                      =   0.0001;

  using TempTempProfiles  = std::array<std::unique_ptr<PolynomApprox>, csTempProfilePointCount>;

  class Static final {
  private:
    std::unique_ptr<TempTempProfiles> mStuff;

  public:
    Static();

    double eval(uint32_t const aIndex, double const mInternalTemp) const { return mStuff->at(aIndex)->eval(mInternalTemp); }
  };

  static constexpr std::array<double, csTempProfilePointCount> csTempProfileHeights = { csMinimalHeight, 0.01, 0.03, 0.13, 0.55, 0.89 }; // meters
  static constexpr std::array<std::array<double, csTempProfileCount>, csTempProfilePointCount> csReferenceTempProfiles = {{ // Absolute temperature based on csReferenceTempAmbient.
    {{ 32.0f, 48.0f, 63.0f }},
    {{ 25.0f, 33.0f, 42.0f }},
    {{ 21.0f, 27.0f, 30.0f }},
    {{ 20.5f, 23.0f, 25.0f }},
    {{ 20.0f, 20.5f, 21.0f }},
    {{ 20.0f, 20.0f, 20.0f }}
  }};

  double mTempAmbient;
  double mTempDiffSurface;

  std::unique_ptr<PolynomApprox> mTempProfile;
  std::unique_ptr<PolynomApprox> mInclinationProfile;

public:
  Angle2apparentMirrorDepth(double const mTempAmbient, double const aTempDiffSurface);

  double getTempAtHeight(double const aHeight)       const { return mTempProfile->eval(aHeight) - csReferenceTempAmbient + mTempAmbient; }
  double getRefractionAtTemp(double const aTemp)     const { return 1.0f + 7.86e-4f * csAtmosphericPressureKpa / (273.15f + aTemp) - 1.5e-11f * csRelativeHumidityPercent * (aTemp * aTemp + 160); }
  double getRefractionAtHeight(double const aHeight) const { return getRefractionAtTemp(getTempAtHeight(aHeight)); }

private:
  void initReflection();
  std::optional<double> getReflectionDepth(double const aInclination) const;
};

#endif
