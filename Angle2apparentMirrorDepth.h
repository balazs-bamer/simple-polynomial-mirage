#ifndef ANGLE2APPARENTMIRRORDEPTH
#define ANGLE2APPARENTMIRRORDEPTH

#include "3dGeomUtil.h"
#include "mathUtil.h"
#include <memory>
#include <optional>

// These calculations do not take relative humidity in account, since it has less, than 0.5% the effect on air refractive index as temperature and pressure.
class Angle2apparentMirrorDepth final {
public:
  static constexpr double   csCelsius2kelvin                = 273.15;

private:
  static constexpr double   csTempAmbient                   = 297.5;                                        // Kelvin
  static constexpr uint32_t csTempProfilePointCount         =   8u;
  static constexpr uint32_t csTempProfileDegree             =   4u;
  static constexpr double   csTplate[csTempProfilePointCount]      = { 336.7, 331.7, 326.7, 321.7, 316.7, 311.7, 306.7, 301.7 };  // Kelvin
  static constexpr double   csHeightLimit[csTempProfilePointCount] = { 0.25,  0.25,  0.25,  0.3,   0.35,  0.35,  0.4,   0.4   };  // cm
  static constexpr double   csB[csTempProfilePointCount]           = { 0.048, 0.041, 0.035, 0.030, 0.024, 0.018, 0.012, 0.006 };
  static constexpr double   csDelta[csTempProfilePointCount]       = { 1.4,   1.4,   1.5,   1.5,   1.6,   1.6,   1.6,   1.6   };
  static constexpr double   csDeltaFallback                 = 1.2;


  static constexpr uint32_t csInclDepthProfilePointCount    = 197u;
  static constexpr uint32_t csInclinationProfileDegree      =  17u;

  static constexpr double   csRelativeHumidityPercent       =  50.0;
  static constexpr double   csAtmosphericPressureKpa        = 101.0;
  static constexpr double   csLayerDeltaTemp                =   0.00005;
  static constexpr double   csInitialHeight                 =   1.0;      // TODO consider if this can be constant for larger plate temperatures.
  static constexpr double   csMinimalHeight                 =   0.001;
  static constexpr double   csAlmostVertical                =   0.01;
  static constexpr double   csAlmostHorizontal              =   (cgPi / 2.0) * 0.9999;
  static constexpr double   csEpsilon                       =   0.00001;

  struct Static final {
    std::unique_ptr<PolynomApprox> mHeightLimit;
    std::unique_ptr<PolynomApprox> mB;
    std::unique_ptr<PolynomApprox> mDelta;

    Static();
  };

  double mTempDiffSurface;
  double mHeightLimit;
  double mB;
  double mDelta;
  double mCriticalInclination;

  std::unique_ptr<PolynomApprox> mInclinationProfile;

public:
  Angle2apparentMirrorDepth(double const aTempDiffSurface); // TODO this should accept ambient temperature in the final version.

  static double getWorkingHeight()                                   { return csInitialHeight; }
  static double getAmbientTemp()    /* Kelvin */                     { return csTempAmbient; }
  double getCriticalInclination()                              const { return mCriticalInclination; }
  double getTempRiseAtHeight(double const aHeight)             const;  // delta Celsius and meter
  double getHeightAtTempRise(double const aTempRise)           const;  // delta Celsius and meter
  double getRefractionAtTempRise(double const aTempRise)       const;  // delta Celsius
  double getRefractionAtHeight(double const aHeight)           const { return getRefractionAtTempRise(getTempRiseAtHeight(aHeight)); }
  double approximateReflectionDepth(double const aInclination) const { return mInclinationProfile->eval(aInclination); }

private:
  void initReflection();
  std::optional<double> getReflectionDepth(double const aInclination) const;
};

#endif
