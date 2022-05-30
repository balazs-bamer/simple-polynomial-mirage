#ifndef EIKONAL
#define EIKONAL

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <cmath>
#include <array>
#include <cstdint>


// These calculations do not take relative humidity in account, since it has less, than 0.5% the effect on air refractive index as temperature and pressure.
class Eikonal final {
public:
  enum class Temperature : uint8_t {
    cAmbient = 0u,
    cBase    = 1u,
    cMinimum = 2u,
    cMaximum = 3u
  };

  enum class Model : uint8_t {
    cConventional = 0u,
    cPorous       = 1u,
    cWater        = 2u
  };

  enum class EarthForm : uint8_t {
    cFlat  = 0u,
    cRound = 1u
  };

  static constexpr double   csCelsius2kelvin                = 273.15;
  static constexpr double   csC                             = 299792458.0; // m/s

private:
  static constexpr uint32_t csTempProfilePointCount         =   8u;
  static constexpr uint32_t csTempProfileDegree             =   4u;
  static constexpr double   csTplate[csTempProfilePointCount]      = { 336.7, 331.7, 326.7, 321.7, 316.7, 311.7, 306.7, 301.7 };  // Kelvin
  static constexpr double   csHeightLimit[csTempProfilePointCount] = { 0.25,  0.25,  0.25,  0.3,   0.35,  0.35,  0.4,   0.4   };  // cm
  static constexpr double   csB[csTempProfilePointCount]           = { 0.048, 0.041, 0.035, 0.030, 0.024, 0.018, 0.012, 0.006 };
  static constexpr double   csDelta[csTempProfilePointCount]       = { 1.4,   1.4,   1.5,   1.5,   1.6,   1.6,   1.6,   1.6   };
  static constexpr double   csDeltaFallback                 = 1.2;

  static constexpr double   csRelativeHumidityPercent       =  50.0;
  static constexpr double   csAtmosphericPressureKpa        = 101.0;

  EarthForm const mEarthForm;
  double    const mEarthRadius;
  Model     const mModel;
  double          mTempAmbient;    // Celsius
  double    const mTempAmbOrig; // Celsius
  double    const mTempAmbMin;
  double    const mTempAmbMax;
  double    const mTempBase;       // Celsius

public:
  static constexpr uint32_t csNvar = 6u;
  using Real                       = double;
  using Variables                  = std::array<Real, csNvar>;
  // For flat Earth, surface has v[4] == 0, the light travels mostly in v[3] direction, v[5] is depth.
  // The light starts close to the origin.
  // 
  // For round Earth, the origin is in the Earth center. We are on the sea on the radius of csRadius.
  // The light starts around v[3] and v[5] == 0, and travels mostly in v[3] direction.

  Eikonal(EarthForm const aEarthForm, double const aEarthRadius, Model const aModel, double const aTempAmbient)
  : mEarthForm(aEarthForm)
  , mEarthRadius(aEarthRadius)
  , mModel(aModel)
  , mTempAmbient(aTempAmbient)
  , mTempAmbOrig(aTempAmbient)
  , mTempAmbMin(aTempAmbient)
  , mTempAmbMax(aTempAmbient)
  , mTempBase(aTempAmbient) {}

  Eikonal(EarthForm const aEarthForm, double const aEarthRadius, Model const aModel, double const aTempAmbient, double const aTempAmbMin, double const aTempAmbMax, double const aTempBase)
  : mEarthForm(aEarthForm)
  , mEarthRadius(aEarthRadius)
  , mModel(aModel)
  , mTempAmbient(aTempAmbient)
  , mTempAmbOrig(aTempAmbient)
  , mTempAmbMin(aTempAmbMin)
  , mTempAmbMax(aTempAmbMax)
  , mTempBase(aTempBase) {}

  Eikonal(Eikonal const&) = default;
  Eikonal(Eikonal &&) = default;
  Eikonal& operator=(Eikonal const&) = delete;
  Eikonal& operator=(Eikonal &&) = delete;

  void setWaterTempAmb(Temperature const aWhich) {
    mTempAmbient = (aWhich == Temperature::cBase ? mTempBase :
                   (aWhich == Temperature::cMinimum ? mTempAmbMin :
                   (aWhich == Temperature::cMaximum ? mTempAmbMax : mTempAmbOrig)));
  }

  EarthForm getEarthForm()   const { return mEarthForm; }
  double    getEarthRadius() const { return mEarthRadius; }

  int differentials(double, const double aY[], double aDydt[]) const {
    int result;
    double n    = getRefract(aY[1]);
    double v    = csC / n;
    double elevation;
    std::array<double, 3u> zenith;
    if(mEarthForm == EarthForm::cFlat) {
      elevation = aY[1];
      zenith[0] = zenith[2] = 0.0;
      zenith[1] = 1.0;
    }
    else {
      double fromCenter = std::sqrt(aY[0] * aY[0] + aY[1] * aY[1] + aY[2] * aY[2]);
      elevation = fromCenter - mEarthRadius;
      if(elevation > 0.0) {
        zenith[0] = aY[0] / fromCenter;
        zenith[1] = aY[1] / fromCenter;
        zenith[2] = aY[2] / fromCenter;
      }
      else {} // nothing to do
    }
    if(elevation > 0.0) {
      result = GSL_SUCCESS;
      aDydt[0] = v * aY[3];
      aDydt[1] = v * aY[4];
      aDydt[2] = v * aY[5];
      double u = getRefractDiff(elevation) / csC;
      aDydt[3] = zenith[0] * u;
      aDydt[4] = zenith[1] * u;
      aDydt[5] = zenith[2] * u;
    }
    else {
      result = GSL_FAILURE;
    }
    return result;
  }

  // Most probably wrong.
  int jacobian(double, const double aY[], double *aDfdy, double aDfdt[]) const {
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (aDfdy, csNvar, csNvar);  // TODO do directly
    gsl_matrix * m = &dfdy_mat.matrix;
//    gsl_matrix_set (m, row, col, 0.0);
    double n    = getRefract(aY[1]);
    double nd1  = getRefractDiff(aY[1]);
    double nd2  = getRefractDiff2(aY[1]);
    double v    = csC / n;
    for(uint32_t i = 0u; i < csNvar; ++i) {
      gsl_matrix_set (m, i, 0u, 0.0);
      gsl_matrix_set (m, i, 2u, 0.0);
      gsl_matrix_set (m, i, 3u, (i == 0u ? v : 0.0));
      gsl_matrix_set (m, i, 4u, (i == 1u ? v : 0.0));
      gsl_matrix_set (m, i, 5u, (i == 2u ? v : 0.0));
      aDfdt[i] = 0.0;
    }
    gsl_matrix_set (m, 0u, 1u, -v * aY[3] * nd1 / n);
    gsl_matrix_set (m, 1u, 1u, -v * aY[4] * nd1 / n);
    gsl_matrix_set (m, 2u, 1u, -v * aY[5] * nd1 / n);
    gsl_matrix_set (m, 3u, 1u, 0.0);
    gsl_matrix_set (m, 4u, 1u, v * (nd2 / n - 2.0 * nd1 * nd1 / n / n));
    gsl_matrix_set (m, 5u, 1u, 0.0);
    return GSL_FAILURE; // because wrong
  }

public:
  double getRefract(double const aH) const {
    return mModel == Model::cConventional ? getConventionalRefract(aH) :
          (mModel == Model::cPorous ? getPorousRefract(aH) : getWaterRefract(aH));
  }

  double getSlowness(double const aH) const {
    return mModel == Model::cConventional ? getConventionalSlowness(aH) :
          (mModel == Model::cPorous ? getPorousSlowness(aH) : getWaterSlowness(aH));
  }

  double getRefractDiff(double const aH) const {
    return mModel == Model::cConventional ? getConventionalRefractDiff(aH) :
          (mModel == Model::cPorous ? getPorousRefractDiff(aH) : getWaterRefractDiff(aH));
  }

  double getRefractDiff2(double const aH) const {
    return mModel == Model::cConventional ? getConventionalRefractDiff2(aH) :
          (mModel == Model::cPorous ? getPorousRefractDiff2(aH) : getWaterRefractDiff2(aH));
  }

private:
  double getConventionalRefract(double const aH) const {
    auto celsius = mTempAmbient + 0.018 + 6.37 * std::exp(-aH * 10.08);
    return 1.0 + 7.86e-4 * 101 / (celsius + 273.15);
  }

  double getConventionalSlowness(double const aH) const {
    return getConventionalRefract(aH) / csC;
  }

  double getConventionalRefractDiff(double const aH) const {
    auto t = mTempAmbient + 6.37 * std::exp(-10.08 * aH) + 273.168;
    return 5.09734 * std::exp(-10.08 * aH) / t / t;
  }

  double getConventionalRefractDiff2(double const aH) const {
    auto t = mTempAmbient + 6.37 * std::exp(-10.08 * aH) + 273.168;
    return 0.079386 / t / t * ((8245.75 * std::exp(-20.16 * aH)) / t
                             - (647.233 * std::exp(-10.08 * aH)));
  }

  double getPorousRefract(double const aH) const {
    auto celsius = mTempAmbient + (-66.8 + 1.9 * mTempAmbient) * (0.002 + 0.994 * std::exp(-aH * 8.35));
    return 1.0 + 7.86e-4 * 101 / (celsius + 273.15);
  }

  double getPorousSlowness(double const aH) const {
    return getPorousRefract(aH) / csC;
  }

  double getPorousRefractDiff(double const aH) const {
    auto t = (1.9 * mTempAmbient - 66.8) * (0.002 + 0.994 * std::exp(-8.35 * aH)) + mTempAmbient + 273.15;
    return 0.658896 * (1.9 * mTempAmbient - 66.8) * std::exp(-8.35 * aH) / t / t;
  }

  double getPorousRefractDiff2(double const aH) const {
    auto t = (1.9 * mTempAmbient - 66.8) * (0.002 + 0.994 * std::exp(-8.35 * aH)) + mTempAmbient + 273.15;
    auto d = 1.9 * mTempAmbient - 66.8;
    return 0.079386 / t / t * d * ((137.777 * d * std::exp(-16.7 * aH)) / t
                                 - (69.3042 * std::exp(-8.35 * aH)));
  }

  double getWaterRefract(double const aH) const {
    auto celsius = mTempAmbient + (mTempBase - mTempAmbient)*(0.011 + 1.05 * std::exp(-20.1 * aH));
    return 1.0 + 7.86e-4 * 101 / (celsius + 273.15);
  }

  double getWaterSlowness(double const aH) const {
    return getWaterRefract(aH) / csC;
  }

  double getWaterRefractDiff(double const aH) const {
    auto t = (mTempBase - mTempAmbient) * (0.011 + 1.05 * std::exp(-20.1 * aH)) + mTempAmbient + 273.15;
    return 1.67544 * std::exp(-20.1 * aH) * (mTempBase - mTempAmbient) / t / t;
  }

  double getWaterRefractDiff2(double const aH) const {
    auto t = (mTempBase - mTempAmbient) * (0.011 + 1.05 * std::exp(-20.1 * aH)) + mTempAmbient + 273.15;
    auto d = mTempBase - mTempAmbient;
    return 0.079386 / t / t * d * ((890.842 * d * std::exp(-40.2 * aH)) / t
                                -  (424.211 * std::exp(-20.1 * aH)));
  }
};

#endif
