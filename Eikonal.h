#ifndef EIKONAL
#define EIKONAL

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <cmath>


// These calculations do not take relative humidity in account, since it has less, than 0.5% the effect on air refractive index as temperature and pressure.
class Eikonal final {
public:
  enum class Mode : uint8_t {
    cConventional = 0u,
    cPorous       = 1u
  };

private:
  static constexpr double   csCelsius2kelvin                = 273.15;
  static constexpr double   csC                             = 299792458.0; // m/s

  static constexpr uint32_t csTempProfilePointCount         =   8u;
  static constexpr uint32_t csTempProfileDegree             =   4u;
  static constexpr double   csTplate[csTempProfilePointCount]      = { 336.7, 331.7, 326.7, 321.7, 316.7, 311.7, 306.7, 301.7 };  // Kelvin
  static constexpr double   csHeightLimit[csTempProfilePointCount] = { 0.25,  0.25,  0.25,  0.3,   0.35,  0.35,  0.4,   0.4   };  // cm
  static constexpr double   csB[csTempProfilePointCount]           = { 0.048, 0.041, 0.035, 0.030, 0.024, 0.018, 0.012, 0.006 };
  static constexpr double   csDelta[csTempProfilePointCount]       = { 1.4,   1.4,   1.5,   1.5,   1.6,   1.6,   1.6,   1.6   };
  static constexpr double   csDeltaFallback                 = 1.2;

  static constexpr double   csRelativeHumidityPercent       =  50.0;
  static constexpr double   csAtmosphericPressureKpa        = 101.0;

/*private:
  struct Static final {
    std::unique_ptr<PolynomApprox> mHeightLimit;
    std::unique_ptr<PolynomApprox> mB;
    std::unique_ptr<PolynomApprox> mDelta;

    Static();
  };*/

  double mTempAmbient;      // Celsius
  Mode   mMode;
  double mHeightLimit;
  double mB;
  double mDelta;

public:
  static constexpr uint32_t csNvar = 6u;
  using Real                       = double;
  using Variables                  = std::array<Real, csNvar>;

  Eikonal(double const aTempAmbient, Mode const aMode) : mTempAmbient(aTempAmbient), mMode(aMode) {}

  void operator() (const double/* aS*/, std::array<Real, csNvar> const &aY, std::array<Real, csNvar> &aDyds) const {
    double n    = getRefract(aY[1]);
    double v    = csC / n;
//    double dvdz = -v / n * refractDiff(aY[2]);

    aDyds[0] = v * aY[3];
    aDyds[1] = v * aY[4];
    aDyds[2] = v * aY[5];
    aDyds[3] = 0.0;
    aDyds[4] = getRefractDiff(aY[1]) / v / n;
    aDyds[5] = 0.0;
    //aDyds[5] = -1.0 / v / v * dvdz;
  }

  int differentials(double, const double aY[], double aDydt[]) const {
    double n    = getRefract(aY[1]);
    double v    = csC / n;

    aDydt[0] = v * aY[3];
    aDydt[1] = v * aY[4];
    aDydt[2] = v * aY[5];
    aDydt[3] = 0.0;
    aDydt[4] = getRefractDiff(aY[1]) / v / n;
    aDydt[5] = 0.0;
    return GSL_SUCCESS;
  }

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
    return GSL_SUCCESS;
  }

public:
  double getRefract(double const aH) const {
    return mMode == Mode::cConventional ? getConventionalRefract(aH) : getPorousRefract(aH);
  }

  double getSlowness(double const aH) const {
    return mMode == Mode::cConventional ? getConventionalSlowness(aH) : getPorousSlowness(aH);
  }

  double getRefractDiff(double const aH) const {
    return mMode == Mode::cConventional ? getConventionalRefractDiff(aH) : getPorousRefractDiff(aH);
  }

  double getRefractDiff2(double const aH) const {
    return mMode == Mode::cConventional ? getConventionalRefractDiff2(aH) : getPorousRefractDiff2(aH);
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
    return 0.079386 * ((8245.75 * std::exp(-20.16 * aH)) / t / t / t
                     - (647.233 * std::exp(-10.08 * aH)) / t / t);
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
    return 0.079386 * ((137.777 * d * d * std::exp(-16.7 * aH)) / t / t / t
                     - (69.3042 * d * std::exp(-8.35 * aH)) / t / t);
  }
};

#endif
