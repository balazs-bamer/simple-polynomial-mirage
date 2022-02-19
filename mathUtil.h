#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <array>

double signum(double const aValue);
double binarySearch(double const aLower, double const aUpper, double const aEpsilon, std::function<double(double)> aLambda);

class PolynomApprox final {
public:
  struct Var {
    double const * const mSamples;
    uint32_t const       mDegree;

    Var(double const * const aSamples, uint32_t const aDegree) : mSamples(aSamples), mDegree(aDegree) {}
    Var(std::vector<double> const& aSamples, uint32_t const aDegree) : mSamples(aSamples.data()), mDegree(aDegree) {}
  };

private:
  uint32_t                         mTotalCoeffCount;
  uint32_t                         mVariableCount;
  std::vector<uint32_t>            mCumulativeCoeffCounts;
  std::vector<uint32_t>            mDegrees;               // For all the outer vectors, the variables follow each other as in the API parameter list.
  std::vector<double>              mXmins;                 // The first variable will be expanded directly using Horner's rule.
  std::vector<double>              mSpanOriginals;
  double                           mSpanFactor;
  double                           mSpanStart;
  Eigen::VectorXd                  mCoefficients;
  double                           mRrmsError;    // https://stats.stackexchange.com/questions/413209/is-there-something-like-a-root-mean-square-relative-error-rmsre-or-what-is-t

  mutable std::vector<uint32_t>            mActualExponents;        // Just to avoid re-allocating it for each evaluation.
  mutable std::vector<std::vector<double>> mActualPowers;           // Just to avoid re-allocating it for each evaluation.

public:
  PolynomApprox(uint32_t const aSampleCount, double const * const aSamplesY, std::initializer_list<Var> aVarsX);

  PolynomApprox(std::vector<double> const& aSamplesY, std::initializer_list<Var> aVarsX) : PolynomApprox(aSamplesY.size(), aSamplesY.data(), aVarsX) {}

  double getRrmsError() const { return mRrmsError; }

  template<template<typename> typename tContainer>         // This allows it to be called with initializer_list and vector
  double eval(tContainer<double> const aVariables) const;

  double eval(double const aX) const { return eval(normalize(aX, 0u), 0u); }

private:
  double normalize(double const aX, uint32_t const aIndex) const { return mSpanStart + mSpanFactor * (aX - mXmins[aIndex]) / mSpanOriginals[aIndex]; }

  double eval(double const aX, uint32_t const aOffset) const;

  // based on Table 1 of Condition number of Vandermonde matrix in least-squares polynomial fitting problems
  static constexpr double getXspan(double const aDegree) { // Optimal for big sample counts
    return 1.90313131 + aDegree * (-0.23114312 + aDegree * (0.02573205 - aDegree * 0.00098032));
  }
};

template<template<typename> typename tContainer>
double PolynomApprox::eval(tContainer<double> const aVariables) const {
  if(aVariables.size() != mVariableCount) {
    throw std::invalid_argument("eval: variable count mismatch.");
  }
  else {} // nothing to do
  double actualNormalized0 = 0.0;
  auto arg = aVariables.begin();
  for(uint32_t v = 0u; v < mVariableCount; ++v) {
    auto normalized = normalize(*arg, v);
    auto& actualPowers = mActualPowers[v];
    if(v == 0u) {
      actualNormalized0 = normalized;
    }
    else {
      double actual = 1.0;
      for(uint32_t e = 0u; e <= mDegrees[v]; ++e) {
        actualPowers[e] = actual;
        actual *= normalized;
      }
    }
    ++arg;
  }
  double result = 0.0;
  uint32_t degree0 = mDegrees[0];
  std::fill(mActualExponents.begin(), mActualExponents.end(), 0u);
  for(uint32_t i = 0u; i < mTotalCoeffCount; i += degree0 + 1u) {
    double rest = 1.0;
    for(uint32_t v = 1u; v < mVariableCount; ++v) {
      rest *= mActualPowers[v][mActualExponents[v]];
    }
    result += rest * eval(actualNormalized0, i);
    ++mActualExponents[mVariableCount - 1u];
    for(uint32_t v = mVariableCount - 1u; v > 0u; --v) {
      if(mActualExponents[v] > mDegrees[v]) {
        mActualExponents[v] = 0u;
        ++mActualExponents[v - 1u];
      }
    }
  }
  return result;
}

#endif // MATHUTIL_H
