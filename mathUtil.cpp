#include "mathUtil.h"
#include <numeric>
#include <stdexcept>


double signum(double const aValue) {
  return std::copysign(1.0f, aValue);
}

double binarySearch(double const aLower, double const aUpper, double const aEpsilon, std::function<bool(double)> aLambda) {
  bool signLower = aLambda(aLower);
  double lower = aLower;
  double upper = aUpper;
  double result = (lower + upper) / 2.0f;
  while(upper - lower > aEpsilon) {
    bool now = aLambda(result);
    if(now == signLower) {
      lower = result;
    }
    else {
      upper = result;
    }
    result = (lower + upper) / 2.0f;
  }
  return result;
}

PolynomApprox::PolynomApprox(uint32_t const aSampleCount, double const * const aSamplesY, std::initializer_list<PolynomApprox::Var> aVarsX) {
  mTotalCoeffCount = 1u;
  for(auto &var : aVarsX) {
    mTotalCoeffCount *= var.mDegree + 1u;
  }
  mVariableCount = aVarsX.size();
  Eigen::MatrixXd vandermonde = Eigen::MatrixXd::Ones(aSampleCount, 1u);
  for(auto &var : aVarsX) {
    mDegrees.push_back(var.mDegree);
    if(mCumulativeCoeffCounts.empty()) {
      mCumulativeCoeffCounts.push_back(var.mDegree + 1u);
    }
    else {
      mCumulativeCoeffCounts.push_back((var.mDegree + 1u) * mCumulativeCoeffCounts.back());
    }
    auto [itMin, itMax] = std::minmax_element(var.mSamples, var.mSamples + aSampleCount);
    mXmins.push_back(*itMin);
    mSpanOriginals.push_back(*itMax - *itMin);
    auto span = getXspan(mTotalCoeffCount);
    mSpanFactor = span * 2.0;
    mSpanStart  = -span;

    Eigen::MatrixXd increment(aSampleCount, var.mDegree + 1u);
    for(uint32_t i = 0u; i < aSampleCount; ++i) {
      for(uint32_t j = 0u; j <= var.mDegree; ++j) {
        increment(i, j) = ::pow(normalize(var.mSamples[i], mDegrees.size() - 1u), j);
      }
    }
    Eigen::MatrixXd previous = vandermonde;
    uint32_t colsPrev = previous.cols();
    uint32_t colsIncr = increment.cols();
    vandermonde = Eigen::MatrixXd(aSampleCount, colsPrev * colsIncr);
    for(uint32_t r = 0u; r < aSampleCount; ++r) {
      for(uint32_t i = 0u; i < colsIncr; ++i) {
        for(uint32_t p = 0u; p < colsPrev; ++p) {
          vandermonde(r, i * colsPrev + p) = previous(r, p) * increment(r, i);
        }
      }
    }
  }
  Eigen::VectorXd y = Eigen::VectorXd::Map(aSamplesY, aSampleCount);
  mCoefficients = vandermonde.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);  // most accurate method according to https://eigen.tuxfamily.org/dox/group__LeastSquares.html

  for(uint32_t i = 0u; i < mVariableCount; ++i) {
    mActualExponents.push_back(0u);
    mActualPowers.emplace_back(mDegrees[i] + 1u, 0.0);
  }

  auto diffs = 0.0f;
  auto desireds = 0.0f;
  std::vector<double> x;
  for(uint32_t i = 0u; i < aSampleCount; ++i) {
    x.clear();
    for(auto &var : aVarsX) {
      x.push_back(var.mSamples[i]);
    }
    auto desired = aSamplesY[i];
    auto diff = desired - eval(x);  // TODO find out why does not compile for Anna
    diffs += diff * diff;
    desireds += desired * desired;
  }
  mRrmsError = (desireds > 0.0f ? ::sqrt(diffs / desireds / aSampleCount) : 0.0f);
}

double PolynomApprox::eval(std::vector<double> const& aVariables) const {
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

double PolynomApprox::eval(double const aX, uint32_t const aOffset) const {
  auto view = mCoefficients.data() + aOffset;
  double result = view[mDegrees[0]] * aX;
  for(uint32_t i = mDegrees[0] - 1u; i > 0u; --i) {
    result = (result + view[i]) * aX;
  }
  result += view[0u];
  return result;
}
