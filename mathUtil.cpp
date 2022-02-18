#include "mathUtil.h"
#include <numeric>
#include <stdexcept>


double signum(double const aValue) {
  return std::copysign(1.0f, aValue);
}

double binarySearch(double const aLower, double const aUpper, double const aEpsilon, std::function<double(double)> aLambda) {
  double signLower = signum(aLambda(aLower));
  double lower = aLower;
  double upper = aUpper;
  double result = (lower + upper) / 2.0f;
  while(upper - lower > aEpsilon) {
    double now = signum(aLambda(result));
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

void PolynomApprox::initActuals() {
  std::fill_n(mActualExponents.begin(), mVariableCount, 0u);
  for(uint32_t i = 0u; i < mVariableCount; ++i) {
    mActualPowers.emplace_back(mDegrees[i] + 1u, 0.0);
  }
}

void PolynomApprox::append(uint32_t const aSampleCount, double const * const aSamplesX, uint32_t const aDegree) {
  ++mVariableCount;
  mDegrees.push_back(aDegree);
  if(mCumulativeCoeffCounts.empty()) {
    mCumulativeCoeffCounts.push_back(aDegree + 1u);
  }
  else {
    mCumulativeCoeffCounts.push_back((aDegree + 1u) * mCumulativeCoeffCounts.back());
  }
  auto [itMin, itMax] = std::minmax_element(aSamplesX, aSamplesX + aSampleCount);
  mXmins.push_back(*itMin);
  mSpanOriginals.push_back(*itMax - *itMin);
  auto span = getXspan(mTotalCoeffCount);
  mSpanFactor = span * 2.0;
  mSpanStart  = -span;
}

Eigen::MatrixXd PolynomApprox::doCalculateVandermonde(uint32_t const aVariableIndex, uint32_t const aSampleCount, double const * const aSamplesX, uint32_t const aDegree) {
  Eigen::MatrixXd result(aSampleCount, aDegree + 1u);

  for(uint32_t i = 0u; i < aSampleCount; ++i) {
    for(uint32_t j = 0u; j <= aDegree; ++j) {
       result(i, j) = ::pow(normalize(aSamplesX[i], aVariableIndex), j);
    }
  }
  return result;
}

Eigen::MatrixXd PolynomApprox::merge(Eigen::MatrixXd const& aOne, Eigen::MatrixXd const& aOther) {
  uint32_t sampleCount = aOne.rows();
  uint32_t colsOne = aOne.cols();
  uint32_t colsOther = aOther.cols();
  Eigen::MatrixXd result(sampleCount, colsOne * colsOther);

  for(uint32_t r = 0u; r < sampleCount; ++r) {
    for(uint32_t o = 0u; o < colsOther; ++o) {
      for(uint32_t t = 0u; t < colsOne; ++t) {
        result(r, o * colsOne + t) = aOne(r, t) * aOther(r, o);
      }
    }
  }
  return result;
}

double PolynomApprox::eval(std::initializer_list<double> const aVariables) {
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
