#include "mathUtil.h"

#include <iostream>


double signum(double const aValue) {
  return std::copysign(1.0f, aValue);
}

double binarySearch(double const aLower, double const aUpper, double const aEpsilon, std::function<double(double)> aLambda) {
//std::cout << "--- " << aLower << "    " << aUpper << '\n';
  double signLower = signum(aLambda(aLower));
  double lower = aLower;
  double upper = aUpper;
  double result = (lower + upper) / 2.0f;
  while(upper - lower > aEpsilon) {
    double now = signum(aLambda(result));
    if(now == signLower) {
//std::cout << "++ " << now << "   " << result << '\n';
      lower = result;
    }
    else {
//std::cout << "-- " << now << "   " << result << '\n';
      upper = result;
    }
    result = (lower + upper) / 2.0f;
  }
//std::cout << "===========================\n";
  return result;
}


PolynomApprox::PolynomApprox(double const* const aSamplesX, double const* const aSamplesY, uint32_t const aSampleCount, uint32_t const aDegreeMinus, uint32_t const aDegreePlus)
  : mDegreeMinus(aDegreeMinus)
  , mCoeffCount(aDegreeMinus + 1u + aDegreePlus)
  , mCoefficients(mCoeffCount) {
  Eigen::MatrixXd fitA(aSampleCount, mCoeffCount);
  Eigen::VectorXd fitB = Eigen::VectorXd::Map(aSamplesY, aSampleCount);

  if(aSampleCount >= mCoeffCount) {
    for(uint32_t i = 0u; i < aSampleCount; ++i) {
      for(uint32_t j = 0u; j < mCoeffCount; ++j) {
        fitA(i, j) = ::pow(aSamplesX[i], static_cast<int32_t>(j) - static_cast<int32_t>(aDegreeMinus));
      }
    }
    mCoefficients = fitA.colPivHouseholderQr().solve(fitB);
    auto diffs = 0.0f;
    auto desireds = 0.0f;
    for(uint32_t i = 0u; i < aSampleCount; ++i) {
      auto desired = aSamplesY[i];
      auto diff = (desired - eval(aSamplesX[i]));
      diffs += diff * diff;
      desireds += desired * desired;
    }
    mRrmsError = (desireds > 0.0f ? ::sqrt(diffs / desireds / aSampleCount) : 0.0f);
    }
  else {
    for(uint32_t j = 0u; j < mCoeffCount; ++j) {
      mCoefficients(j) = 0.0f;
    }
    mRrmsError = std::numeric_limits<double>::max();
  }
}

double PolynomApprox::eval(double const aX) const {
  double result = mCoefficients(mCoeffCount - 1u) * aX;
  for(uint32_t i = mCoeffCount - 2u; i > 0u; --i) {
    result = (result + mCoefficients(i)) * aX;
  }
  result += mCoefficients(0);
  return result / ::pow(aX, mDegreeMinus);
}


/* TODO create unit tests with these and more:

TEST(polynomApprox, x2few) {
  PolynomApprox poly({1.0f, 2.0f}, {1.0f, 4.0f, 9.0f}, 0u, 2u);
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(1.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f),  0.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), std::numeric_limits<double>::max()));
}

TEST(polynomApprox, x20) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f}, {0.0f, 0.0f, 0.0f}, 0u, 2u);
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(5.0f),  0.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_21_27_30) {
  PolynomApprox poly({32.0f, 48.0f, 63.0f}, {21.0f, 27.0f, 30.0f}, 0u, 2u);
  EXPECT_TRUE(eq(poly.eval(32.0f),  21.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(48.0f),  27.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(63.0f),  30.0f, 0.1f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 9.0f}, 0u, 2u);
  EXPECT_TRUE(poly.size() == 3u);
  EXPECT_TRUE(eq(poly[0u],  0.0f));
  EXPECT_TRUE(eq(poly[1u],  0.0f));
  EXPECT_TRUE(eq(poly[2u],  1.0f));
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f), 16.0f));
  EXPECT_TRUE(eq(poly.eval(5.0f), 25.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_5) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f}, {6.0f, 9.0f, 14.0f}, 0u, 2u);
  EXPECT_TRUE(poly.size() == 3u);
  EXPECT_TRUE(eq(poly[0u],  5.0f));
  EXPECT_TRUE(eq(poly[1u],  0.0f));
  EXPECT_TRUE(eq(poly[2u],  1.0f));
  EXPECT_TRUE(eq(poly.eval(0.0f),  5.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(4.0f), 21.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(5.0f), 30.0f, 0.1f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2off) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f, 4.0f}, {1.0f, 4.0f, 9.0f, 15.9f}, 0u, 2u);
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(4.0f), 16.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(5.0f), 25.0f, 0.3f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.00059f));
}

TEST(polynomApprox, x2veryoff) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f, 4.0f}, {1.0f, 3.0f, 13.0f, 13.0f}, 0u, 2u);
  EXPECT_FALSE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_FALSE(eq(poly.eval(4.0f), 16.0f));
  EXPECT_FALSE(eq(poly.eval(5.0f), 25.0f));
  EXPECT_TRUE(poly.getRrmsError() > 0.1f);
}

TEST(polynomApprox, 2x3_3x2_4x_5) {
  PolynomApprox poly({0.0f, 1.0f, 2.0f, 3.0f}, {5.0f, 14.0f, 41.0f, 98.0f}, 0u, 3u);
  EXPECT_TRUE(poly.size() == 4u);
  EXPECT_TRUE(eq(poly.eval(0.5f),  8.0f));
  EXPECT_TRUE(eq(poly.eval(1.5f), 24.5f));
  EXPECT_TRUE(eq(poly.eval(2.5f), 65.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}
*/
