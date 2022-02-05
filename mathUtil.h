#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <Eigen/Dense>
#include <functional>

double signum(double const aValue);
double binarySearch(double const aLower, double const aUpper, double const aEpsilon, std::function<double(double)> aLambda);

class PolynomApprox final {
private:
  uint32_t const  mDegreeMinus;  // I wnt for Laurent-polynomial. Apparently the normal least squares algorithm can't cope with it.
  uint32_t const  mCoeffCount;   // Further mathematical investigations needed.

  double          mXmin;
  double          mSpanOriginal;
  double          mSpanFactor;
  double          mSpanStart;
  Eigen::VectorXd mCoefficients;
  double          mRrmsError;    // https://stats.stackexchange.com/questions/413209/is-there-something-like-a-root-mean-square-relative-error-rmsre-or-what-is-t

public:
  PolynomApprox(double const* const aSamplesX, double const* const aSamplesY, uint32_t const aSampleCount, uint32_t const aDegreeMinus, uint32_t const aDegreePlus);
  PolynomApprox(std::vector<double> const& aSamplesX, std::vector<double> const& aSamplesY, uint32_t const aDegreeMinus, uint32_t const aDegreePlus)
    : PolynomApprox(aSamplesX.data(), aSamplesY.data(), std::min(aSamplesX.size(), aSamplesY.size()), aDegreeMinus, aDegreePlus) {}

  uint32_t size()                          const { return mCoeffCount; }
  double operator[](uint32_t const aIndex) const { return mCoefficients(aIndex); }
  double getRrmsError()                    const { return mRrmsError; }

  double eval(double const aX) const;

private:
  double normalize(double const aX) const { return mSpanStart + mSpanFactor * (aX - mXmin) / mSpanOriginal; }

  // based on Table 1 of Condition number of Vandermonde matrix in least-squares polynomial fitting problems
  static constexpr double getXspan(double const aDegree) { // Optimal for big sample counts
    return 1.90313131 + aDegree * (-0.23114312 + aDegree * (0.02573205 - aDegree * 0.00098032));
  }
};

#endif // MATHUTIL_H
