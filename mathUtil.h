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

  double eval(std::initializer_list<double> const& aVariables) const {
    std::vector<double> variables;  // Hope they fit in the default reserve. Probabaly yes.
    for(auto d : aVariables) {
      variables.push_back(d);
    }
    return eval(variables);
  }

  double eval(std::vector<double> const& aVariables) const;

  double eval(double const aX) const { return eval(normalize(aX, 0u), 0u); }

private:
  double normalize(double const aX, uint32_t const aIndex) const { return mSpanStart + mSpanFactor * (aX - mXmins[aIndex]) / mSpanOriginals[aIndex]; }

  double eval(double const aX, uint32_t const aOffset) const;

  // based on Table 1 of Condition number of Vandermonde matrix in least-squares polynomial fitting problems
  static constexpr double getXspan(double const aDegree) { // Optimal for big sample counts
    return 1.90313131 + aDegree * (-0.23114312 + aDegree * (0.02573205 - aDegree * 0.00098032));
  }
};

#endif // MATHUTIL_H
