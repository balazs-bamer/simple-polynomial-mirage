#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <array>

double signum(double const aValue);
double binarySearch(double const aLower, double const aUpper, double const aEpsilon, std::function<double(double)> aLambda);

class PolynomApprox final {
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
  std::vector<uint32_t>            mActualExponents;        // Just to avoid re-allocating it for each evaluation.
  std::vector<std::vector<double>> mActualPowers;           // Just to avoid re-allocating it for each evaluation.
  // TODO consider double          mRrmsError;    // https://stats.stackexchange.com/questions/413209/is-there-something-like-a-root-mean-square-relative-error-rmsre-or-what-is-t

public:
  PolynomApprox() : mTotalCoeffCount(1u), mVariableCount(0u) {}

  template<typename... tArgs>
  void init(uint32_t const aSampleCount, double const * const aSamplesY, tArgs... aRest) {
    mTotalCoeffCount = gatherCoeffCount(aRest...);
    Eigen::VectorXd y = Eigen::VectorXd::Map(aSamplesY, aSampleCount);
    auto vandermonde = calculateVandermonde(aSampleCount, aRest...);
    mCoefficients = (vandermonde.transpose() * vandermonde).inverse() * vandermonde.transpose() * y;
    initActuals();
  }

  template<typename... tArgs>
  void init(std::vector<double> const& aSamplesY, tArgs... aRest) {
    mTotalCoeffCount = gatherCoeffCount(aRest...);
    Eigen::VectorXd y = Eigen::VectorXd::Map(aSamplesY.data(), aSamplesY.size());
    auto vandermonde = calculateVandermonde(aSamplesY.size(), aRest...);
    mCoefficients = (vandermonde.transpose() * vandermonde).inverse() * vandermonde.transpose() * y;
    initActuals();
  }

  double eval(std::initializer_list<double> const aVariables);
  double eval(double const aX) { return eval(normalize(aX, 0u), 0u); }

private:
  uint32_t gatherCoeffCount() {
    return 1u;
  }

  template<typename... tArgs>
  uint32_t gatherCoeffCount(double const * const, uint32_t const aDegree, tArgs... aRest) {
    return (aDegree + 1u) * gatherCoeffCount(aRest...);
  }

  template<typename... tArgs>
  uint32_t gatherCoeffCount(std::vector<double> const&, uint32_t const aDegree, tArgs... aRest) {
    return (aDegree + 1u) * gatherCoeffCount(aRest...);
  }

  Eigen::MatrixXd calculateVandermonde(uint32_t const aSampleCount) {
    return Eigen::MatrixXd::Ones(aSampleCount, 1u);
  }

  template<typename... tArgs>
  Eigen::MatrixXd calculateVandermonde(uint32_t const aSampleCount, double const * const aSamplesX, uint32_t const aDegree, tArgs... aRest) {
    append(aSampleCount, aSamplesX, aDegree);
    return merge(doCalculateVandermonde(mVariableCount - 1u, aSampleCount, aSamplesX, aDegree), calculateVandermonde(aSampleCount, aRest...));
  }

  template<typename... tArgs>
  Eigen::MatrixXd calculateVandermonde(uint32_t const aSampleCount, std::vector<double> const& aSamplesX, uint32_t const aDegree, tArgs... aRest) {
    append(aSampleCount, aSamplesX.data(), aDegree);
    return merge(doCalculateVandermonde(mVariableCount - 1u, aSamplesX.size(), aSamplesX.data(), aDegree), calculateVandermonde(aSampleCount, aRest...));
  }

  void initActuals();
  void append(uint32_t const aSampleCount, double const * const aSamplesX, uint32_t const aDegree);
  Eigen::MatrixXd doCalculateVandermonde(uint32_t const aVariableIndex, uint32_t const aSampleCount, double const * const aSamplesX, uint32_t const aDegree);
  Eigen::MatrixXd merge(Eigen::MatrixXd const& aOne, Eigen::MatrixXd const& aOther);

  double normalize(double const aX, uint32_t const aIndex) const { return mSpanStart + mSpanFactor * (aX - mXmins[aIndex]) / mSpanOriginals[aIndex]; }

  double eval(double const aX, uint32_t const aOffset) const;

  // based on Table 1 of Condition number of Vandermonde matrix in least-squares polynomial fitting problems
  static constexpr double getXspan(double const aDegree) { // Optimal for big sample counts
    return 1.90313131 + aDegree * (-0.23114312 + aDegree * (0.02573205 - aDegree * 0.00098032));
  }
};

#endif // MATHUTIL_H
