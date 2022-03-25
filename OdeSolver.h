#include <array>
#include <cmath>
#include <functional>


template<typename tStepper>
class OdeSolver {
private:
  using Real          = typename tStepper::Real;
  using OdeDefinition = typename tStepper::OdeDefinition;

	static constexpr uint32_t csMaxStep = 50000u;
  static constexpr uint32_t csNvar    = tStepper::csNvar;
	std::array<Real, csNvar>  mYstart;
	Real                      mXstart;
	Real                      mXend;
	Real                      mStepStart;
	Real                      mStepMin;
	OdeDefinition&            mOdeDef;
	tStepper                  mStepper;

public:
	OdeSolver(std::array<Real, csNvar> const &ayStart, const Real aXstart, const Real aXend,
		const Real aAtol, const Real aRtol, const Real aStepStart,
		const Real aStepMin, OdeDefinition &aOdeDefs);

	std::array<Real, csNvar> solve(std::function<bool(std::array<double, csNvar> const&)> aJudge);
};

template<typename tStepper>
OdeSolver<tStepper>::OdeSolver(std::array<Real, csNvar> const &aYstart, const Real aXstart, const Real aXend,
		const Real aAtol, const Real aRtol, const Real aStepStart,
		const Real aStepMin, OdeDefinition &aOdeDef)
: mYstart(aYstart)
, mXstart(aXstart)
, mXend(aXend)
, mStepStart(aStepStart)
, mStepMin(aStepMin)
, mOdeDef(aOdeDef)
,	mStepper(aAtol, aRtol, aOdeDef) {
}

template<typename tStepper>
std::array<typename OdeSolver<tStepper>::Real, OdeSolver<tStepper>::csNvar> OdeSolver<tStepper>::solve(std::function<bool(std::array<double, csNvar> const&)> aJudge) {
  std::array<Real, csNvar> result;
  result.fill(0.0);   // TODO remove
	auto x = mXstart;
  auto h = std::copysign(mStepStart, mXend - mXstart);
  mStepper.init(mXstart, mYstart);
	for (uint32_t steps = 0u; steps < csMaxStep; steps++) {
		if ((x + h * 1.0001 - mXend) * (mXend - mXstart) > 0.0) {
			h = mXend - x;
    }
		std::tie(x, h) = mStepper.step(x, h);
std::cout << x << '\n';
		if (aJudge(mStepper.getY())) {
			return mStepper.getY();
		}
		if (abs(h) <= mStepMin) throw("Step size too small in OdeSolver");
	}
	throw("Too many steps in routine OdeSolver");
}
