#include "mathUtil.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <functional>


template<typename tStepper>
class OdeSolver {
private:
  using Real          = typename tStepper::Real;
  using Variables     = typename tStepper::Variables;
  using OdeDefinition = typename tStepper::OdeDefinition;

  static constexpr uint32_t csMaxStep = 1000u;
  static constexpr uint32_t csNvar    = tStepper::csNvar;
	Real                      mXstart;
	Real                      mXend;
	Real                      mStepStart;
	Real                      mStepMin;
	OdeDefinition&            mOdeDef;
	tStepper                  mStepper;

public:
	OdeSolver(const Real aXstart, const Real aXend, const Real aAtol, const Real aRtol, const Real aStepStart, const Real aStepMin, OdeDefinition &aOdeDefs);

	Variables solve(Variables const &aYstart, std::function<bool(std::array<double, csNvar> const&)> aJudge, Real const aEpsilon);
};

template<typename tStepper>
OdeSolver<tStepper>::OdeSolver(const Real aXstart, const Real aXend,
		const Real aAtol, const Real aRtol, const Real aStepStart,
		const Real aStepMin, OdeDefinition &aOdeDef)
: mXstart(aXstart)
, mXend(aXend)
, mStepStart(aStepStart)
, mStepMin(aStepMin)
, mOdeDef(aOdeDef)
,	mStepper(aAtol, aRtol, aOdeDef) {
}

template<typename tStepper>
typename OdeSolver<tStepper>::Variables OdeSolver<tStepper>::solve(Variables const &aYstart, std::function<bool(std::array<double, csNvar> const&)> aJudge, Real const aEpsilon) {
  Variables result;
  bool ready = false;
	auto x = mXstart;
  auto h = std::copysign(mStepStart, mXend - mXstart);
  mStepper.init(mXstart, aYstart);
	for (uint32_t steps = 0u; steps < csMaxStep; steps++) {
		if ((x + h * 1.0001 - mXend) * (mXend - mXstart) > 0.0) {
			h = mXend - x;
    }
		auto stepData = mStepper.step(x, h);
		if (aJudge(mStepper.getY())) {
      mStepper.prepareDense(stepData.xNow, stepData.hNow);
      binarySearch(stepData.xNow, stepData.xNext, aEpsilon, [&result, &stepData, this, &aJudge](auto const aWhere){
        result = mStepper.interpolate(stepData.xNow, aWhere, stepData.hNow);
        return aJudge(result);
      });
			ready = true;
      break;
		}
		if (abs(h) <= mStepMin) throw std::out_of_range("Step size too small in OdeSolver");
    x = stepData.xNext;
    h = stepData.hNext;
	}
  if(!ready) {
  	throw std::out_of_range("Too many steps in routine OdeSolver");
  }
  return result;
}
