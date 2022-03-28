#include "mathUtil.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <functional>

#include<iostream>


template<typename tStepper>
class OdeSolver {
private:
  using Real          = typename tStepper::Real;
  using Variables     = typename tStepper::Variables;
  using OdeDefinition = typename tStepper::OdeDefinition;

	static constexpr uint32_t csMaxStep = 50000u;
  static constexpr uint32_t csNvar    = tStepper::csNvar;
	Variables                 mYstart;
	Real                      mXstart;
	Real                      mXend;
	Real                      mStepStart;
	Real                      mStepMin;
	OdeDefinition&            mOdeDef;
	tStepper                  mStepper;

public:
	OdeSolver(Variables const &ayStart, const Real aXstart, const Real aXend,
		const Real aAtol, const Real aRtol, const Real aStepStart,
		const Real aStepMin, OdeDefinition &aOdeDefs);

	Variables solve(std::function<bool(std::array<double, csNvar> const&)> aJudge, Real const aEpsilon);
};

template<typename tStepper>
OdeSolver<tStepper>::OdeSolver(Variables const &aYstart, const Real aXstart, const Real aXend,
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
typename OdeSolver<tStepper>::Variables OdeSolver<tStepper>::solve(std::function<bool(std::array<double, csNvar> const&)> aJudge, Real const aEpsilon) {
  Variables result;
  bool ready = false;
	auto x = mXstart;
  auto h = std::copysign(mStepStart, mXend - mXstart);
  mStepper.init(mXstart, mYstart);
	for (uint32_t steps = 0u; steps < csMaxStep; steps++) {
		if ((x + h * 1.0001 - mXend) * (mXend - mXstart) > 0.0) {
			h = mXend - x;
    }
		auto stepData = mStepper.step(x, h);
std::cout << "xNow: " << stepData.xNow << " y: " << mStepper.getY()[0] << '\n';
		if (aJudge(mStepper.getY())) {
      mStepper.prepareDense(stepData.xNow, stepData.hNow);
std::cout << "ready: " << stepData.xNow << '-' << stepData.xNext << '\n';
      binarySearch(stepData.xNow, stepData.xNext, aEpsilon, [&result, &stepData, this, &aJudge](auto const aWhere){
std::cout << "s: " << aWhere << '\n';
        result = mStepper.interpolate(stepData.xNow, aWhere, stepData.hNow);
        return aJudge(result);
      });
			ready = true;
      break;
		}
		if (abs(h) <= mStepMin) throw std::out_of_range("Step size too small in OdeSolver");
    x = stepData.xNext;
    h = stepData.hNext;
std::cout << "-------------------------------------\n";
	}
  if(!ready) {
  	throw std::out_of_range("Too many steps in routine OdeSolver");
  }
  return result;
}
