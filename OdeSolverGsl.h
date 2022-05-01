#ifndef ODESOLVERGSL_H
#define ODESOLVERGSL_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <array>
#include <stdexcept>
#include <functional>


enum class StepperType : uint8_t {
  cRungeKutta23                = 0u,
  cRungeKuttaClass4            = 1u,
  cRungeKuttaFehlberg45        = 2u,
  cRungeKuttaCashKarp45        = 3u,
  cRungeKuttaPrinceDormand89   = 4u,
  cBulirschStoerBaderDeuflhard = 5u
};

template <typename tOdeDefinition>
class OdeSolverGsl final {
public:

  static constexpr uint32_t csNvar = tOdeDefinition::csNvar;
  using Variables                  = std::array<double, csNvar>;

private:
  using OdeDefinition              = tOdeDefinition;

  static constexpr uint32_t csMaxStep = 1000u;
  double                    mTstart;
  double                    mTend;
  double                    mStepStart;
  double                    mStepMax;
  OdeDefinition const&      mOdeDef;
  gsl_odeiv2_step          *mStepper = nullptr;
  gsl_odeiv2_control       *mController;
  gsl_odeiv2_evolve        *mEvolver;
  gsl_odeiv2_system         mSystem;

public:
  OdeSolverGsl(StepperType const aStepper, const double aTstart, const double aTend, const double aAtol, const double aRtol, const double aStepStart, double const aStepMax, OdeDefinition const& aOdeDef);
  ~OdeSolverGsl();

  std::pair<double, Variables> solve(Variables const &aYstart, std::function<bool(double const, std::array<double, csNvar> const&)> aJudge);
};

template <typename tOdeDefinition>
OdeSolverGsl<tOdeDefinition>::OdeSolverGsl(StepperType const aStepper, const double aTstart, const double aTend, const double aAtol, const double aRtol, const double aStepStart, double const aStepMax, OdeDefinition const& aOdeDef)
  : mTstart(aTstart)
  , mTend(aTend)
  , mStepStart(aStepStart)
  , mStepMax(aStepMax)
  , mOdeDef(aOdeDef)
  , mController(gsl_odeiv2_control_y_new(aAtol, aRtol))
  , mEvolver(gsl_odeiv2_evolve_alloc(csNvar)) {

  if(aStepper == StepperType::cRungeKutta23) {
    mStepper = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk2, csNvar);
  }
  else if(aStepper == StepperType::cRungeKuttaClass4) {
    mStepper = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk4, csNvar);
  }
  if(aStepper == StepperType::cRungeKuttaFehlberg45) {
    mStepper = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, csNvar);
  }
  if(aStepper == StepperType::cRungeKuttaCashKarp45) {
    mStepper = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck, csNvar);
  }
  else if(aStepper == StepperType::cRungeKuttaPrinceDormand89) {
    mStepper = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, csNvar);
  }
  else if(aStepper == StepperType::cBulirschStoerBaderDeuflhard) {
    mStepper = gsl_odeiv2_step_alloc(gsl_odeiv2_step_bsimp, csNvar);
  }
  else {} // nothing to do

  mSystem.function = [](double aT, double const aY[], double aDydt[], void *aObject)->int {
    auto object = reinterpret_cast<OdeDefinition*>(aObject);
    return object->differentials(aT, aY, aDydt);
  };
  mSystem.jacobian = [](double aT, double const aY[], double *aDfdy, double aDfdt[], void *aObject)->int {
    auto object = reinterpret_cast<OdeDefinition*>(aObject);
    return object->jacobian(aT, aY, aDfdy, aDfdt);
  };
  mSystem.dimension = csNvar;
  mSystem.params = const_cast<OdeDefinition*>(&mOdeDef);
}

template <typename tOdeDefinition>
OdeSolverGsl<tOdeDefinition>::~OdeSolverGsl() {
  gsl_odeiv2_evolve_free(mEvolver);
  gsl_odeiv2_control_free (mController);
  gsl_odeiv2_step_free(mStepper);
}

template <typename tOdeDefinition>
std::pair<double, typename OdeSolverGsl<tOdeDefinition>::Variables> OdeSolverGsl<tOdeDefinition>::solve(Variables const &aYstart, std::function<bool(double const, Variables const&)> aJudge) {
  std::pair<double, Variables> result;
  double start = mTstart;
  double end = mTend;
  Variables y = aYstart;
  uint32_t stepsAll = 0;
  while(true) {
    double h = mStepStart;
    double t = start;
    bool verdictPrev = aJudge(t, y);
    Variables yPrev;
    double tPrev;
    uint32_t stepsNow = 0;
    bool wasBigH = false;
    while (t < end && stepsAll < csMaxStep) {
      yPrev = y;
      tPrev = t;
      int status;
      do {
        status = gsl_odeiv2_evolve_apply (mEvolver, mController, mStepper,
                                         &mSystem,
                                         &t, end,
                                         &h, y.data());
        h /= 2.0;
      } while(status == GSL_FAILURE);
      if (status != GSL_SUCCESS) {
        gsl_odeiv2_evolve_reset(mEvolver);
        gsl_odeiv2_step_reset(mStepper);
        throw std::out_of_range("OdeSolverGsl: Can't apply step in evolver.");
      }
      else {} // Nothing to do
      ++stepsAll;
      ++stepsNow;
      if(verdictPrev != aJudge(t, y)) {
        break;
      }
      else {} // Nothing to do
      if(h > mStepMax) {              // If h is too big, it may make a too big step yielding false results GSL unable to detect.
        wasBigH = true;
        break;
      }
      else {} // Nothing to do
    }
    gsl_odeiv2_evolve_reset(mEvolver);
    gsl_odeiv2_step_reset(mStepper);
    if(!wasBigH && stepsNow == 1u) {
      result.first = t;
      result.second = y;
      break;
    }
    else {
      y = yPrev;
      start = tPrev;
      if(!wasBigH) {
        end = t;
      }
      else{} // nothing to do
    }
    if(stepsAll == csMaxStep) {
      throw std::out_of_range("OdeSolverGsl: Too many steps.");
    }
    else {} // Nothing to do
  }
  return result;
}

#endif // ODESOLVERGSL_H
