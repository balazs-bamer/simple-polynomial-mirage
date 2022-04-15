#ifndef ODESOLVERGSL_H
#define ODESOLVERGSL_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <array>
#include <stdexcept>
#include <functional>


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
  OdeDefinition const&      mOdeDef;
  gsl_odeiv2_step          *mStepper;
  gsl_odeiv2_control       *mController;
  gsl_odeiv2_evolve        *mEvolver;
  gsl_odeiv2_system         mSystem;

public:
  OdeSolverGsl(const double aTstart, const double aTend, const double aAtol, const double aRtol, const double aStepStart, OdeDefinition const& aOdeDef);
  ~OdeSolverGsl();

  std::pair<double, Variables> solve(Variables const &aYstart, std::function<bool(std::array<double, csNvar> const&)> aJudge);
};

template <typename tOdeDefinition>
OdeSolverGsl<tOdeDefinition>::OdeSolverGsl(const double aTstart, const double aTend, const double aAtol, const double aRtol, const double aStepStart, OdeDefinition const& aOdeDef)
  : mTstart(aTstart)
  , mTend(aTend)
  , mStepStart(aStepStart)
  , mOdeDef(aOdeDef)
  , mStepper(gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, csNvar))
  , mController(gsl_odeiv2_control_y_new(aAtol, aRtol))
  , mEvolver(gsl_odeiv2_evolve_alloc(csNvar)) {
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
std::pair<double, typename OdeSolverGsl<tOdeDefinition>::Variables> OdeSolverGsl<tOdeDefinition>::solve(Variables const &aYstart, std::function<bool(Variables const&)> aJudge) {
  std::pair<double, Variables> result;
  double start = mTstart;
  double end = mTend;
  Variables y = aYstart;
  uint32_t stepsAll = 0;
  while(stepsAll < csMaxStep) {
    double h = mStepStart;
    double t = start;
    bool verdictPrev = aJudge(y);
    Variables yPrev;
    double tPrev;
    uint32_t stepsNow = 0;
    while (t < end && stepsAll < csMaxStep) {
      yPrev = y;
      tPrev = t;
      int status = gsl_odeiv2_evolve_apply (mEvolver, mController, mStepper,
                                           &mSystem,
                                           &t, end,
                                           &h, y.data());
      if (status != GSL_SUCCESS) {
        gsl_odeiv2_evolve_reset(mEvolver);
        gsl_odeiv2_step_reset(mStepper);
        throw std::out_of_range("OdeSolverGsl: Can't apply step in evolver.");
      }
      else {} // Nothing to do
      ++stepsAll;
      ++stepsNow;
      if(verdictPrev != aJudge(y)) {
        break;
      }
      else {} // Nothing to do
    }
    gsl_odeiv2_evolve_reset(mEvolver);
    gsl_odeiv2_step_reset(mStepper);
    if(stepsNow == 1u) {
      result.first = t;
      result.second = y;
      break;
    }
    else {
      y = yPrev;
      start = tPrev;
      end = t;
    }
  }
  if(stepsAll == csMaxStep) {
    throw std::out_of_range("OdeSolverGsl: Too many steps.");
  }
  else {} // Nothing to do
  return result;
}

#endif // ODESOLVERGSL_H
