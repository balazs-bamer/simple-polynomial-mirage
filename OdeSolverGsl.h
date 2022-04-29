#ifndef ODESOLVERGSL_H
#define ODESOLVERGSL_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <array>
#include <stdexcept>
#include <functional>

#include <vector>
#include <iostream>
#include <iomanip>

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
  OdeDefinition const&      mOdeDef;
  gsl_odeiv2_step          *mStepper = nullptr;
  gsl_odeiv2_control       *mController;
  gsl_odeiv2_evolve        *mEvolver;
  gsl_odeiv2_system         mSystem;

public:
  OdeSolverGsl(StepperType const aStepper, const double aTstart, const double aTend, const double aAtol, const double aRtol, const double aStepStart, OdeDefinition const& aOdeDef);
  ~OdeSolverGsl();

  std::pair<double, Variables> solve(Variables const &aYstart, std::function<bool(double const, std::array<double, csNvar> const&)> aJudge);
};

template <typename tOdeDefinition>
OdeSolverGsl<tOdeDefinition>::OdeSolverGsl(StepperType const aStepper, const double aTstart, const double aTend, const double aAtol, const double aRtol, const double aStepStart, OdeDefinition const& aOdeDef)
  : mTstart(aTstart)
  , mTend(aTend)
  , mStepStart(aStepStart)
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
std::vector<Variables> stuff;
  while(stepsAll < csMaxStep) {
    double h = mStepStart;
    double t = start;
    bool verdictPrev = aJudge(t, y);
    Variables yPrev;
    double tPrev;
    uint32_t stepsNow = 0;
if(stuff.empty() || y[0] > stuff.back()[0]) {
  std::cout << "pre while x: " << y[0] << " y: " << y[1] << " z: " << y[2] << '\n';
  stuff.push_back(y);
}
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
if(y[0] < mTend) {
  std::cout << " in while x: " << y[0] << " y: " << y[1] << " z: " << y[2] << '\n';
  stuff.push_back(y);
}
      if(verdictPrev != aJudge(t, y)) {
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

std::cout << "x=[";
for (int i = 0; i < stuff.size(); ++i) {
  std::cout << std::setprecision(10) << stuff[i][0] << (i < stuff.size() - 1 ? ", " : "];\n");
}
std::cout << "y=[";
for (int i = 0; i < stuff.size(); ++i) {
  std::cout << std::setprecision(10) << stuff[i][1] - Eikonal::csRadius << (i < stuff.size() - 1 ? ", " : "];\n");
}

  return result;
}

/*
./eikonal --dir -0.1 --target 1000 --dist 3333 --earthForm round
target = 980

0, 0.0099999845, 0.059999907, 0.3099995195, 1.559997582, 7.809987895, 39.05993946, 195.3096973, 976.5594881, [3333]
976.5694882, 976.6194885, 976.8694902, 978.1194985, 984.3695399, 
978.1294985, 978.1794989, 978.4295005, 979.6795088, 984.3695399, 
979.6895089, 979.7395092, 979.9895109, 981.2395192,
979.9995109, 980.0495113,
980.009511

1, 0.999982547, 0.9998952802, 0.9994589482, 0.9972772906, 0.9863690641, 0.9318298856, 0.6593213715, 10.19546658, [?]
10.19609492, 10.19923658, 10.21494488, 10.29348637, 10.68619384, 
10.2941147, 10.29725636, 10.31296466, 10.39150616, 10.68619384,
10.39213449, 10.39527615, 10.41098445, 10.48952594,
10.41161278, 10.41475444,
10.41224111
*/

#endif // ODESOLVERGSL_H
