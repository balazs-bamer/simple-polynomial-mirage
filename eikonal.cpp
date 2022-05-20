#include "Eikonal.h"
#include "3dGeomUtil.h"
#include "OdeSolverGsl.h"
#include "RungeKuttaRayBending.h"
#include "mathUtil.h"
#include "CLI11.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>


// g++ -I/usr/include/eigen3 -I../repos/eigen-initializer_list/src -DEIGEN_MATRIX_PLUGIN=\"Matrix_initializer_list.h\" -DEIGEN_ARRAY_PLUGIN=\"Array_initializer_list.h\" -std=c++17 eikonal.cpp RungeKuttaRayBending.cpp mathUtil.cpp -o eikonal -O2 -lgsl

RungeKuttaRayBending::Result comp1(RungeKuttaRayBending::Parameters const& aParameters, Eikonal::EarthForm const aEarthForm, double const aEarthRadius, Eikonal::Model aMode, double aTempAmb, double aTempBase,
double aDir, double aHeight, double aTarget) {

  Eikonal eikonal(aEarthForm, aEarthRadius, aMode, aTempAmb, aTempBase);
  RungeKuttaRayBending rk(aParameters, eikonal);
  Vertex start(0.0, aHeight, 0.0);
  Vector dir(std::cos(aDir / 180.0 * cgPi), std::sin(aDir / 180.0 * cgPi), 0.0);
  RungeKuttaRayBending::Result solution;
  if(aTarget == 0.0) {
    solution.mValid = true;
    solution.mValue[0] = 0.0;
    solution.mValue[1] = aHeight;
    solution.mValue[2] = 0.0;
  }
  else {
    solution = rk.solve4x(start, dir, aTarget);
  }
  return solution;
}

void comp(RungeKuttaRayBending::Parameters const& aParameters, Eikonal::EarthForm const aEarthForm, double const aEarthRadius, Eikonal::Model aMode, double aTempAmb, double aTempBase,
double aDir, double aHeight, double aTarget, uint32_t aSamples) {
  std::vector<Vertex> stuff;
  std::ofstream out("values.txt");
  auto end = aTarget * (1.0 + 0.5 / aSamples);
  for(double t = 0.0; t <= end; t += aTarget / aSamples) {
    RungeKuttaRayBending::Result solution = comp1(aParameters, aEarthForm, aEarthRadius, aMode, aTempAmb, aTempBase, aDir, aHeight, t);
    if(solution.mValid) {
      stuff.push_back(solution.mValue);
      out << std::setprecision(10) << solution.mValue[0] << '\t' << std::setprecision(10) << solution.mValue[1] << '\n';
    }
  }
  std::cout << "x=[";
  for (int i = 0; i < stuff.size(); ++i) {
    std::cout << std::setprecision(10) << stuff[i][0] << (i < stuff.size() - 1 ? ", " : "];\n");
  }
  std::cout << "y=[";
  for (int i = 0; i < stuff.size(); ++i) {
    std::cout << std::setprecision(10) << stuff[i][1] << (i < stuff.size() - 1 ? ", " : "];\n");
  }
  if(aEarthForm == Eikonal::EarthForm::cRound) {
    std::cout << "d=[";
    for (int i = 0; i < stuff.size(); ++i) {
      std::cout << std::setprecision(10) << std::sqrt(aEarthRadius * aEarthRadius - stuff[i][0] * stuff[i][0]) - aEarthRadius << (i < stuff.size() - 1 ? ", " : "];\n");
    }
  } 
  else {} // nothing to do
  std::cout << "\n";
}

int main(int aArgc, char **aArgv) {
  CLI::App opt{"Usage"};
  std::string nameBase = "water";
  opt.add_option("--base", nameBase, "base type (conventional / porous / water) [water]");
  double dir = std::nan("");
  opt.add_option("--dir", dir, "start direction, neg downwards (degrees) [computed to touch surface]");
  bool findTouch = false;
  opt.add_option("--findTouch", findTouch, "whether to find direction to touch the surface (true, false) [false]");
  double dist = 2000.0;
  opt.add_option("--dist", dist, "distance along the ray to track (m) [2000]");
  std::string nameForm = "round";
  opt.add_option("--earthForm", nameForm, "Earth form (flat / round) [round]");
  double rawRadius = 6371.0;
  opt.add_option("--earthRadius", rawRadius, "Earth radius (km) [6371.0]");
  double height = 1.1;
  opt.add_option("--height", height, "start height (m) [1.1]");
  double maxCosDirChange = 0.99999999999;
  opt.add_option("--maxCosDirChange", maxCosDirChange, "Maximum of cos of direction change to reset big step [0.99999999999]");
  double samples = 100;
  opt.add_option("--samples", samples, "number of samples on ray [100]");
  bool silent = false;
  opt.add_option("--silent", silent, "surpress parameter echo (true, false) [false]");
  double step1 = 0.01;
  opt.add_option("--step1", step1, "initial step size (m) [0.01]");
  double stepMin = 1e-7;
  opt.add_option("--stepMin", stepMin, "maximal step size (m) [1e-7]");
  double stepMax = 55.5;
  opt.add_option("--stepMax", stepMax, "maximal step size (m) [55.5]");
  std::string nameStepper = "RungeKuttaFehlberg45";
  opt.add_option("--stepper", nameStepper, "stepper type (RungeKutta23 / RungeKuttaClass4 / RungeKuttaFehlberg45 / RungeKuttaCashKarp45 / RungeKuttaPrinceDormand89 / BulirschStoerBaderDeuflhard) [RungeKuttaFehlberg45]");
  double target = 1000.0;
  opt.add_option("--target", target, "horizontal distance to travel (m) [1000]");
  double tempAmb = std::nan("");
  opt.add_option("--tempAmb", tempAmb, "ambient temperature (Celsius) [20 for conventional, 38.5 for porous, 10 for water]");
  double tempBase = 13.0;
  opt.add_option("--tempBase", tempBase, "base temperature, only for water (Celsius) [13]");
  double tolAbs = 0.001;
  opt.add_option("--tolAbs", tolAbs, "absolute tolerance (m) [1e-3]");
  double tolRel = 0.001;
  opt.add_option("--tolRel", tolRel, "relative tolerance (m) [1e-3]");
  CLI11_PARSE(opt, aArgc, aArgv);

  Eikonal::Model base;
  if(nameBase == "conventional") {
    base = Eikonal::Model::cConventional;
  }
  else if(nameBase == "porous") {
    base = Eikonal::Model::cPorous;
  }
  else if(nameBase == "water") {
    base = Eikonal::Model::cWater;
  }
  else {
    std::cerr << "Illegal base value: " << nameBase << '\n';
    return 1;
  }

  Eikonal::EarthForm earthForm;
  if(nameForm == "flat") {
    earthForm = Eikonal::EarthForm::cFlat;
  }
  else if(nameForm == "round") {
    earthForm = Eikonal::EarthForm::cRound;
  }
  else {
    std::cerr << "Illegal Earth form value: " << nameForm << '\n';
    return 1;
  }

  double earthRadius = rawRadius * 1000.0;

  StepperType stepper;
  if(nameStepper == "RungeKutta23") {
    stepper = StepperType::cRungeKutta23;
  }
  else if(nameStepper == "RungeKuttaClass4") {
    stepper = StepperType::cRungeKuttaClass4;
  }
  else if(nameStepper == "RungeKuttaFehlberg45") {
    stepper = StepperType::cRungeKuttaFehlberg45;
  }
  else if(nameStepper == "RungeKuttaCashKarp45") {
    stepper = StepperType::cRungeKuttaCashKarp45;
  }
  else if(nameStepper == "RungeKuttaPrinceDormand89") {
    stepper = StepperType::cRungeKuttaPrinceDormand89;
  }
  else if(nameStepper == "BulirschStoerBaderDeuflhard") {
    stepper = StepperType::cBulirschStoerBaderDeuflhard;
  }
  else {
    std::cerr << "Illegal stepper value: " << nameStepper << '\n';
    return 1;
  }

  if(std::isnan(tempAmb)) {
    tempAmb = (base == Eikonal::Model::cConventional ? 20.0 :
              (base == Eikonal::Model::cPorous ? 38.5 : 10.0));
  }
  else {} // nothing to do

  double const cTolerance = 1e-6;

  RungeKuttaRayBending::Parameters parameters;
  parameters.mStepper         = stepper;
  parameters.mDistAlongRay    = dist;
  parameters.mTolAbs          = cTolerance;
  parameters.mTolRel          = cTolerance;
  parameters.mStep1           = step1;
  parameters.mStepMin         = stepMin;
  parameters.mStepMax         = stepMax;
  parameters.mMaxCosDirChange = maxCosDirChange;

  if(std::isnan(dir)) {
    dir = binarySearch(-45, 0.0, cTolerance, [&parameters, earthForm, earthRadius, base, tempAmb, tempBase, height, target](auto const angle){
      auto solution = comp1(parameters, earthForm, earthRadius, base, tempAmb, tempBase, angle, height, target);
      return solution.mValid;
    }) + cTolerance;
  }

  double horizonDirection = 0.0;
  bool   downwardsPrev = true;
  bool   validHorizon = true;
  double delta = dir;
  while(std::abs(delta) > cTolerance) {
    horizonDirection += delta;
std::cout << horizonDirection << '\n';
    auto solution = comp1(parameters, earthForm, earthRadius, base, tempAmb, tempBase, horizonDirection, height, target);
    if(!solution.mValid) {
      validHorizon = false;
      break;
    }
    else {} // nothing to do
    bool downwardsNow = (solution.mDirection(1) < 0);
    if(downwardsPrev != downwardsNow) {
      delta /= -3.0;
      downwardsPrev = downwardsNow;
    }
    else {} // nothing to do
  }

  parameters.mTolAbs          = tolAbs;
  parameters.mTolRel          = tolRel;

  if(!silent) {
    std::cout << "base type:                               .  .  .  " << nameBase << ' ' << static_cast<int>(base) << '\n';
    std::cout << "start direction, neg downwards (degrees):         " << dir << '\n';
    std::cout << "distance along the ray to track (m):              " << dist << '\n';
    std::cout << "Earth form:         .  .  .  .  .  .  .  .  .  .  " << nameForm << ' ' << static_cast<int>(earthForm) << '\n';
    std::cout << "Earth radius (km):                                " << earthRadius / 1000.0 << '\n';
    std::cout << "start height (m):                                 " << height << '\n';
    std::cout << "max of cos of direction change to reset big step: " << std::setprecision(17) << maxCosDirChange << '\n';
    std::cout << "number of samples on ray:                         " << samples << '\n';
    std::cout << "initial step size (m):                            " << step1 << '\n';
    std::cout << "minimal step size (m): .  .  .  .  .  .  .  .  .  " << stepMin << '\n';
    std::cout << "maximal step size (m):                            " << stepMax << '\n';
    std::cout << "stepper type:                                     " << nameStepper << ' ' << static_cast<int>(stepper) << '\n';
    std::cout << "horizontal distance to travel (m): .  .  .  .  .  " << target << '\n';
    std::cout << "ambient temperature (Celsius):                    " << tempAmb << '\n';
    std::cout << "base temperature, only for water (Celsius):       " << tempBase << '\n';
    std::cout << "absolute tolerance (m):   .  .  .  .  .  .  .  .  " << tolAbs << '\n';
    std::cout << "relative tolerance (m):                           " << tolRel << '\n';
    if(validHorizon) {
      std::cout << "horizon (computed) (degrees):                     " << horizonDirection << '\n';
    }
    else {
      std::cout << "horizon not calculable, inrease target\n";
    }
  }
  else {} // nothing to do
 
  comp(parameters, earthForm, earthRadius, base, tempAmb, tempBase, dir, height, target, samples);
  return 0;
}
