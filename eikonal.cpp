#include  <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "CLI11.hpp"
#include "Eikonal.h"
#include "3dGeomUtil.h"
#include "OdeSolverGsl.h"
#include "RungeKuttaRayBending.h"


//clang++ -I/usr/include/eigen3 -I../repos/eigen-initializer_list/src -DEIGEN_MATRIX_PLUGIN=\"Matrix_initializer_list.h\" -DEIGEN_ARRAY_PLUGIN=\"Array_initializer_list.h\" -std=c++17 eikonal.cpp RungeKuttaRayBending.cpp -o eikonal -ggdb -lgsl

void comp(StepperType aStepper,
Eikonal::EarthForm const aEarthForm, Eikonal::Model aMode, double aTempAmb, double aTempBase,
double aDir, double aDist, double aHeight, double aStep1, double aStepMax, double aTolAbs, double aTolRel, double aTarget, uint32_t aSamples) {

  Eikonal eikonal(aEarthForm, aMode, aTempAmb, aTempBase);
  RungeKuttaRayBending rk(aStepper, aDist, aTolAbs, aTolRel, aStep1, aStepMax, eikonal);
  Vertex start(0.0, aHeight, 0.0);
  Vector dir(std::cos(aDir / 180.0 * 3.1415926539), std::sin(aDir / 180.0 * 3.1415926539), 0.0);

  std::vector<Vertex> stuff;
  std::ofstream out("values.txt");
  for(double t = 0.0; t < aTarget; t += aTarget / aSamples) {
    RungeKuttaRayBending::Result solution;
    if(t == 0.0) {
      solution.mValid = true;
      solution.mValue[0] = 0.0;
      solution.mValue[1] = aHeight;
      solution.mValue[2] = 0.0;
    }
    else {
      solution = rk.solve4x(start, dir, t);
    }
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
      std::cout << std::setprecision(10) << std::sqrt(Eikonal::csRadius * Eikonal::csRadius - stuff[i][0] * stuff[i][0]) - Eikonal::csRadius << (i < stuff.size() - 1 ? ", " : "];\n");
    }
  } 
  else {} // nothing to do
  std::cout << "\n";
/*  std::cout << "h=[";
  for (auto i = 0.01; i < 1; i += .01) {
    std::cout << std::setprecision(10) << i << ", ";
  }
  std::cout << "n=[";
  for (auto i = 0.01; i < 1; i += .01) {
    std::cout << std::setprecision(10) << eikonal.getRefract(i) << ", ";
  }
  std::cout << "\n";*/
}

int main(int aArgc, char **aArgv) {
  CLI::App opt{"Usage"};
  std::string nameBase = "water";
  opt.add_option("--base", nameBase, "base type (conventional / porous / water) [water]");
  double dir = 0.0;
  opt.add_option("--dir", dir, "start direction, neg downwards (degrees) [0]");
  double dist = 2000.0;
  opt.add_option("--dist", dist, "distance along the ray to track (m) [2000]");
  std::string nameForm = "round";
  opt.add_option("--earthForm", nameForm, "Earth form (flat / round) [round]");
  double height = 1.0;
  opt.add_option("--height", height, "start height (m) [1]");
  double samples = 100;
  opt.add_option("--samples", samples, "number of samples on ray [100]");
  bool silent = false;
  opt.add_option("--silent", silent, "surpress parameter echo (true, false) [false]");
  double step1 = 0.01;
  opt.add_option("--step1", step1, "initial step size (m) [0.01]");
  double stepMax = 111.1;
  opt.add_option("--stepMax", stepMax, "maximal step size (m) [111.1]");
  std::string nameStepper = "RungeKutta23";
  opt.add_option("--stepper", nameStepper, "stepper type (RungeKutta23 / RungeKuttaClass4 / RungeKuttaFehlberg45 / RungeKuttaCashKarp45 / RungeKuttaPrinceDormand89 / BulirschStoerBaderDeuflhard) [RungeKutta23]");
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

  if(!silent) {
    std::cout << "base type:                                  " << nameBase << ' ' << static_cast<int>(base) << '\n';
    std::cout << "start direction, neg downwards (degrees):   " << dir << '\n';
    std::cout << "distance along the ray to track (m):        " << dist << '\n';
    std::cout << "Earth form:                                 " << nameForm << ' ' << static_cast<int>(earthForm) << '\n';
    std::cout << "start height (m):         .  .  .  .  .  .  " << height << '\n';
    std::cout << "number of samples on ray:                   " << samples << '\n';
    std::cout << "initial step size (m):                      " << step1 << '\n';
    std::cout << "maximal step size (m): .  .  .  .  .  .  .  " << stepMax << '\n';
    std::cout << "stepper type:                               " << nameStepper << ' ' << static_cast<int>(stepper) << '\n';
    std::cout << "horizontal distance to travel (m):          " << target << '\n';
    std::cout << "ambient temperature (Celsius):              " << tempAmb << '\n';
    std::cout << "base temperature, only for water (Celsius): " << tempBase << '\n';
    std::cout << "absolute tolerance (m):                     " << tolAbs << '\n';
    std::cout << "relative tolerance (m):                     " << tolRel << '\n';
  }
  else {} // nothing to do
  
  comp(stepper, earthForm, base, tempAmb, tempBase, dir, dist, height, step1, stepMax, tolAbs, tolRel, target, samples);
  return 0;
}
