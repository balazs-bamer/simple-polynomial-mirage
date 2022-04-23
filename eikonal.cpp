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

void comp(StepperType aStepper, Eikonal::EarthForm const aEarthForm, double aDir, double aDist, double aHeight, double aStep1, double aTempAmb, Eikonal::Model aMode, double aTolAbs, double aTolRel, double aTarget, uint32_t aSamples) {

  Eikonal eikonal(aTempAmb, aMode, aEarthForm);
  RungeKuttaRayBending rk(aStepper, aDist, aTolAbs, aTolRel, aStep1, eikonal);
  Vertex start(0.0, aHeight, 0.0);
  Vector dir(std::cos(aDir / 180.0 * 3.1415926539), std::sin(aDir / 180.0 * 3.1415926539), 0.0);

  std::vector<Vertex> stuff;
  std::ofstream out("values.txt");
  for(double t = 0.0; t < aTarget; t += aTarget / aSamples) {
    Vertex y;
    if(t == 0.0) {
      y[0] = 0.0;
      y[1] = aHeight;
      y[2] = 0.0;
    }
    else {
      y = rk.solve4x(start, dir, t);
    }
    stuff.push_back(y);
    out << std::setprecision(10) << y[0] << '\t' << std::setprecision(10) << y[1] << '\n';
  }
  std::cout << "x=[";
  for (int i = 0; i < stuff.size(); ++i) {
    std::cout << std::setprecision(10) << stuff[i][0] << (i < stuff.size() - 1 ? ", " : "];\n");
  }
  std::cout << "y=[";
  for (int i = 0; i < stuff.size(); ++i) {
    std::cout << std::setprecision(10) << stuff[i][1] << (i < stuff.size() - 1 ? ", " : "];\n");
  }
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
  std::string nameAsphalt = "conventional";
  opt.add_option("--asphalt", nameAsphalt, "asphalt type (conventional / porous) [conventional]");
  double dir = 0.0;
  opt.add_option("--dir", dir, "start direction (degrees), neg is down, 0 horizontal [0]");
  double dist = 2000.0;
  opt.add_option("--dist", dist, "distance along the ray to track [2000]");
  std::string nameForm = "round";
  opt.add_option("--earthForm", nameForm, "Earth form (flat / round) [round]");
  double height = 1.0;
  opt.add_option("--height", height, "start height (m) [1]");
  double samples = 100;
  opt.add_option("--samples", samples, "number of samples on ray [100]");
  double step1 = 0.01;
  opt.add_option("--step1", step1, "initial step size [0.01]");
  std::string nameStepper = "RungeKutta23";
  opt.add_option("--stepper", nameStepper, "stepper type (RungeKutta23 / RungeKuttaClass4 / RungeKuttaFehlberg45 / RungeKuttaCashKarp45 / RungeKuttaPrinceDormand89 / BulirschStoerBaderDeuflhard) [RungeKutta23]");
  double target = 1000.0;
  opt.add_option("--target", target, "horizontal distance to travel [1000]");
  double tempAmb = std::nan("");
  opt.add_option("--tempAmb", tempAmb, "ambient temperature (Celsius) [20 for conventional, 38.5 for porous]");
  double tolAbs = 0.001;
  opt.add_option("--tolAbs", tolAbs, "absolute tolerance [1e-3]");
  double tolRel = 0.001;
  opt.add_option("--tolRel", tolRel, "relative tolerance [1e-3]");
  CLI11_PARSE(opt, aArgc, aArgv);

  Eikonal::Model asphalt;
  if(nameAsphalt == "conventional") {
    asphalt = Eikonal::Model::cConventional;
  }
  else if(nameAsphalt == "porous") {
    asphalt = Eikonal::Model::cPorous;
  }
  else {
    std::cerr << "Illegal asphalt value: " << nameAsphalt << '\n';
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
    tempAmb = (asphalt == Eikonal::Model::cConventional ? 20.0 : 38.5);
  }
  else {} // nothing to do
  
  comp(stepper, earthForm, dir, dist, height, step1, tempAmb, asphalt, tolAbs, tolRel, target, samples);
  return 0;
}
