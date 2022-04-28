#include "simpleRaytracer.h"
#include "CLI11.hpp"
#include <iostream>


int main(int aArgc, char **aArgv) {
  CLI::App opt{"Usage"};
  std::string nameBase = "water";
  opt.add_option("--base", nameBase, "base type (conventional / porous / water) [water]");
  double bullLift = 0.5;
  opt.add_option("--bullLift", bullLift, "lift of bulletin from ground (m) [0.5]");
  double camCenter = 1.1;
  opt.add_option("--camCenter", camCenter, "height of camera center (m) [1.1]");
  double dist = 1000.0;
  opt.add_option("--dist", dist, "distance of bulletin and camera [1000]");
  std::string nameForm = "round";
  opt.add_option("--earthForm", nameForm, "Earth form (flat / round) [round]");
  double height = 9.0;
  opt.add_option("--height", height, "height of bulletin (m) [9.0]  its width will be calculated");
  std::string nameIn = "monoscopeRca.png";
  opt.add_option("--nameIn", nameIn, "input filename [monoscopeRca.png]");
  std::string nameOut = "result.png";
  opt.add_option("--nameOut", nameOut, "output filename [result.png]");
  double pinholeDist = 4.0;
  opt.add_option("--pinholeDist", pinholeDist, "pinhole distance from film (m) [20.0]");
  uint32_t resolution = 1000u;
  opt.add_option("--resolution", resolution, "film resulution in both directions (pixel) [1000]");
  double step1 = 0.01;
  opt.add_option("--step1", step1, "initial step size [0.01]");
  std::string nameStepper = "RungeKutta23";
  opt.add_option("--stepper", nameStepper, "stepper type (RungeKutta23 / RungeKuttaClass4 / RungeKuttaFehlberg45 / RungeKuttaCashKarp45 / RungeKuttaPrinceDormand89 / BulirschStoerBaderDeuflhard) [RungeKutta23]");
  uint32_t subsample = 2u;
  opt.add_option("--subsample", subsample, "subsampling each pixel in both directions (count) [2]");
  double tempAmb = std::nan("");
  opt.add_option("--tempAmb", tempAmb, "ambient temperature (Celsius) [20 for conventional, 38.5 for porous, 10 for water]");
  double tempBase = 13.0;
  opt.add_option("--tempBase", tempAmb, "base temperature, only for water (Celsius) [13]");
  double tilt = 0.0;
  opt.add_option("--tilt", tilt, "camera tilt in degrees, neg downwards [0.0]");
  double tolAbs = 0.001;
  opt.add_option("--tolAbs", tolAbs, "absolute tolerance [1e-3]");
  double tolRel = 0.001;
  opt.add_option("--tolRel", tolRel, "relative tolerance [1e-3]");
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

  Object object(nameIn.c_str(), dist, bullLift, height);
  Medium medium(stepper, earthForm, base, tempAmb, tempBase, dist * 2.0, tolAbs, tolRel, step1, object);
  Image image(true, camCenter, tilt, pinholeDist, 0.1 / resolution, resolution, resolution, subsample, medium);
  image.process(nameOut.c_str());
  return 0;
}
