#include "simpleRaytracer.h"
#include "CLI11.hpp"
#include <iostream>
#include <thread>


int main(int aArgc, char **aArgv) {
  CLI::App opt{"Usage"};
  std::string nameBase = "water";
  opt.add_option("--base", nameBase, "base type (conventional / porous / water) [water]");
  double bullLift = 0.0;
  opt.add_option("--bullLift", bullLift, "lift of bulletin from ground (m) [0.0]");
  double camCenter = 1.1;
  opt.add_option("--camCenter", camCenter, "height of camera center (m) [1.1]");
  double dist = 1000.0;
  opt.add_option("--dist", dist, "distance of bulletin and camera [1000]");
  std::string nameForm = "round";
  opt.add_option("--earthForm", nameForm, "Earth form (flat / round) [round]");
  double rawRadius = 6371.0;
  opt.add_option("--earthRadius", rawRadius, "Earth radius (km) [6371.0]");
  double height = 9.0;
  opt.add_option("--height", height, "height of bulletin (m) [9.0]  its width will be calculated");
  std::string nameIn = "monoscopeRca.png";
  opt.add_option("--nameIn", nameIn, "input filename [monoscopeRca.png]");
  std::string nameOut = "result.png";
  opt.add_option("--nameOut", nameOut, "output filename [result.png]");
  double pinholeDist = 4.0;
  opt.add_option("--pinholeDist", pinholeDist, "pinhole distance from film (m) [4.0]");
  uint32_t resolution = 1000u;
  opt.add_option("--resolution", resolution, "film resulution in both directions (pixel) [1000]");
  uint32_t saveCpus = 0u;
  opt.add_option("--saveCpus", saveCpus, "amount of CPUs to save to keep the system responsive (natural integer) [0]");
  bool silent = false;
  opt.add_option("--silent", silent, "surpress parameter echo (true, false) [false]");
  double step1 = 0.01;
  opt.add_option("--step1", step1, "initial step size (m) [0.01]");
  double stepMin = 1e-4;
  opt.add_option("--stepMin", stepMin, "maximal step size (m) [1e-4]");
  double stepMax = 55.5;
  opt.add_option("--stepMax", stepMax, "maximal step size (m) [55.5]");
  std::string nameStepper = "RungeKuttaFehlberg45";
  opt.add_option("--stepper", nameStepper, "stepper type (RungeKutta23 / RungeKuttaClass4 / RungeKuttaFehlberg45 / RungeKuttaCashKarp45 / RungeKuttaPrinceDormand89 / BulirschStoerBaderDeuflhard) [RungeKuttaFehlberg45]");
  uint32_t subsample = 2u;
  opt.add_option("--subsample", subsample, "subsampling each pixel in both directions (count) [2]");
  double tempAmb = std::nan("");
  opt.add_option("--tempAmb", tempAmb, "ambient temperature (Celsius) [20 for conventional, 38.5 for porous, 10 for water]");
  double tempBase = 13.0;
  opt.add_option("--tempBase", tempBase, "base temperature, only for water (Celsius) [13]");
  double tilt = 0.0;
  opt.add_option("--tilt", tilt, "camera tilt, neg downwards (degrees) [0.0]");
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

  if(!silent) {
    std::cout << "base type:                        .  .  .  .  .  . " << nameBase << ' ' << static_cast<int>(base) << '\n';
    std::cout << "lift of bulletin from ground (m):                  " << bullLift << '\n';
    std::cout << "height of camera center (m):                       " << camCenter << '\n';
    std::cout << "distance of bulletin and camera (m): .  .  .  .  . " << dist << '\n';
    std::cout << "Earth form:                                        " << nameForm << ' ' << static_cast<int>(earthForm) << '\n';
    std::cout << "Earth radius (km):                                 " << earthRadius / 1000.0 << '\n';
    std::cout << "height of bulletin (m):    .   .  .  .  .  .  .  . " << height << '\n';
    std::cout << "input filename:                                    " << nameIn << '\n';
    std::cout << "output filename:                                   " << nameOut << '\n';
    std::cout << "pinhole distance from film (m):            .  .  . " << pinholeDist << '\n';
    std::cout << "film resolution in both directions (pixel):        " << resolution << '\n';
    std::cout << "initial step size (m):                             " << step1 << '\n';
    std::cout << "minimal step size (m):   .  .  .  .  .  .  .  .  . " << stepMin << '\n';
    std::cout << "maximal step size (m):                             " << stepMax << '\n';
    std::cout << "stepper type:                                      " << nameStepper << ' ' << static_cast<int>(stepper) << '\n';
    std::cout << "subsampling each pixel in both directions (count): " << subsample << '\n';
    std::cout << "ambient temperature (Celsius):                     " << tempAmb << '\n';
    std::cout << "base temperature, only for water (Celsius):        " << tempBase << '\n';
    std::cout << "camera tilt, neg downwards (degrees):      .  .  . " << tilt << '\n';
    std::cout << "absolute tolerance (m):                            " << tolAbs << '\n';
    std::cout << "relative tolerance (m):                            " << tolRel << '\n';
    uint32_t nCpus = std::thread::hardware_concurrency();
    nCpus -= (nCpus <= saveCpus ? nCpus - 1u : saveCpus);
    std::cout << "Using " << nCpus << " thread(s)" << std::endl;
  }
  else {} // nothing to do

  RungeKuttaRayBending::Parameters parameters;
  parameters.mStepper      = stepper;
  parameters.mDistAlongRay = dist * 2.0;
  parameters.mTolAbs       = tolAbs;
  parameters.mTolRel       = tolRel;
  parameters.mStep1        = step1;
  parameters.mStepMin      = stepMin;
  parameters.mStepMax      = stepMax;
  Object object(nameIn.c_str(), dist, bullLift, height);
  Medium medium(parameters, earthForm, earthRadius, base, tempAmb, tempBase, object);
  Image image(saveCpus, camCenter, tilt, pinholeDist, 0.1 / resolution, resolution, resolution, subsample, medium);
  image.process(nameOut.c_str());
  return 0;
}
