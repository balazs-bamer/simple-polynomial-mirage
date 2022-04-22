#include "simpleRaytracer.h"
#include "CLI11.hpp"
#include <iostream>


int main(int aArgc, char **aArgv) {
  CLI::App opt{"Usage"};
  std::string name = "conventional";
  opt.add_option("--asphalt", name, "asphalt type (conventional / porous) [conventional]");
  double bullCenter = 5.0;
  opt.add_option("--bullCenter", bullCenter, "height of bulletin center (m) [5.0]");
  double camCenter = 1.1;
  opt.add_option("--camCenter", camCenter, "height of camera center (m) [1.1]");
  double dist = 1000.0;
  opt.add_option("--dist", dist, "distance of bulletin and camera [1000]");
  double height = 9.0;
  opt.add_option("--height", height, "height of bulletin (m) [9.0]");
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
  std::string name2 = "RungeKutta23";
  opt.add_option("--stepper", name2, "stepper type (RungeKutta23 / RungeKuttaClass4 / RungeKuttaFehlberg45 / RungeKuttaCashKarp45 / RungeKuttaPrinceDormand89 / BulirschStoerBaderDeuflhard) [RungeKutta23]");
  uint32_t subsample = 2u;
  opt.add_option("--subsample", subsample, "subsampling each pixel in both directions (count) [2]");
  double tempAmb = std::nan("");
  opt.add_option("--tempAmb", tempAmb, "ambient temperature (Celsius) [20 for conventional, 38.5 for porous]");
  double tilt = 0.0;
  opt.add_option("--tilt", tilt, "camera tilt in degrees, neg downwards [0.0]");
  double tolAbs = 0.001;
  opt.add_option("--tolAbs", tolAbs, "absolute tolerance [1e-3]");
  double tolRel = 0.001;
  opt.add_option("--tolRel", tolRel, "relative tolerance [1e-3]");
  double width = 12.0;
  opt.add_option("--width", width, "width of bulletin [12.0]");
  CLI11_PARSE(opt, aArgc, aArgv);

  Eikonal::Mode asphalt;
  if(name == "conventional") {
    asphalt = Eikonal::Mode::cConventional;
  }
  else if(name == "porous") {
    asphalt = Eikonal::Mode::cPorous;
  }
  else {
    std::cerr << "Illegal asphalt value: " << name << '\n';
    return 1;
  }

  StepperType stepper;
  if(name2 == "RungeKutta23") {
    stepper = StepperType::cRungeKutta23;
  }
  else if(name2 == "RungeKuttaClass4") {
    stepper = StepperType::cRungeKuttaClass4;
  }
  else if(name2 == "RungeKuttaFehlberg45") {
    stepper = StepperType::cRungeKuttaFehlberg45;
  }
  else if(name2 == "RungeKuttaCashKarp45") {
    stepper = StepperType::cRungeKuttaCashKarp45;
  }
  else if(name2 == "RungeKuttaPrinceDormand89") {
    stepper = StepperType::cRungeKuttaPrinceDormand89;
  }
  else if(name2 == "BulirschStoerBaderDeuflhard") {
    stepper = StepperType::cBulirschStoerBaderDeuflhard;
  }
  else {
    std::cerr << "Illegal stepper value: " << name2 << '\n';
    return 1;
  }

  if(std::isnan(tempAmb)) {
    tempAmb = (asphalt == Eikonal::Mode::cConventional ? 20.0 : 38.5);
  }
  else {} // nothing to do

  Object object(nameIn.c_str(), dist, bullCenter, width, height);
  Medium medium(stepper, tempAmb, asphalt, dist * 2.0, tolAbs, tolRel, step1, object);
  Image image(true, camCenter, tilt, pinholeDist, 0.1 / resolution, resolution, resolution, subsample, medium);
  image.process(nameOut.c_str());
  return 0;
}
