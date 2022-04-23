#ifndef SIMPLERAYTRACER_H
#define SIMPLERAYTRACER_H

#include "RungeKuttaRayBending.h"
#include "3dGeomUtil.h"
#include <png++/png.hpp>

class Object final {
private:
  png::image<png::gray_pixel> mImage;
  double const mDy;
  double const mDz;
  double const mMinY;
  double const mMaxY;
  double const mMinZ;
  double const mMaxZ;
  double const mX;

public:
  Object(char const * const aName, double const aDispX, double const aLiftY, double const aHeight);
  double  getX() const { return mX; }
  uint8_t getPixel(Vertex const &aHit) const;
};

class Medium final {
private:
  Eikonal              mEikonal;
  RungeKuttaRayBending mSolver;
  Object const&        mObject;

public:
  Medium(StepperType const aStepper, double const aTempAmb, Eikonal::Model const aModel, Eikonal::EarthForm aEarthForm,
         double const aDistAlongRay, double const aTolAbs, double const aTolRel, double const aStep1, Object const& aObject)
  : mEikonal(aTempAmb, aModel, aEarthForm)
  , mSolver(aStepper, aDistAlongRay, aTolAbs, aTolRel, aStep1, mEikonal)
  , mObject(aObject) {}

  uint8_t trace(Ray const& aRay);
};

class Image final {
private:
  bool     const  mRestrictCpu;
  png::image<png::gray_pixel> mImage;
  uint32_t const  mSubSample;
  double   const  mSsFactor;
  double   const  mPixelSize;
  Vertex   const  mCenter;
  Vector   const  mNormal;
  Vector   const  mInPlaneZ;
  Vector   const  mInPlaneY;
  Vertex   const  mPinhole;
  double   const  mBiasZ;
  double   const  mBiasY;
  double   const  mBiasSub;
  Medium         &mMedium;

public:
  Image(bool const aRestrictCpu, double const aCenterY,
        double const aTilt, double const aPinholeDist,
        double const aPixelSize,
        uint32_t const aResZ, uint32_t const aResY,
        uint32_t const aSubSample, Medium &aMedium);

  void process(char const * const aName);
};

#endif
