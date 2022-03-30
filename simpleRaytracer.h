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
  Plane  mPlane;
  double const mX;

public:
  Object(char const * const aName, double const aCenterX, double const aCenterY, double const aWidth, double const aHeight);
  double  getX() const { return mX; }
  uint8_t getPixel(Vertex const &aHit) const;
};

class Medium final {
private:
  Eikonal mEikonal(aTempAmb, aTempDiff);
  RungeKuttaRayBending mSolver(aDist, aTolAbs, aTolRel, aStep1, aStepMin, eikonal);
  Object const&     mObject;

public:
  Medium(double const aTempAmb, double const aTempDiff,
         double const aDistAlongRay, double const aTolAbs, double const aTolRel, double const aStep1, double const aStepMin)
  : Eikonal mEikonal(aTempAmb, aTempDiff)
  , mSolver(aDist, aTolAbs, aTolRel, aStep1, aStepMin, mEikonal) {}

  uint8_t trace(Ray const& aRay);
};

class Image final {
private:
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
  Image(double const aCenterX, double const aCenterY,
        double const aTilt, double const aPinholeDist,
        double const aPixelSize,
        uint32_t const aResZ, uint32_t const aResY,
        uint32_t const aSubSample, Medium &aMedium);

  void process(char const * const aName);
};

#endif
