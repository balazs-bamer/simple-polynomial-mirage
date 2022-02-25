#include "simpleRaytracer.h"


#include <iostream> // TODO remove
#include <iomanip> // TODO remove

Object::Object(char const * const aName, double const aCenterX, double const aCenterY, double const aWidth, double const aHeight)
  : mImage(aName)
  , mDy(aHeight / mImage.get_height())
  , mDz(aWidth / mImage.get_width())
  , mMinY(aCenterY - aHeight / 2.0)
  , mMaxY(aCenterY + aHeight / 2.0)
  , mMinZ(-aWidth / 2.0)
  , mMaxZ(aWidth / 2.0) {
  mPlane = Plane::createFrom2vectors1point({0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}, {aCenterX, 0.0f, 0.0f});
}

bool gNeed;
double gIncl;
double gDepth;
double gBottom;
double gNewStart;

uint8_t Object::getPixel(Ray const &aRay) const {
  uint8_t result = 0u;
  auto intersection = mPlane.intersect(aRay);
  if(intersection.mValid && intersection.mPoint(1) > mMinY && intersection.mPoint(1) < mMaxY && intersection.mPoint(2) > mMinZ && intersection.mPoint(2) < mMaxZ) {
    uint32_t x = /*mImage.get_width() -*/ (intersection.mPoint(2) - mMinZ) / mDz/* - 1u*/;
    uint32_t y = mImage.get_height() - (intersection.mPoint(1) - mMinY) / mDy - 1u;
    result = mImage.get_pixel(x, y);
if(gNeed) {
/*  std::cout << " inc:  " << std::setw(17) << gIncl
            << " dpt:  " << std::setw(17) << gDepth
            << " bot:  " << std::setw(17) << gBottom
            << " new: " << std::setw(17) << gNewStart
            << " x: " << std::setw(5) << x
            << " y: " << std::setw(5) << y
            << " p: " << (int)(result) << '\n';*/
}
  }
  else {} // Nothing to do
  return result;
}

uint8_t Medium::trace(Ray const& aRay) {
  uint8_t result = 0u;
  Plane horizontal;
  horizontal.mNormal = Vector{0.0f, 1.0f, 0.0f};
  horizontal.mConstant = mHotPlate.getWorkingHeight();
  auto intersection = horizontal.intersect(aRay);
  double inclination;
gNeed = false;
  if(intersection.mValid) {
     if(aRay.mDirection(1) < 0.0f   // Simplest case: must point downwards.
        && (inclination = ::acos(static_cast<double>(-aRay.mDirection.dot(horizontal.mNormal)))) >= mHotPlate.getCriticalInclination()) {
gIncl = inclination;
       Ray inner(intersection.mPoint, aRay.mDirection);  // Virtual ray downwards inside the bending air.
       horizontal.mConstant = 0; // TODO mHotPlate.approximateReflectionDepth(inclination);
gDepth = horizontal.mConstant;
       intersection = horizontal.intersect(inner);
gBottom = intersection.mPoint(1);
       inner.mStart = intersection.mPoint;
       inner.mDirection(1) = -inner.mDirection(1);       // Virtual reflection, ray still virtual, still inside the bending air.
       horizontal.mConstant = mHotPlate.getWorkingHeight();
       intersection = horizontal.intersect(inner);
gNewStart = intersection.mPoint(1);
       inner.mStart = intersection.mPoint;               // Normal ray in homogeneous air.
gNeed = true;
       result = mObject.getPixel(inner);
     }
     else {} // Nothing to do
  }
  else {
    result = mObject.getPixel(aRay);
  }
  return result;
}


Image::Image(double const aCenterX, double const aCenterY,
        double const aTilt, double const aPinholeDist,
        double const aPixelSize,
        uint32_t const aResZ, uint32_t const aResY,
        uint32_t const aSubSample, Medium &aMedium)
  : mImage(aResZ, aResY)
  , mSubSample(aSubSample)
  , mSsFactor(1.0 / aSubSample)
  , mPixelSize(aPixelSize)
  , mCenter(aCenterX, aCenterY, 0.0)
  , mNormal(::cos(aTilt), -::sin(aTilt), 0.0)
  , mInPlaneZ(0.0, 0.0, 1.0)
  , mInPlaneY(mNormal.cross(mInPlaneZ))
  , mPinhole(mCenter + aPinholeDist * mNormal)
  , mBiasZ((aResZ - 1.0) / 2.0)
  , mBiasY((aResY - 1.0) / 2.0)
  , mBiasSub((aSubSample - 1.0) / 2.0)
  , mMedium(aMedium) {}

void Image::process(char const * const aName) {
  Ray ray;
  ray.mStart = mPinhole;
  for(int y = 0; y < mImage.get_height(); ++y) {
    for(int z = 0; z < mImage.get_width(); ++z) {
      double sum = 0.0;
      for(uint32_t i = 0; i < mSubSample; ++i) {
        for(uint32_t j = 0; j < mSubSample; ++j) {
          Vertex subpixel = mCenter + mPixelSize * (
                (z - mBiasZ + mSsFactor * (i - mBiasSub)) * mInPlaneZ +
                (y - mBiasY + mSsFactor * (j - mBiasSub)) * mInPlaneY);
          ray.mDirection = (mPinhole - subpixel).normalized();
          sum += mMedium.trace(ray);
        }
      }
      mImage.set_pixel(mImage.get_width() - z - 1u, mImage.get_height() - y - 1u, ::round(sum / static_cast<double>(mSubSample * mSubSample)));
    }
  }
  mImage.write(aName);
}
