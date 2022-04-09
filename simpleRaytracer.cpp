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
  , mMaxZ(aWidth / 2.0)
  , mX(aCenterX) {
std::cout << mMinZ << ' ' << mMaxZ << ' ' << mMinY << ' ' << mMaxY << '\n';
}

uint8_t Object::getPixel(Vertex const &aHit) const {
  uint8_t result = 0u;
  if(aHit(1) > mMinY && aHit(1)  < mMaxY && aHit(2) > mMinZ && aHit(2) < mMaxZ) {
    uint32_t x = std::max(0.0, (aHit(2) - mMinZ) / mDz);
    uint32_t y = mImage.get_height() - (aHit(1) - mMinY) / mDy - 1u;
    result = mImage.get_pixel(x, y);
  }
  else {} // Nothing to do
  return result;
}


bool eq(double a, double b) { return std::abs(a-b)<1e-7; }

uint8_t Medium::trace(Ray const& aRay) {
  if(eq(aRay.mStart(0),4.0) && eq(aRay.mStart(1),1.1) && eq(aRay.mStart(2),0)
  && eq(aRay.mDirection(0),0.999891) && eq(aRay.mDirection(1),-0.00783665) && eq(aRay.mDirection(2),0.0124861)) {
std::cout<<"itt\n";
}
  try {
    auto hit = mSolver.solve4x(aRay.mStart, aRay.mDirection, mObject.getX());
    return mObject.getPixel(hit);
  }
  catch(...) {
std::cout << aRay.mStart(0) << ' ' << aRay.mStart(1) << ' ' << aRay.mStart(2) << ' '
          << aRay.mDirection(0) << ' ' << aRay.mDirection(1) << ' ' << aRay.mDirection(2) << '\n';
    throw 0u;
  }
}


Image::Image(double const aCenterY,
        double const aTilt, double const aPinholeDist,
        double const aPixelSize,
        uint32_t const aResZ, uint32_t const aResY,
        uint32_t const aSubSample, Medium &aMedium)
  : mImage(aResZ, aResY)
  , mSubSample(aSubSample)
  , mSsFactor(1.0 / aSubSample)
  , mPixelSize(aPixelSize)
  , mCenter(0.0, aCenterY, 0.0)
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
