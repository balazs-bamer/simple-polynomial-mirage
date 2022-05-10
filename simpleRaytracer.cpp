#include "simpleRaytracer.h"
#include <iomanip>
#include <iostream>
#include <thread>


Object::Object(char const * const aName, double const aDispX, double const aLiftY, double const aHeight)
  : mImage(aName)
  , mDy(aHeight / mImage.get_height())
  , mDz(mDy)
  , mMinY(aLiftY)
  , mMaxY(aLiftY + aHeight)
  , mMinZ(-static_cast<double>(mImage.get_width()) * aHeight / static_cast<double>(mImage.get_height()) / 2.0)
  , mMaxZ(-mMinZ)
  , mX(aDispX) {
std::cout << mMinZ << ' ' << mMaxZ << ' ' << mMinY << ' ' << mMaxY << '\n';
}

bool Object::hasPixel(Vertex const &aHit) const {
  return aHit(1) > mMinY && aHit(1)  < mMaxY && aHit(2) > mMinZ && aHit(2) < mMaxZ;
}

uint8_t Object::getPixel(Vertex const &aHit) const {
  uint8_t result = 0u;
  if(hasPixel(aHit)) {
    uint32_t x = std::max(0.0, (aHit(2) - mMinZ) / mDz);
    uint32_t y = mImage.get_height() - (aHit(1) - mMinY) / mDy - 1u;
    result = mImage.get_pixel(x, y);
  }
  else {} // Nothing to do
  return result;
}


bool eq(double a, double b) { return std::abs(a-b)<1e-7; }

uint8_t Medium::trace(Ray const& aRay) {
  try {
    auto hit = mSolver.solve4x(aRay.mStart, aRay.mDirection, mObject.getX());
    if(hit.mValid) {
      return mObject.getPixel(hit.mValue);
    }
    else {
      return 0;
    }
  }
  catch(...) {
std::cout << aRay.mStart(0) << ' ' << aRay.mStart(1) << ' ' << aRay.mStart(2) << ' '
          << aRay.mDirection(0) << ' ' << aRay.mDirection(1) << ' ' << aRay.mDirection(2) << '\n';
    throw 0;
  }
}

bool Medium::hits(Ray const& aRay) {
  try {
    auto hit = mSolver.solve4x(aRay.mStart, aRay.mDirection, mObject.getX());
    return hit.mValid && mObject.hasPixel(hit.mValue);
  }
  catch(...) {
std::cout << aRay.mStart(0) << ' ' << aRay.mStart(1) << ' ' << aRay.mStart(2) << ' '
          << aRay.mDirection(0) << ' ' << aRay.mDirection(1) << ' ' << aRay.mDirection(2) << '\n';
    return false;
  }
}


Image::Image(uint32_t const aRestrictCpu, double const aCenterY,
        double const aTilt, double const aPinholeDist,
        double const aPixelSize,
        uint32_t const aResZ, uint32_t const aResY,
        uint32_t const aSubSample, Medium &aMedium)
  : mRestrictCpu(aRestrictCpu)
  , mBuffer(aResZ * aResY)
  , mImage(aResZ, aResY)
  , mSubSample(aSubSample)
  , mSsFactor(1.0 / aSubSample)
  , mPixelSize(aPixelSize)
  , mCenter(0.0, aCenterY, 0.0)
  , mNormal(::cos(aTilt * cgPi / 180.0), ::sin(aTilt * cgPi / 180.0), 0.0)
  , mInPlaneZ(0.0, 0.0, 1.0)
  , mInPlaneY(mNormal.cross(mInPlaneZ))
  , mPinhole(mCenter + aPinholeDist * mNormal)
  , mBiasZ((aResZ - 1.0) / 2.0)
  , mBiasY((aResY - 1.0) / 2.0)
  , mBiasSub((aSubSample - 1.0) / 2.0)
  , mMedium(aMedium) {}

void Image::dumpLimits() const {
  Ray ray;
  ray.mStart = mPinhole;
  ray.mDirection = getDirectionInXy(csLimitLow - csLimitDelta);
  auto lastHit = mMedium.hits(ray);
  for(auto angle = csLimitLow; angle <= csLimitHigh; angle += csLimitDelta) {
    ray.mDirection = getDirectionInXy(angle);
    auto thisHit = mMedium.hits(ray);
    if(lastHit != thisHit) {
      auto critical = binarySearch(angle - csLimitDelta, angle, csLimitEpsilon, [this, &ray](auto const search){
        ray.mDirection = getDirectionInXy(search);
        return mMedium.hits(ray);
      });
      critical += (thisHit ? csLimitEpsilon : 0.0);
      std::cout << (thisHit ? "enter: " : "leave: ") << std::setprecision(10) << (critical * 180.0 / cgPi) << '\n';
    }
    else {} // nothing to do
    lastHit = thisHit;
  }
}

void Image::process(char const * const aName) {
  uint32_t nCpus = std::thread::hardware_concurrency();
  nCpus -= (nCpus <= mRestrictCpu ? nCpus - 1u : mRestrictCpu);
  std::vector<std::thread> threads(nCpus);
  for (uint32_t i = 0u; i < nCpus; ++i) {
    threads[i] = std::thread([this, nCpus, i] {
      Ray ray;
      ray.mStart = mPinhole;
      Medium localMedium(mMedium);
      auto yBegin = i * mImage.get_height() / nCpus;
      auto yEnd = (i + 1u) * mImage.get_height() / nCpus;
      for(int y = yBegin; y < yEnd; ++y) {
        for(int z = 0; z < mImage.get_width(); ++z) {
          double sum = 0.0;
          for(uint32_t i = 0; i < mSubSample; ++i) {
            for(uint32_t j = 0; j < mSubSample; ++j) {
              Vertex subpixel = mCenter + mPixelSize * (
                    (z - mBiasZ + mSsFactor * (i - mBiasSub)) * mInPlaneZ +
                    (y - mBiasY + mSsFactor * (j - mBiasSub)) * mInPlaneY);
              ray.mDirection = (mPinhole - subpixel).normalized();
              sum += localMedium.trace(ray);
            }
          }
          mBuffer[(mImage.get_width() - z - 1u) + mImage.get_width() * (mImage.get_height() - y - 1u)] = ::round(sum / static_cast<double>(mSubSample * mSubSample));
        }
      }
    });
  }
  for (auto& t : threads) {
    t.join();
  }
  for(int y = 0; y < mImage.get_height(); ++y) {
    for(int z = 0; z < mImage.get_width(); ++z) {
      mImage.set_pixel(z, y, mBuffer[y * mImage.get_width() + z]);
    }
  }
  mImage.write(aName);
}
