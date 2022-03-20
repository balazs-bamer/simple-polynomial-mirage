#ifndef SHEPARDINTERPOLATION
#define SHEPARDINTERPOLATION

#include <cmath>
#include <array>
#include <deque>
#include <limits>
#include <memory>
#include <vector>
#include <variant>
#include <stdexcept>
#include <algorithm>
#include <type_traits>

template<typename tValue, size_t tArrayCount, size_t tArraySize>
class FixedStack final {
private:
  struct Iterator final {
    uint32_t    mIndexArray;
    uint32_t    mIndexValue; // both point to valid
    FixedStack &mHost;
    Iterator(FixedStack &aHost, uint32_t const aIndexArray, uint32_t const aIndexValue) : mHost(aHost), mIndexArray(aIndexArray), mIndexValue(aIndexValue) {}
    Iterator(FixedStack const &aHost, uint32_t const aIndexArray, uint32_t const aIndexValue) : mHost(const_cast<FixedStack&>(aHost)), mIndexArray(aIndexArray), mIndexValue(aIndexValue) {}

    void operator++() {
      if(++mIndexValue == tArraySize) {
        ++mIndexArray;
        mIndexValue = 0u;
      }
      else {} // Nothing to do
    }

    void operator--() {
      if(mIndexValue == 0u) {
        --mIndexArray;
        mIndexValue = tArraySize - 1u;
      }
      else {
        --mIndexValue;
      }
    }

    tValue const& operator*() const { return mHost.get(mIndexArray, mIndexValue); }
    tValue const& operator->() const { return mHost.get(mIndexArray, mIndexValue); }
    tValue&       operator*() { return mHost.get(mIndexArray, mIndexValue); }
    tValue&       operator->() { return mHost.get(mIndexArray, mIndexValue); }
    bool          operator==(Iterator const &aOther) const { return &mHost == &(aOther.mHost) && mIndexArray == aOther.mIndexArray && mIndexValue == aOther.mIndexValue; }
    bool          operator!=(Iterator const &aOther) const { return !(*this == aOther); }
  };

  using Array = std::array<tValue, tArraySize>;

  std::array<std::unique_ptr<Array>, tArrayCount> mArrays;
  uint32_t                                        mIndexArray = 0u; // Points to full unused but valid one if mIndexValue == 0u
  uint32_t                                        mIndexValue = 0u; // Can't be tArraySize. Both point to next unused.

public:
  FixedStack() { mArrays[0u] = std::move(std::make_unique<Array>()); }

  Iterator begin()       { return Iterator(*this, 0u, 0u); }
  Iterator end()         { return Iterator(*this, mIndexArray, mIndexValue); }
  Iterator begin() const { return Iterator(*this, 0u, 0u); }
  Iterator end()   const { return Iterator(*this, mIndexArray, mIndexValue); }

  size_t size() const { return tArraySize * mIndexArray + mIndexValue; }

  void push_back(tValue&& aValue) {
    (*mArrays[mIndexArray])[mIndexValue] = std::move(aValue);
    if(++mIndexValue == tArraySize) {
      mArrays[++mIndexArray] = std::move(std::make_unique<Array>());
      mIndexValue = 0u;
    }
    else {} // Nothing to do
  }

  void push_back(tValue const& aValue) {
    (*mArrays[mIndexArray])[mIndexValue] = aValue;
    if(++mIndexValue == tArraySize) {
      mArrays[++mIndexArray] = std::move(std::make_unique<Array>());
      mIndexValue = 0u;
    }
    else {} // Nothing to do
  }

  tValue& back() {
    uint32_t ia = (mIndexValue == 0u ? mIndexArray - 1u : mIndexArray);
    uint32_t iv = (mIndexValue == 0u ? tArraySize - 1u : mIndexValue - 1u);
    return (*mArrays[ia])[iv];
  }

  void pop_back() {
    if(mIndexValue == 0u) {
      mArrays[mIndexArray].reset(nullptr);
      mIndexValue = tArraySize - 1u;
      --mIndexArray;
    }
    else {
      --mIndexValue;
    }
  }

  void clear() {
    std::fill(mArrays.begin() + 1u, mArrays.end(), nullptr);
    mIndexValue = mIndexArray = 0u;
  }

  bool isValid() {
    return std::all_of(mArrays.begin() + mIndexArray + 1u, mArrays.end(), [](auto const &item){ return !item; });
  }

private:
  tValue const& get(uint32_t const aIndexArray, uint32_t const& aIndexValue) const { return (*mArrays[aIndexArray])[aIndexValue]; }
  tValue&       get(uint32_t const aIndexArray, uint32_t const& aIndexValue)       { return (*mArrays[aIndexArray])[aIndexValue]; }
};

template<typename tCoordinate, size_t tSize>
class CoefficientWise final {
private:
  std::array<tCoordinate, tSize> mArray;

public:
  using Scalar = tCoordinate;

  CoefficientWise() = default;
  CoefficientWise(tCoordinate const aInit) { std::fill(mArray.begin(), mArray.end(), aInit); }

  CoefficientWise(std::initializer_list<tCoordinate> aInit) {
    if(aInit.size() != tSize) {
      throw std::invalid_argument("CoefficientWise: invalid initializer_list.");
    }
    else {} // nothing to do
    std::copy(aInit.begin(), aInit.end(), mArray.begin());
  }

  CoefficientWise(CoefficientWise const &) = default;
  CoefficientWise(CoefficientWise &&) = default;
  CoefficientWise& operator=(CoefficientWise const &) = default;
  CoefficientWise& operator=(CoefficientWise &&) = default;

  tCoordinate& operator[](uint32_t const aIndex) { return mArray[aIndex]; }
  tCoordinate const& operator[](uint32_t const aIndex) const { return mArray[aIndex]; }

  static constexpr CoefficientWise Zero() { return CoefficientWise(0.0); }

  CoefficientWise operator+(CoefficientWise const& aOther) const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = mArray[i] + aOther.mArray[i];
    }
    return result;
  }

  CoefficientWise operator-(CoefficientWise const& aOther) const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = mArray[i] - aOther.mArray[i];
    }
    return result;
  }

  CoefficientWise operator*(CoefficientWise const& aOther) const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = mArray[i] * aOther.mArray[i];
    }
    return result;
  }

  CoefficientWise operator/(CoefficientWise const& aOther) const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = mArray[i] / aOther.mArray[i];
    }
    return result;
  }

  CoefficientWise operator+(tCoordinate const aOther) const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = mArray[i] + aOther;
    }
    return result;
  }

  CoefficientWise operator-(tCoordinate const aOther) const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = mArray[i] - aOther;
    }
    return result;
  }

  CoefficientWise operator*(tCoordinate const aOther) const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = mArray[i] * aOther;
    }
    return result;
  }

  CoefficientWise operator/(tCoordinate const aOther) const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = mArray[i] / aOther;
    }
    return result;
  }

  CoefficientWise& operator+=(tCoordinate const aOther) {
    for(uint32_t i = 0u; i < tSize; ++i) {
      mArray[i] += aOther;
    }
    return *this;
  }

  CoefficientWise& operator-=(tCoordinate const aOther) {
    for(uint32_t i = 0u; i < tSize; ++i) {
      mArray[i] -= aOther;
    }
    return *this;
  }

  CoefficientWise& operator*=(tCoordinate const aOther) {
    for(uint32_t i = 0u; i < tSize; ++i) {
      mArray[i] *= aOther;
    }
    return *this;
  }

  CoefficientWise& operator/=(tCoordinate const aOther) {
    for(uint32_t i = 0u; i < tSize; ++i) {
      mArray[i] /= aOther;
    }
    return *this;
  }

  CoefficientWise& operator+=(CoefficientWise const& aOther) {
    for(uint32_t i = 0u; i < tSize; ++i) {
      mArray[i] += aOther[i];
    }
    return *this;
  }

  CoefficientWise& operator-=(CoefficientWise const& aOther) {
    for(uint32_t i = 0u; i < tSize; ++i) {
      mArray[i] -= aOther[i];
    }
    return *this;
  }

  CoefficientWise& operator*=(CoefficientWise const& aOther) {
    for(uint32_t i = 0u; i < tSize; ++i) {
      mArray[i] *= aOther[i];
    }
    return *this;
  }

  CoefficientWise& operator/=(CoefficientWise const& aOther) {
    for(uint32_t i = 0u; i < tSize; ++i) {
      mArray[i] /= aOther[i];
    }
    return *this;
  }

  tCoordinate norm() const {
    tCoordinate result = 0.0;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result += mArray[i] * mArray[i];
    }
    return ::sqrt(result);
  }

  CoefficientWise floor() const {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = ::floor(mArray[i]);
    }
    return result;
  }

  static CoefficientWise min(CoefficientWise const& aOne, CoefficientWise const& aOther) {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = std::min(aOne.mArray[i], aOther.mArray[i]);
    }
    return result;
  }

  static CoefficientWise max(CoefficientWise const& aOne, CoefficientWise const& aOther) {
    CoefficientWise result;
    for(uint32_t i = 0u; i < tSize; ++i) {
      result.mArray[i] = std::max(aOne.mArray[i], aOther.mArray[i]);
    }
    return result;
  }
};

template<class> struct sfinaeTrue : std::true_type {};

template<class tType>
static auto testMin( int) -> sfinaeTrue<decltype(std::declval<tType>().min(std::declval<tType>(), std::declval<tType>()))>;

template<class tType>
static auto testMin(long) -> std::false_type;

template<class tType>
struct hasMin : decltype(testMin<tType>(0)){};

template<class tType>
static auto testMax( int) -> sfinaeTrue<decltype(std::declval<tType>().max(std::declval<tType>(), std::declval<tType>()))>;

template<class tType>
static auto testMax(long) -> std::false_type;

template<class tType>
struct hasMax : decltype(testMax<tType>(0)){};

// TODO store the locations in one template type (possibly float) to save place and use an other one (possibly long double) for calculations

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
class ShepardInterpolation final {
public:
  using Location = CoefficientWise<tCoordinate, tDimensions>;
  struct Data {
    Location mLocation;
    tPayload mPayload;
  };
  using DataTransfer = FixedStack<Data, 16384u, 65536u>;

  static constexpr size_t getAverageCount() {
    size_t result = 1u;
    for(size_t d = 0; d < tDimensions; ++d) {
      result *= tAverageCount1d;
    }
    return result;
  }

  static constexpr size_t      csAverageCount               = getAverageCount();
  static constexpr uint32_t    csChildCount                 =  1u << tDimensions;
  static constexpr uint32_t    csMaxLevels                  =  50u;
  static constexpr tCoordinate csInflateBounds              =  1.01;
  static constexpr tCoordinate csTargetEpsilonFactor        =  0.01;
  static constexpr tCoordinate csDefaultAverageRelativeSize =  0.5;
  static constexpr tCoordinate csDefaultShepardExponent     =  3.0;
  static constexpr tCoordinate csDefaultBiasSize            = 13.0;

  static_assert(tDimensions > 0u && tDimensions <= 10u);
  static_assert(tInPlace > 1u && tInPlace <= 1024u);
  static_assert(tAverageCount1d > 0u && tAverageCount1d <= 5u && csAverageCount <= 1024u);

private:
  struct Node final {
  public:
    using Children = std::array<std::unique_ptr<Node>, csChildCount>;
    using Payload  = std::array<Data, tInPlace>;

    uint32_t                                        mCountTotal;
    uint32_t                                        mCountHere;
    std::variant<std::monostate, Children, Payload> mContents;
    Location                                        mCenter;
    Location                                        mSizeDiv4;

    Node(Location const &aCenter, Location const& aSize) : mCountTotal(0u), mCountHere(0u), mCenter(aCenter), mSizeDiv4(aSize / 4.0) {}

    std::pair<uint32_t, Location> getChildIndexCenter(Location const &aTarget) const {
      uint32_t index = 0u;
      Location location = mCenter;
      for(uint32_t i = 0u; i < tDimensions; ++i) {
        auto diff = mSizeDiv4[i];
        if(aTarget[i] >= mCenter[i]) {
          index += 1u << i;
          location[i] += diff;
        }
        else {
          location[i] -= diff;
        }
      }
      return std::pair(index, location);
    }
  };

  std::array<std::unique_ptr<Node>, csChildCount> mRoots;
  Location                                        mBoundsMin;
  Location                                        mBoundsMax;
  uint32_t const                                 cmSamplesToConsider;
  tCoordinate const                              cmShepardExponentMod;
  Location                                        mLocationScale;
  tPayload                                        mBias;
  uint32_t                                        mTargetLevelInChild0;  // Where average of total count >= cmSamplesToConsider
  Location                                        mTargetSize;
  Location                                        mTargetSizeDiv2;
  Location                                        mTargetSizeDiv4;
  tCoordinate                                     mTargetEpsilonSquared;
  std::vector<uint32_t>                           mNodesPerLevel;
  std::vector<uint32_t>                           mItemsPerLevel;
  mutable std::vector<Node*>                      mInterpolatorQueue;  // Just to avoid frequent allocations
  std::array<Location, csAverageCount>            mAverageLocations;

public:
  ShepardInterpolation() = delete;
  ShepardInterpolation(ShepardInterpolation const &) = delete;
  ShepardInterpolation(ShepardInterpolation &&) = delete;
  ShepardInterpolation& operator=(ShepardInterpolation const &) = delete;
  ShepardInterpolation& operator=(ShepardInterpolation &&) = delete;

  ShepardInterpolation(DataTransfer const &aData, uint32_t const aSamplesToConsider, tCoordinate const aAverageRelativeSize, tCoordinate const aShepardExponent, tCoordinate const aBiasSize);
  ShepardInterpolation(DataTransfer const &aData, uint32_t const aSamplesToConsider, tCoordinate const aAverageRelativeSize, tCoordinate const aShepardExponent)
    : ShepardInterpolation(aData, aSamplesToConsider, aAverageRelativeSize, aShepardExponent, csDefaultBiasSize) {}
  ShepardInterpolation(DataTransfer const &aData, uint32_t const aSamplesToConsider, tCoordinate const aAverageRelativeSize)
    : ShepardInterpolation(aData, aSamplesToConsider, aAverageRelativeSize, csDefaultShepardExponent, csDefaultBiasSize) {}
  ShepardInterpolation(DataTransfer const &aData, uint32_t const aSamplesToConsider)
    : ShepardInterpolation(aData, aSamplesToConsider, csDefaultAverageRelativeSize, csDefaultShepardExponent, csDefaultBiasSize) {}
  uint32_t getTargetLevel()                    const { return mTargetLevelInChild0; }
  uint32_t getLevelCount()                     const { return mNodesPerLevel.size(); }
  uint32_t getNodeCount(uint32_t const aLevel) const { return mNodesPerLevel[aLevel]; }
  uint32_t getItemCount(uint32_t const aLevel) const { return mItemsPerLevel[aLevel]; }

  // TODO Octave output of interpolated function.
  // TODO RRMSE for vector of test points.
  tPayload interpolate(Location const& aLocation) const;
  tCoordinate getDistanceFromTargetCenter(Location const& aLocation) const;

private:
  void                       buildTree(size_t const aWhichRoot, Location const aCenter, Location const aBoundsMax, DataTransfer const& aData);
  void                       addLeaf(Node * aBranch, Location const& aCenter, Location const& aSize, Data const& aItem, uint32_t const aLevel);
  void                       calculateTargetLevelFromChild0();
  std::pair<Node*, uint32_t> getTargetNodeLevelDiff(uint32_t const aWhichRoot, Location const& aLoc) const;
  std::pair<Node*, uint32_t> getTargetNodeLevelDiff(Location const& aLoc) const;
};

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::ShepardInterpolation(DataTransfer const &aData, uint32_t const aSamplesToConsider, tCoordinate const aAverageRelativeSize, tCoordinate const aShepardExponent, tCoordinate const aBiasSize)
  : mBoundsMin ( std::numeric_limits<tCoordinate>::max())
  , mBoundsMax (-std::numeric_limits<tCoordinate>::max())
  , cmSamplesToConsider(aSamplesToConsider)
  , cmShepardExponentMod(-0.5 * aShepardExponent)
  , mTargetLevelInChild0(0u) {
  if(aData.size() < aSamplesToConsider || aSamplesToConsider == 0u) {
    throw std::invalid_argument("ShepardInterpolation: invalid constructor arguments.");
  }
  else {} // nothing to do

  tPayload valueMin;
  tPayload valueMax;
  if constexpr(std::is_fundamental_v<tPayload>) {
    valueMin = std::numeric_limits<tPayload>::max();
    valueMax = -std::numeric_limits<tPayload>::max();
  }
  else {
    valueMin = std::numeric_limits<typename tPayload::Scalar>::max();
    valueMax = -std::numeric_limits<typename tPayload::Scalar>::max();
  }
  Location means = Location::Zero();
  for(auto const &item : aData) {
    if constexpr(hasMin<tPayload>() && hasMax<tPayload>()) {
      valueMin = tPayload::min(valueMin, item.mPayload);
      valueMax = tPayload::max(valueMax, item.mPayload);
    }
    else {
      valueMin = std::min(valueMin, item.mPayload);
      valueMax = std::max(valueMax, item.mPayload);
    }
    // TODO do for 1d?
    means += item.mLocation;   // TODO perhaps compensating sum
  }

  means /= aData.size();
  Location variances = Location::Zero();
  for(auto const &item : aData) {
    auto diff = item.mLocation - means;
    variances += diff * diff;            // TODO perhaps compensating sum
  }
  variances /= aData.size() - 1u;
  for(uint32_t d = 0u; d < tDimensions; ++d) {
    mLocationScale[d] = 1.0 / std::sqrt(variances[d]);
  }

  for(auto const &item : aData) {
    mBoundsMin = Location::min(mBoundsMin, item.mLocation * mLocationScale);
    mBoundsMax = Location::max(mBoundsMax, item.mLocation * mLocationScale);
  }
  auto size = (mBoundsMax - mBoundsMin) * csInflateBounds;
  auto center = (mBoundsMax + mBoundsMin) / 2.0;
  mBoundsMin = center - size / 2.0;
  mBoundsMax = center + size / 2.0;
  mBias = (valueMax - valueMin) * aBiasSize - valueMin;

  buildTree(0u, (mBoundsMin + mBoundsMax) / 2u, size, aData);
  calculateTargetLevelFromChild0();
  mTargetSize = size / ::pow(2.0, mTargetLevelInChild0);
  mTargetSizeDiv2 = mTargetSize / 2.0;
  mTargetSizeDiv4 = mTargetSize / 4.0;
  auto tmp = mTargetSize.norm() * csTargetEpsilonFactor;
  mTargetEpsilonSquared = tmp * tmp;

  for(uint32_t i = 1u; i < csChildCount; ++i) {
    auto newBoundsMax = mBoundsMax;
    auto newBoundsMin = mBoundsMin;
    for(uint32_t dim = 0u; dim < tDimensions; ++dim) {
      if(i && (1u << dim)) {
        newBoundsMin[dim] -= mTargetSizeDiv2[dim];
        newBoundsMax[dim] += size[dim] - mTargetSizeDiv2[dim];
      }
      else {
        newBoundsMax[dim] += size[dim];
      }
    }
    buildTree(i, (newBoundsMin + newBoundsMax) / 2u, newBoundsMax - newBoundsMin, aData);
  }

  if constexpr(tAverageCount1d > 1u) {
    auto increment = mTargetSizeDiv2 * (aAverageRelativeSize / (tAverageCount1d - 1u));
    auto start = mTargetSizeDiv4 * -aAverageRelativeSize;
    for(size_t i = 0u; i < csAverageCount; ++i) {
      auto soFar = i;
      for(uint32_t d = 0u; d < tDimensions; ++d) {
        auto where = soFar % tAverageCount1d;
        mAverageLocations[i][d] = start[d] + increment[d] * where;
        soFar /= tAverageCount1d;
      }
    }
  }
  else {
    mAverageLocations[0] = Location::Zero();
  }
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
tPayload ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::interpolate(Location const &aLocation) const {
  auto location = aLocation * mLocationScale;
  auto[branch, levelDiff] = getTargetNodeLevelDiff(location);
  mInterpolatorQueue.emplace_back(branch);
  std::array<bool,        csAverageCount> readys;
  std::array<tPayload,    csAverageCount> results;
  std::array<tPayload,    csAverageCount> sampleSums;
  std::array<tCoordinate, csAverageCount> weightSums;
  std::fill(sampleSums.begin(), sampleSums.end(), tPayload::Zero());
  std::fill(weightSums.begin(), weightSums.end(), 0.0);
  std::fill(readys.begin(),     readys.end(),    false);
  while(mInterpolatorQueue.size() > 0u) {
    branch = mInterpolatorQueue.back();
    mInterpolatorQueue.pop_back();
    if(std::holds_alternative<typename Node::Children>(branch->mContents)) {
      auto &children = std::get<typename Node::Children>(branch->mContents);
      for(auto const &child : children) {
        if(child) {
          mInterpolatorQueue.emplace_back(child.get());
        }
        else {} // Nothing to do
      }
    }
    else if(std::holds_alternative<typename Node::Payload>(branch->mContents)) {
      auto &payload = std::get<typename Node::Payload>(branch->mContents);
      for(uint32_t s = 0u; s < branch->mCountHere; ++s) {
        auto const& sample = payload[s];
        for(uint32_t a = 0u; a < csAverageCount; ++a) {
          tCoordinate diff = 0.0;
          auto const& disp = mAverageLocations[a];
          for(uint32_t d = 0u; d < tDimensions; ++d) {
            auto diffScal = sample.mLocation[d] - (location[d] + disp[d]);
            diff += diffScal * diffScal;
          }
          if(diff < mTargetEpsilonSquared) {
            results[a] = sample.mPayload;
            readys[a] = true;
          }
          else {
            if(!readys[a]) {
              auto weight = ::pow(diff, cmShepardExponentMod);
              weightSums[a] += weight;
              sampleSums[a] += sample.mPayload * weight;
            }
            else {} // nothing to do
          }
        }
      }
    }
  }
  tPayload result = tPayload::Zero();
  for(uint32_t a = 0u; a < csAverageCount; ++a) {
    if(readys[a]) {
      result += results[a];
    }
    else {
      result += sampleSums[a] / weightSums[a];
    }
  }
  return result / csAverageCount - mBias;
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
tCoordinate ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::getDistanceFromTargetCenter(Location const& aLocation) const {
  tCoordinate result;
  auto[branch, levelDiff] = getTargetNodeLevelDiff(aLocation);
  if(levelDiff > 0u) {
    result = 0.0;
  }
  else {
    result = (branch->mCenter - aLocation).norm() / mTargetSize.norm();
  }
  return result;
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
void ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::buildTree(size_t const aWhichRoot, Location const aCenter, Location const aSize, DataTransfer const& aData) {
  mRoots[aWhichRoot] = std::move(std::make_unique<Node>(aCenter, aSize));
  auto root = mRoots[aWhichRoot].get();
  for(auto const &item : aData) {
    auto biased = item;
    biased.mLocation *= mLocationScale;
    biased.mPayload += mBias;
    addLeaf(root, aCenter, aSize, biased, 0u);
  }
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
void ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::addLeaf(Node * aBranch, Location const& aCenter, Location const& aSize, Data const& aItem, uint32_t const aLevel) {
  Node * branch   = aBranch;
  Location center = aCenter;
  Location size   = aSize;
  auto level      = aLevel;
  while(true) {
    if(level > csMaxLevels) {
      throw std::invalid_argument("ShepardInterpolation: maximum tree levels reached.");
    }
    else {} // nothing to do
    ++branch->mCountTotal;
    if(branch->mCountTotal <= tInPlace) {
      if(branch->mCountHere == 0u) {
        branch->mContents = typename Node::Payload();
      }
      else {} // nothing to do
      std::get<typename Node::Payload>(branch->mContents)[branch->mCountHere] = aItem;
      ++branch->mCountHere;
      break;
    }
    else {
      size /= 2.0;
      if(branch->mCountHere == tInPlace) {
        branch->mCountHere = 0u;
        typename Node::Payload contentsSoFar = std::get<typename Node::Payload>(branch->mContents);
        branch->mContents = typename Node::Children();
        auto& children = std::get<typename Node::Children>(branch->mContents);
        for(auto const &item : contentsSoFar) {
          auto [childIndex, childCenter] = branch->getChildIndexCenter(item.mLocation);
          if(!children[childIndex]) {
            children[childIndex] = std::move(std::make_unique<Node>(childCenter, size));
          }
          else {} // nothing to do
          addLeaf(children[childIndex].get(), childCenter, size, item, level + 1u);  // Recursive, but has only a depth of 1
        }
      }
      else {} // Nothing to do
      uint32_t childIndex;
      std::tie(childIndex, center) = branch->getChildIndexCenter(aItem.mLocation);
      auto& children = std::get<typename Node::Children>(branch->mContents);
      if(!children[childIndex]) {
        children[childIndex] = std::move(std::make_unique<Node>(center, size));
      }
      else {} // nothing to do
      branch = children[childIndex].get();
      ++level;
    }
  }
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
void ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::calculateTargetLevelFromChild0() {
  struct Item {
    Node*    mNode;
    uint32_t mLevel;
    Item(Node * const aNode, uint32_t const aLevel) : mNode(aNode), mLevel(aLevel) {}
  };
  std::deque<Item> queue;
  uint32_t level = 0u;
  mTargetLevelInChild0 = std::numeric_limits<uint32_t>::max();
  auto root = mRoots[0].get();
  queue.emplace_back(root, level);
  while(queue.size() > 0u) {
    tCoordinate average = 0.0;
    int32_t count = 0u;
    auto item = queue.front();
    mNodesPerLevel.push_back(0u);
    mItemsPerLevel.push_back(0u);
    while(item.mLevel == level) {
      ++mNodesPerLevel.back();
      mItemsPerLevel.back() += item.mNode->mCountHere;
      average += item.mNode->mCountTotal;
      ++count;
      if(std::holds_alternative<typename Node::Children>(item.mNode->mContents)) {
        auto &children = std::get<typename Node::Children>(item.mNode->mContents);
        for(auto const &child : children) {
          if(child) {
            queue.emplace_back(child.get(), level + 1u);
          }
          else {} // Nothing to do
        }
      }
      else {} // Nothing to do
      queue.pop_front();
      if(queue.size() > 0u) {
        item = queue.front();
      }
      else {
        break;
      }
    }
    average /= count;
    if(average >= cmSamplesToConsider) {
      mTargetLevelInChild0 = level;
    }
    else {} // Nothing to do
    ++level;
  }
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
std::pair<typename ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::Node*, uint32_t>
ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::getTargetNodeLevelDiff(uint32_t const aWhichRoot, Location const& aLoc) const {
  uint32_t level = (aWhichRoot == 0u ? 0u : 1u);
  uint32_t actualTargetLevel = mTargetLevelInChild0 + level;
  auto branch = mRoots[aWhichRoot].get();
  uint32_t levelDiff = actualTargetLevel;
  Node* result = branch;
  while(branch != nullptr && level < actualTargetLevel) {
    if(branch->mCountTotal >= cmSamplesToConsider) {
      result = branch;
      levelDiff = actualTargetLevel - level;
    }
    else{} // nothing to do
    if(std::holds_alternative<typename Node::Children>(branch->mContents)) {
      auto const& children = std::get<typename Node::Children>(branch->mContents);
      auto [childIndex, childCenter] = branch->getChildIndexCenter(aLoc);
      branch = children[childIndex].get();
      ++level;
    }
    else {
      branch = nullptr;
    }
  }
  return std::pair(result, levelDiff);
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace, size_t tAverageCount1d>
std::pair<typename ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::Node*, uint32_t>
ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace, tAverageCount1d>::getTargetNodeLevelDiff(Location const& aLocation) const {
  Location fromCenter;
  for(uint32_t i = 0u; i < tDimensions; ++i) {
    auto fromMin = aLocation[i] - mBoundsMin[i];
    fromCenter[i] = fromMin - ::floor(fromMin / mTargetSize[i]) * mTargetSize[i] - mTargetSizeDiv2[i]; // TODO fails if space size is 0 in some dimension
  }
  uint32_t index = 0u;
  for(uint32_t i = 0u; i < tDimensions; ++i) {
    if(std::abs(fromCenter[i]) >= mTargetSizeDiv4[i]) {
      index += 1u << i;
    }
    else {} // nothing to do
  }
  return getTargetNodeLevelDiff(index, aLocation);
}

#endif  // SHEPARD_H
