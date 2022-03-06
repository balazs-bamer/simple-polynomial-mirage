#ifndef SHEPARDINTERPOLATION
#define SHEPARDINTERPOLATION

#include <array>
#include <deque>
#include <limits>
#include <memory>
#include <vector>
#include <variant>
#include <stdexcept>

#include <iostream> // TODO remove


template<typename tCoordinate, size_t tSize>
class CoefficientWise final {
private:
  std::array<tCoordinate, tSize> mArray;

public:
  CoefficientWise() = default;
  CoefficientWise(tCoordinate const aInit) { std::fill(mArray.begin(), mArray.end(), aInit); }

  CoefficientWise(CoefficientWise const &) = default;
  CoefficientWise(CoefficientWise &&) = default;
  CoefficientWise& operator=(CoefficientWise const &) = default;
  CoefficientWise& operator=(CoefficientWise &&) = default;

  tCoordinate& operator[](uint32_t const aIndex) { return mArray[aIndex]; }
  tCoordinate const& operator[](uint32_t const aIndex) const { return mArray[aIndex]; }

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

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
class ShepardInterpolation final {
public:
  using Location = CoefficientWise<tCoordinate, tDimensions>;
  struct Data {
    Location mLocation;
    tPayload mPayload;
  };

  static constexpr uint32_t csChildCount = 1u << tDimensions;
  static constexpr uint32_t csMaxLevels  = 50u;

  static_assert(tDimensions > 0u && tDimensions <= 10u);
  static_assert(tInPlace > 1u && tInPlace <= 1024u);

private:
  struct Node final {
  public:
    using Children = std::array<std::unique_ptr<Node>, csChildCount>;
    using Payload  = std::array<Data, tInPlace>;

    uint32_t                        mCountTotal;
    uint32_t                        mCountHere;
    std::variant<Children, Payload> mContents;
    Location                        mCenter;
    Location                        mSize;

    Node(Location const &aCenter, Location const& aSize) : mCountTotal(0u), mCountHere(0u), mCenter(aCenter), mSize(aSize) {}

    std::pair<uint32_t, Location> getChildIndexCenter(Data const &aItem) const {
      uint32_t index = 0u;
      Location location = mCenter;
      for(uint32_t i = 0u; i < tDimensions; ++i) {
        auto diff = mSize[i] / 4.0;
        if(aItem.mLocation[i] >= mCenter[i]) {
          index += 1u << i;
          location[i] += diff;
        }
        else {
          location[i] -= diff;
        }
      }
      return std::pair(index, location);
    }

void print(uint32_t aLevel) const {
  std::cout << "Level(" << aLevel <<") [CT:" << mCountTotal << " CH:" << mCountHere << " C:" << mCenter[0] << " S:" << mSize[0] << ' ';
  if(mCountHere > 0u) {
    auto& payload = std::get<Payload>(mContents);
    for(uint32_t i = 0u; i < mCountHere; ++i) {
      std::cout << " (l:" << payload[i].mLocation[0] << " p:" << payload[i].mPayload << ")";
    }
  }
  else if(mCountTotal > 0u) {
    auto& children = std::get<Children>(mContents);
    for(uint32_t i = 0u; i < csChildCount; ++i) {
      std::cout << (children[i] ? " (C)" : " (-)");
    }
  }
  std::cout << "]\n";
}
  };

  std::array<std::unique_ptr<Node>, csChildCount> mRoots;
  Location                                        mBoundsMin;
  Location                                        mBoundsMax;
  uint32_t const                                  cmSamplesToConsider;
  uint32_t                                        mTargetLevelInChild0;
  std::vector<uint32_t>                           mNodesPerLevel;
  std::vector<uint32_t>                           mItemsPerLevel;

public:
  ShepardInterpolation() = delete;
  ShepardInterpolation(ShepardInterpolation const &) = delete;
  ShepardInterpolation(ShepardInterpolation &&) = delete;
  ShepardInterpolation& operator=(ShepardInterpolation const &) = delete;
  ShepardInterpolation& operator=(ShepardInterpolation &&) = delete;

  ShepardInterpolation(std::vector<Data> const &aData, uint32_t const aSamplesToConsider);
  uint32_t getTargetLevel()                    const { return mTargetLevelInChild0; }
  uint32_t getLevelCount()                     const { return mNodesPerLevel.size(); }
  uint32_t getNodeCount(uint32_t const aLevel) const { return mNodesPerLevel[aLevel]; }
  uint32_t getItemCount(uint32_t const aLevel) const { return mItemsPerLevel[aLevel]; }

  tPayload interpolate(tCoordinate const aX, tCoordinate const aY, tCoordinate const aZ) const;
  tPayload interpolate(Location const& aLocation) const { return interpolate(aLocation[0], aLocation[1], aLocation[2]); }

private:
  void buildTree(size_t const aWhich, Location const aCenter, Location const aBoundsMax, std::vector<Data> const& aData);
  void addLeaf(Node * aBranch, Location const& aCenter, Location const& aSize, Data const& aItem, uint32_t const aLevel);
  void calculateTargetLevelFromChild0();
};

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace>::ShepardInterpolation(std::vector<Data> const &aData, uint32_t const aSamplesToConsider)
  : mBoundsMin {  std::numeric_limits<tCoordinate>::max()}
  , mBoundsMax { -std::numeric_limits<tCoordinate>::max()}
  , cmSamplesToConsider(aSamplesToConsider)
  , mTargetLevelInChild0(0u) {
  for(auto const &item : aData) {
    for(size_t j = 0u; j < tDimensions; ++j) {
      mBoundsMin = Location::min(mBoundsMin, item.mLocation);
      mBoundsMax = Location::max(mBoundsMax, item.mLocation);
    }
  }

  buildTree(0u, (mBoundsMin + mBoundsMax) / 2u, mBoundsMax - mBoundsMin, aData);
// TODO restore calculateTargetLevelFromChild0();
if(mNodesPerLevel.size() == 0u) calculateTargetLevelFromChild0();
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
void ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace>::buildTree(size_t const aWhich, Location const aCenter, Location const aSize, std::vector<Data> const& aData) {
  mRoots[aWhich] = std::move(std::make_unique<Node>(aCenter, aSize));
  auto root = mRoots[aWhich].get();
  for(auto const &item : aData) {
    addLeaf(root, aCenter, aSize, item, 0u);
mNodesPerLevel.clear();
mItemsPerLevel.clear();
calculateTargetLevelFromChild0();
std::cout << '\n';
  }
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
void ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace>::addLeaf(Node * aBranch, Location const& aCenter, Location const& aSize, Data const& aItem, uint32_t const aLevel) {
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
          auto [childIndex, childCenter] = branch->getChildIndexCenter(item);
          if(!children[childIndex]) {
            children[childIndex] = std::move(std::make_unique<Node>(childCenter, size));
          }
          else {} // nothing to do
          addLeaf(children[childIndex].get(), childCenter, size, item, level + 1u);  // Recursive, but has only a depth of 1
        }
      }
      else {} // Nothing to do
      uint32_t childIndex;
      std::tie(childIndex, center) = branch->getChildIndexCenter(aItem);
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

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
void ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace>::calculateTargetLevelFromChild0() {
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
item.mNode->print(item.mLevel);
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

#endif // SHEPARD_H
