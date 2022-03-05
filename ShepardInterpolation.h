#ifndef SHEPARDINTERPOLATION
#define SHEPARDINTERPOLATION

#include <array>
#include <deque>
#include <limits>
#include <memory>
#include <vector>
#include <variant>
#include <Eigen/Dense>


template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
class ShepardInterpolation final {
public:
  using Location = Eigen::Array<tCoordinate, 1, tDimensions>;
  struct Data {
    Location mLocation;
    tPayload mPayload;
  };

  static constexpr uint32_t csChildCount = 1u << tDimensions;

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
  };

  std::array<std::unique_ptr<Node>, csChildCount> mRoots;
  Location                                        mBoundsMin;
  Location                                        mBoundsMax;
  uint32_t const                                  cmSamplesToConsider;
  uint32_t                                        mTargetLevelInChild0;

public:
  ShepardInterpolation(std::vector<Data> const &aData, uint32_t const aSamplesToConsider);
  uint32_t getTargetLevel() const { return mTargetLevelInChild0; }

  tPayload interpolate(tCoordinate const aX, tCoordinate const aY, tCoordinate const aZ) const;
  tPayload interpolate(Location const& aLocation) const { return interpolate(aLocation[0], aLocation[1], aLocation[2]); }

private:
  void buildTree(size_t const aWhich, Location const aCenter, Location const aBoundsMax, std::vector<Data> const& aData);
  void addLeaf(Node * aBranch, Location const& aCenter, Location const& aSize, Data const& aItem);
  uint32_t calculateTargetLevelFromChild0() const;
};

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace>::ShepardInterpolation(std::vector<Data> const &aData, uint32_t const aSamplesToConsider)
  : mBoundsMin {  std::numeric_limits<tCoordinate>::max(),  std::numeric_limits<tCoordinate>::max(),  std::numeric_limits<tCoordinate>::max()}
  , mBoundsMax { -std::numeric_limits<tCoordinate>::max(), -std::numeric_limits<tCoordinate>::max(), -std::numeric_limits<tCoordinate>::max()}
  , cmSamplesToConsider(aSamplesToConsider)
  , mTargetLevelInChild0(0) {
  for(auto const &item : aData) {
    for(size_t j = 0u; j < tDimensions; ++j) {
      mBoundsMin[j] = std::min(mBoundsMin[j], item.mLocation[j]);
      mBoundsMax[j] = std::max(mBoundsMax[j], item.mLocation[j]);
    }
  }

  buildTree(0u, (mBoundsMin + mBoundsMax) / 2u, mBoundsMax - mBoundsMin, aData);
  mTargetLevelInChild0 = calculateTargetLevelFromChild0();
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
void ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace>::buildTree(size_t const aWhich, Location const aCenter, Location const aSize, std::vector<Data> const& aData) {
  mRoots[aWhich] = std::move(std::make_unique<Node>(aCenter, aSize));
  auto root = mRoots[aWhich].get();
  for(auto const &item : aData) {
    addLeaf(root, aCenter, aSize, item);
  }
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
void ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace>::addLeaf(Node * aBranch, Location const& aCenter, Location const& aSize, Data const& aItem) {
  Node * branch    = aBranch;
  Location center = aCenter;
  Location size   = aSize;
  while(true) {
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
          addLeaf(children[childIndex].get(), childCenter, size, item);
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
    }
  }
}

template<typename tCoordinate, uint32_t tDimensions, typename tPayload, size_t tInPlace>
uint32_t ShepardInterpolation<tCoordinate, tDimensions, tPayload, tInPlace>::calculateTargetLevelFromChild0() const {
  struct Item {
    Node const * const cmNode;
    uint32_t const     cmLevel;
    Item(Node const * const aNode, uint32_t const aLevel) : cmNode(aNode), cmLevel(aLevel) {}
  };
  std::deque<Item> queue;
  uint32_t level = 0u;
  auto root = mRoots[0].get();
  queue.emplace_back(root, level);
  while(queue.size() > 0u) {
    double average = 0.0;
    int32_t count = 0u;
    auto item = queue.front();
    while(item.cmLevel == level) {
      average += item.cmNode->mCountTotal;
      ++count;
      if(std::holds_alternative<typename Node::Children>(item.cmNode->mContents)) {
        auto &children = std::get<typename Node::Children>(item.cmNode->mContents);
        for(auto const &child : children) {
          if(child) {
            queue.emplace_back(child.get(), level + 1u);
          }
          else {} // Nothing to do
        }
      }
      else {} // Nothing to do
      queue.pop_front();
    }
    average /= count;
    if(average < cmSamplesToConsider) {
      break;
    }
    else {} // Nothing to do
    ++level;
  }
  return level;
}

#endif // SHEPARD_H
