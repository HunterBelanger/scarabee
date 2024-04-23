#ifndef FLAT_SOURCE_REGION_H
#define FLAT_SOURCE_REGION_H

#include <moc/surface.hpp>
#include <moc/vector.hpp>
#include <moc/direction.hpp>
#include <transport_xs.hpp>
#include <utils/constants.hpp>

#include <htl/static_vector.hpp>

#include <xtensor/xtensor.hpp>

#include <memory>
#include <vector>

namespace scarabee {

struct RegionToken {
  std::shared_ptr<Surface> surface;
  Surface::Side side;

  bool inside(const Vector& r, const Direction& u) const {
    const auto current_side = surface->side(r, u);
    return current_side == side;
  }
};

class FlatSourceRegion {
 public:
  FlatSourceRegion() : tokens_(), xs_(), volume_(), id_(id_counter++) {}

  bool inside(const Vector& r, const Direction& u) const {
    for (const auto& t : tokens_) {
      if (t.inside(r, u) == false) return false;
    }
    return true;
  }

  double distance(const Vector& r, const Direction& u) const {
    double min_dist = INF;
    for (const auto& t : tokens_) {
      double token_dist = t.surface->distance(r, u);
      if (token_dist < min_dist) min_dist = token_dist;
    }
    return min_dist;
  }

  std::size_t id() const { return id_; }

  std::shared_ptr<TransportXS>& xs() { return xs_; }
  const std::shared_ptr<TransportXS>& xs() const { return xs_; }

  htl::static_vector<RegionToken, MAX_SURFS>& tokens() { return tokens_; }
  const htl::static_vector<RegionToken, MAX_SURFS>& tokens() const {
    return tokens_;
  }

  double& volume() { return volume_; }
  const double& volume() const { return volume_; }

 private:
  htl::static_vector<RegionToken, MAX_SURFS> tokens_;
  std::shared_ptr<TransportXS> xs_;
  double volume_;
  std::size_t id_;

  static std::size_t id_counter;
};

struct UniqueFSR {
  const FlatSourceRegion* fsr{nullptr};
  std::size_t instance{0};
};

}  // namespace scarabee

#endif
