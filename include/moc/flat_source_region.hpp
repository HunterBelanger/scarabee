#ifndef FLAT_SOURCE_REGION_H
#define FLAT_SOURCE_REGION_H

#include <moc/surface.hpp>
#include <moc/vector.hpp>
#include <moc/direction.hpp>
#include <transport_xs.hpp>
#include <utils/constants.hpp>

#include <htl/static_vector.hpp>

#include <memory>
#include <vector>

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
  FlatSourceRegion() : tokens_(), flux_(), source_(), xs_(), volume_() {}

  bool inside(const Vector& r, const Direction& u) const {
    for (const auto& t : tokens_) {
      if (t.inside(r, u) == false) return false;
    }
    return true;
  }

  void initialize();

  double distance(const Vector& r, const Direction& u) const {
    double min_dist = INF;
    for (const auto& t : tokens_) {
      double token_dist = t.surface->distance(r, u);
      if (token_dist < min_dist) min_dist = token_dist;
    }
    return min_dist;
  }

  std::shared_ptr<TransportXS>& xs() { return xs_; }
  const std::shared_ptr<TransportXS>& xs() const { return xs_; }

  htl::static_vector<RegionToken, MAX_SURFS>& tokens() { return tokens_; }
  const htl::static_vector<RegionToken, MAX_SURFS>& tokens() const {
    return tokens_;
  }

  double& volume() { return volume_; }
  const double& volume() const { return volume_; }

  std::vector<double>& flux() { return flux_; }
  const std::vector<double>& flux() const { return flux_; }

  std::vector<double>& source() { return source_; }
  const std::vector<double>& source() const { return source_; }

 private:
  htl::static_vector<RegionToken, MAX_SURFS> tokens_;
  std::vector<double> flux_;
  std::vector<double> source_;
  std::shared_ptr<TransportXS> xs_;
  double volume_;
};

#endif