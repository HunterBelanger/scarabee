#ifndef SCARABEE_GEOMETRY_H
#define SCARABEE_GEOMETRY_H

#include <material.hpp>

#include <cstdint>
#include <vector>

#include <yaml-cpp/yaml.h>

struct MatTile {
  enum class Type {Vacuum, Reflection, Material};
  Type type = Type::Vacuum;
  Material* material = nullptr;
};

class Geom1D {
  public:
    Geom1D(const YAML::Node& node);

    const MatTile& operator()(std::size_t i) const { return tiles_[i]; }

    float dx(std::size_t i) const { return bounds_[i+1] - bounds_[i]; }

    std::size_t size() const { return tiles_.size(); }
    
    const std::vector<float>& bounds() const { return bounds_; }

    std::uint32_t ngroups() const { return ngroups_; }

  private:
    std::uint32_t ngroups_;
    std::vector<float> bounds_;
    std::vector<MatTile> tiles_;
};

#endif