#include <geometry.hpp>
#include <material.hpp>

#include <algorithm>
#include <exception>

Geom1D::Geom1D(const YAML::Node& node): ngroups_(0), bounds_(), tiles_() {
  // Get the tile boundaries
  if (!node["bounds"]) {
    throw std::runtime_error("No bounds entry in geometry.");
  } else if (node["bounds"].IsSequence() == false) {
    throw std::runtime_error("bounds entry in geometry must be a sequence of floats.");
  }
  
  if (node["bounds"].size() < 4) {
    throw std::runtime_error("Must have at least 3 tiles in geometry definition.");
  }

  // Make sure we have the tiles
  if (!node["tiles"]) {
    throw std::runtime_error("No tiles entry in geometry.");
  } else if (node["tiles"].IsSequence() == false) {
    throw std::runtime_error("tiles entry in geometry must be a sequence.");
  }

  if (node["tiles"].size() != node["bounds"].size()-1) {
    throw std::runtime_error("Length of tiles should be one less than length of bounds.");
  }

  // Check if subdivisions are given
  bool has_subdivs = false;
  if (node["subdivs"] && node["subdivs"].IsSequence() &&
      (node["subdivs"].size() == node["tiles"].size())) {
    has_subdivs = true;
  } else if (node["subdivs"]) {
    throw std::runtime_error("Invalid subdivs entry in geometry.");
  }

  // Go through all tiles an either initialize BC or get material
  for (std::size_t t = 0; t < node["tiles"].size(); t++) {
    std::string tile_val = node["tiles"][t].as<std::string>();
    float lower_bound = node["bounds"][t].as<float>();
    float upper_bound = node["bounds"][t+1].as<float>();
    std::size_t subdivs = 1;

    if (has_subdivs) {
      subdivs = node["subdivs"][t].as<std::size_t>();
    }
    const float dx = (upper_bound - lower_bound) / static_cast<float>(subdivs);

    // Make the tile
    MatTile tile;
    if (tile_val == "v" || tile_val == "V") {
      tile.type = MatTile::Type::Vacuum;
    } else if (tile_val == "r" || tile_val == "R") {
      tile.type = MatTile::Type::Reflection;
    } else {
      std::uint32_t mat_id = 0;

      try {
        mat_id = static_cast<std::uint32_t>(std::stoul(tile_val));
      } catch (...) {
        throw std::runtime_error("Could not interpret geometry tile " + tile_val + ".");
      }

      if (materials.find(mat_id) == materials.end()) {
        throw std::runtime_error("Unknown material id tile " + tile_val + ".");
      }

      tile.type = MatTile::Type::Material;
      tile.material = &materials.at(mat_id);

      if (ngroups_ == 0) {
        // Get number of groups from first material
        ngroups_ = tile.material->ngroups();
      }
    }

    // Add the tile subdiv times, also adding the bounds
    for (std::size_t s = 0; s < subdivs; s++) {
      tiles_.push_back(tile);
      bounds_.push_back(lower_bound);
      lower_bound += dx;
    }
  }

  if (ngroups_ == 0) {
    // No materials were given !
    throw std::runtime_error("No material tiles were given in geometry.");
  }

  // Make sure the bounds are sorted !
  if (std::is_sorted(bounds_.begin(), bounds_.end()) == false) {
    throw std::runtime_error("Bounds in geometry are not sorted.");
  }
}