#include <diffusion/diffusion_geometry.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <xtensor/xstrides.hpp>

#include <sstream>

namespace scarabee {

DiffusionGeometry::DiffusionGeometry(const std::vector<double>& dx,
                                     double albedo_xn, double albedo_xp)
    : tiles_(),
      xn_(),
      xp_(),
      yn_(),
      yp_(),
      zn_(),
      zp_(),
      x_bounds_(),
      y_bounds_(),
      z_bounds_(),
      n_mat_tiles_(0),
      mat_indx_to_flat_tile_indx_() {
  // Make sure we have at least 3 bins in each direction
  if (dx.size() == 0) {
    auto mssg = "Must provide at least 3 x widths.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure all widths > 0
  double dx_tot = 0.;
  for (std::size_t i = 0; i < dx.size(); i++) {
    if (dx[i] <= 0.) {
      std::stringstream mssg;
      mssg << "dx at index " << i << " is <= 0.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    dx_tot += dx[i];
  }
  
  // Create bounds
  x_bounds_.reserve(dx.size() + 1);
  x_bounds_.push_back(-0.5 * dx_tot);
  for (const auto& d : dx) {
    const double new_x0 = x_bounds_.back() + d;
    x_bounds_.push_back(new_x0);
  }

  // Assign boundary conditions
  if (albedo_xn < 0. || albedo_xn > 1.) {
    std::stringstream mssg;
    mssg << "Negative x albedo is invalid. Must be in range [0, 1].";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }
  xn_.albedo = albedo_xn;
  xn_.xs = nullptr;

  if (albedo_xp < 0. || albedo_xp > 1.) {
    std::stringstream mssg;
    mssg << "Positive x albedo is invalid. Must be in range [0, 1].";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }
  xp_.albedo = albedo_xp;
  xp_.xs = nullptr;

  // Reshape tiles array
  tiles_.resize({nx()});
}

void DiffusionGeometry::set_tiles(const std::vector<TileFill>& tiles) {
  if (ndims() == 1) {
    this->set_tiles_1d(tiles);
  } else if (ndims() == 2) {
    throw ScarabeeException("2D not yet implemented.");
  } else {
    // ndims == 3
    throw ScarabeeException("3D not yet implemented.");
  }
}

std::size_t DiffusionGeometry::ngroups() const {
  if (nmats() == 0) return 0;

  return mat(0)->ngroups();
}

std::pair<DiffusionGeometry::Tile, std::optional<std::size_t>> DiffusionGeometry::neighbor(std::size_t m, Neighbor n) const {
  std::vector<std::size_t> array_indx;
  try {
    array_indx = this->mat_indxs(m);
  } catch (ScarabeeException& err) {
    auto mssg = "Could not obtain material array index.";
    spdlog::error(mssg);
    err.add_to_exception(mssg);
    throw err;
  }

  // Based on our dimensions and what we got, we return the new tile
  if (ndims() == 1) {
    if (n != Neighbor::XN && n != Neighbor::XP) {
      std::stringstream mssg;
      mssg << "A 1D DiffusionGeometry was passed a non x-direction neighbor "
              "request.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (array_indx[0] == 0 && n == Neighbor::XN) {
      return {xn_, std::nullopt};
    } else if (array_indx[0] == nx() - 1 && n == Neighbor::XP) {
      return {xp_, std::nullopt};
    }

    // Now we change the index to get the neighbor
    if (n == Neighbor::XN) {
      array_indx[0] -= 1;
    } else {
      array_indx[0] += 1;
    }
  } else if (ndims() == 2) {
    if (n == Neighbor::ZN || n == Neighbor::ZP) {
      std::stringstream mssg;
      mssg << "A 2D DiffusionGeometry was passed a z-direction neighbor "
              "request.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    if (array_indx[0] == 0 && n == Neighbor::XN) {
      return {xn_, std::nullopt};
    } else if (array_indx[0] == nx() - 1 && n == Neighbor::XP) {
      return {xp_, std::nullopt};
    } else if (array_indx[1] == 0 && n == Neighbor::YN) {
      return {yn_, std::nullopt};
    } else if (array_indx[1] == ny() - 1 && n == Neighbor::YP) {
      return {yp_, std::nullopt};
    }

    // Now we change the index to get the neighbor
    if (n == Neighbor::XN) {
      array_indx[0] -= 1;
    } else if (n == Neighbor::XP) {
      array_indx[0] += 1;
    } else if (n == Neighbor::YN) {
      array_indx[1] -= 1;
    } else {
      array_indx[1] += 1;
    }
  } else {  // ndims() == 3
    if (array_indx[0] == 0 && n == Neighbor::XN) {
      return {xn_, std::nullopt};
    } else if (array_indx[0] == nx() - 1 && n == Neighbor::XP) {
      return {xp_, std::nullopt};
    } else if (array_indx[1] == 0 && n == Neighbor::YN) {
      return {yn_, std::nullopt};
    } else if (array_indx[1] == ny() - 1 && n == Neighbor::YP) {
      return {yp_, std::nullopt};
    } else if (array_indx[2] == 0 && n == Neighbor::ZN) {
      return {zn_, std::nullopt};
    } else if (array_indx[2] == nz() - 1 && n == Neighbor::ZP) {
      return {zp_, std::nullopt};
    }

    if (n == Neighbor::XN) {
      array_indx[0] -= 1;
    } else if (n == Neighbor::XP) {
      array_indx[0] += 1;
    } else if (n == Neighbor::YN) {
      array_indx[1] -= 1;
    } else if (n == Neighbor::YP) {
      array_indx[1] += 1;
    } else if (n == Neighbor::ZN) {
      array_indx[2] -= 1;
    } else if (n == Neighbor::ZP) {
      array_indx[2] += 1;
    }
  }

  Tile out_tile = tiles_.element(array_indx.begin(), array_indx.end());
  std::optional<std::size_t> out_mat_indx {std::nullopt};

  if (out_tile.xs) {
    std::vector<std::vector<std::size_t>> vec_of_indx {array_indx};
    out_mat_indx = xt::ravel_indices(vec_of_indx, tiles_.shape(), tiles_.layout())[0];
  }

  return {out_tile, out_mat_indx};
}

const std::shared_ptr<DiffusionCrossSection>& DiffusionGeometry::mat(std::size_t m) const {
  if (m >= nmats()) {
    std::stringstream mssg;
    mssg << "Requested material index out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  return tiles_.flat(mat_indx_to_flat_tile_indx_[m]).xs;
}

double DiffusionGeometry::volume(std::size_t m ) const {
  if (m >= nmats()) {
    std::stringstream mssg;
    mssg << "Requested material index out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  auto inds = this->mat_indxs(m);

  if (ndims() == 1) {
    return dx(inds[0]);
  } else if (ndims() == 2) {
    return dx(inds[0]) * dy(inds[1]);
  } else {
    return dx(inds[0]) * dy(inds[1]) * dz(inds[2]);
  }
}

std::vector<std::size_t> DiffusionGeometry::mat_indxs(std::size_t m) const {
  if (m >= nmats()) {
    std::stringstream mssg;
    mssg << "Requested material index out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Get flat index for material
  const std::size_t i = mat_indx_to_flat_tile_indx_[m];

  // Now we get the array_indx
  std::array<std::size_t, 1> indx{i};
  auto array_indx =
      xt::unravel_indices(indx, tiles_.shape(), tiles_.layout())[0];

  return std::vector<std::size_t>(array_indx.begin(), array_indx.end());
}

void DiffusionGeometry::set_tiles_1d(const std::vector<TileFill>& tiles) {
  // Check for the right dimensions
  if (tiles.size() != nx()) {
    std::stringstream mssg;
    mssg << "The number of provided tiles disagrees with the problem geometry.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Go through all tiles. Should only have xs tiles, as a 1D problem shouldn't
  // have BCs in the middle of the problem domain.
  std::size_t i = 0;
  for (const auto& tile : tiles) {
    if (std::holds_alternative<double>(tile)) {
      std::stringstream mssg;
      mssg << "A 1D diffusion problem cannot have albedo tiles.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    tiles_(i).albedo = std::nullopt;
    tiles_(i).xs = std::get<std::shared_ptr<DiffusionCrossSection>>(tile);
    n_mat_tiles_++;

    mat_indx_to_flat_tile_indx_.push_back(i);
    
    i++;
  }
}

}  // namespace scarabee
