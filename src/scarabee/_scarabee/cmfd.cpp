#include <moc/cmfd.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <cmath>

namespace scarabee {

CMFD::CMFD(const std::vector<double>& dx, const std::vector<double>& dy, const std::vector<std::pair<std::size_t, std::size_t>>& groups)
    : x_bounds_(),
      y_bounds_(),
      moc_to_cmfd_group_map_(),
      nx_(),
      ny_(),
      nx_surfs_(),
      ny_surfs_(),
      temp_fsrs_(),
      fsrs_() {
  // Make sure we have at least 1 bin in each direction
  if (dx.size() == 0) {
    auto mssg = "Must provide at least 1 x width.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (dy.size() == 0) {
    auto mssg = "Must provide at least 1 y width.";
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

  double dy_tot = 0.;
  for (std::size_t j = 0; j < dy.size(); j++) {
    if (dy[j] <= 0.) {
      std::stringstream mssg;
      mssg << "dy at index " << j << " is <= 0.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    dy_tot += dy[j];
  }

  // Create surfaces
  x_bounds_.reserve(dx.size() + 1);
  x_bounds_.emplace_back();
  x_bounds_.back().type() = Surface::Type::XPlane;
  x_bounds_.back().x0() = -0.5 * dx_tot;
  for (const auto& d : dx) {
    const double new_x0 = x_bounds_.back().x0() + d;

    x_bounds_.emplace_back();
    x_bounds_.back().type() = Surface::Type::XPlane;
    x_bounds_.back().x0() = new_x0;
  }

  y_bounds_.reserve(dy.size() + 1);
  y_bounds_.emplace_back();
  y_bounds_.back().type() = Surface::Type::YPlane;
  y_bounds_.back().y0() = -0.5 * dy_tot;
  for (const auto& d : dy) {
    const double new_y0 = y_bounds_.back().y0() + d;

    y_bounds_.emplace_back();
    y_bounds_.back().type() = Surface::Type::YPlane;
    y_bounds_.back().y0() = new_y0;
  }

  nx_ = dx.size();
  ny_ = dy.size();

  nx_surfs_ = x_bounds_.size() * (y_bounds_.size()-1);
  ny_surfs_ = y_bounds_.size() * (x_bounds_.size()-1);

  // Allocate temp_fsrs_ 
  temp_fsrs_.resize(nx_*ny_);

  // Check group condensation scheme
  if (groups.size() == 0) {
    auto mssg = "Empty energy condensation scheme provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (groups.front().first != 0) {
    auto mssg = "The energy condensation scheme does not start with 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (std::size_t i = 0; i < groups.size() - 1; i++) {
    if (groups[i].second + 1 != groups[i + 1].first) {
      std::stringstream mssg;
      mssg << "Condensed groups " << i << " and " << i + 1
           << " are not continuous.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
  }

  moc_to_cmfd_group_map_.resize(groups.back().second+1);
  for (std::size_t g = 0; g < moc_to_cmfd_group_map_.size(); g++) {
    for (std::size_t G = 0; G < groups.size(); G++) {
      if (groups[G].first <= g && g <= groups[G].second) {
        moc_to_cmfd_group_map_[g] = G;
        break;
      }
    }
  }

  // Allocate surfaces array
  surface_currents_ = xt::zeros<double>({groups.size(), nx_surfs_+ny_surfs_});
}

std::optional<std::array<std::size_t, 2>> CMFD::get_tile(
    const Vector& r, const Direction& u) const {
  for (std::size_t i = 0; i < nx(); i++) {
    for (std::size_t j = 0; j < ny(); j++) {
      // Get the surfaces that make up our tile
      const auto& xl = x_bounds_[i];
      const auto& xh = x_bounds_[i + 1];
      const auto& yl = y_bounds_[j];
      const auto& yh = y_bounds_[j + 1];

      // Check if we are in the tile. If so, return index
      if (xl.side(r, u) == Surface::Side::Positive &&
          xh.side(r, u) == Surface::Side::Negative &&
          yl.side(r, u) == Surface::Side::Positive &&
          yh.side(r, u) == Surface::Side::Negative) {
        return std::array<std::size_t, 2>{i, j};
      }
    }
  }

  // If we get here, we weren't in any tile.
  // Return nullopt
  return std::nullopt;
}

std::optional<std::size_t> CMFD::get_surface(const Vector& r,
                                             const Direction& u) const {
  // First, we check all the x surfaces, and see if we are on one
  for (std::size_t xsi = 0; xsi < x_bounds_.size(); xsi++) {
    if (std::abs(x_bounds_[xsi].x0() - r.x()) < SURFACE_COINCIDENT) {
      // we are on this surface ! get our tile (just used for the y bin)
      auto otile = this->get_tile(r, u);
      if (otile.has_value() == false) otile = this->get_tile(r, -u);
      if (otile.has_value() == false) return std::nullopt;

      const auto& tile = otile.value();

      return tile[1] * x_bounds_.size() + xsi;
    }
  }

  // We apparently weren't on an x surface, so try all y surfaces.
  for (std::size_t ysi = 0; ysi < y_bounds_.size(); ysi++) {
    if (std::abs(y_bounds_[ysi].y0() - r.y()) < SURFACE_COINCIDENT) {
      // we are on this surface ! get our tile (just used for the x bin)
      auto otile = this->get_tile(r, u);
      if (otile.has_value() == false) otile = this->get_tile(r, -u);
      if (otile.has_value() == false) return std::nullopt;

      const auto& tile = otile.value();
      
      return nx_surfs_ + tile[0] * y_bounds_.size() + ysi;
    }
  }

  // If we get here, we aren't on a surface
  return std::nullopt;
}

std::size_t CMFD::tile_to_indx(const std::array<std::size_t, 2>& tile) const {
  return tile[1] * nx_ + tile[0];
}

void CMFD::insert_fsr(const std::array<std::size_t, 2>& tile, std::size_t fsr) {
  // Compute linear index
  const std::size_t i = this->tile_to_indx(tile);
  temp_fsrs_.at(i).insert(fsr);
}

void CMFD::pack_fsr_lists() {
  fsrs_.resize(nx_ * ny_, std::vector<std::size_t>());

  for (std::size_t i = 0; i < fsrs_.size(); i++) {
    fsrs_[i].insert(fsrs_[i].begin(), temp_fsrs_[i].begin(),
                    temp_fsrs_[i].end());
  }

  temp_fsrs_.clear();
  temp_fsrs_.shrink_to_fit();
}

std::size_t CMFD::moc_to_cmfd_group(std::size_t g) const {
  if (g >= moc_to_cmfd_group_map_.size()) {
    auto mssg = "Group index is out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return moc_to_cmfd_group_map_[g];
}

double& CMFD::current(const std::size_t G, const std::size_t surface) {
  if (G >= surface_currents_.shape()[0]) {
    auto mssg = "Group index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (surface >= surface_currents_.shape()[1]) {
    auto mssg = "Surface index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return surface_currents_(G, surface);
}

const double& CMFD::current(const std::size_t G,
                            const std::size_t surface) const {
  if (G >= surface_currents_.shape()[0]) {
    auto mssg = "Group index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (surface >= surface_currents_.shape()[1]) {
    auto mssg = "Surface index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return surface_currents_(G, surface);
}

void CMFD::tally_current(double aflx, const Direction& u, std::size_t G, const std::size_t surf) {
  if (G >= surface_currents_.shape()[0]) {
    auto mssg = "Group index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (surf >= surface_currents_.shape()[1]) {
    auto mssg = "Surface index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (surf < nx_surfs_) {
#pragma omp atomic
    surface_currents_(G, surf) += aflx * u.x();
  } else {
#pragma omp atomic
    surface_currents_(G, surf) += aflx * u.y();
  }
}

}  // namespace scarabee