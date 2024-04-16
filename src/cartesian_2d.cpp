#include <moc/cartesian_2d.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>
#include <utils/constants.hpp>

#include <algorithm>
#include <optional>

namespace scarabee {

Cartesian2D::Cartesian2D(const std::vector<std::shared_ptr<Surface>>& x_bounds,
                         const std::vector<std::shared_ptr<Surface>>& y_bounds)
    : x_bounds_(x_bounds), y_bounds_(y_bounds), tiles_(), nx_(), ny_() {
  // First, make sure we have at least 2 bounds in each direction
  if (x_bounds_.size() < 2) {
    auto mssg = "Must provide at least 2 x-bounds.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (y_bounds_.size() < 2) {
    auto mssg = "Must provide at least 2 y-bounds.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure bounds are sorted
  if (std::is_sorted(x_bounds_.begin(), x_bounds_.end(),
                     [](const auto& a, const auto& b) {
                       return a->x0() < b->x0();
                     }) == false) {
    auto mssg = "Provided x-bounds are not sorted.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (std::is_sorted(y_bounds_.begin(), y_bounds_.end(),
                     [](const auto& a, const auto& b) {
                       return a->y0() < b->y0();
                     }) == false) {
    auto mssg = "Provided y-bounds are not sorted.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Get sizes and allocate tiles array
  nx_ = x_bounds_.size() - 1;
  ny_ = y_bounds_.size() - 1;

  // All tiles start as uninitialized
  tiles_.resize({nx_, ny_});
  tiles_.fill(Tile{nullptr, nullptr});
}

Cartesian2D::Cartesian2D(const std::vector<double>& dx,
                         const std::vector<double>& dy)
    : x_bounds_(), y_bounds_(), tiles_(), nx_(), ny_() {
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
  x_bounds_.push_back(std::make_shared<Surface>());
  x_bounds_.back()->type() = Surface::Type::XPlane;
  x_bounds_.back()->x0() = -0.5 * dx_tot;
  for (const auto& d : dx) {
    const double new_x0 = x_bounds_.back()->x0() + d;

    x_bounds_.push_back(std::make_shared<Surface>());
    x_bounds_.back()->type() = Surface::Type::XPlane;
    x_bounds_.back()->x0() = new_x0;
  }

  y_bounds_.reserve(dy.size() + 1);
  y_bounds_.push_back(std::make_shared<Surface>());
  y_bounds_.back()->type() = Surface::Type::YPlane;
  y_bounds_.back()->y0() = -0.5 * dy_tot;
  for (const auto& d : dy) {
    const double new_y0 = y_bounds_.back()->y0() + d;

    y_bounds_.push_back(std::make_shared<Surface>());
    y_bounds_.back()->type() = Surface::Type::YPlane;
    y_bounds_.back()->y0() = new_y0;
  }

  nx_ = dx.size();
  ny_ = dy.size();

  // All tiles start as uninitialized
  tiles_.resize({nx_, ny_});
  tiles_.fill(Tile{nullptr, nullptr});
}

std::optional<Cartesian2D::TileIndex> Cartesian2D::get_tile_index(
    const Vector& r, const Direction& u) const {
  for (std::size_t i = 0; i < nx(); i++) {
    for (std::size_t j = 0; j < ny(); j++) {
      // Get the surfaces that make up our tile
      const auto& xl = x_bounds_[i];
      const auto& xh = x_bounds_[i + 1];
      const auto& yl = y_bounds_[j];
      const auto& yh = y_bounds_[j + 1];

      // Check if we are in the tile. If so, return index
      if (xl->side(r, u) == Surface::Side::Positive &&
          xh->side(r, u) == Surface::Side::Negative &&
          yl->side(r, u) == Surface::Side::Positive &&
          yh->side(r, u) == Surface::Side::Negative) {
        return TileIndex{i, j};
      }
    }
  }

  // If we get here, we weren't in any tile.
  // Return nullopt
  return std::nullopt;
}

Vector Cartesian2D::get_tile_center(const TileIndex& ti) const {
  // Get the lower and upper x and y bounds
  const auto& xl = x_bounds_[ti.i];
  const auto& xh = x_bounds_[ti.i + 1];
  const auto& yl = y_bounds_[ti.j];
  const auto& yh = y_bounds_[ti.j + 1];

  return Vector(0.5 * (xl->x0() + xh->x0()), 0.5 * (yl->y0() + yh->y0()));
}

double Cartesian2D::trace_segments(Vector& r, const Direction& u,
                                   std::vector<Segment>& segments) {
  if (this->tiles_valid() == false) {
    auto mssg =
        "Cannot trace segments on a Cartesian2D geometry with invalid tiles.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  double dist_total = 0.;

  // Get our current tile index
  auto ti = this->get_tile_index(r, u);

  // While we have a tile index
  while (ti) {
    // Get tile center
    auto r_tile_center = this->get_tile_center(*ti);

    // Get reference to the tile
    auto& tile = tiles_(ti->i, ti->j);

    // Trace segments on tile
    auto r_tile = r - r_tile_center;
    double dist = tile.trace_segments(r_tile, u, segments);
    r = r + dist * u;
    dist_total += dist;

    // Get new tile index
    ti = this->get_tile_index(r, u);
  }

  return dist_total;
}

const Cartesian2D::Tile& Cartesian2D::tile(
    const Cartesian2D::TileIndex& ti) const {
  check_tile_index(ti);
  return tiles_(ti.i, ti.j);
}

void Cartesian2D::set_tile(const TileIndex& ti,
                           const std::shared_ptr<Cartesian2D>& c2d) {
  check_tile_index(ti);
  const auto dxdy = tile_dx_dy(ti);

  if (std::abs(dxdy.first - c2d->dx()) > VEC_FP_TOL) {
    auto mssg = "x width does not agree with tile x width";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (std::abs(dxdy.second - c2d->dy()) > VEC_FP_TOL) {
    auto mssg = "y width does not agree with tile y width";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  tiles_(ti.i, ti.j).c2d = std::make_shared<Cartesian2D>(*c2d);
}

void Cartesian2D::set_tile(const TileIndex& ti,
                           const std::shared_ptr<Cell>& cell) {
  check_tile_index(ti);
  const auto dxdy = tile_dx_dy(ti);

  std::pair<double, double> cdxdy = {cell->dx(), cell->dy()};

  if (std::abs(dxdy.first - cdxdy.first) > VEC_FP_TOL) {
    auto mssg = "x width does not agree with tile x width";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (std::abs(dxdy.second - cdxdy.second) > VEC_FP_TOL) {
    auto mssg = "y width does not agree with tile y width";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  tiles_(ti.i, ti.j).cell = cell->clone();
}

void Cartesian2D::set_tiles(const std::vector<TileFill>& fills) {
  // Make sure dimensions match
  if (fills.size() != tiles_.size()) {
    auto mssg = "Number of provided file fills does not match number of files.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  std::size_t indx = 0;
  for (std::size_t j = ny(); j > 0; j--) {
    for (std::size_t i = 0; i < nx(); i++) {
      // Check which fill type we have
      const TileFill& fill = fills[indx];
      if (std::holds_alternative<std::shared_ptr<Cartesian2D>>(fill)) {
        this->set_tile({i, j - 1},
                       std::get<std::shared_ptr<Cartesian2D>>(fill));
      } else {
        this->set_tile({i, j - 1}, std::get<std::shared_ptr<Cell>>(fill));
      }

      indx++;
    }
  }
}

std::pair<double, double> Cartesian2D::tile_dx_dy(const TileIndex& ti) const {
  // Get the lower and upper x and y bounds
  const double xl = x_bounds_[ti.i]->x0();
  const double xh = x_bounds_[ti.i + 1]->x0();
  const double yl = y_bounds_[ti.j]->y0();
  const double yh = y_bounds_[ti.j + 1]->y0();

  return {xh - xl, yh - yl};
}

void Cartesian2D::check_tile_index(const TileIndex& ti) const {
  if (ti.i >= nx_) {
    auto mssg = "TileIndex i out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ti.j >= ny_) {
    auto mssg = "TileIndex j out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

bool Cartesian2D::tiles_valid() const {
  for (const auto& t : tiles_) {
    if (t.valid() == false) return false;
  }

  return true;
}

std::shared_ptr<TransportXS> Cartesian2D::get_xs(const Vector& r,
                                                 const Direction& u) const {
  try {
    const auto& fsr = this->get_fsr(r, u);
    return fsr.xs();
  } catch (ScarabeeException& err) {
    err.add_to_exception("Could not find flat source region.");
    throw err;
  }
}

FlatSourceRegion& Cartesian2D::get_fsr(const Vector& r, const Direction& u) {
  auto ti = this->get_tile_index(r, u);

  if (ti.has_value() == false) {
    // We apparently are not in a tile
    std::stringstream mssg;
    mssg << "Position r = " << r << ", and Direction u = " << u
         << ", are not in Cartesian2D geometry.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  check_tile_index(*ti);
  auto& t = tiles_(ti->i, ti->j);
  const Vector tc = get_tile_center(*ti);
  Vector r_tile = r - tc;

  if (t.valid() == false) {
    std::stringstream mssg;
    mssg << "Tile for Position r = " << r << ", and Direction u = " << u
         << " is empty.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (t.c2d) {
    return t.c2d->get_fsr(r_tile, u);
  } else {
    return t.cell->get_fsr(r_tile, u);
  }
}

const FlatSourceRegion& Cartesian2D::get_fsr(const Vector& r,
                                             const Direction& u) const {
  auto ti = this->get_tile_index(r, u);

  if (ti.has_value() == false) {
    // We apparently are not in a tile
    std::stringstream mssg;
    mssg << "Position r = " << r << ", and Direction u = " << u
         << ", are not in Cartesian2D geometry.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  const auto& t = this->tile(*ti);
  const Vector tc = get_tile_center(*ti);
  Vector r_tile = r - tc;

  if (t.valid() == false) {
    std::stringstream mssg;
    mssg << "Tile for Position r = " << r << ", and Direction u = " << u
         << " is empty.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  if (t.c2d) {
    return t.c2d->get_fsr(r_tile, u);
  } else {
    return t.cell->get_fsr(r_tile, u);
  }
}

std::size_t Cartesian2D::num_fsrs() const {
  std::size_t n = 0;

  for (const auto& tile : tiles_) n += tile.num_fsrs();

  return n;
}

void Cartesian2D::append_fsrs(std::vector<FlatSourceRegion*>& fsrs) {
  for (auto& tile : tiles_) tile.append_fsrs(fsrs);
}

double Cartesian2D::Tile::trace_segments(Vector& r, const Direction& u,
                                         std::vector<Segment>& segments) {
  if (this->valid() == false) {
    auto mssg = "Cannot trace a Tile which is empty.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  double dist = 0.;

  if (c2d) {
    dist = c2d->trace_segments(r, u, segments);
  } else {
    dist = cell->trace_segments(r, u, segments);
  }

  return dist;
}

std::size_t Cartesian2D::Tile::num_fsrs() const {
  if (this->valid() == false) {
    return 0;
  }

  if (c2d) {
    return c2d->num_fsrs();
  } else {
    return cell->num_fsrs();
  }
}

void Cartesian2D::Tile::append_fsrs(std::vector<FlatSourceRegion*>& fsrs) {
  if (this->valid() == false) {
    return;
  }

  if (c2d) {
    c2d->append_fsrs(fsrs);
  } else {
    cell->append_fsrs(fsrs);
  }
}

}  // namespace scarabee
