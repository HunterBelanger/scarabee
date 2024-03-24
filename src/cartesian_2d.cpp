#include <moc/cartesian_2d.hpp>
#include <utils/scarabee_exception.hpp>

#include <algorithm>
#include <optional>

Cartesian2D::Cartesian2D(const std::vector<std::shared_ptr<Surface>>& x_bounds,
                         const std::vector<std::shared_ptr<Surface>>& y_bounds)
    : x_bounds_(x_bounds), y_bounds_(y_bounds), tiles_(), nx_(), ny_() {
  // First, make sure we have at least 2 bounds in each direction
  if (x_bounds_.size() < 2) {
    throw ScarabeeException("Must provide at least 2 x-bounds");
  } else if (y_bounds_.size() < 2) {
    throw ScarabeeException("Must provide at least 2 y-bounds");
  }

  // Make sure bounds are sorted
  if (std::is_sorted(x_bounds_.begin(), x_bounds_.end(),
                     [](const auto& a, const auto& b) {
                       return a->x0() < b->x0();
                     }) == false) {
    throw ScarabeeException("Provided x-bounds are not sorted.");
  }

  if (std::is_sorted(y_bounds_.begin(), y_bounds_.end(),
                     [](const auto& a, const auto& b) {
                       return a->y0() < b->y0();
                     }) == false) {
    throw ScarabeeException("Provided y-bounds are not sorted.");
  }

  // Get sizes and allocate tiles array
  nx_ = x_bounds_.size() - 1;
  ny_ = y_bounds_.size() - 1;

  // All tiles start as uninitialized
  tiles_.resize({nx_, ny_});
  tiles_.fill(Tile{nullptr, std::nullopt});
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

void Cartesian2D::trace_segments(Vector& r, const Direction& u,
                                 std::vector<Segment>& segments) {
  if (this->tiles_valid() == false) {
    throw ScarabeeException(
        "Cannot trace segments on a Cartesian2D geometry with invalid tiles.");
  }

  // Get our current tile index
  auto ti = this->get_tile_index(r, u);

  // While we have a tile index
  while (ti) {
    // Get reference to the tile
    auto& tile = this->tile(*ti);

    // Trace segments on tile
    tile.trace_segments(r, u, segments);

    // Get new tile index
    ti = this->get_tile_index(r, u);
  }
}

Cartesian2D::Tile& Cartesian2D::tile(const Cartesian2D::TileIndex& ti) {
  if (ti.i >= nx_) {
    throw ScarabeeException("TileIndex i out of range.");
  }

  if (ti.j >= ny_) {
    throw ScarabeeException("TileIndex j out of range.");
  }

  return tiles_(ti.i, ti.j);
}

const Cartesian2D::Tile& Cartesian2D::tile(
    const Cartesian2D::TileIndex& ti) const {
  if (ti.i >= nx_) {
    throw ScarabeeException("TileIndex i out of range.");
  }

  if (ti.j >= ny_) {
    throw ScarabeeException("TileIndex j out of range.");
  }

  return tiles_(ti.i, ti.j);
}

bool Cartesian2D::tiles_valid() const {
  for (const auto& t : tiles_) {
    if (t.valid() == false) return false;
  }

  return true;
}

FlatSourceRegion& Cartesian2D::get_fsr(const Vector& r, const Direction& u) {
  auto ti = this->get_tile_index(r, u);

  if (ti.has_value() == false) {
    // We apparently are not in a tile
    std::stringstream mssg;
    mssg << "Position r = " << r << ", and Direction u = " << u
         << ", are not in Cartesian2D geometry.";
    throw ScarabeeException(mssg.str());
  }

  auto& t = this->tile(*ti);

  if (t.valid() == false) {
    std::stringstream mssg;
    mssg << "Tile for Position r = " << r << ", and Direction u = " << u
         << " is empty.";
    throw ScarabeeException(mssg.str());
  }

  if (t.c2d) {
    return t.c2d->get_fsr(r, u);
  } else {
    auto cell_get_fsr = [&r, &u ](auto& cell) -> auto& {
      return cell.get_fsr(r, u);
    };
    return std::visit(cell_get_fsr, t.cell.value());
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
    throw ScarabeeException(mssg.str());
  }

  const auto& t = this->tile(*ti);

  if (t.valid() == false) {
    std::stringstream mssg;
    mssg << "Tile for Position r = " << r << ", and Direction u = " << u
         << " is empty.";
    throw ScarabeeException(mssg.str());
  }

  if (t.c2d) {
    return t.c2d->get_fsr(r, u);
  } else {
    auto cell_get_fsr = [&r, &u ](const auto& cell) -> const auto& {
      return cell.get_fsr(r, u);
    };
    return std::visit(cell_get_fsr, t.cell.value());
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

void Cartesian2D::Tile::trace_segments(Vector& r, const Direction& u,
                                       std::vector<Segment>& segments) {
  if (this->valid() == false) {
    throw ScarabeeException("Cannot trace a Tile which is empty.");
  }

  if (c2d) {
    c2d->trace_segments(r, u, segments);
  } else {
    auto trace_lmbda = [&r, &u, &segments](auto& cell) {
      cell.trace_segments(r, u, segments);
    };
    std::visit(trace_lmbda, *cell);
  }
}

std::size_t Cartesian2D::Tile::num_fsrs() const {
  if (this->valid() == false) {
    return 0;
  }

  if (c2d) {
    return c2d->num_fsrs();
  } else {
    auto num_fsrs = [](auto& cell) { return cell.num_fsrs(); };
    return std::visit(num_fsrs, *cell);
  }
}

void Cartesian2D::Tile::append_fsrs(std::vector<FlatSourceRegion*>& fsrs) {
  if (this->valid() == false) {
    return;
  }

  if (c2d) {
    c2d->append_fsrs(fsrs);
  } else {
    auto append_fsrs = [&fsrs](auto& cell) { cell.append_fsrs(fsrs); };
    std::visit(append_fsrs, *cell);
  }
}