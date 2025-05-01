#include <moc/cmfd.hpp>
#include <moc/moc_driver.hpp>
#include <utils/logging.hpp>
#include <utils/timer.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/constants.hpp>

#include <cmath>

namespace scarabee {

CMFD::CMFD(const std::vector<double>& dx, const std::vector<double>& dy,
           const std::vector<std::pair<std::size_t, std::size_t>>& groups)
    : dx_(dx),
      dy_(dy),
      x_bounds_(),
      y_bounds_(),
      moc_to_cmfd_group_map_(),
      group_condensation_(groups),
      nx_(dx_.size()),
      ny_(dy_.size()),
      ng_(groups.size()),
      nx_surfs_(),
      ny_surfs_(),
      temp_fsrs_(),
      fsrs_() {
  // Make sure we have at least 1 bin in each direction
  if (dx_.size() == 0) {
    auto mssg = "Must provide at least 1 x width.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (dy_.size() == 0) {
    auto mssg = "Must provide at least 1 y width.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure all widths > 0
  double dx_tot = 0.;
  for (std::size_t i = 0; i < dx_.size(); i++) {
    if (dx_[i] <= 0.) {
      std::stringstream mssg;
      mssg << "dx at index " << i << " is <= 0.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    dx_tot += dx_[i];
  }

  double dy_tot = 0.;
  for (std::size_t j = 0; j < dy_.size(); j++) {
    if (dy_[j] <= 0.) {
      std::stringstream mssg;
      mssg << "dy at index " << j << " is <= 0.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }

    dy_tot += dy_[j];
  }

  // Create surfaces
  x_bounds_.reserve(dx_.size() + 1);
  x_bounds_.emplace_back();
  x_bounds_.back().type() = Surface::Type::XPlane;
  x_bounds_.back().x0() = -0.5 * dx_tot;
  for (const auto& d : dx_) {
    const double new_x0 = x_bounds_.back().x0() + d;

    x_bounds_.emplace_back();
    x_bounds_.back().type() = Surface::Type::XPlane;
    x_bounds_.back().x0() = new_x0;
  }

  y_bounds_.reserve(dy_.size() + 1);
  y_bounds_.emplace_back();
  y_bounds_.back().type() = Surface::Type::YPlane;
  y_bounds_.back().y0() = -0.5 * dy_tot;
  for (const auto& d : dy_) {
    const double new_y0 = y_bounds_.back().y0() + d;

    y_bounds_.emplace_back();
    y_bounds_.back().type() = Surface::Type::YPlane;
    y_bounds_.back().y0() = new_y0;
  }

  nx_surfs_ = x_bounds_.size() * (y_bounds_.size() - 1);
  ny_surfs_ = y_bounds_.size() * (x_bounds_.size() - 1);

  // Allocate temp_fsrs_
  temp_fsrs_.resize(nx_ * ny_);

  // Check group condensation scheme
  if (group_condensation_.size() == 0) {
    auto mssg = "Empty energy condensation scheme provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (group_condensation_.front().first != 0) {
    auto mssg = "The energy condensation scheme does not start with 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (std::size_t i = 0; i < group_condensation_.size() - 1; i++) {
    if (group_condensation_[i].second + 1 != group_condensation_[i + 1].first) {
      std::stringstream mssg;
      mssg << "Condensed groups " << i << " and " << i + 1
           << " are not continuous.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
  }

  moc_to_cmfd_group_map_.resize(group_condensation_.back().second + 1);
  for (std::size_t g = 0; g < moc_to_cmfd_group_map_.size(); g++) {
    for (std::size_t G = 0; G < group_condensation_.size(); G++) {
      if (group_condensation_[G].first <= g &&
          g <= group_condensation_[G].second) {
        moc_to_cmfd_group_map_[g] = G;
        break;
      }
    }
  }

  // Allocate surfaces array
  surface_currents_ =
      xt::zeros<double>({group_condensation_.size(), nx_surfs_ + ny_surfs_});

  // Allocate the xs array
  xs_.resize({nx_, ny_});
  xs_.fill(nullptr);

  volumes_.resize(ng_ * nx_ * ny_);

  // Allocate the flux, Et, and D_transp_corr arrays
  flux_ = xt::zeros<double>({ng_, nx_, ny_});
  Et_ = xt::zeros<double>({ng_, nx_, ny_});
  D_transp_corr_ = xt::zeros<double>({ng_, nx_, ny_});
  Er_ = xt::zeros<double>({ng_, nx_, ny_});
  Chi_ = xt::zeros<double>({ng_, nx_, ny_});
  vEf_ = xt::zeros<double>({ng_, nx_, ny_});
  Es_ = xt::zeros<double>({ng_, ng_, nx_, ny_});
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

CMFDSurfaceCrossing CMFD::get_surface(const Vector& r, const Direction& u) const {
  CMFDSurfaceCrossing surface;

  // First, we try to get a tile
  auto otile = this->get_tile(r, u);
  if (otile.has_value() == false){
    surface.is_valid = false;
    return surface;
  } 

  const std::size_t i = (*otile)[0];
  const std::size_t j = (*otile)[1];

  surface.cell_index = tile_to_indx(*otile);

  // Now we get our surfaces for this tile 
  const auto& x_n = x_bounds_[i];
  const auto& x_p = x_bounds_[i+1];
  const auto& y_n = y_bounds_[j];
  const auto& y_p = y_bounds_[j+1];

  // Get the distances
  const double dx_n = std::abs(x_n.x0() - r.x());
  const double dx_p = std::abs(x_p.x0() - r.x());
  const double dy_n = std::abs(y_n.y0() - r.y());
  const double dy_p = std::abs(y_p.y0() - r.y());

  //start with corners, then sides
  //top right
  if (dx_p < SURFACE_COINCIDENT && dy_p < SURFACE_COINCIDENT) {
    surface.crossing = CMFDSurfaceCrossing::Type::TR;
  } else if (dx_p < SURFACE_COINCIDENT && dy_n < SURFACE_COINCIDENT) {
    surface.crossing = CMFDSurfaceCrossing::Type::BR;
  } else if (dx_n < SURFACE_COINCIDENT && dy_n < SURFACE_COINCIDENT) {
    surface.crossing = CMFDSurfaceCrossing::Type::BL;
  } else if (dx_n < SURFACE_COINCIDENT && dy_p < SURFACE_COINCIDENT) {
    surface.crossing = CMFDSurfaceCrossing::Type::TL;
  } else if (dx_p < SURFACE_COINCIDENT) {
    surface.crossing = CMFDSurfaceCrossing::Type::XP;
  } else if (dx_n < SURFACE_COINCIDENT) {
    surface.crossing = CMFDSurfaceCrossing::Type::XN;
  } else if (dy_p < SURFACE_COINCIDENT) {
    surface.crossing = CMFDSurfaceCrossing::Type::YP;
  } else if (dy_n < SURFACE_COINCIDENT) {
    surface.crossing = CMFDSurfaceCrossing::Type::YN;
  }
  
  else {
    surface.is_valid = false;
    return surface;
  }
  surface.is_valid = true;

  return surface;
}

std::size_t CMFD::tile_to_indx(const std::array<std::size_t, 2>& tile) const {
  return this->tile_to_indx(tile[0], tile[1]);
}

std::size_t CMFD::tile_to_indx(const std::size_t& i,
                               const std::size_t& j) const {
  return j * nx_ + i;
}

std::array<std::size_t, 2> CMFD::indx_to_tile(std::size_t cell_index){
  std::array<std::size_t, 2> tile;
  tile[0] = cell_index % nx_;
  tile[1] = (cell_index - tile[0] ) / nx_;
  return tile;
}

void CMFD::insert_fsr(const std::array<std::size_t, 2>& tile, std::size_t fsr) {
  // Compute linear index
  const std::size_t i = this->tile_to_indx(tile);
  temp_fsrs_.at(i).insert(fsr);
}

void CMFD::insert_fsr(std::size_t tile_indx, std::size_t fsr) {
  temp_fsrs_.at(tile_indx).insert(fsr);
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

void CMFD::tally_current(double aflx, const Direction& u, std::size_t G,
                         const CMFDSurfaceCrossing& surf) {
  if (G >= surface_currents_.shape()[0]) {
    auto mssg = "Group index out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
  
  const auto tile = indx_to_tile(surf.cell_index);
  std::size_t i = tile[0];
  std::size_t j = tile[1];


  //Get surface index(s) from CMFDSurfaceCrossing
  //need cell indexes
  htl::static_vector<std::size_t,4> surf_indexes;
  bool is_corner = false;
  if (surf.crossing == CMFDSurfaceCrossing::Type::XN){
    surf_indexes.push_back(get_x_neg_surf(i,j));
  } else if (surf.crossing == CMFDSurfaceCrossing::Type::XP){
    surf_indexes.push_back(get_x_pos_surf(i,j));
  } else if (surf.crossing == CMFDSurfaceCrossing::Type::YN){
    surf_indexes.push_back(get_y_neg_surf(i,j));
  } else if (surf.crossing == CMFDSurfaceCrossing::Type::YP){
    surf_indexes.push_back(get_y_pos_surf(i,j));
  } else if (surf.crossing == CMFDSurfaceCrossing::Type::TR){
    is_corner = true;
    surf_indexes.push_back(get_x_pos_surf(i,j));
    surf_indexes.push_back(get_y_pos_surf(i,j));
    if (i+1 < nx_){
      surf_indexes.push_back(get_y_pos_surf(i+1,j));
    }
    if (j+1 < ny_){
      surf_indexes.push_back(get_x_pos_surf(i,j+1));
    }
  } else if (surf.crossing == CMFDSurfaceCrossing::Type::BR){
    is_corner = true;
    surf_indexes.push_back(get_x_pos_surf(i,j));
    surf_indexes.push_back(get_y_neg_surf(i,j));
    if (i+1 < nx_){
      surf_indexes.push_back(get_y_neg_surf(i+1,j));
    }
    if (j != 0){
      surf_indexes.push_back(get_x_pos_surf(i,j-1));
    }
  } else if (surf.crossing == CMFDSurfaceCrossing::Type::BL){
    is_corner = true;
    surf_indexes.push_back(get_x_neg_surf(i,j));
    surf_indexes.push_back(get_y_neg_surf(i,j));
    if (i != 0){
      surf_indexes.push_back(get_y_neg_surf(i-1,j));
    }
    if (j != 0){
      surf_indexes.push_back(get_x_neg_surf(i,j-1));
    }
  } else if (surf.crossing == CMFDSurfaceCrossing::Type::TL){
    is_corner = true;
    surf_indexes.push_back(get_x_neg_surf(i,j));
    surf_indexes.push_back(get_y_pos_surf(i,j));
    if (i != 0){
      surf_indexes.push_back(get_y_pos_surf(i-1,j));
    }
    if (j+1 < ny_){
      surf_indexes.push_back(get_x_neg_surf(i,j+1));
    }
  }
  
  //Check surface indexes are not out of range
  for (std::size_t k=0; k < surf_indexes.size(); ++k){
    if (surf_indexes[k] >= surface_currents_.shape()[1]) {
      auto mssg = "Surface index out of range.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  if (is_corner){
    //Split flux between two surfaces evenly
    aflx *= 0.5;
  }

  for (auto si: surf_indexes){
    if (si < nx_surfs_) {
      #pragma omp atomic
          surface_currents_(G, si) += std::copysign(aflx, u.x());
        } else {
      #pragma omp atomic
          surface_currents_(G, si) += std::copysign(aflx, u.y());
        }
  }
}

void CMFD::normalize_currents() {
  // We must normalize the currents by the lengths of each surface

  // First, go through all y-levels, and normalize the x surfaces
  std::size_t s = 0;
  for (std::size_t j = 0; j < ny_; j++) {
    const double dy = y_bounds_[j + 1].y0() - y_bounds_[j].y0();
    const double invs_dy = 1. / dy;

    for (std::size_t xs = 0; xs < x_bounds_.size(); xs++) {
      for (std::size_t g = 0; g < ng_; g++)
        surface_currents_.at(g, s) *= invs_dy;
      s++;
    }
  }

  // Now go through all x-levels, and normalize the y surfaces
  for (std::size_t i = 0; i < nx_; i++) {
    const double dx = x_bounds_[i + 1].x0() - x_bounds_[i].x0();
    const double invs_dx = 1. / dx;

    for (std::size_t ys = 0; ys < y_bounds_.size(); ys++) {
      for (std::size_t g = 0; g < ng_; g++)
        surface_currents_.at(g, s) *= invs_dx;
      s++;
    }
  }

  surface_currents_normalized_ = true;
}

void CMFD::compute_homogenized_xs_and_flux(const MOCDriver& moc) {
  for (std::size_t i = 0; i < nx_; i++) {
    for (std::size_t j = 0; j < ny_; j++) {
      const auto indx = this->tile_to_indx(i, j);
      const auto fg_xs = moc.homogenize(fsrs_[indx]);
      const auto fg_dxs = fg_xs->diffusion_xs();
      const auto flux_spec = moc.homogenize_flux_spectrum(fsrs_[indx]);
      auto& xs = xs_(i, j);
      if (xs)
        *xs = *(fg_dxs->condense(group_condensation_, flux_spec));
      else
        xs = fg_dxs->condense(group_condensation_, flux_spec);

      //hopefully this is getting the correct fsr volumes
      double cell_volume = 0.;
      for (auto fsr: fsrs_[indx]){
        cell_volume += moc.volume(fsr);
      }

      // Generate the flux and Et values for the tile
      xt::view(flux_, xt::all(), i, j) = 0.;
      xt::view(Et_, xt::all(), i, j) = 0.;
      xt::view(D_transp_corr_,xt::all(),i,j) = 0.;
      xt::view(Er_, xt::all(),i,j) = 0.;
      xt::view(Chi_, xt::all(), i, j) = 0.;
      xt::view(vEf_, xt::all(), i, j) = 0.;
      xt::view(Es_, xt::all(), xt::all(), i, j) = 0.;
      for (std::size_t g = 0; g < moc_to_cmfd_group_map_.size(); g++) {
        flux_(moc_to_cmfd_group_map_[g], i, j) += flux_spec(g);
        Et_(moc_to_cmfd_group_map_[g], i, j) += fg_xs->Etr(g) * flux_spec(g);
      }
      xt::view(Et_, xt::all(), i, j) /= xt::view(flux_, xt::all(), i, j);

      //assign collapsed diffusion XS, g is CMFD group 
      for (std::size_t g = 0; g < ng_; g++){
        std::size_t linear_idx = g * nx_ * ny_ + indx;
        volumes_[linear_idx] = cell_volume;
        D_transp_corr_(g,i,j) = fg_dxs->D(g);
        Er_(g,i,j) = fg_dxs -> Er(g);
        Chi_(g,i,j) = fg_dxs -> chi(g);
        vEf_(g,i,j) = fg_dxs -> vEf(g);
        for (std::size_t gg = 0; gg < ng_; gg++){
          Es_(g,gg,i,j) = fg_dxs -> Es(g,gg);
        }
      }
    }
  }
}

void CMFD::check_neutron_balance(const std::size_t i, const std::size_t j, std::size_t g, const double keff) const {
  // First, we get the spacing of the tile
  const double dx = dx_[i];
  const double dy = dy_[j];

  // Next, get the surfaces for out tile
  const double J_xn = surface_currents_.at(g, j*x_bounds_.size() + i);
  const double J_xp = surface_currents_.at(g, j*x_bounds_.size() + i + 1);
  const double J_yn = surface_currents_.at(g, nx_surfs_ + i*y_bounds_.size() + j);
  const double J_yp = surface_currents_.at(g, nx_surfs_ + i*y_bounds_.size() + j + 1);

  // Get the xs for our tile
  const auto& xs = *xs_(i, j);
  const double chi_g_keff = xs.chi(g) / keff;

  // Compute the fission and scatter sources
  double fiss_source = 0.;
  double scat_source = 0.;
  for (std::size_t gg = 0; gg < ng_; gg++) {
    const double flx_gg = flux_.at(gg, i, j);
    fiss_source += chi_g_keff * xs.vEf(gg) * flx_gg;
    scat_source += xs.Es(gg, g) * flx_gg;
  }

  // Compute the total leakage from the currents
  const double leak_rate = ((J_xp - J_xn) / dx) + ((J_yp - J_yn) / dy);

  // Now compute the removal reaction rate
  const double tot_reac_rate = Et_(g, i, j) * flux_.at(g, i, j);

  // Compute the residual of the balance equation
  const double residual = leak_rate + tot_reac_rate - (scat_source + fiss_source);
  //const double req_leak = scat_source + fiss_source - tot_reac_rate;

  //if (std::abs(residual) >= 1.E-5) {
    spdlog::error("CMFD tile ({:d}, {:d}) in group {:d} has a neutron balance residual of {:.5E}.", i, j, g, residual);
    //spdlog::error("CMFD tile ({:d}, {:d}) in group {:d} requires leakage rate of {:.5E}.", i, j, g, req_leak);
  //}
  //debugging for testing
  //spdlog::error("CMFD tile ({:d}, {:d}) in group {:d} has a true leakage rate of {:.5E}", i, j, g, leak_rate);
  //spdlog::error("CMFD tile ({:d}, {:d}) in group {:d} has a leakage ratio of {:.5E}.", i, j, g, req_leak/leak_rate);
}

double CMFD::calc_surf_diffusion_coef(std::size_t i, std::size_t j, std::size_t g, std::size_t surf,const MOCDriver& moc) const {
  //Calculates the surface diffusion coefficent from cell average diffusion coefficient
  //surf is int: 0 for right, 1 for up, 2 for left, 3 for down
  //TODO: the order of checks can likely be rearranged to reduce repeated statements
  //TODO: May need to correct for optical thickness using Larsens effective diffusion coefficient
  double Ds;
  double Dij = D_transp_corr_(g,i,j);
  double dx_ij;
  double dx_s;
  if (surf==0){
    //x positive
    dx_ij = dx_[i];
    if (i+1 == nx_){
      auto x_max_bc = moc.x_max_bc();
      //If we are on x positive boundary
      if (x_max_bc == BoundaryCondition::Reflective){
        Ds = Dij;
        dx_s = dx_ij;
      } else if (x_max_bc == BoundaryCondition::Vacuum){
        //TODO Implement
        auto mssg = "Vacuum BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      } else if (x_max_bc == BoundaryCondition::Periodic){
        //TODO Implement
        auto mssg = "Periodic BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    } else {
      Ds = D_transp_corr_(g,i+1,j);
      dx_s = dx_[i+1];
    }
    
  } else if (surf == 1){
    //y positive
    dx_ij = dy_[j];
    if (j+1 == ny_){
      auto y_max_bc = moc.y_max_bc();
      //If we are on y positive boundary
      if (y_max_bc == BoundaryCondition::Reflective){
        Ds = Dij;
        dx_s = dx_ij;
      } else if (y_max_bc == BoundaryCondition::Vacuum){
        auto mssg = "Vacuum BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      } else if (y_max_bc == BoundaryCondition::Periodic){
        auto mssg = "Periodic BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    } else {
      Ds = D_transp_corr_(g,i,j+1);
      dx_s = dy_[j+1];
    }

  } else if (surf == 2){
    //x negative
    dx_ij = dx_[i];
    if (i == 0){
      auto x_min_bc = moc.x_min_bc();
      //If we are on x negative boundary
      if (x_min_bc == BoundaryCondition::Reflective){
        Ds= Dij;
        dx_s = dx_ij;
      } else if (x_min_bc == BoundaryCondition::Vacuum){
        auto mssg = "Vacuum BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      } else if (x_min_bc == BoundaryCondition::Periodic){
        auto mssg = "Periodic BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    } else {
      Ds = D_transp_corr_(g,i-1,j);
      dx_s = dx_[i-1];
    }

  } else if (surf == 3){
    //y negative
    dx_ij = dy_[j];
    if (j == 0){
      auto y_min_bc = moc.y_min_bc();
      //If we are on y negative boundary
      if (y_min_bc == BoundaryCondition::Reflective){
        Ds = Dij;
        dx_s = dx_ij;
      } else if (y_min_bc == BoundaryCondition::Vacuum){
        auto mssg = "Vacuum BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      } else if (y_min_bc == BoundaryCondition::Periodic){
        auto mssg = "Periodic BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    } else {
      Ds = D_transp_corr_(g,i,j-1);
      dx_s = dy_[j-1];
    }
  }

  //spdlog::info("calc_surf_diffusion_coef: Dij = {}, Ds = {}, dx_s = {}, dx_ij = {}", Dij, Ds, dx_s, dx_ij);

  //from OpenMOC documentation
  double D_surf = (2*Dij*Ds)/(Dij*dx_s + Ds*dx_ij);
  return D_surf;
}

double CMFD::calc_nonlinear_diffusion_coef(std::size_t i, std::size_t j, std::size_t g, std::size_t surf,
   double D_surf,const MOCDriver& moc) const {

  //Calculates nonlinear diffusion coefficients for a given surface 
  
  std::size_t cell_index = tile_to_indx(i,j);
  double flx_ij = flux_(g,i,j);
  double flx_s;
  double dx_ij;
  std::size_t surf_index;
  if (surf==0){
    dx_ij = dx_[i];
    surf_index = get_x_pos_surf(i,j);
    if (i+1 == nx_){
      auto x_max_bc = moc.x_max_bc();
      if (x_max_bc == BoundaryCondition::Reflective){
        flx_s = flx_ij;
      } else if (x_max_bc == BoundaryCondition::Vacuum){
        auto mssg = "Vacuum BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      } else if (x_max_bc == BoundaryCondition::Periodic){
        auto mssg = "Periodic BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    } else {
      flx_s = flux_(g,i+1,j);
    }
    
  } else if (surf == 1){
    dx_ij = dy_[j];
    surf_index = get_y_pos_surf(i,j);
    if (j+1 == ny_){
      auto y_max_bc = moc.y_max_bc();
      if (y_max_bc == BoundaryCondition::Reflective){
        flx_s = flx_ij;
      } else if (y_max_bc == BoundaryCondition::Vacuum){
        auto mssg = "Vacuum BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      } else if (y_max_bc == BoundaryCondition::Periodic){
        auto mssg = "Periodic BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    } else {
      flx_s = flux_(g,i,j+1);
    }
  } else if (surf == 2){
    dx_ij = dx_[i];
    surf_index = get_x_neg_surf(i,j);
    if (i == 0){
      auto x_min_bc = moc.x_min_bc();
      if (x_min_bc == BoundaryCondition::Reflective){
        flx_s = flx_ij;
      } else if (x_min_bc == BoundaryCondition::Vacuum){
        auto mssg = "Vacuum BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      } else if (x_min_bc == BoundaryCondition::Periodic){
        auto mssg = "Periodic BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    } else {
      //flx_s = moc.flux(next_cell,g);
      flx_s = flux_(g,i-1,j);
    }
    
  } else if (surf == 3){
    dx_ij = dy_[j];
    surf_index = get_y_neg_surf(i,j);
    if (j == 0){
      auto y_min_bc = moc.y_min_bc();
      if (y_min_bc == BoundaryCondition::Reflective){
        flx_s = flx_ij;
      } else if (y_min_bc == BoundaryCondition::Vacuum){
        auto mssg = "Vacuum BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      } else if (y_min_bc == BoundaryCondition::Periodic){
        auto mssg = "Periodic BC not implemented for CMFD";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    } else {
      flx_s = flux_(g,i,j-1);
    }
  }
  double current = surface_currents_(g, surf_index);
  
  //this might not be right for every direction, should check
  double D_nl = (-D_surf*(flx_s - flx_ij)-(current/dx_ij))/(flx_s + flx_ij);
  //spdlog::info("calc_surf_diffusion_coef: Dsurf = {}, Flx_next = {}, Flx_ij = {}, current = {}, dx_ij = {}", D_surf, flx_s, flx_ij, current, dx_ij);

  return D_nl;
}

void CMFD::create_loss_matrix(const MOCDriver& moc){

  //have to implement different boundary conditions
  //Reflective is simple
  std::size_t tot_cells = nx_ * ny_ ;
  M_.resize(ng_*tot_cells, ng_*tot_cells);
  std::vector<Eigen::Triplet<double>> global_triplets;

  //each equation is independent so it can be parallized, might be better to do so over cells rather than groups
#pragma omp parallel for
  for (std::size_t g=0; g < ng_; ++g){
    std::vector<Eigen::Triplet<double>> groupwise_vals;
    double Dxp;
    double Dyp;
    double Dxn;
    double Dyn;
    double Dnl_xp;
    double Dnl_yp;
    double Dnl_xn;
    double Dnl_yn;
    double flx_ij;
    double flx_yp;
    double flx_yn;
    double flx_xp;
    double flx_xn;

    for (std::size_t l=0; l < nx_ * ny_; ++l){
      auto [i, j] = indx_to_tile(l);
      std::size_t row_indx = g *tot_cells + l;

      double loss_xp;
      double loss_xn;
      double loss_yp;
      double loss_yn;

      //calculate surface diffusion coefficients
      Dxp = calc_surf_diffusion_coef(i,j,g,0,moc);
      Dyp = calc_surf_diffusion_coef(i,j,g,1,moc);
      Dxn = calc_surf_diffusion_coef(i,j,g,2,moc);
      Dyn = calc_surf_diffusion_coef(i,j,g,3,moc);

      //calculate nonlinear diffusion coefficients
      Dnl_xp = calc_nonlinear_diffusion_coef(i,j,g,0,Dxp,moc);
      Dnl_yp = calc_nonlinear_diffusion_coef(i,j,g,1,Dyp,moc);
      Dnl_xn = calc_nonlinear_diffusion_coef(i,j,g,2,Dxn,moc);
      Dnl_yn = calc_nonlinear_diffusion_coef(i,j,g,3,Dyn,moc);

      //spdlog::info("Test Dnl_xp {:f}",Dnl_xp);
      //spdlog::info("Test Dnl_yp {:f}",Dnl_yp);
      //spdlog::info("Test Dnl_xn {:f}",Dnl_xn);
      //spdlog::info("Test Dnl_yn {:f}",Dnl_yn);
      const double dx = dx_[i];
      const double dy = dy_[j];

      //not efficient to do this again after calculating is it
      flx_ij = flux_(g,i,j);

      double loss_ij = dy*(Dxn + Dnl_xn + Dxp - Dnl_xp) +
                dx*(Dyn + Dnl_yn + Dyp - Dnl_yp) +
                dx*dy*Er_(g,i,j)*flx_ij;
      groupwise_vals.emplace_back(row_indx,row_indx,loss_ij);

      //set x current diff for current cell
      //X pos boundary
      if (i+1 == nx_ ){
        if (moc.x_max_bc() == BoundaryCondition::Reflective){
          flx_xp = flx_ij;
          flx_xn = flux_(g,i-1,j);
          loss_xn = dy*(-Dxn + Dnl_xn)*flx_xn;
          groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i-1,j), loss_xn);
          loss_xp = dy*(Dxp + Dnl_xp)*flx_xp;
          groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i,j), loss_xp);
        }
      //X negative boundary
      } else if (i == 0){
          if (moc.x_min_bc() == BoundaryCondition::Reflective){
            flx_xn = flx_ij;
            flx_xp = flux_(g,i+1,j);
            loss_xn = dy*(-Dxn + Dnl_xn)*flx_xn;
            groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i,j), loss_xn);
            loss_xp = dy*(Dxp + Dnl_xp)*flx_xp;
            groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i+1,j), loss_xp);
          }
      } else {
        flx_xp = flux_(g,i+1,j);
        flx_xn = flux_(g,i-1,j);
        loss_xn = dy*(-Dxn + Dnl_xn)*flx_xn;
        groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i-1,j), loss_xn);
        loss_xp = dy*(Dxp + Dnl_xp)*flx_xp;
        groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i+1,j), loss_xp);
      }

      //set y current diff for current cell
      if (j + 1 == ny_){
        if (moc.y_max_bc() == BoundaryCondition::Reflective){
          flx_yp = flx_ij;
          flx_yn = flux_(g,i,j-1);
          loss_yn = dx*(-Dyn + Dnl_yn)*flx_yn;
          groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i,j-1), loss_yn);
          loss_yp = dx*(Dyp + Dnl_yp)*flx_yp;
          groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i,j), loss_yp);
        }
      } else if (j == 0){
        if (moc.y_min_bc() == BoundaryCondition::Reflective){
          flx_yn = flx_ij;
          flx_yp = flux_(g,i,j+1);
          loss_yn = dx*(-Dyn + Dnl_yn)*flx_yn;
          groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i,j), loss_yn);
          loss_yp = dx*(Dyp + Dnl_yp)*flx_yp;
          groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i,j+1), loss_yp);
        }
      } else {
        flx_yp = flux_(g,i,j+1);
        flx_yn = flux_(g,i,j-1);
        loss_yn = dx*(-Dyn + Dnl_yn)*flx_yn;
        groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i,j-1), loss_yn);
        loss_yp = dx*(Dyp + Dnl_yp)*flx_yp;
        groupwise_vals.emplace_back(row_indx,g*tot_cells + tile_to_indx(i,j+1), loss_yp);
      }

      for (std::size_t gg = 0; gg < ng_; ++gg){
        //subtract scattering source Es_(g_in, g_out)
        if (gg != g){
          groupwise_vals.emplace_back(row_indx, gg*tot_cells + l, -dx*dy*Es_(gg,g));
        }
      }
    }
    #pragma omp critical
        {
            global_triplets.insert(global_triplets.end(), groupwise_vals.begin(), groupwise_vals.end());
        }
  }
  M_.setFromTriplets(global_triplets.begin(), global_triplets.end());
  M_.makeCompressed();
}

void CMFD::create_source_matrix(){

  const std::size_t tot_cells = nx_ * ny_ ;
  QM_.resize(ng_*tot_cells, ng_*tot_cells);
  std::vector<Eigen::Triplet<double>> global_triplets;

  //not entirely sure this order can be used, but the source matrix needs to be ordered the same as the loss matrix
  //It makes sense from first glance because if we are in group g for cell l, we still want to sum over g' for that group g
#pragma omp parallel for
  for (std::size_t g=0; g<ng_; ++g){
    std::vector<Eigen::Triplet<double>> groupwise_vals;

    for (std::size_t l=0; l < nx_ * ny_; ++l){
      auto [i, j] = indx_to_tile(l);
      std::size_t row_indx = g *tot_cells + l;
      const double dx = dx_[i];
      const double dy = dy_[j];
      for (std::size_t gg=0; gg < ng_; ++gg){
        //Fission source
        groupwise_vals.emplace_back(row_indx,gg*tot_cells + l, dx*dy*Chi_(g,i,j)*vEf_(gg,i,j));
      }
    }
    #pragma omp critical
    {
        global_triplets.insert(global_triplets.end(), groupwise_vals.begin(), groupwise_vals.end());
    }
  }
  
  QM_.setFromTriplets(global_triplets.begin(), global_triplets.end());
  QM_.makeCompressed();
}

Eigen::VectorXd CMFD::flatten_flux() const {
  //Turns g, i, j flux to linearly indexed flux
  const std::size_t tot_cells = nx_ * ny_;
  Eigen::VectorXd flx_flat(ng_ * tot_cells);

  for (std::size_t g = 0; g < ng_; ++g) {
    //for each cell in row-major order (i fastest, then j)
    for (std::size_t j = 0; j < ny_; ++j) {
      for (std::size_t i = 0; i < nx_; ++i) {
        std::size_t cell_idx = tile_to_indx(i,j);
        std::size_t linear_idx = g * tot_cells + cell_idx;
        flx_flat[linear_idx] = flux_(g, i, j);
      }
    }
  }
  return flx_flat;
}

void CMFD::solve(MOCDriver& moc, double keff) {
  spdlog::info("Starting CMFD");
  spdlog::info("Normalizing current");
  this->normalize_currents();
  spdlog::info("Computing homogenized xs");
  this->compute_homogenized_xs_and_flux(moc);
  spdlog::info("Creating loss and source matrix");
  this->create_loss_matrix(moc);
  this->create_source_matrix();
  spdlog::info("Starting CMFD solve");
  spdlog::info("Keff passed to CMFD: {}", keff);

  //Debug info
  spdlog::info("M_ nonzeros: {}", M_.nonZeros());
  spdlog::info("QM_ nonzeros: {}", QM_.nonZeros());
  spdlog::info("QM_ sum: {}", QM_.sum());
  
  //Power Iteration to solve for Keff
  if (solve_ == 0){
    std::size_t max_iter = 1000; //just for initial testing, shouldn't use more than this.

    Eigen::VectorXd flux_moc = flatten_flux();

    // Initialize flux and source vectors
    Eigen::VectorXd flux(ng_*nx_*ny_);
    Eigen::VectorXd new_flux(ng_*nx_*ny_);
    Eigen::VectorXd Q(ng_*nx_*ny_);
    Eigen::VectorXd Q_new(ng_*nx_*ny_);

    // Initialize a vector for computing keff faster
    Eigen::VectorXd VvEf(ng_ * nx_ * ny_);
    for (std::size_t l = 0; l < nx_*ny_; l++) {
      auto [i, j] = indx_to_tile(l);
      for (std::size_t g = 0; g < ng_; g++) {
        VvEf(l + g * ny_*nx_) = volumes_[l + g*nx_*ny_]*vEf_(g,i,j);
      }
    }

    //Use moc homogenized flux as initial guess?
    //Maybe should just start with 1
    //flux = flux_moc;
    flux.fill(1.);
    flux.normalize();

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(M_);
    solver.factorize(M_);

    // Begin power iteration
    double keff_diff = 100.;
    double flux_diff = 100.;
    std::size_t iteration = 0;
    Timer iteration_timer;
    while (keff_diff > keff_tol_ || flux_diff > flux_tol_) {
      iteration_timer.reset();
      iteration_timer.start();
      iteration++;

      // Compute source vector
      Q = (1. / keff) * QM_ * flux;

      new_flux = solver.solve(Q);
      if (solver.info() != Eigen::Success) {
        spdlog::error("Solution impossible.");
        throw ScarabeeException("Solution impossible");
      }

      spdlog::info("Sum of source vector: {}", Q.sum());

      // Estiamte keff - not sure about this part, might need to be volume-weighed?  Q.dot(volumes)
      double prev_keff = keff;
      keff = prev_keff * VvEf.dot(new_flux) / VvEf.dot(flux);
      keff_diff = std::abs(keff - prev_keff) / keff;

      spdlog::info("Sum of old flux: {}", flux.sum());
      spdlog::info("Sum of new flux: {}", new_flux.sum());


      // Normalize our new flux
      new_flux *= prev_keff / keff;

      // Find the max flux error
      flux_diff = 0.;
      for (std::size_t i = 0; i < ng_ * nx_ * ny_; i++) {
        double flux_diff_i = std::abs(new_flux(i) - flux(i)) / new_flux(i);
        if (flux_diff_i > flux_diff) flux_diff = flux_diff_i;
      }
      flux = new_flux;

      // Write information
      spdlog::info("-----------------CMFD-----------------");
      spdlog::info("Iteration {:>4d}          keff: {:.5f}", iteration, keff);
      spdlog::info("     keff difference:     {:.5E}", keff_diff);
      spdlog::info("     max flux difference: {:.5E}", flux_diff);
      spdlog::info("     iteration time: {:.5E} s",
                  iteration_timer.elapsed_time());
      
      if (iteration > max_iter){
        auto mssg = "Max iterations exceeded, maybe problem can't converge";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
    }
  }
  
  
  /** 
  for (std::size_t i = 0; i < nx_; i++) {
    for (std::size_t j = 0; j < ny_; j++) {
      for (std::size_t g = 0; g < ng_; g++) {
        this->check_neutron_balance(i, j, g, keff);
      }
    }
  }
  */
}

}  // namespace scarabee