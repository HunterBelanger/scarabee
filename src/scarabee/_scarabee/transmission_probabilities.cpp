#include <transmission_probabilities.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>
#include <utils/math.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/timer.hpp>

#include <cmath>

namespace scarabee {

TransmissionProbabilities::TransmissionProbabilities(
    const std::vector<std::shared_ptr<CrossSection>>& xs,
    const std::vector<double>& dx, const std::vector<double>& dy,
    BoundaryCondition x_min_bc, BoundaryCondition x_max_bc,
    BoundaryCondition y_min_bc, BoundaryCondition y_max_bc)
    : dx_(),
      dy_(),
      xs_(),
      flux_incurrent_(),
      T_(),
      x_min_bc_(x_min_bc),
      x_max_bc_(x_max_bc),
      y_min_bc_(y_min_bc),
      y_max_bc_(y_max_bc),
      ngroups_(0),
      k_(1.),
      k_tol_(1.E-5),
      flux_tol_(1.E-5),
      solved_(false) {
  // Make sure we have valide widths
  if (dx.size() == 0) {
    auto mssg = "Must have at least one x width.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (const auto& delt_x : dx) {
    if (delt_x <= 0.) {
      const auto mssg = "All x widths must be > 0.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  dx_ = xt::zeros<double>({dx.size()});
  for (std::size_t i = 0; i < dx.size(); i++) dx_(i) = dx[i];

  if (dy.size() == 0) {
    auto mssg = "Must have at least one y width.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  for (const auto& delt_y : dy) {
    if (delt_y <= 0.) {
      const auto mssg = "All y widths must be > 0.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  dy_ = xt::zeros<double>({dy.size()});
  for (std::size_t j = 0; j < dy.size(); j++) dy_(j) = dy[j];

  if (xs.size() != (nx() * ny())) {
    auto mssg =
        "The number of cross sections does not agree with the number of x and "
        "y widths.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Check all xs objects are defined
  for (const auto& x : xs) {
    if (x == nullptr) {
      auto mssg = "All cross sections must be defined.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Make sure all xs objects have same number of energy groups
  ngroups_ = xs.front()->ngroups();
  for (const auto& x : xs) {
    if (x->ngroups() != this->ngroups()) {
      auto mssg = "Not all cross sections have the same number of groups.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
  }

  // Fill flux/incurrent array with zeros
  flux_incurrent_ = xt::zeros<double>({ngroups(), nx(), ny(), Jyp + 1});

  // Now we need to allocate and fill the xs array
  xs_.resize({nx(), ny()});
  xs_.fill(nullptr);

  std::size_t xs_ind = 0;
  for (int ij = static_cast<int>(ny()) - 1; ij >= 0; ij--) {
    const std::size_t j = static_cast<std::size_t>(ij);
    for (std::size_t i = 0; i < nx(); i++) {
      xs_(i, j) = xs[xs_ind];
      xs_ind++;
    }
  }
}

void TransmissionProbabilities::set_flux_tolerance(double ftol) {
  if (ftol <= 0.) {
    auto mssg = "Tolerance for flux must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ftol >= 0.1) {
    auto mssg = "Tolerance for flux must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  flux_tol_ = ftol;
}

void TransmissionProbabilities::set_keff_tolerance(double ktol) {
  if (ktol <= 0. || ktol >= 0.1) {
    auto mssg = "Tolerance for keff must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  k_tol_ = ktol;
}

double TransmissionProbabilities::flux(std::size_t i, std::size_t j,
                                       std::size_t g) const {
  if (i >= this->nx()) {
    auto mssg = "Index along x is out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (j >= this->nx()) {
    auto mssg = "Index along y is out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (g >= this->ngroups()) {
    auto mssg = "Energy group index is out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return flux_incurrent_(g, i, j, Flux);
}

xt::xtensor<double, 1> TransmissionProbabilities::flux_spectrum(
    std::size_t i, std::size_t j) const {
  if (i >= this->nx()) {
    auto mssg = "Index along x is out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (j >= this->nx()) {
    auto mssg = "Index along y is out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return xt::view(flux_incurrent_, xt::all(), i, j, Flux);
}

void TransmissionProbabilities::trace_tile(const std::size_t i,
                                           const std::size_t j,
                                           const std::uint32_t NA,
                                           const double delta_phi,
                                           const double d) {
  // Initialize volume and area variables
  double V = 0.;
  double Axm = 0.;
  double Axp = 0.;
  double Aym = 0.;
  double Ayp = 0.;

  // Get our cross section !
  const CrossSection& xs = *this->xs_(i, j);

  // Get cell dimensions !
  const double X = dx_(i);
  const double Y = dy_(j);

  const double invs_NA = 1. / static_cast<double>(NA);
  constexpr double PI_4 = PI / 4.;

  std::size_t Ntot = 0;

  // Do all angles
  for (std::uint32_t a = 0; a < NA/2; a++) {
    // Compute angle-dependent constants
    const double phi = delta_phi * (static_cast<double>(a) + 0.5);

    const double ux = phi < 0.5*PI ? std::cos(phi) : -std::cos(phi);
    const double uy = phi < 0.5*PI ? std::sin(phi) : -std::sin(phi);

    const double tmp_nx = (X * std::abs(uy)) / d;
    const double tmp_ny = (Y * std::abs(ux)) / d;
    const std::size_t ntot = static_cast<std::size_t>(std::round(tmp_nx + tmp_ny));
    Ntot += 2 * ntot; // Multiply by 2 for forwards and backwards directions

    const double dx = d / std::abs(uy);
    const double dy = d / std::abs(ux);

    double x0 = 0.;
    double y0 = (phi < 0.5 * PI) ? Y + 0.5 * dy : -0.5 * dy;

    double ddx = 0.;
    double ddy = std::abs(dy);

    // Do all tracks for the given angle
    for (std::size_t t = 0; t < ntot; t++) {
      double s = 0.;  // Will hold the track length
      bool xm = false, xp = false, ym = false,
           yp = false;  // Surfaces implicated

      if (phi < 0.5 * PI) {
        x0 = x0 + ddx;
        y0 = y0 - ddy;

        if (y0 < 0.) {
          const double f = std::abs(y0 / ddy);
          y0 = 0.;
          ddy = 0.;
          ddx = dx;
          x0 = f * ddx;
        }

        // Indicate starting surface (one or other)
        if (y0 == 0.)
          ym = true;
        else
          xm = true;

        const double sx = (X - x0) / ux;
        const double sy = (Y - y0) / uy;

        // Get track length and exit surface
        if (sx < sy) {
          xp = true;
          s = sx;
        } else {
          yp = true;
          s = sy;
        }
      } else { // phi >= 0.5 * PI
        x0 = x0 + ddx;
        y0 = y0 + ddy;

        if (y0 > Y) {
          const double f = std::abs((y0 - Y) / ddy);
          y0 = Y;
          ddy = 0.;
          ddx = dx;
          x0 = f * ddx;
        }

        // Indicate starting surface (one or other)
        if (y0 == Y)
          yp = true;
        else
          xm = true;

        const double sx = (X - x0) / ux;
        const double sy = (0. - y0) / uy;

        // Get track length and exit surface
        if (sx < sy) {
          xp = true;
          s = sx;
        } else {
          ym = true;
          s = sy;
        }
      }

      if (s < 0.) {
        const auto mssg = "Track length must be >= 0.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }

      // At this point, we have the track length, starting surface, and exit
      // surface. First, lets accumulate volume and areas.
      V += 2. * d * s * invs_NA; // Multiply by 2 because we only do half the angles [0,pi]
      if (xm) Axm += 2. * d * invs_NA / std::abs(ux);
      if (xp) Axp += 2. * d * invs_NA / std::abs(ux);
      if (ym) Aym += 2. * d * invs_NA / std::abs(uy);
      if (yp) Ayp += 2. * d * invs_NA / std::abs(uy);

      // Now, we need to iterate through all groups
      for (std::size_t g = 0; g < this->ngroups(); g++) {
        // Get the optical path
        const double tau = s * xs.Etr(g);
        const double Ki3_tau = Ki3(tau);

        // The i->i and s->i probabilities sum (pi/4 - Ki4(t))
        const double val = PI_4 - Ki3_tau;
        T_(g, i, j, i_i) += -val;

        if (xm) T_(g, i, j, xm_i) += val;
        if (xp) T_(g, i, j, xp_i) += val;
        if (ym) T_(g, i, j, ym_i) += val;
        if (yp) T_(g, i, j, yp_i) += val;

        if (xm && xp) T_(g, i, j, xm_xp) += Ki3_tau;
        if (xm && ym) T_(g, i, j, xm_ym) += Ki3_tau;
        if (xm && yp) T_(g, i, j, xm_yp) += Ki3_tau;

        if (xp && ym) T_(g, i, j, xp_ym) += Ki3_tau;
        if (xp && yp) T_(g, i, j, xp_yp) += Ki3_tau;

        if (ym && yp) T_(g, i, j, ym_yp) += Ki3_tau;
      }
    }
  }

  // We can now correct the values
  for (std::size_t g = 0; g < this->ngroups(); g++) {
    const double Etr = xs.Etr(g);

    T_(g, i, j, i_i) += static_cast<double>(Ntot) * Etr * V;
    T_(g, i, j, i_i) *= d * invs_NA / ((Etr * V) * (Etr * V));

    T_(g, i, j, xm_i) *= d * invs_NA / (2. * Etr * V * Axm);
    T_(g, i, j, xp_i) *= d * invs_NA / (2. * Etr * V * Axp);
    T_(g, i, j, ym_i) *= d * invs_NA / (2. * Etr * V * Aym);
    T_(g, i, j, yp_i) *= d * invs_NA / (2. * Etr * V * Ayp);

    T_(g, i, j, xm_xp) *= invs_NA / (Axm * Axp);
    T_(g, i, j, xm_ym) *= invs_NA / (Axm * Aym);
    T_(g, i, j, xm_yp) *= invs_NA / (Axm * Ayp);

    T_(g, i, j, xp_ym) *= invs_NA / (Axp * Aym);
    T_(g, i, j, xp_yp) *= invs_NA / (Axp * Ayp);

    T_(g, i, j, ym_yp) *= invs_NA / (Aym * Ayp);

    // Check values
    if (T_(g, i, j, i_i) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission ij -> ij.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, xm_i) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission xm -> ij.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, xp_i) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission xp -> ij.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, ym_i) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission ym -> ij.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, yp_i) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission yp -> ij.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, xm_xp) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission xm -> xp.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, xm_ym) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission xm -> ym.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, xm_yp) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission xm -> yp.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, xp_ym) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission xp -> ym.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, xp_yp) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission xp -> yp.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    } else if (T_(g, i, j, ym_yp) <= 0.) {
      std::stringstream mssg;
      mssg << "Negative transmission probability in group "  << g << ", tile (" << i << ", " << j << "), transmission ym -> yp.";
      spdlog::error(mssg.str());
      throw ScarabeeException(mssg.str());
    }
  }
}

void TransmissionProbabilities::solve_tile(
    const std::size_t g, const std::size_t i, const std::size_t j,
    xt::xtensor<double, 4>& flux_incurrent,
    const xt::xtensor<double, 3>& Q) const {
  // Get constants for cell (i,j)
  const auto T = xt::view(T_, g, i, j, xt::all());
  auto FlxJ = xt::view(flux_incurrent, g, i, j, xt::all());
  const double Q_g = Q(g, i, j);
  const double Ax = dx_(i);
  const double Ay = dy_(j);
  const double V = Ax * Ay;

  // First, get the new flux
  FlxJ(Flux) = T(i_i) * Q_g * V +
               4. * (T(xm_i) * FlxJ(Jxm) * Ax + T(xp_i) * FlxJ(Jxp) * Ax +
                     T(ym_i) * FlxJ(Jym) * Ay + T(yp_i) * FlxJ(Jyp) * Ay);

  // Now we compute each outgoing current
  const double j_xm_out = Q_g * V * T(xm_i) + 4. * (FlxJ(Jxp) * T(xm_xp) * Ax +
                                                    FlxJ(Jym) * T(xm_ym) * Ay +
                                                    FlxJ(Jyp) * T(xm_yp) * Ay);
  const double j_xp_out = Q_g * V * T(xp_i) + 4. * (FlxJ(Jxm) * T(xm_xp) * Ax +
                                                    FlxJ(Jym) * T(xp_ym) * Ay +
                                                    FlxJ(Jyp) * T(xp_yp) * Ay);
  const double j_ym_out = Q_g * V * T(ym_i) + 4. * (FlxJ(Jxm) * T(xm_ym) * Ax +
                                                    FlxJ(Jxp) * T(xp_ym) * Ax +
                                                    FlxJ(Jyp) * T(ym_yp) * Ay);
  const double j_yp_out = Q_g * V * T(yp_i) + 4. * (FlxJ(Jxm) * T(xm_yp) * Ax +
                                                    FlxJ(Jxp) * T(xp_yp) * Ax +
                                                    FlxJ(Jym) * T(ym_yp) * Ay);

  // Finally, assign outgoing current to in-currents for other tiles, keeping in
  // mind the boundary conditions
  const std::size_t I = nx() - 1;  // Max x index
  const std::size_t J = ny() - 1;  // Max y index

  // Jxm Outgoing Current
  if (i == 0) {
    switch (x_min_bc_) {
      case BoundaryCondition::Reflective:
        FlxJ(Jxm) = j_xm_out;  // In-coming flux is the out-going
        break;

      case BoundaryCondition::Vacuum:
        FlxJ(Jxm) = 0.;  // Zero incoming flux
        break;

      case BoundaryCondition::Periodic:
        flux_incurrent(g, I, j, Jxp) =
            j_xm_out;  // In-coming flux on the OTHER SIDE of the geometry is
                       // our ourgoing flux
        break;
    }
  } else {
    flux_incurrent(g, i - 1, j, Jxp) = j_xm_out;
  }

  // Jxp Outgoing Current
  if (i == I) {
    switch (x_max_bc_) {
      case BoundaryCondition::Reflective:
        FlxJ(Jxp) = j_xp_out;  // In-coming flux is the out-going
        break;

      case BoundaryCondition::Vacuum:
        FlxJ(Jxp) = 0.;  // Zero incoming flux
        break;

      case BoundaryCondition::Periodic:
        flux_incurrent(g, 0, j, Jxm) =
            j_xp_out;  // In-coming flux on the OTHER SIDE of the geometry is
                       // our ourgoing flux
        break;
    }
  } else {
    flux_incurrent(g, i + 1, j, Jxm) = j_xp_out;
  }

  // Jym Outgoing Current
  if (j == 0) {
    switch (y_min_bc_) {
      case BoundaryCondition::Reflective:
        FlxJ(Jym) = j_ym_out;  // In-coming flux is the out-going
        break;

      case BoundaryCondition::Vacuum:
        FlxJ(Jym) = 0.;  // Zero incoming flux
        break;

      case BoundaryCondition::Periodic:
        flux_incurrent(g, i, J, Jyp) =
            j_ym_out;  // In-coming flux on the OTHER SIDE of the geometry is
                       // our ourgoing flux
        break;
    }
  } else {
    flux_incurrent(g, i, j - 1, Jyp) = j_ym_out;
  }

  // Jyp Outgoing Current
  if (j == J) {
    switch (y_max_bc_) {
      case BoundaryCondition::Reflective:
        FlxJ(Jyp) = j_yp_out;  // In-coming flux is the out-going
        break;

      case BoundaryCondition::Vacuum:
        FlxJ(Jyp) = 0.;  // Zero incoming flux
        break;

      case BoundaryCondition::Periodic:
        flux_incurrent(g, i, 0, Jym) =
            j_yp_out;  // In-coming flux on the OTHER SIDE of the geometry is
                       // our ourgoing flux
        break;
    }
  } else {
    flux_incurrent(g, i, j + 1, Jym) = j_yp_out;
  }
}

void TransmissionProbabilities::generate_tracks(const std::uint32_t n_angles,
                                                const double d) {
  if (d <= 0. || d >= 0.1) {
    auto mssg = "Track spacing must be in the interval (0, 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (n_angles < 4 || (n_angles % 4) != 0) {
    auto mssg =
        "The number of angles must be >= 4 and also be a multiple of 4.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  const double delta_phi = 2. * PI / static_cast<double>(n_angles);

  solved_ = false;

  // Now we can allocate the probabilities array
  T_ = xt::zeros<double>({this->ngroups(), this->nx(), this->ny(), ym_yp + 1});

// We can now calculate all the proabilities in parallel
#pragma omp parallel for
  for (int ii = 0; ii < static_cast<int>(nx()); ii++) {
    const std::size_t i = static_cast<std::size_t>(ii);
    for (std::size_t j = 0; j < ny(); j++) {
      trace_tile(i, j, n_angles, delta_phi, d);
    }
  }
}

void TransmissionProbabilities::fill_src(
    xt::xtensor<double, 3>& Q, const xt::xtensor<double, 4>& flux_incurrent) {
  const double invs_k = 1. / k_;

#pragma omp parallel for
  for (int ig = 0; ig < static_cast<int>(this->ngroups()); ig++) {
    const std::size_t g = static_cast<std::size_t>(ig);

    for (std::size_t i = 0; i < this->nx(); i++) {
      for (std::size_t j = 0; j < this->ny(); j++) {
        const auto& xs = *xs_(i, j);

        const double chi_g = xs.chi(g);
        double Qout = 0.;

        for (std::size_t gg = 0; gg < this->ngroups(); gg++) {
          // Scatter source
          const double flux_gg_i = flux_incurrent(gg, i, j, Flux);
          const double Es_gg_to_g = xs.Es_tr(gg, g);
          Qout += Es_gg_to_g * flux_gg_i;

          // Fission Source
          const double vEf_gg = xs.vEf(gg);
          Qout += invs_k * chi_g * vEf_gg * flux_gg_i;
        }

        Q(g, i, j) = Qout;
      }
    }
  }
}

void TransmissionProbabilities::solve() {
  Timer sim_timer;
  sim_timer.start();

  // Make sure tiles have been traced.
  if (T_.size() == 0) {
    auto mssg =
        "Cannot solve Transmission Probabilities problem. Geometry has not "
        "been traced.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Create array for the initial flux / in-currents
  flux_incurrent_ = xt::zeros<double>({this->ngroups(), this->nx(), this->ny(), Jyp+1});
  // Initial flux in all tiles is 1
  xt::view(flux_incurrent_, xt::all(), xt::all(), xt::all(), Flux) = 1.;
  // Initial inward currents are all 1/4 on all OUTSIDE boundaries 
  xt::view(flux_incurrent_, xt::all(), 0, xt::all(), Jxm) = 0.25;
  xt::view(flux_incurrent_, xt::all(), this->nx()-1, xt::all(), Jxp) = 0.25;
  xt::view(flux_incurrent_, xt::all(), xt::all(), 0, Jym) = 0.25;
  xt::view(flux_incurrent_, xt::all(), xt::all(), this->ny()-1, Jym) = 0.25;

  // Create a copy for next iteration
  auto next_flux_incurrent = flux_incurrent_;

  // Initialize source array
  xt::xtensor<double, 3> src =
      xt::zeros<double>({this->ngroups(), this->nx(), this->ny()});

  double max_flx_diff = 100.;
  double rel_diff_keff = 100.;
  k_ = 1.;
  double prev_k = k_;

  std::size_t iteration = 0;
  Timer iteration_timer;
  while ((rel_diff_keff > k_tol_ || max_flx_diff > flux_tol_) && iteration < 100) {
    iteration_timer.reset();
    iteration_timer.start();
    iteration++;

    // At this point in loop, flux_incurrent_ should always be identical to next_flux_incurrent
    fill_src(src, flux_incurrent_);
    sweep(src, next_flux_incurrent);
    // Now next_flux_incurrent should be filled with new values

    prev_k = k_;
    k_ = calc_keff(next_flux_incurrent, flux_incurrent_);
    rel_diff_keff = std::abs(k_ - prev_k) / k_;

    // Get flux difference
    const auto next_flx =
        xt::view(next_flux_incurrent, xt::all(), xt::all(), xt::all(), Flux);
    const auto flx =
        xt::view(flux_incurrent_, xt::all(), xt::all(), xt::all(), Flux);
    max_flx_diff = xt::amax(xt::abs(next_flx - flx) / next_flx)();

    flux_incurrent_ = next_flux_incurrent;

    spdlog::info("-------------------------------------");
    spdlog::info("Iteration {:>4d}          keff: {:.5f}", iteration, k_);
    spdlog::info("     keff difference:     {:.5E}", rel_diff_keff);
    spdlog::info("     max flux difference: {:.5E}", max_flx_diff);
    spdlog::info("     iteration time: {:.5E} s",
                 iteration_timer.elapsed_time());
  }

  solved_ = true;
  sim_timer.stop();
  spdlog::info("");
  spdlog::info("Simulation Time: {:.5E} s", sim_timer.elapsed_time());
}

double TransmissionProbabilities::calc_keff(
    const xt::xtensor<double, 4>& new_flux_incurrent,
    const xt::xtensor<double, 4>& flux_incurrent) const {
  double num = 0.;
  double denom = 0.;

#pragma omp parallel
  {
    double num_thrd = 0.;
    double denom_thrd = 0.;

#pragma omp for
    for (int ig = 0; ig < static_cast<int>(this->ngroups()); ig++) {
      const std::size_t g = static_cast<std::size_t>(ig);

      for (std::size_t i = 0; i < this->nx(); i++) {
        for (std::size_t j = 0; j < this->ny(); j++) {
          const auto& xs = *xs_(i, j);
          const double VvEf = dx_(i) * dy_(j) * xs.vEf(g);
          const double new_flx = new_flux_incurrent(g, i, j, Flux);
          const double old_flx = flux_incurrent(g, i, j, Flux);

          num_thrd += VvEf * new_flx;
          denom_thrd += VvEf * old_flx;
        }
      }
    }

#pragma omp atomic
    num += num_thrd;

#pragma omp atomic
    denom += denom_thrd;
  }

  return k_ * num / denom;
}

void TransmissionProbabilities::sweep(
    const xt::xtensor<double, 3>& Q,
    xt::xtensor<double, 4>& flux_incurrent) const {
// We can do energy groups in parallel
#pragma omp parallel for
  for (int ig = 0; ig < static_cast<int>(this->ngroups()); ig++) {
    const std::size_t g = static_cast<std::size_t>(ig);

    // Do all red tiles
    for (std::size_t i = 0; i < this->nx(); i++) {
      for (std::size_t j = 0; j < this->ny(); j++) {
        if (((i + j) % 2) == 0) {
          this->solve_tile(g, i, j, flux_incurrent, Q);
        }
      }
    }

    // Do all black tiles
    for (std::size_t i = 0; i < this->nx(); i++) {
      for (std::size_t j = 0; j < this->ny(); j++) {
        if (((i + j) % 2) != 0) {
          this->solve_tile(g, i, j, flux_incurrent, Q);
        }
      }
    }
  }
}

}  // namespace scarabee