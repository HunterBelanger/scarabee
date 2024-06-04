#include <diffusion/fd_diffusion_driver.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/timer.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <cmath>
#include <sstream>

namespace scarabee {

void set_x_current_diff(DiffusionGeometry& geom, Eigen::SparseMatrix<double>& M, const std::size_t g, const std::size_t m) {
  // Get material index
  const auto indxs = geom.mat_indxs(m);

  // Get our left and right neighbors info
  DiffusionGeometry::Tile tile_mm1, tile_mp1;
  std::optional<std::size_t> op_mm1, op_mp1;

  std::tie(tile_mm1, op_mm1) = geom.neighbor(m, DiffusionGeometry::Neighbor::XN);
  std::tie(tile_mp1, op_mp1) = geom.neighbor(m, DiffusionGeometry::Neighbor::XP);

  // Get parameters for the tile in question
  const double D_m   = geom.mat(m)->D(g);
  const double dx_m   = geom.dx(indxs[0]);

  // Our fill method depends on if we have boundaries or not.
  if (op_mm1 && op_mp1) {
    // Both the left and right are materials. We can fill it like a normal tile
    const std::size_t mm1 = op_mm1.value();
    const std::size_t mp1 = op_mp1.value();

    // Get diffusion coefficients
    const double D_mm1 = geom.mat(mm1)->D(g);
    const double D_mp1 = geom.mat(mp1)->D(g);

    // Get widths
    const double dx_mm1 = geom.dx(indxs[0]-1);
    const double dx_mp1 = geom.dx(indxs[0]+1);

    // Now we compute the coefficient for each flux term
    const double a = (2. / dx_m) * (D_m*D_mp1 / (dx_mp1*D_m + dx_m*D_mp1)); // for flux m+1
    const double c = (2. / dx_m) * (D_m*D_mm1 / (dx_mm1*D_m + dx_m*D_mm1)); // for flux m-1
    const double b = -(a + c); // for flux m

    // Set matrix components
    M.coeffRef(m + g*geom.nmats(), mp1 + g*geom.nmats()) += a;
    M.coeffRef(m + g*geom.nmats(), m   + g*geom.nmats()) += b;
    M.coeffRef(m + g*geom.nmats(), mm1 + g*geom.nmats()) += c;
  } else if (op_mm1) {
    // Boundary condition on the right
    const std::size_t mm1 = op_mm1.value();
    const double D_mm1 = geom.mat(mm1)->D(g);
    const double dx_mm1 = geom.dx(indxs[0]-1);
    
    const double c = (2. / dx_m) * (D_m*D_mm1 / (dx_mm1*D_m + dx_m*D_mm1)); // for flux m-1
    const double b = -c - (2.*D_m/(dx_m*dx_m)); // for flux m

    M.coeffRef(m + g*geom.nmats(), m   + g*geom.nmats()) += b;
    M.coeffRef(m + g*geom.nmats(), mm1 + g*geom.nmats()) += c;
  } else if (op_mp1) {
    // Boundary condition on the left
    const std::size_t mp1 = op_mp1.value();
    const double D_mp1 = geom.mat(mp1)->D(g);
    const double dx_mp1 = geom.dx(indxs[0]+1);
    
    const double a = (2. / dx_m) * (D_m*D_mp1 / (dx_mp1*D_m + dx_m*D_mp1)); // for flux m+1
    const double b = -a - (2.*D_m/(dx_m*dx_m)); // for flux m

    M.coeffRef(m + g*geom.nmats(), mp1 + g*geom.nmats()) += a;
    M.coeffRef(m + g*geom.nmats(), m   + g*geom.nmats()) += b;
  } else {
    // Both sides are boundary conditions. Nothing we can do here.
    std::stringstream mssg;
    mssg << "Boundary condition on left and right of tile in x-direciton.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }
}

void load_1d_matrix(DiffusionGeometry& geom, Eigen::SparseMatrix<double>& M) {
  // Allocate M
  M.resize(geom.ngroups()*geom.nmats(), geom.ngroups()*geom.nmats());

  for (std::size_t g = 0; g < geom.ngroups(); g++) {
    for (std::size_t m = 0; m < geom.nmats(); m++) {
      set_x_current_diff(geom, M, g, m);
      
      // Get removal xs for m
      const double Er = geom.mat(m)->Er(g);
      
      // Add removal xs along the diagonal
      M.coeffRef(m + g*geom.nmats(), m + g*geom.nmats()) -= Er;
    }
  }
  M = -M; 
  M.makeCompressed();
}

void fill_source(const DiffusionGeometry& geom, const Eigen::VectorXd& flux, const double keff, Eigen::VectorXd& Q) {
  const double invs_keff = 1. / keff;

  for (std::size_t m = 0; m < geom.nmats(); m++) {
    const DiffusionCrossSection& xs = *geom.mat(m);

    for (std::size_t g = 0; g < geom.ngroups(); g++) {
      const double chi_g_keff = xs.chi(g) * invs_keff;
      double src = 0.;

      for (std::size_t gg = 0; gg < geom.ngroups(); gg++) {
        src += chi_g_keff * xs.vEf(gg) * flux(m + gg*geom.nmats());

        if (gg != g) {
          src += xs.Es(gg, g) * flux(m + gg*geom.nmats());
        }
      }

      Q(m + g*geom.nmats()) = src;
    }
  }
}

double calc_keff(const DiffusionGeometry& geom, const Eigen::VectorXd& flux, const Eigen::VectorXd& old_flux, const double keff) {
  double num = 0.;
  double denom = 0.;

  for (std::size_t m = 0; m < geom.nmats(); m++) {
    const double Vr = geom.volume(m);
    const auto& mat = *geom.mat(m);
    for (std::uint32_t g = 0; g < geom.ngroups(); g++) {
      const double VvEf = Vr * mat.vEf(g);
      const double flx = flux(m + g*geom.nmats());
      const double oflx = old_flux(m + g*geom.nmats());

      num += VvEf * flx;
      denom += VvEf * oflx;
    }
  }

  return keff * num / denom;
}

FDDiffusionDriver::FDDiffusionDriver(std::shared_ptr<DiffusionGeometry> geom) : geom_(geom), flux_() {
  if (geom_ == nullptr) {
    auto mssg = "FDDiffusionDriver provided with nullptr geometry.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

void FDDiffusionDriver::set_flux_tolerance(double ftol) {
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

void FDDiffusionDriver::set_keff_tolerance(double ktol) {
  if (ktol <= 0.) {
    auto mssg = "Tolerance for keff must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ktol >= 0.1) {
    auto mssg = "Tolerance for keff must be in the interval (0., 0.1).";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  keff_tol_ = ktol;
}

void FDDiffusionDriver::solve() {
  Timer sim_timer;
  sim_timer.start();

  spdlog::info("Solving for keff.");
  spdlog::info("keff tolerance: {:.5E}", keff_tol_);
  spdlog::info("Flux tolerance: {:.5E}", flux_tol_);


  // First, we create our loss matrix
  Eigen::SparseMatrix<double> M;

  // Load the loss matrix
  if (geom_->ndims() == 1) {
    load_1d_matrix(*geom_, M);
  } else if (geom_->ndims() == 2) {
    spdlog::error("2D Diffusion not yet supported.");
    throw ScarabeeException("2D Diffusion not yet supported.");
  } else {
    spdlog::error("3D Diffusion not yet supported.");
    throw ScarabeeException("3D Diffusion not yet supported.");
  }

  std::cout << "\n" << M << "\n";

  // Create a solver for the problem
  spdlog::info("Performing LU decomposition.");
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(M);
  if (solver.info() != Eigen::Success) {
    std::stringstream mssg;
    mssg << "Could not perform LU factorization on loss matrix.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Initialize flux and source vectors
  Eigen::VectorXd flux(geom_->ngroups()*geom_->nmats());
  Eigen::VectorXd new_flux(geom_->ngroups()*geom_->nmats());
  Eigen::VectorXd Q(geom_->ngroups()*geom_->nmats());

  flux.fill(1.);
  flux.normalize();
  //std::cout << "\n" << flux << "\n";

  // Begin power iteration
  double keff_diff = 100.;
  double flux_diff = 100.;
  std::size_t iteration = 0;
  Timer iteration_timer;
  while ((keff_diff > keff_tol_ || flux_diff > flux_tol_) && iteration < 100) {
    iteration_timer.reset();
    iteration_timer.start();
    iteration++;

    // Compute source vector
    fill_source(*geom_, flux, keff_, Q);
    std::cout << "\n" << Q << "\n";

    // Get new flux
    new_flux = solver.solve(Q);
    if (solver.info() != Eigen::Success) {
      spdlog::error("Solution impossible.");
      throw ScarabeeException("Solution impossible");
    }

    // Estiamte keff
    double prev_keff = keff_;
    keff_ = calc_keff(*geom_, new_flux, flux, prev_keff);
    keff_diff = std::abs(keff_ - prev_keff) / keff_;

    // Normalize our new flux
    new_flux *= prev_keff / keff_;

    // Find the max flux error
    flux_diff = 0.;
    for (std::size_t i = 0; i < geom_->ngroups()*geom_->nmats(); i++) {
      double flux_diff_i = std::abs(new_flux(i) - flux(i)) / new_flux(i);
      if (flux_diff_i > flux_diff)
        flux_diff = flux_diff_i;
    }
    flux = new_flux;
    
    //std::cout << "\n" << flux << "\n";

    // Write information
    spdlog::info("-------------------------------------");
    spdlog::info("Iteration {:>4d}          keff: {:.5f}", iteration, keff_);
    spdlog::info("     keff difference:     {:.5E}", keff_diff);
    spdlog::info("     max flux difference: {:.5E}", flux_diff);
    spdlog::info("     Iteration time: {:.5E} s",
                 iteration_timer.elapsed_time());
  }

  // Copy flux into the permanent xtensor array
  flux_.resize({geom_->ngroups()*geom_->nmats()});
  for (std::size_t i = 0; i < geom_->ngroups()*geom_->nmats(); i++) {
    flux_(i) = flux(i);
  }

  solved_ = true;

  sim_timer.stop();
  spdlog::info("");
  spdlog::info("Simulation Time: {:.5E} s", sim_timer.elapsed_time());
}

}