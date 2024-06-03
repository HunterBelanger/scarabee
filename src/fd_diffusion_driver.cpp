#include <diffusion/fd_diffusion_driver.hpp>

#include <Eigen/Sparse>

namespace scarabee {

void load_1d_matrix(DiffusionGeometry& geom, Eigen::SparseMatrix<double>& M) {
  // Allocate M
  M.resize(geom.ngroups()*geom.nmats(), geom.ngroups()*geom.nmats());

  for (std::size_t g = 0; g < geom.ngroups(); g++) {
    for (std::size_t m = 0; m < geom.nmats(); m++) {
      // Get material index
      const auto indxs = geom.mat_indxs(m);

      if (m != 0 && m != geom.nmats()-1) {
        // Get diffusion coefficients
        const double Dg_mm1 = geom.mat(m-1)->D(g);
        const double Dg_m   = geom.mat(m)->D(g);
        const double Dg_mp1 = geom.mat(m+1)->D(g);

        // Get widths
        const double dx_mm1 = geom.dx(indxs[0]-1);
        const double dx_m   = geom.dx(indxs[0]);
        const double dx_mp1 = geom.dx(indxs[0]+1);

        // Compute reduced diffusion coefficients
        const double dg_mm1 = Dg_mm1 / dx_mm1;
        const double dg_m   = Dg_m   / dx_m;
        const double dg_mp1 = Dg_mp1 / dx_mp1;

        // Get removal xs for m
        const double Er = geom.mat(m)->Er(g);

        // Now we compute the coefficient for each flux term
        const double a = (1. / dx_m) * (-2.*dg_m*dg_mp1 / (dg_m + dg_mp1)); // for flux m+1
        const double b = (1. / dx_m) * 2. * (((dg_m*dg_mp1)/(dg_m+dg_mp1)) + ((dg_m*dg_mm1)/(dg_m+dg_mm1))) + Er; // for flux m
        const double c = (1. / dx_m) * (-2.*dg_m*dg_mm1 / (dg_m + dg_mm1)); // for flux m-1

        // Set matrix components
        M.coeffRef(g*m,g*(m+1)) = a;
        M.coeffRef(g*m,g*m) = b;
        M.coeffRef(g*m,g*(m-1)) = c;
      } else if (m == 0) {
        // Left boundary
        DiffusionGeometry::Tile bc;
        std::optional<std::size_t> n;
        std::tie(bc, n) = geom.neighbor(0, DiffusionGeometry::Neighbor::XN);

        const double Dg_m  = geom.mat(m)->D(g);
        const double dx_m  = geom.dx(indxs[0]);
        const double dg_m  = Dg_m / dx_m;

        const double Dg_mp1  = geom.mat(m+1)->D(g);
        const double dx_mp1  = geom.dx(indxs[0]+1);
        const double dg_mp1  = Dg_mp1 / dx_mp1;

        // Get the albedo



      } else {
        // Right boundary
      }
    }

    M.makeCompressed();
  }
}

}