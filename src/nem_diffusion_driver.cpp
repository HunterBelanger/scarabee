#include <diffusion/nem_diffusion_driver.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/timer.hpp>
#include <utils/constants.hpp>

#include <array>
#include <cmath>
#include <cstdarg>

namespace scarabee {

NEMDiffusionDriver::NEMDiffusionDriver(std::shared_ptr<DiffusionGeometry> geom)
    : geom_(geom), NG_(geom_->ngroups()), NM_(geom_->nmats()) {
  if (geom_ == nullptr) {
    auto mssg = "NEMDiffusionDriver provided with nullptr geometry.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (geom_->ndims() != 3) {
    auto mssg = "NEMDiffusionDriver requires a 3D diffusion geometry.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

void NEMDiffusionDriver::set_flux_tolerance(double ftol) {
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

void NEMDiffusionDriver::set_keff_tolerance(double ktol) {
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

double NEMDiffusionDriver::calc_keff(
    double keff, const xt::xtensor<double, 2>& old_flux,
    const xt::xtensor<double, 2>& new_flux) const {
  double num = 0.;
  double denom = 0.;

  for (std::size_t m = 0; m < NM_; m++) {
    const double Vr = geom_->volume(m);
    const auto& mat = geom_->mat(m);
    for (std::size_t g = 0; g < NG_; g++) {
      const double VvEf = Vr * mat->vEf(g);
      const double nflx = new_flux(g, m);
      const double oflx = old_flux(g, m);

      num += VvEf * nflx;
      denom += VvEf * oflx;
    }
  }

  return keff_ * num / denom;
}

void NEMDiffusionDriver::fill_coupling_matrices() {
  spdlog::info("Loading coupling matrices");

  for (std::size_t m = 0; m < NM_; m++) {
    const auto geom_indx = geom_inds_(m);
    const double del_x = geom_->dx(geom_indx[0]);
    const double del_y = geom_->dy(geom_indx[1]);
    const double del_z = geom_->dz(geom_indx[2]);
    const auto& xs = *geom_->mat(m);

    for (std::size_t g = 0; g < NG_; g++) {
      const double D = xs.D(g);
      const double Er = xs.Er(g);
      const double dx = D / del_x;
      const double dy = D / del_y;
      const double dz = D / del_z;
      const double lx = 1. / (Er * del_x);
      const double ly = 1. / (Er * del_y);
      const double lz = 1. / (Er * del_z);

      // Create matrices
      const double ax = 1. + 32. * dx + 120. * dx * lx + 960. * dx * dx * lx +
                        840. * dx * dx * lx * lx;
      const double ay = 1. + 32. * dy + 120. * dy * ly + 960. * dy * dy * ly +
                        840. * dy * dy * ly * ly;
      const double az = 1. + 32. * dz + 120. * dz * lz + 960. * dz * dz * lz +
                        840. * dz * dz * lz * lz;
      const double ax1 = 8. * dx + 60. * dx * lx + 720. * dx * dx * lx +
                         840. * dx * dx * lx * lx;
      const double ay1 = 8. * dy + 60. * dy * ly + 720. * dy * dy * ly +
                         840. * dy * dy * ly * ly;
      const double az1 = 8. * dz + 60. * dz * lz + 720. * dz * dz * lz +
                         840. * dz * dz * lz * lz;
      const double xy = 20. * dx * ly + 840. * dx * dx * lx * ly;
      const double xz = 20. * dx * lz + 840. * dx * dx * lx * lz;
      const double yx = 20. * dy * lx + 840. * dy * dy * ly * lx;
      const double yz = 20. * dy * lz + 840. * dy * dy * ly * lz;
      const double zx = 20. * dz * lx + 840. * dz * dz * lz * lx;
      const double zy = 20. * dz * ly + 840. * dz * dz * lz * ly;
      Eigen::Matrix<double, 6, 6> A{
          {ax, ax1, xy, xy, xz, xz}, {ax1, ax, xy, xy, xz, xz},
          {yx, yx, ay, ay1, yz, yz}, {yx, yx, ay1, ay, yz, yz},
          {zx, zx, zy, zy, az, az1}, {zx, zx, zy, zy, az1, az}};

      const double bx = 1. - 32. * dx + 120. * dx * lx - 960. * dx * dx * lx +
                        840. * dx * dx * lx * lx;
      const double by = 1. - 32. * dy + 120. * dy * ly - 960. * dy * dy * ly +
                        840. * dy * dy * ly * ly;
      const double bz = 1. - 32. * dz + 120. * dz * lz - 960. * dz * dz * lz +
                        840. * dz * dz * lz * lz;
      const double bx1 = -8. * dx + 60. * dx * lx - 720. * dx * dx * lx +
                         840. * dx * dx * lx * lx;
      const double by1 = -8. * dy + 60. * dy * ly - 720. * dy * dy * ly +
                         840. * dy * dy * ly * ly;
      const double bz1 = -8. * dz + 60. * dz * lz - 720. * dz * dz * lz +
                         840. * dz * dz * lz * lz;
      Eigen::Matrix<double, 6, 6> B{
          {bx, bx1, xy, xy, xz, xz}, {bx1, bx, xy, xy, xz, xz},
          {yx, yx, by, by1, yz, yz}, {yx, yx, by1, by, yz, yz},
          {zx, zx, zy, zy, bz, bz1}, {zx, zx, zy, zy, bz1, bz}};

      const double cx =
          20. * dx * lx * del_x + 840. * dx * dx * lx * lx * del_x;
      const double cy =
          20. * dy * ly * del_y + 840. * dy * dy * ly * ly * del_y;
      const double cz =
          20. * dz * lz * del_z + 840. * dz * dz * lz * lz * del_z;
      const double cx1 = 60. * dx * lx * del_x;
      const double cy1 = 60. * dy * ly * del_y;
      const double cz1 = 60. * dz * lz * del_z;
      const double cx2 = 140. * dx * lx * del_x;
      const double cy2 = 140. * dy * ly * del_y;
      const double cz2 = 140. * dz * lz * del_z;
      Eigen::Matrix<double, 6, 7> C{
          {cx, cx1, 0., 0., cx2, 0., 0.}, {cx, -cx1, 0., 0., cx2, 0., 0.},
          {cy, 0., cy1, 0., 0., cy2, 0.}, {cy, 0., -cy1, 0., 0., cy2, 0.},
          {cz, 0., 0., cz1, 0., 0., cz2}, {cz, 0., 0., -cz1, 0., 0., cz2}};

      auto Ainvs = A.inverse();
      Rmats_(g, m) = Ainvs * B;
      Pmats_(g, m) = Ainvs * C;
    }
  }
}

void NEMDiffusionDriver::fill_source() {
  for (std::size_t m = 0; m < NM_; m++) {
    const auto& mat = geom_->mat(m);
    for (std::size_t g = 0; g < NG_; g++) {
      auto& Q = Q_(g, m);
      Q.fill(0.);

      const double chi_g = mat->chi(g);
      const double invs_keff = 1. / keff_;

      for (std::size_t gg = 0; gg < NG_; gg++) {
        const double vEf_gg = mat->vEf(gg);
        const double Es_gg_g = mat->Es(gg, g);
        const double flx_avg = flux_avg_(gg, m);
        const double flx_x1 = flux_x1_(gg, m);
        const double flx_x2 = flux_x2_(gg, m);
        const double flx_y1 = flux_y1_(gg, m);
        const double flx_y2 = flux_y2_(gg, m);
        const double flx_z1 = flux_z1_(gg, m);
        const double flx_z2 = flux_z2_(gg, m);

        Q(MomentIndx::AVG) += invs_keff * chi_g * vEf_gg * flx_avg;
        if (gg != g) Q(MomentIndx::AVG) += Es_gg_g * flx_avg;

        Q(MomentIndx::X1) += invs_keff * chi_g * vEf_gg * flx_x1;
        if (gg != g) Q(MomentIndx::X1) += Es_gg_g * flx_x1;

        Q(MomentIndx::X2) += invs_keff * chi_g * vEf_gg * flx_x2;
        if (gg != g) Q(MomentIndx::X2) += Es_gg_g * flx_x2;

        Q(MomentIndx::Y1) += invs_keff * chi_g * vEf_gg * flx_y1;
        if (gg != g) Q(MomentIndx::Y1) += Es_gg_g * flx_y1;

        Q(MomentIndx::Y2) += invs_keff * chi_g * vEf_gg * flx_y2;
        if (gg != g) Q(MomentIndx::Y2) += Es_gg_g * flx_y2;

        Q(MomentIndx::Z1) += invs_keff * chi_g * vEf_gg * flx_z1;
        if (gg != g) Q(MomentIndx::Z1) += Es_gg_g * flx_z1;

        Q(MomentIndx::Z2) += invs_keff * chi_g * vEf_gg * flx_z2;
        if (gg != g) Q(MomentIndx::Z2) += Es_gg_g * flx_z2;
      }
    }
  }
}

void NEMDiffusionDriver::fill_neighbors_and_geom_inds() {
  neighbors_.resize({NM_, 6});
  geom_inds_.resize({NM_});

  // Go through all mats
  for (std::size_t m = 0; m < NM_; m++) {
    geom_inds_(m) = geom_->geom_indx(m);

    neighbors_(m, 0) = geom_->neighbor(m, DiffusionGeometry::Neighbor::XP);
    neighbors_(m, 1) = geom_->neighbor(m, DiffusionGeometry::Neighbor::XN);
    neighbors_(m, 2) = geom_->neighbor(m, DiffusionGeometry::Neighbor::YP);
    neighbors_(m, 3) = geom_->neighbor(m, DiffusionGeometry::Neighbor::YN);
    neighbors_(m, 4) = geom_->neighbor(m, DiffusionGeometry::Neighbor::ZP);
    neighbors_(m, 5) = geom_->neighbor(m, DiffusionGeometry::Neighbor::ZN);
  }
}

void NEMDiffusionDriver::update_Jin_from_Jout(std::size_t g, std::size_t m) {
  // Get neighbor info
  const auto& n_xp = neighbors_(m, 0);
  const auto& n_xm = neighbors_(m, 1);
  const auto& n_yp = neighbors_(m, 2);
  const auto& n_ym = neighbors_(m, 3);
  const auto& n_zp = neighbors_(m, 4);
  const auto& n_zm = neighbors_(m, 5);

  // UPDATE INCOMING CURRENTS IN NEIGHBORING NODES / B.C.
  // x+ surface
  if (n_xp.second) {
    j_ins_(g, n_xp.second.value())(CurrentIndx::XM) =
        j_outs_(g, m)(CurrentIndx::XP);
  } else {
    const double albedo = n_xp.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::XP) = albedo * j_outs_(g, m)(CurrentIndx::XP);
  }

  // x- surface
  if (n_xm.second) {
    j_ins_(g, n_xm.second.value())(CurrentIndx::XP) =
        j_outs_(g, m)(CurrentIndx::XM);
  } else {
    const double albedo = n_xm.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::XM) = albedo * j_outs_(g, m)(CurrentIndx::XM);
  }

  // y+ surface
  if (n_yp.second) {
    j_ins_(g, n_yp.second.value())(CurrentIndx::YM) =
        j_outs_(g, m)(CurrentIndx::YP);
  } else {
    const double albedo = n_yp.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::YP) = albedo * j_outs_(g, m)(CurrentIndx::YP);
  }

  // y- surface
  if (n_ym.second) {
    j_ins_(g, n_ym.second.value())(CurrentIndx::YP) =
        j_outs_(g, m)(CurrentIndx::YM);
  } else {
    const double albedo = n_ym.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::YM) = albedo * j_outs_(g, m)(CurrentIndx::YM);
  }

  // z+ surface
  if (n_zp.second) {
    j_ins_(g, n_zp.second.value())(CurrentIndx::ZM) =
        j_outs_(g, m)(CurrentIndx::ZP);
  } else {
    const double albedo = n_zp.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::ZP) = albedo * j_outs_(g, m)(CurrentIndx::ZP);
  }

  // z- surface
  if (n_zm.second) {
    j_ins_(g, n_zm.second.value())(CurrentIndx::ZP) =
        j_outs_(g, m)(CurrentIndx::ZM);
  } else {
    const double albedo = n_zm.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::ZM) = albedo * j_outs_(g, m)(CurrentIndx::ZM);
  }
}

NEMDiffusionDriver::MomentsVector NEMDiffusionDriver::calc_leakage_moments(std::size_t g, std::size_t m) const {
  const auto& Jout = j_outs_(g, m);
  const auto& Jin = j_ins_(g, m);

  // Get the currents along each axis at the positive and negative bounds
  const double Jxp = calc_net_current(Jin, Jout, CurrentIndx::XP);
  const double Jxm = calc_net_current(Jin, Jout, CurrentIndx::XM);
  const double Jyp = calc_net_current(Jin, Jout, CurrentIndx::YP);
  const double Jym = calc_net_current(Jin, Jout, CurrentIndx::YM);
  const double Jzp = calc_net_current(Jin, Jout, CurrentIndx::ZP);
  const double Jzm = calc_net_current(Jin, Jout, CurrentIndx::ZM);

  // Compute average transverse leakages in each direction
  const double Lx = Jxp - Jxm;
  const double Ly = Jyp - Jym;
  const double Lz = Jzp - Jzm;

  // Get neighbor info
  const auto& n_xp = neighbors_(m, 0);
  const auto& n_xm = neighbors_(m, 1);
  const auto& n_yp = neighbors_(m, 2);
  const auto& n_ym = neighbors_(m, 3);
  const auto& n_zp = neighbors_(m, 4);
  const auto& n_zm = neighbors_(m, 5);

  // Obtain geometry spacings for node
  const auto geom_indxs = geom_inds_(m);
  const double dx = geom_->dx(geom_indxs[0]);
  const double dy = geom_->dy(geom_indxs[1]);
  const double dz = geom_->dz(geom_indxs[2]);
  const double invs_dx = 1. / dx;
  const double invs_dy = 1. / dy;
  const double invs_dz = 1. / dz;

  // This returns the average transverse leakage moments for a given node
  auto comp_avg_trans_lks = [this](std::size_t g, std::size_t m) {
    const auto& Jin = j_ins_(g, m);
    const auto& Jout = j_outs_(g, m);

    const double Jxp = calc_net_current(Jin, Jout, CurrentIndx::XP);
    const double Jxm = calc_net_current(Jin, Jout, CurrentIndx::XM);
    const double Jyp = calc_net_current(Jin, Jout, CurrentIndx::YP);
    const double Jym = calc_net_current(Jin, Jout, CurrentIndx::YM);
    const double Jzp = calc_net_current(Jin, Jout, CurrentIndx::ZP);
    const double Jzm = calc_net_current(Jin, Jout, CurrentIndx::ZM);

    const double Lx = Jxp - Jxm;
    const double Ly = Jyp - Jym;
    const double Lz = Jzp - Jzm;

    return std::array<double, 3>{Lx, Ly, Lz};
  };

  // Compute net moments
  MomentsVector L;
  L.fill(0.);

  std::array<double, 3> tmp;
  // x-axis
  if (n_xp.second && n_xm.second) {
    const double dx_xp = geom_->dx(geom_indxs[0] + 1);
    const double dx_xm = geom_->dx(geom_indxs[0] - 1);
    const double eta_xp = dx_xp * invs_dx;
    const double eta_xm = dx_xm * invs_dx;
    const double p1xm = eta_xm + 1.;
    const double p2xm = 2. * eta_xm + 1.;
    const double p1xp = eta_xp + 1.;
    const double p2xp = 2. * eta_xp + 1.;

    tmp = comp_avg_trans_lks(g, n_xp.second.value());
    const double Ly_xp = n_xp.second ? tmp[1] : 0.;
    const double Lz_xp = n_xp.second ? tmp[2] : 0.;

    tmp = comp_avg_trans_lks(g, n_xm.second.value());
    const double Ly_xm = n_xm.second ? tmp[1] : 0.;
    const double Lz_xm = n_xm.second ? tmp[2] : 0.;

    const double rho1yx = (p1xm * p2xm * Ly_xp - p1xp * p2xp * Ly_xm +
                           (p1xp * p2xp - p1xm * p2xm) * Ly) /
                          (p1xp * p1xm * (eta_xp + eta_xm + 1.));
    const double rho2yx =
        (p1xm * Ly_xp + p1xp * Ly_xm - (eta_xp + eta_xm + 2.) * Ly) /
        (p1xp * p1xm * (eta_xp + eta_xm + 1.));
    const double rho1zx = (p1xm * p2xm * Lz_xp - p1xp * p2xp * Lz_xm +
                           (p1xp * p2xp - p1xm * p2xm) * Lz) /
                          (p1xp * p1xm * (eta_xp + eta_xm + 1.));
    const double rho2zx =
        (p1xm * Lz_xp + p1xp * Lz_xm - (eta_xp + eta_xm + 2.) * Lz) /
        (p1xp * p1xm * (eta_xp + eta_xm + 1.));

    const double Lyx1 = rho1yx / 12.;
    const double Lyx2 = rho2yx / 20.;
    const double Lzx1 = rho1zx / 12.;
    const double Lzx2 = rho2zx / 20.;

    L(MomentIndx::X1) = invs_dy * Lyx1 + invs_dz * Lzx1;
    L(MomentIndx::X2) = invs_dy * Lyx2 + invs_dz * Lzx2;
  }

  // y-axis
  if (n_yp.second && n_ym.second) {
    const double dy_yp = geom_->dy(geom_indxs[1] + 1);
    const double dy_ym = geom_->dy(geom_indxs[1] - 1);
    const double eta_yp = dy_yp * invs_dy;
    const double eta_ym = dy_ym * invs_dy;
    const double p1ym = eta_ym + 1.;
    const double p2ym = 2. * eta_ym + 1.;
    const double p1yp = eta_yp + 1.;
    const double p2yp = 2. * eta_yp + 1.;

    tmp = comp_avg_trans_lks(g, n_yp.second.value());
    const double Lx_yp = n_yp.second ? tmp[0] : 0.;
    const double Lz_yp = n_yp.second ? tmp[2] : 0.;

    tmp = comp_avg_trans_lks(g, n_ym.second.value());
    const double Lx_ym = n_ym.second ? tmp[0] : 0.;
    const double Lz_ym = n_ym.second ? tmp[2] : 0.;

    const double rho1xy = (p1ym * p2ym * Lx_yp - p1yp * p2yp * Lx_ym +
                           (p1yp * p2yp - p1ym * p2ym) * Lx) /
                          (p1yp * p1ym * (eta_yp + eta_ym + 1.));
    const double rho2xy =
        (p1ym * Lx_yp + p1yp * Lx_ym - (eta_yp + eta_ym + 2.) * Lx) /
        (p1yp * p1ym * (eta_yp + eta_ym + 1.));
    const double rho1zy = (p1ym * p2ym * Lz_yp - p1yp * p2yp * Lz_ym +
                           (p1yp * p2yp - p1ym * p2ym) * Lz) /
                          (p1yp * p1ym * (eta_yp + eta_ym + 1.));
    const double rho2zy =
        (p1ym * Lz_yp + p1yp * Lz_ym - (eta_yp + eta_ym + 2.) * Lz) /
        (p1yp * p1ym * (eta_yp + eta_ym + 1.));

    const double Lxy1 = rho1xy / 12.;
    const double Lxy2 = rho2xy / 20.;
    const double Lzy1 = rho1zy / 12.;
    const double Lzy2 = rho2zy / 20.;

    L(MomentIndx::Y1) = invs_dx * Lxy1 + invs_dz * Lzy1;
    L(MomentIndx::Y2) = invs_dx * Lxy2 + invs_dz * Lzy2;
  }

  // z-axis
  if (n_zp.second && n_zm.second) {
    const double dz_zp = geom_->dz(geom_indxs[2] + 1);
    const double dz_zm = geom_->dz(geom_indxs[2] - 1);
    const double eta_zp = dz_zp * invs_dz;
    const double eta_zm = dz_zm * invs_dz;
    const double p1zm = eta_zm + 1.;
    const double p2zm = 2. * eta_zm + 1.;
    const double p1zp = eta_zp + 1.;
    const double p2zp = 2. * eta_zp + 1.;

    tmp = comp_avg_trans_lks(g, n_zp.second.value());
    const double Lx_zp = n_zp.second ? tmp[0] : 0.;
    const double Ly_zp = n_zp.second ? tmp[1] : 0.;

    tmp = comp_avg_trans_lks(g, n_zm.second.value());
    const double Lx_zm = n_zm.second ? tmp[0] : 0.;
    const double Ly_zm = n_zm.second ? tmp[1] : 0.;

    const double rho1xz = (p1zm * p2zm * Lx_zp - p1zp * p2zp * Lx_zm +
                           (p1zp * p2zp - p1zm * p2zm) * Lx) /
                          (p1zp * p1zm * (eta_zp + eta_zm + 1.));
    const double rho2xz =
        (p1zm * Lx_zp + p1zp * Lx_zm - (eta_zp + eta_zm + 2.) * Lx) /
        (p1zp * p1zm * (eta_zp + eta_zm + 1.));
    const double rho1yz = (p1zm * p2zm * Ly_zp - p1zp * p2zp * Ly_zm +
                           (p1zp * p2zp - p1zm * p2zm) * Ly) /
                          (p1zp * p1zm * (eta_zp + eta_zm + 1.));
    const double rho2yz =
        (p1zm * Ly_zp + p1zp * Ly_zm - (eta_zp + eta_zm + 2.) * Ly) /
        (p1zp * p1zm * (eta_zp + eta_zm + 1.));

    const double Lxz1 = rho1xz / 12.;
    const double Lxz2 = rho2xz / 20.;
    const double Lyz1 = rho1yz / 12.;
    const double Lyz2 = rho2yz / 20.;

    L(MomentIndx::Z1) = invs_dx * Lxz1 + invs_dy * Lxz1;
    L(MomentIndx::Z2) = invs_dx * Lxz2 + invs_dy * Lxz2;
  }

  // Return the transverse leakage moments vector
  return L;
}

void NEMDiffusionDriver::calc_node(const std::size_t g, const std::size_t m,
                                   const double invs_dx, const double invs_dy,
                                   const double invs_dz,
                                   const DiffusionCrossSection& xs) {
  //----------------------------------------------------------------------------
  // OBTAIN NECESSARY ARRAYS AND DATA
  const auto& Q = Q_(g, m);
  const auto& R = Rmats_(g, m);
  const auto& P = Pmats_(g, m);
  auto& Jout = j_outs_(g, m);  // Not const as we update this here
  auto& Jin = j_ins_(g, m);
  const double Er = xs.Er(g);  // Removal cross section
  const double D = xs.D(g);    // Diffusion coefficient

  //----------------------------------------------------------------------------
  // CALCULATE TRANSVERSE LEAKAGES AND THEIR MOMENTS

  // Get the currents along each axis at the positive and negative bounds
  double Jxp = calc_net_current(Jin, Jout, CurrentIndx::XP);
  double Jxm = calc_net_current(Jin, Jout, CurrentIndx::XM);
  double Jyp = calc_net_current(Jin, Jout, CurrentIndx::YP);
  double Jym = calc_net_current(Jin, Jout, CurrentIndx::YM);
  double Jzp = calc_net_current(Jin, Jout, CurrentIndx::ZP);
  double Jzm = calc_net_current(Jin, Jout, CurrentIndx::ZM);

  // Compute the transverse leakage moments vector
  MomentsVector L = calc_leakage_moments(g, m);

  //----------------------------------------------------------------------------
  // UPDATE OUTGOING PARTIAL CURRENTS
  Jout = R * Jin + P * (Q - L);

  //----------------------------------------------------------------------------
  // UPDATE FLUX AVERAGE AND MOMENTS
  // Get the currents along each axis at the positive and negative bounds
  Jxp = calc_net_current(Jin, Jout, CurrentIndx::XP);
  Jxm = calc_net_current(Jin, Jout, CurrentIndx::XM);
  Jyp = calc_net_current(Jin, Jout, CurrentIndx::YP);
  Jym = calc_net_current(Jin, Jout, CurrentIndx::YM);
  Jzp = calc_net_current(Jin, Jout, CurrentIndx::ZP);
  Jzm = calc_net_current(Jin, Jout, CurrentIndx::ZM);

  // Compute average transverse leakages in each direction
  const double Lx = Jxp - Jxm;
  const double Ly = Jyp - Jym;
  const double Lz = Jzp - Jzm;

  // Compute average T in each direction
  const double Tx = Jxp + Jxm;
  const double Ty = Jyp + Jym;
  const double Tz = Jzp + Jzm;

  // From the partial currents, compute the surface fluxes
  const double flx_xp = 2. * (Jout(CurrentIndx::XP) + Jin(CurrentIndx::XP));
  const double flx_xm = 2. * (Jout(CurrentIndx::XM) + Jin(CurrentIndx::XM));
  const double flx_yp = 2. * (Jout(CurrentIndx::YP) + Jin(CurrentIndx::YP));
  const double flx_ym = 2. * (Jout(CurrentIndx::YM) + Jin(CurrentIndx::YM));
  const double flx_zp = 2. * (Jout(CurrentIndx::ZP) + Jin(CurrentIndx::ZP));
  const double flx_zm = 2. * (Jout(CurrentIndx::ZM) + Jin(CurrentIndx::ZM));

  // Compute the new average flux
  const double flx_avg =
      (Q(MomentIndx::AVG) - invs_dx * Lx - invs_dy * Ly - invs_dz * Lz) / Er;
  flux_avg_(g, m) = flx_avg;

  // Calculate the first two polynomial coefficients along each direction
  const double ax1 = flx_xp - flx_xm;
  const double ax2 = flx_xp + flx_xm - 2. * flx_avg;
  const double ay1 = flx_yp - flx_ym;
  const double ay2 = flx_yp + flx_ym - 2. * flx_avg;
  const double az1 = flx_zp - flx_zm;
  const double az2 = flx_zp + flx_zm - 2. * flx_avg;

  // Calculate the first flux moments
  const double flx_x1 = (Q(MomentIndx::X1) - L(MomentIndx::X1) -
                         0.5 * invs_dx * Tx - invs_dx * invs_dx * D * ax1) /
                        Er;
  flux_x1_(g, m) = flx_x1;

  const double flx_y1 = (Q(MomentIndx::Y1) - L(MomentIndx::Y1) -
                         0.5 * invs_dy * Ty - invs_dy * invs_dy * D * ay1) /
                        Er;
  flux_y1_(g, m) = flx_y1;

  const double flx_z1 = (Q(MomentIndx::Z1) - L(MomentIndx::Z1) -
                         0.5 * invs_dz * Tz - invs_dz * invs_dz * D * az1) /
                        Er;
  flux_z1_(g, m) = flx_z1;

  // Calculate the second flux moments
  const double flx_x2 =
      (Q(MomentIndx::X2) - L(MomentIndx::X2) - 0.5 * invs_dx * Lx -
       3. * invs_dx * invs_dx * D * ax2) /
      Er;
  flux_x2_(g, m) = flx_x2;

  const double flx_y2 =
      (Q(MomentIndx::Y2) - L(MomentIndx::Y2) - 0.5 * invs_dy * Ly -
       3. * invs_dy * invs_dy * D * ay2) /
      Er;
  flux_y2_(g, m) = flx_y2;

  const double flx_z2 =
      (Q(MomentIndx::Z2) - L(MomentIndx::Z2) - 0.5 * invs_dz * Lz -
       3. * invs_dz * invs_dz * D * az2) /
      Er;
  flux_z2_(g, m) = flx_z2;

  //----------------------------------------------------------------------------
  // UPDATE INCOMING CURRENTS IN NEIGHBORING NODES / B.C.
  update_Jin_from_Jout(g, m);
}

void NEMDiffusionDriver::inner_iteration() {
  // Iterate through all nodes
  for (std::size_t m = 0; m < NM_; m++) {
    const auto geom_indx = geom_inds_(m);
    const double dx = geom_->dx(geom_indx[0]);
    const double dy = geom_->dy(geom_indx[1]);
    const double dz = geom_->dz(geom_indx[2]);
    const double invs_dx = 1. / dx;
    const double invs_dy = 1. / dy;
    const double invs_dz = 1. / dz;
    const auto& xs = *geom_->mat(m);

    for (std::size_t g = 0; g < NG_; g++) {
      calc_node(g, m, invs_dx, invs_dy, invs_dz, xs);
    }
  }
}

void NEMDiffusionDriver::solve() {
  Timer sim_timer;
  sim_timer.start();

  spdlog::info("Solving for keff.");
  spdlog::info("keff tolerance: {:.5E}", keff_tol_);
  spdlog::info("Flux tolerance: {:.5E}", flux_tol_);

  // Allocate all arrays
  flux_avg_.resize({NG_, NM_});
  xt::xtensor<double, 2> old_flux_avg = flux_avg_;
  flux_x1_.resize({NG_, NM_});
  flux_x2_.resize({NG_, NM_});
  flux_y1_.resize({NG_, NM_});
  flux_y2_.resize({NG_, NM_});
  flux_z1_.resize({NG_, NM_});
  flux_z2_.resize({NG_, NM_});
  j_outs_.resize({NG_, NM_});
  j_ins_.resize({NG_, NM_});
  Rmats_.resize({NG_, NM_});
  Pmats_.resize({NG_, NM_});
  Q_.resize({NG_, NM_});

  // Load the flux and current values with an initial guess
  flux_avg_ = xt::ones<double>({NG_, NM_});
  old_flux_avg = xt::ones<double>({NG_, NM_});
  flux_x1_ = xt::ones<double>({NG_, NM_});
  flux_x2_ = xt::ones<double>({NG_, NM_});
  flux_y1_ = xt::ones<double>({NG_, NM_});
  flux_y2_ = xt::ones<double>({NG_, NM_});
  flux_z1_ = xt::ones<double>({NG_, NM_});
  flux_z2_ = xt::ones<double>({NG_, NM_});
  for (std::size_t g = 0; g < NG_; g++) {
    for (std::size_t m = 0; m < NM_; m++) {
      j_outs_(g, m).fill(1.);
      j_ins_(g, m).fill(1.);
    }
  }

  // Fill the coupling matrices
  fill_neighbors_and_geom_inds();
  fill_coupling_matrices();

  // Begin power iteration
  double keff_diff = 100.;
  double flux_diff = 100.;
  std::size_t iteration = 0;
  Timer iteration_timer;
  while (keff_diff > keff_tol_ || flux_diff > flux_tol_) {
    iteration_timer.reset();
    iteration_timer.start();
    iteration++;

    // Copy current flux into old flux
    old_flux_avg = flux_avg_;

    // Calculating the source
    fill_source();

    // Perform 2 inner iterations per outer generation
    inner_iteration();
    inner_iteration();

    // Compute new keff
    double prev_keff = keff_;
    keff_ = calc_keff(prev_keff, old_flux_avg, flux_avg_);
    keff_diff = std::abs(keff_ - prev_keff) / keff_;

    // Find the max flux error
    flux_diff = xt::amax(xt::abs(flux_avg_ - old_flux_avg) / flux_avg_)();

    // Write information
    spdlog::info("-------------------------------------");
    spdlog::info("Iteration {:>4d}          keff: {:.5f}", iteration, keff_);
    spdlog::info("     keff difference:     {:.5E}", keff_diff);
    spdlog::info("     max flux difference: {:.5E}", flux_diff);
    spdlog::info("     iteration time: {:.5E} s",
                 iteration_timer.elapsed_time());
  }

  solved_ = true;

  sim_timer.stop();
  spdlog::info("");
  spdlog::info("Simulation Time: {:.5E} s", sim_timer.elapsed_time());

  // Deallocate arrays that are not needed for reconstruction
  Rmats_.resize({0, 0});
  Pmats_.resize({0, 0});
  Q_.resize({0, 0});

  // Calculate flux reconstruction parameters for each node
  spdlog::info("Fitting flux reconstruction parameters");
  recon_params.resize({NG_, NM_});
  for (std::size_t m = 0; m < NM_; m++) {
    for (std::size_t g = 0; g < NG_; g++) {
      recon_params(g, m) = fit_node_recon_params(g, m);
    }
  }
}

double NEMDiffusionDriver::flux(double x, double y, double z,
                                std::size_t g) const {
  // Check group index
  if (g >= ngroups()) {
    std::stringstream mssg;
    mssg << "Group index g = " << g << " is out of range.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  // Get geometry index
  const auto oi = geom_->x_to_i(x);
  const auto oj = geom_->y_to_j(y);
  const auto ok = geom_->z_to_k(z);
  if (oi.has_value() == false || oj.has_value() == false ||
      ok.has_value() == false)
    return 0.;
  const std::size_t i = oi.value();
  const std::size_t j = oj.value();
  const std::size_t k = ok.value();
  const xt::svector<std::size_t> geom_inds{i, j, k};

  // Get material index
  const auto om = geom_->geom_to_mat_indx(geom_inds);
  if (om.has_value() == false) return 0.;
  const std::size_t m = om.value();

  //return recon_params(g, m)(x, y, z);
  // Get the bounds and coordinates along each direction
  const double x_low = geom_->x_bounds()[i];
  const double x_hi = geom_->x_bounds()[i + 1];
  const double x_mid = 0.5 * (x_low + x_hi);
  const double xi_x = (x - x_mid) / (x_hi - x_low);

  const double y_low = geom_->y_bounds()[j];
  const double y_hi = geom_->y_bounds()[j + 1];
  const double y_mid = 0.5 * (y_low + y_hi);
  const double xi_y = (y - y_mid) / (y_hi - y_low);

  const double z_low = geom_->z_bounds()[k];
  const double z_hi = geom_->z_bounds()[k + 1];
  const double z_mid = 0.5 * (z_low + z_hi);
  const double xi_z = (z - z_mid) / (z_hi - z_low);

  // Construct flux coefficients along each direction
  const auto& Jout = j_outs_(g, m);
  const auto& Jin = j_ins_(g, m);
  const double flx_avg = flux_avg_(g, m);

  const double flx_xp = 2. * (Jout(CurrentIndx::XP) + Jin(CurrentIndx::XP));
  const double flx_xm = 2. * (Jout(CurrentIndx::XM) + Jin(CurrentIndx::XM));
  const double ax1 = flx_xp - flx_xm;
  const double ax2 = flx_xp + flx_xm - 2. * flx_avg;
  const double ax3 = -120. * flux_x1_(g, m) + 10. * ax1;
  const double ax4 = -700. * flux_x2_(g, m) + 35. * ax2;

  const double flx_yp = 2. * (Jout(CurrentIndx::YP) + Jin(CurrentIndx::YP));
  const double flx_ym = 2. * (Jout(CurrentIndx::YM) + Jin(CurrentIndx::YM));
  const double ay1 = flx_yp - flx_ym;
  const double ay2 = flx_yp + flx_ym - 2. * flx_avg;
  const double ay3 = -120. * flux_y1_(g, m) + 10. * ay1;
  const double ay4 = -700. * flux_y2_(g, m) + 35. * ay2;

  const double flx_zp = 2. * (Jout(CurrentIndx::ZP) + Jin(CurrentIndx::ZP));
  const double flx_zm = 2. * (Jout(CurrentIndx::ZM) + Jin(CurrentIndx::ZM));
  const double az1 = flx_zp - flx_zm;
  const double az2 = flx_zp + flx_zm - 2. * flx_avg;
  const double az3 = -120. * flux_z1_(g, m) + 10. * az1;
  const double az4 = -700. * flux_z2_(g, m) + 35. * az2;

  // Calculate flux
  const double flx_x = flx_avg + ax1 * f1(xi_x) + ax2 * f2(xi_x) +
                       ax3 * f3(xi_x) + ax4 * f4(xi_x);
  const double flx_y = flx_avg + ay1 * f1(xi_y) + ay2 * f2(xi_y) +
                       ay3 * f3(xi_y) + ay4 * f4(xi_y);
  const double flx_z = flx_avg + az1 * f1(xi_z) + az2 * f2(xi_z) +
                       az3 * f3(xi_z) + az4 * f4(xi_z);

  return flx_x * flx_y * flx_z / (flx_avg * flx_avg); 
  //return flx_x + flx_y + flx_z - 2.*flx_avg; 
}

xt::xtensor<double, 4> NEMDiffusionDriver::flux(
    const xt::xtensor<double, 1>& x, const xt::xtensor<double, 1>& y,
    const xt::xtensor<double, 1>& z) const {
  // Make sure x, y, and z have at least 1 coordinate
  if (x.size() == 0) {
    auto mssg = "Array of x coordinates must have at least one entry.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
  if (y.size() == 0) {
    auto mssg = "Array of y coordinates must have at least one entry.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
  if (z.size() == 0) {
    auto mssg = "Array of z coordinates must have at least one entry.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  xt::xtensor<double, 4> flux_out;
  flux_out.resize({ngroups(), x.size(), y.size(), z.size()});

  for (std::size_t g = 0; g < ngroups(); g++) {
    for (std::size_t i = 0; i < x.size(); i++) {
      for (std::size_t j = 0; j < y.size(); j++) {
        for (std::size_t k = 0; k < z.size(); k++) {
          flux_out(g, i, j, k) = this->flux(x(i), y(j), z(k), g);
        }
      }
    }
  }

  return flux_out;
}

xt::xtensor<double, 4> NEMDiffusionDriver::avg_flux() const {
  const std::size_t nx = geom_->nx();
  const std::size_t ny = geom_->ny();
  const std::size_t nz = geom_->nz();

  xt::xtensor<double, 4> flux_out;
  flux_out.resize({ngroups(), nx, ny, nz});

  for (std::size_t g = 0; g < ngroups(); g++) {
    for (std::size_t i = 0; i < nx; i++) {
      for (std::size_t j = 0; j < ny; j++) {
        for (std::size_t k = 0; k < nz; k++) {
          const auto om = geom_->geom_to_mat_indx({i, j, k});

          if (om.has_value() == false)
            flux_out(g, i, j, k) = 0.;
          else
            flux_out(g, i, j, k) = flux_avg_(g, om.value());
        }
      }
    }
  }

  return flux_out;
}

NEMDiffusionDriver::FluxRecon NEMDiffusionDriver::fit_node_recon_params(std::size_t g, std::size_t m) const {
  // Get node parameters
  const auto geom_indx = geom_->geom_indx(m); 
  const double dx = geom_->dx(geom_indx[0]);
  const double dy = geom_->dy(geom_indx[1]);
  const double dz = geom_->dz(geom_indx[2]);
  const auto& xs = *geom_->mat(m);
  const double D = xs.D(g);
  const double Er = xs.Er(g);

  const double kx = std::sqrt(dx*dx*Er/D);
  const double ky = std::sqrt(dy*dy*Er/D);

  const double x_low = geom_->x_bounds()[geom_indx[0]];
  const double x_hi = geom_->x_bounds()[geom_indx[0]+1];
  const double y_low = geom_->y_bounds()[geom_indx[1]];
  const double y_hi = geom_->y_bounds()[geom_indx[1]+1];
  const double z_low = geom_->z_bounds()[geom_indx[2]];
  const double z_hi = geom_->z_bounds()[geom_indx[2]+1];

  const auto& Jout = j_outs_(g, m);
  const auto& Jin = j_ins_(g, m);
  const double flx_avg = flux_avg_(g, m);
  const double Jxp = calc_net_current(Jin, Jout, CurrentIndx::XP);
  const double Jxm = calc_net_current(Jin, Jout, CurrentIndx::XM);
  const double Jyp = calc_net_current(Jin, Jout, CurrentIndx::YP);
  const double Jym = calc_net_current(Jin, Jout, CurrentIndx::YM);
  const double Jzp = calc_net_current(Jin, Jout, CurrentIndx::ZP);
  const double Jzm = calc_net_current(Jin, Jout, CurrentIndx::ZM);

  //const double xp = 1.; // Cooridnates in in scaled space
  //const double xm = -1.;
  //const double yp = 1.;
  //const double ym = -1.;
  const double xp = 0.5*dx;
  const double xm = -xp;
  const double yp = 0.5*dy;
  const double ym = -yp;

  const double flx_xp = 2. * (Jout(CurrentIndx::XP) + Jin(CurrentIndx::XP));
  const double flx_xm = 2. * (Jout(CurrentIndx::XM) + Jin(CurrentIndx::XM));
  const double flx_yp = 2. * (Jout(CurrentIndx::YP) + Jin(CurrentIndx::YP));
  const double flx_ym = 2. * (Jout(CurrentIndx::YM) + Jin(CurrentIndx::YM));

  const double flx_pp = avg_corner_flux(g, m, Corner::PP);
  const double flx_pm = avg_corner_flux(g, m, Corner::PM);
  const double flx_mp = avg_corner_flux(g, m, Corner::MP);
  const double flx_mm = avg_corner_flux(g, m, Corner::MM);

  // Define all functions
  F00 f00; F01 f01; F02 f02;
  F10 f10; F11 f11; F12 f12;
  F20 f20; F21 f21; F22 f22;

  //f01.fy.k = ky;
  //f02.fy.k = ky;

  //f10.fx.k = kx;
  //f11.fx.k = kx; f11.fy.k = ky;
  //f12.fx.k = kx; f12.fy.k = ky;

  //f20.fx.k = kx;
  //f21.fx.k = kx; f21.fy.k = ky;
  //f22.fx.k = kx; f22.fy.k = ky;

  Eigen::Matrix<double, 9, 1> b;
  b(0) = dx*dy*flx_avg; // Average flux

  //b(1) = (-dy/D) * Jxp; // Current on +x
  //b(2) = (-dy/D) * Jxm; // Current on -x
  //b(3) = (-dx/D) * Jyp; // Current on +y
  //b(4) = (-dx/D) * Jym; // Current on -y

  b(1) = dy * flx_xp;
  b(2) = dy * flx_xm;
  b(3) = dx * flx_yp;
  b(4) = dx * flx_ym;

  b(5) = flx_pp; // Corner flux at (xp, yp)
  b(6) = flx_pm; // Corner flux at (xp, ym)
  b(7) = flx_mp; // Corner flux at (xm, yp)
  b(8) = flx_mm; // Corner flux at (xm, ym)

  Eigen::Matrix<double, 9, 9> A 
      {{f00.ingr(dx, dy), f01.ingr(dx, dy), f02.ingr(dx, dy), f10.ingr(dx, dy), f11.ingr(dx, dy), f12.ingr(dx, dy), f20.ingr(dx, dy), f21.ingr(dx, dy), f22.ingr(dx, dy)},

       //{f00.dxiy(dy, xp), f01.dxiy(dy, xp), f02.dxiy(dy, xp), f10.dxiy(dy, xp), f11.dxiy(dy, xp), f12.dxiy(dy, xp), f20.dxiy(dy, xp), f21.dxiy(dy, xp), f22.dxiy(dy, xp)},
       //{f00.dxiy(dy, xm), f01.dxiy(dy, xm), f02.dxiy(dy, xm), f10.dxiy(dy, xm), f11.dxiy(dy, xm), f12.dxiy(dy, xm), f20.dxiy(dy, xm), f21.dxiy(dy, xm), f22.dxiy(dy, xm)},
       //{f00.dyix(dx, yp), f01.dyix(dx, yp), f02.dyix(dx, yp), f10.dyix(dx, yp), f11.dyix(dx, yp), f12.dyix(dx, yp), f20.dyix(dx, yp), f21.dyix(dx, yp), f22.dyix(dx, yp)},
       //{f00.dyix(dx, ym), f01.dyix(dx, ym), f02.dyix(dx, ym), f10.dyix(dx, ym), f11.dyix(dx, ym), f12.dyix(dx, ym), f20.dyix(dx, ym), f21.dyix(dx, ym), f22.dyix(dx, ym)},

       {f00.iy(dy, xp),   f01.iy(dy, xp),   f02.iy(dy, xp),   f10.iy(dy, xp),   f11.iy(dy, xp),   f12.iy(dy, xp),   f20.iy(dy, xp),   f21.iy(dy, xp),   f22.iy(dy, xp)},
       {f00.iy(dy, xm),   f01.iy(dy, xm),   f02.iy(dy, xm),   f10.iy(dy, xm),   f11.iy(dy, xm),   f12.iy(dy, xm),   f20.iy(dy, xm),   f21.iy(dy, xm),   f22.iy(dy, xm)},
       {f00.ix(dx, yp),   f01.ix(dx, yp),   f02.ix(dx, yp),   f10.ix(dx, yp),   f11.ix(dx, yp),   f12.ix(dx, yp),   f20.ix(dx, yp),   f21.ix(dx, yp),   f22.ix(dx, yp)},
       {f00.ix(dx, ym),   f01.ix(dx, ym),   f02.ix(dx, ym),   f10.ix(dx, ym),   f11.ix(dx, ym),   f12.ix(dx, ym),   f20.ix(dx, ym),   f21.ix(dx, ym),   f22.ix(dx, ym)},

       {f00(xp, yp),      f01(xp, yp),      f02(xp, yp),      f10(xp, yp),      f11(xp, yp),      f12(xp, yp),      f20(xp, yp),      f21(xp, yp),      f22(xp, yp)},
       {f00(xp, ym),      f01(xp, ym),      f02(xp, ym),      f10(xp, ym),      f11(xp, ym),      f12(xp, ym),      f20(xp, ym),      f21(xp, ym),      f22(xp, ym)},
       {f00(xm, yp),      f01(xm, yp),      f02(xm, yp),      f10(xm, yp),      f11(xm, yp),      f12(xm, yp),      f20(xm, yp),      f21(xm, yp),      f22(xm, yp)},
       {f00(xm, ym),      f01(xm, ym),      f02(xm, ym),      f10(xm, ym),      f11(xm, ym),      f12(xm, ym),      f20(xm, ym),      f21(xm, ym),      f22(xm, ym)}
      };

  // Calculate coefficients
  
  Eigen::FullPivLU<Eigen::Matrix<double, 9, 9>> lu(A);
  if (lu.isInvertible() == false) {
    spdlog::error("Could not invert coefficients matrix.");
    std::cout << A << "\n\n";
    std::cout << b << "\n\n";
    std::cout << "det(A) = " << A.determinant() << "\n";
    std::exit(1);
  }
  Eigen::Matrix<double, 9, 1> c = lu.solve(b);

  //Eigen::Matrix<double, 9, 1> c = A.colPivHouseholderQr().solve(b);

  FluxRecon out;
  out.radial[0] = c(0);
  out.radial[1] = c(1);
  out.radial[2] = c(2);
  out.radial[3] = c(3);
  out.radial[4] = c(4);
  out.radial[5] = c(5);
  out.radial[6] = c(6);
  out.radial[7] = c(7);
  out.radial[8] = c(8);

  out.f00 = f00;
  out.f01 = f01;
  out.f02 = f02;
  out.f10 = f10;
  out.f11 = f11;
  out.f12 = f12;
  out.f20 = f20;
  out.f21 = f21;
  out.f22 = f22;

  const double flx_zp = 2. * (Jout(CurrentIndx::ZP) + Jin(CurrentIndx::ZP));
  const double flx_zm = 2. * (Jout(CurrentIndx::ZM) + Jin(CurrentIndx::ZM));
  out.axial[0] = flx_avg;
  out.axial[1] = flx_zp - flx_zm;
  out.axial[2] = flx_zp + flx_zm - 2. * flx_avg;
  out.axial[3] = -120. * flux_z1_(g, m) + 10. * out.axial[1];
  out.axial[4] = -700. * flux_z2_(g, m) + 35. * out.axial[2];

  out.x_low = x_low;
  out.x_hi = x_hi;
  out.y_low = y_low;
  out.y_hi = y_hi;
  out.z_low = z_low;
  out.z_hi = z_hi;

  return out;
}

double NEMDiffusionDriver::eval_corner_flux(std::size_t g, std::size_t m, Corner c) const {
  const auto geom_inds = geom_inds_(m);
  const std::size_t i = geom_inds[0];
  const std::size_t j = geom_inds[1];
  const std::size_t k = geom_inds[2];

  // Get the bounds and coordinates along each direction
  const double xi_x = (c == Corner::PP || c == Corner::PM) ? 0.5 : -0.5;
  const double xi_y = (c == Corner::PP || c == Corner::MP) ? 0.5 : -0.5;
  const double xi_z = 0.;

  // Construct flux coefficients along each direction
  const auto& Jout = j_outs_(g, m);
  const auto& Jin = j_ins_(g, m);
  const double flx_avg = flux_avg_(g, m);

  const double flx_xp = 2. * (Jout(CurrentIndx::XP) + Jin(CurrentIndx::XP));
  const double flx_xm = 2. * (Jout(CurrentIndx::XM) + Jin(CurrentIndx::XM));
  const double ax1 = flx_xp - flx_xm;
  const double ax2 = flx_xp + flx_xm - 2. * flx_avg;
  const double ax3 = -120. * flux_x1_(g, m) + 10. * ax1;
  const double ax4 = -700. * flux_x2_(g, m) + 35. * ax2;

  const double flx_yp = 2. * (Jout(CurrentIndx::YP) + Jin(CurrentIndx::YP));
  const double flx_ym = 2. * (Jout(CurrentIndx::YM) + Jin(CurrentIndx::YM));
  const double ay1 = flx_yp - flx_ym;
  const double ay2 = flx_yp + flx_ym - 2. * flx_avg;
  const double ay3 = -120. * flux_y1_(g, m) + 10. * ay1;
  const double ay4 = -700. * flux_y2_(g, m) + 35. * ay2;

  const double flx_zp = 2. * (Jout(CurrentIndx::ZP) + Jin(CurrentIndx::ZP));
  const double flx_zm = 2. * (Jout(CurrentIndx::ZM) + Jin(CurrentIndx::ZM));
  const double az1 = flx_zp - flx_zm;
  const double az2 = flx_zp + flx_zm - 2. * flx_avg;
  const double az3 = -120. * flux_z1_(g, m) + 10. * az1;
  const double az4 = -700. * flux_z2_(g, m) + 35. * az2;

  // Calculate flux
  const double flx_x = flx_avg + ax1 * f1(xi_x) + ax2 * f2(xi_x) + ax3 * f3(xi_x) + ax4 * f4(xi_x);
  const double flx_y = flx_avg + ay1 * f1(xi_y) + ay2 * f2(xi_y) + ay3 * f3(xi_y) + ay4 * f4(xi_y);
  const double flx_z = flx_avg + az1 * f1(xi_z) + az2 * f2(xi_z) + az3 * f3(xi_z) + az4 * f4(xi_z);
  return flx_x * flx_y * flx_z / (flx_avg * flx_avg);
}

double NEMDiffusionDriver::avg_corner_flux(std::size_t g, std::size_t m, Corner c) const {
  const auto geom_inds = geom_inds_(m);

  double flux_sum = 0.;
  int npts = 0;

  // First, we add our contribution to the corner flux estimation
  flux_sum += eval_corner_flux(g, m, c);
  npts++;

  if (c == Corner::PP) {
    const auto& n_xp = neighbors_(m, 0);
    const auto& n_yp = neighbors_(m, 2);

    if (n_xp.second) {
      flux_sum += eval_corner_flux(g, n_xp.second.value(), Corner::MP);
      npts++;
    }
    if (n_yp.second) {
      flux_sum += eval_corner_flux(g, n_yp.second.value(), Corner::PM);
      npts++;
    }

    const auto om = geom_->geom_to_mat_indx({geom_inds[0]+1, geom_inds[1]+1, geom_inds[2]});
    if (om) {
      flux_sum += eval_corner_flux(g, om.value(), Corner::MM);
      npts++;
    }
  } else if (c == Corner::PM) {
    const auto& n_xp = neighbors_(m, 0);
    const auto& n_ym = neighbors_(m, 3);
    
    if (n_xp.second) {
      flux_sum += eval_corner_flux(g, n_xp.second.value(), Corner::MM);
      npts++;
    }
    if (n_ym.second) {
      flux_sum += eval_corner_flux(g, n_ym.second.value(), Corner::PP);
      npts++;
    }

    const auto om = geom_->geom_to_mat_indx({geom_inds[0]+1, geom_inds[1]-1, geom_inds[2]});
    if (om) {
      flux_sum += eval_corner_flux(g, om.value(), Corner::MP);
      npts++;
    }
  } else if (c == Corner::MM) {
    const auto& n_xm = neighbors_(m, 1);
    const auto& n_ym = neighbors_(m, 3);
    
    if (n_xm.second) {
      flux_sum += eval_corner_flux(g, n_xm.second.value(), Corner::PM);
      npts++;
    }
    if (n_ym.second) {
      flux_sum += eval_corner_flux(g, n_ym.second.value(), Corner::MP);
      npts++;
    }

    const auto om = geom_->geom_to_mat_indx({geom_inds[0]-1, geom_inds[1]-1, geom_inds[2]});
    if (om) {
      flux_sum += eval_corner_flux(g, om.value(), Corner::PP);
      npts++;
    }
  } else { // c = Corner::MP
    const auto& n_xm = neighbors_(m, 1);
    const auto& n_yp = neighbors_(m, 2);
    
    if (n_xm.second) {
      flux_sum += eval_corner_flux(g, n_xm.second.value(), Corner::PP);
      npts++;
    }
    if (n_yp.second) {
      flux_sum += eval_corner_flux(g, n_yp.second.value(), Corner::MM);
      npts++;
    }

    const auto om = geom_->geom_to_mat_indx({geom_inds[0]-1, geom_inds[1]+1, geom_inds[2]});
    if (om) {
      flux_sum += eval_corner_flux(g, om.value(), Corner::PM);
      npts++;
    }
  }

  return flux_sum / static_cast<double>(npts);
}

}  // namespace scarabee