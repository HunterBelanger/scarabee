#include <diffusion/nem_diffusion_driver.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/timer.hpp>

#include <array>
#include <cmath>
#include <cstdarg>

namespace scarabee {

double var_abs_max(int count, ...) {
  double mx = 0.;

  va_list args;
  va_start(args, count);

  for (int i = 0; i < count; i++) {
    const double val = va_arg(args, double);

    if (val > mx) mx = val;
  }

  va_end(args);

  return mx;
}

double NEMDiffusionDriver::calc_net_current(const Current& Jin, const Current& Jout, CurrentIndx indx) {
  if (indx == CurrentIndx::XP ||
      indx == CurrentIndx::YP ||
      indx == CurrentIndx::ZP) {
    return Jout(indx) - Jin(indx);
  } else {
    return Jin(indx) - Jout(indx);
  }
}

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

double NEMDiffusionDriver::calc_keff(double keff, const xt::xtensor<double, 2>& old_flux, const xt::xtensor<double, 2>& new_flux) const {
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
    const auto geom_indx = geom_->geom_indx(m);
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
      const double ax = 1. + 32.*dx + 120.*dx*lx + 960.*dx*dx*lx + 840.*dx*dx*lx*lx;
      const double ay = 1. + 32.*dy + 120.*dy*ly + 960.*dy*dy*ly + 840.*dy*dy*ly*ly;
      const double az = 1. + 32.*dz + 120.*dz*lz + 960.*dz*dz*lz + 840.*dz*dz*lz*lz;
      const double ax1 = 8.*dx + 60.*dx*lx + 720.*dx*dx*lx + 840.*dx*dx*lx*lx;
      const double ay1 = 8.*dy + 60.*dy*ly + 720.*dy*dy*ly + 840.*dy*dy*ly*ly;
      const double az1 = 8.*dz + 60.*dz*lz + 720.*dz*dz*lz + 840.*dz*dz*lz*lz;
      const double xy = 20.*dx*ly + 840.*dx*dx*lx*ly;
      const double xz = 20.*dx*lz + 840.*dx*dx*lx*lz;
      const double yx = 20.*dy*lx + 840.*dy*dy*ly*lx;
      const double yz = 20.*dy*lz + 840.*dy*dy*ly*lz;
      const double zx = 20.*dz*lx + 840.*dz*dz*lz*lx;
      const double zy = 20.*dz*ly + 840.*dz*dz*lz*ly;
      Eigen::Matrix<double, 6, 6> A {{ax,  ax1, xy,  xy,  xz,  xz },
                                     {ax1, ax,  xy,  xy,  xz,  xz },
                                     {yx,  yx,  ay,  ay1, yz,  yz },
                                     {yx,  yx,  ay1, ay,  yz,  yz },
                                     {zx,  zx,  zy,  zy,  az,  az1},
                                     {zx,  zx,  zy,  zy,  az1, az }};

      const double bx = 1. - 32.*dx + 120.*dx*lx - 960.*dx*dx*lx + 840.*dx*dx*lx*lx;
      const double by = 1. - 32.*dy + 120.*dy*ly - 960.*dy*dy*ly + 840.*dy*dy*ly*ly;
      const double bz = 1. - 32.*dz + 120.*dz*lz - 960.*dz*dz*lz + 840.*dz*dz*lz*lz;
      const double bx1 = -8.*dx + 60.*dx*lx - 720.*dx*dx*lx + 840.*dx*dx*lx*lx;
      const double by1 = -8.*dy + 60.*dy*ly - 720.*dy*dy*ly + 840.*dy*dy*ly*ly;
      const double bz1 = -8.*dz + 60.*dz*lz - 720.*dz*dz*lz + 840.*dz*dz*lz*lz;
      Eigen::Matrix<double, 6, 6> B {{bx,  bx1, xy,  xy,  xz,  xz },
                                     {bx1, bx,  xy,  xy,  xz,  xz },
                                     {yx,  yx,  by,  by1, yz,  yz },
                                     {yx,  yx,  by1, by,  yz,  yz },
                                     {zx,  zx,  zy,  zy,  bz,  bz1},
                                     {zx,  zx,  zy,  zy,  bz1, bz }};

      const double cx = 20.*dx*lx*del_x + 840.*dx*dx*lx*lx*del_x;
      const double cy = 20.*dy*ly*del_y + 840.*dy*dy*ly*ly*del_y;
      const double cz = 20.*dz*lz*del_z + 840.*dz*dz*lz*lz*del_z;
      const double cx1 = 60.*dx*lx*del_x;
      const double cy1 = 60.*dy*ly*del_y;
      const double cz1 = 60.*dz*lz*del_z;
      const double cx2 = 140.*dx*lx*del_x;
      const double cy2 = 140.*dy*ly*del_y;
      const double cz2 = 140.*dz*lz*del_z;
      Eigen::Matrix<double, 6, 7> C {{cx,  cx1, 0.,  0.,  cx2, 0.,  0.},
                                     {cx, -cx1, 0.,  0.,  cx2, 0.,  0.},
                                     {cy,  0.,  cy1, 0.,  0.,  cy2, 0.},
                                     {cy,  0., -cy1, 0.,  0.,  cy2, 0.},
                                     {cz,  0.,  0.,  cz1, 0.,  0.,  cz2},
                                     {cz,  0.,  0., -cz1, 0.,  0.,  cz2}};

      auto Ainvs = A.inverse();
      Rmats_(g, m) = Ainvs*B;
      Pmats_(g, m) = Ainvs*C;
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

void NEMDiffusionDriver::update_Jin_from_Jout(std::size_t g, std::size_t m) {
  // Get neighbor info
  const auto n_xp = geom_->neighbor(m, DiffusionGeometry::Neighbor::XP);
  const auto n_xm = geom_->neighbor(m, DiffusionGeometry::Neighbor::XN);
  const auto n_yp = geom_->neighbor(m, DiffusionGeometry::Neighbor::YP);
  const auto n_ym = geom_->neighbor(m, DiffusionGeometry::Neighbor::YN);
  const auto n_zp = geom_->neighbor(m, DiffusionGeometry::Neighbor::ZP);
  const auto n_zm = geom_->neighbor(m, DiffusionGeometry::Neighbor::ZN);

  // UPDATE INCOMING CURRENTS IN NEIGHBORING NODES / B.C.
  // x+ surface
  if (n_xp.second) {
    j_ins_(g, n_xp.second.value())(CurrentIndx::XM) = j_outs_(g, m)(CurrentIndx::XP);
  } else {
    const double albedo = n_xp.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::XP) = albedo * j_outs_(g, m)(CurrentIndx::XP);
  }

  // x- surface
  if (n_xm.second) {
    j_ins_(g, n_xm.second.value())(CurrentIndx::XP) = j_outs_(g, m)(CurrentIndx::XM);
  } else {
    const double albedo = n_xm.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::XM) = albedo * j_outs_(g, m)(CurrentIndx::XM);
  }

  // y+ surface
  if (n_yp.second) {
    j_ins_(g, n_yp.second.value())(CurrentIndx::YM) = j_outs_(g, m)(CurrentIndx::YP);
  } else {
    const double albedo = n_yp.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::YP) = albedo * j_outs_(g, m)(CurrentIndx::YP);
  }

  // y- surface
  if (n_ym.second) {
    j_ins_(g, n_ym.second.value())(CurrentIndx::YP) = j_outs_(g, m)(CurrentIndx::YM);
  } else {
    const double albedo = n_ym.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::YM) = albedo * j_outs_(g, m)(CurrentIndx::YM);
  }

  // z+ surface
  if (n_zp.second) {
    j_ins_(g, n_zp.second.value())(CurrentIndx::ZM) = j_outs_(g, m)(CurrentIndx::ZP);
  } else {
    const double albedo = n_zp.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::ZP) = albedo * j_outs_(g, m)(CurrentIndx::ZP);
  }

  // z- surface
  if (n_zm.second) {
    j_ins_(g, n_zm.second.value())(CurrentIndx::ZP) = j_outs_(g, m)(CurrentIndx::ZM);
  } else {
    const double albedo = n_zm.first.albedo.value();
    j_ins_(g, m)(CurrentIndx::ZM) = albedo * j_outs_(g, m)(CurrentIndx::ZM);
  }
}

NEMDiffusionDriver::MomentsVector NEMDiffusionDriver::calc_leakage_moments(std::size_t g, std::size_t m) const {
  const auto& Jout = j_outs_(g, m);
  const auto& Jin  = j_ins_(g, m);

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
  const auto n_xp = geom_->neighbor(m, DiffusionGeometry::Neighbor::XP);
  const auto n_xm = geom_->neighbor(m, DiffusionGeometry::Neighbor::XN);
  const auto n_yp = geom_->neighbor(m, DiffusionGeometry::Neighbor::YP);
  const auto n_ym = geom_->neighbor(m, DiffusionGeometry::Neighbor::YN);
  const auto n_zp = geom_->neighbor(m, DiffusionGeometry::Neighbor::ZP);
  const auto n_zm = geom_->neighbor(m, DiffusionGeometry::Neighbor::ZN);

  // Obtain geometry spacings for adjacent nodes
  const auto geom_indxs = geom_->geom_indx(m);
  const double dx = geom_->dx(geom_indxs[0]);
  const double dy = geom_->dy(geom_indxs[1]);
  const double dz = geom_->dz(geom_indxs[2]);
  const double invs_dx = 1. / dx;
  const double invs_dy = 1. / dy;
  const double invs_dz = 1. / dz;
  const double dx_xp = n_xp.second ? geom_->dx(geom_indxs[0]+1) : 0.;
  const double dx_xm = n_xm.second ? geom_->dx(geom_indxs[0]-1) : 0.;
  const double dy_yp = n_yp.second ? geom_->dy(geom_indxs[1]+1) : 0.;
  const double dy_ym = n_ym.second ? geom_->dy(geom_indxs[1]-1) : 0.;
  const double dz_zp = n_zp.second ? geom_->dz(geom_indxs[2]+1) : 0.;
  const double dz_zm = n_zm.second ? geom_->dz(geom_indxs[2]-1) : 0.;

  // Constants for computing transverse leakage coefficients
  const double eta_xp = dx_xp * invs_dx;
  const double eta_xm = dx_xm * invs_dx;
  const double eta_yp = dy_yp * invs_dy;
  const double eta_ym = dy_ym * invs_dy;
  const double eta_zp = dz_zp * invs_dz;
  const double eta_zm = dz_zm * invs_dz;
  const double p1xm =    eta_xm + 1.;
  const double p2xm = 2.*eta_xm + 1.;
  const double p1xp =    eta_xp + 1.;
  const double p2xp = 2.*eta_xp + 1.;
  const double p1ym =    eta_ym + 1.;
  const double p2ym = 2.*eta_ym + 1.;
  const double p1yp =    eta_yp + 1.;
  const double p2yp = 2.*eta_yp + 1.;
  const double p1zm =    eta_zm + 1.;
  const double p2zm = 2.*eta_zm + 1.;
  const double p1zp =    eta_zp + 1.;
  const double p2zp = 2.*eta_zp + 1.;

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

    return std::array<double, 3> {Lx, Ly, Lz};
  };

  // Compute average transverse leakages in each node
  std::array<double, 3> tmp;
  if (n_xp.second) tmp = comp_avg_trans_lks(g, n_xp.second.value());
  const double Ly_xp = n_xp.second ? tmp[1] : 0.;
  const double Lz_xp = n_xp.second ? tmp[2] : 0.;

  if (n_xm.second) tmp = comp_avg_trans_lks(g, n_xm.second.value());
  const double Ly_xm = n_xm.second ? tmp[1] : 0.;
  const double Lz_xm = n_xm.second ? tmp[2] : 0.;

  if (n_yp.second) tmp = comp_avg_trans_lks(g, n_yp.second.value());
  const double Lx_yp = n_yp.second ? tmp[0] : 0.;
  const double Lz_yp = n_yp.second ? tmp[2] : 0.;

  if (n_ym.second) tmp = comp_avg_trans_lks(g, n_ym.second.value());
  const double Lx_ym = n_ym.second ? tmp[0] : 0.;
  const double Lz_ym = n_ym.second ? tmp[2] : 0.;

  if (n_zp.second) tmp = comp_avg_trans_lks(g, n_zp.second.value());
  const double Lx_zp = n_zp.second ? tmp[0] : 0.;
  const double Ly_zp = n_zp.second ? tmp[1] : 0.;

  if (n_zm.second) tmp = comp_avg_trans_lks(g, n_zm.second.value());
  const double Lx_zm = n_zm.second ? tmp[0] : 0.;
  const double Ly_zm = n_zm.second ? tmp[1] : 0.;

  // x-axis
  const double rho1yx = (p1xm*p2xm*Ly_xp - p1xp*p2xp*Ly_xm + (p1xp*p2xp - p1xm*p2xm)*Ly) / (p1xp * p1xm * (eta_xp + eta_xm + 1.));
  const double rho2yx = (p1xm*Ly_xp + p1xp*Ly_xm - (eta_xp + eta_xm + 2.)*Ly) / (p1xp*p1xm*(eta_xp + eta_xm + 1.));
  const double rho1zx = (p1xm*p2xm*Lz_xp - p1xp*p2xp*Lz_xm + (p1xp*p2xp - p1xm*p2xm)*Lz) / (p1xp * p1xm * (eta_xp + eta_xm + 1.));
  const double rho2zx = (p1xm*Lz_xp + p1xp*Lz_xm - (eta_xp + eta_xm + 2.)*Lz) / (p1xp*p1xm*(eta_xp + eta_xm + 1.));

  const double Lyx1 = rho1yx / 12.;
  const double Lyx2 = rho2yx / 20.;
  const double Lzx1 = rho1zx / 12.;
  const double Lzx2 = rho2zx / 20.;

  // y-axis
  const double rho1xy = (p1ym*p2ym*Lx_yp - p1yp*p2yp*Lx_ym + (p1yp*p2yp - p1ym*p2ym)*Lx) / (p1yp * p1ym * (eta_yp + eta_ym + 1.));
  const double rho2xy = (p1ym*Lx_yp + p1yp*Lx_ym - (eta_yp + eta_ym + 2.)*Lx) / (p1yp*p1ym*(eta_yp + eta_ym + 1.));
  const double rho1zy = (p1ym*p2ym*Lz_yp - p1yp*p2yp*Lz_ym + (p1yp*p2yp - p1ym*p2ym)*Lz) / (p1yp * p1ym * (eta_yp + eta_ym + 1.));
  const double rho2zy = (p1ym*Lz_yp + p1yp*Lz_ym - (eta_yp + eta_ym + 2.)*Lz) / (p1yp*p1ym*(eta_yp + eta_ym + 1.));

  const double Lxy1 = rho1xy / 12.;
  const double Lxy2 = rho2xy / 20.;
  const double Lzy1 = rho1zy / 12.;
  const double Lzy2 = rho2zy / 20.;

  // z-axis
  const double rho1xz = (p1zm*p2zm*Lx_zp - p1zp*p2zp*Lx_zm + (p1zp*p2zp - p1zm*p2zm)*Lx) / (p1zp * p1zm * (eta_zp + eta_zm + 1.));
  const double rho2xz = (p1zm*Lx_zp + p1zp*Lx_zm - (eta_zp + eta_zm + 2.)*Lx) / (p1zp*p1zm*(eta_zp + eta_zm + 1.));
  const double rho1yz = (p1zm*p2zm*Ly_zp - p1zp*p2zp*Ly_zm + (p1zp*p2zp - p1zm*p2zm)*Ly) / (p1zp * p1zm * (eta_zp + eta_zm + 1.));
  const double rho2yz = (p1zm*Ly_zp + p1zp*Ly_zm - (eta_zp + eta_zm + 2.)*Ly) / (p1zp*p1zm*(eta_zp + eta_zm + 1.));

  const double Lxz1 = rho1xz / 12.;
  const double Lxz2 = rho2xz / 20.;
  const double Lyz1 = rho1yz / 12.;
  const double Lyz2 = rho2yz / 20.;

  // Compute net moments
  MomentsVector L;
  L(MomentIndx::AVG) = 0.;
  L(MomentIndx::X1) = invs_dy * Lyx1 + invs_dz * Lzx1;
  L(MomentIndx::X2) = invs_dy * Lyx2 + invs_dz * Lzx2;

  L(MomentIndx::Y1) = invs_dx * Lxy1 + invs_dz * Lzy1;
  L(MomentIndx::Y2) = invs_dx * Lxy2 + invs_dz * Lzy2;
  
  L(MomentIndx::Z1) = invs_dx * Lxz1 + invs_dy * Lxz1;
  L(MomentIndx::Z2) = invs_dx * Lxz2 + invs_dy * Lxz2;

  // Return the transverse leakage moments vector
  return L;
}

double NEMDiffusionDriver::calc_node(const std::size_t g, const std::size_t m, const xt::svector<std::size_t>& geom_indx,  const double invs_dx, const double invs_dy, const double invs_dz, const DiffusionCrossSection& xs) {
  //----------------------------------------------------------------------------
  // OBTAIN NECESSARY ARRAYS AND DATA
  const auto& Q = Q_(g, m);
  const auto& R = Rmats_(g, m);
  const auto& P = Pmats_(g, m);
  auto& Jout = j_outs_(g, m); // Not const as we update this here
  auto& Jin  = j_ins_(g, m);
  const double Er = xs.Er(g); // Removal cross section
  const double D = xs.D(g);   // Diffusion coefficient
  
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
  Jout = R*Jin + P*(Q - L);

  //update_Jin_from_Jout(g, m);
  
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
  const double old_flx_avg = flux_avg_(g, m);
  const double flx_avg = (Q(MomentIndx::AVG) - invs_dx*Lx - invs_dy*Ly - invs_dz*Lz) / Er;
  flux_avg_(g, m) = flx_avg;
  if (flx_avg <= 0.) {
    std::cout << ">> flux(g = " << g << ", m = " << m << ") = " << flx_avg << "\n";
    std::cout << ">> Q = " << Q(MomentIndx::AVG) << "\n";
    std::cout << ">> dx = " << 1. / invs_dx << "\n";
    std::cout << ">> dy = " << 1. / invs_dy << "\n";
    std::cout << ">> dy = " << 1. / invs_dz << "\n";
    std::cout << ">> Lx = " << Lx << "\n";
    std::cout << ">> Ly = " << Ly << "\n";
    std::cout << ">> Lz = " << Lz << "\n";
    std::cout << ">> Er = " << Er << "\n\n";
  }

  // Calculate the first two polynomial coefficients along each direction
  const double ax1 = flx_xp - flx_xm;
  const double ax2 = flx_xp + flx_xm - 2.*flx_avg;
  const double ay1 = flx_yp - flx_ym;
  const double ay2 = flx_yp + flx_ym - 2.*flx_avg;
  const double az1 = flx_zp - flx_zm;
  const double az2 = flx_zp + flx_zm - 2.*flx_avg;

  // Calculate the first flux moments
  const double old_flx_x1 = flux_x1_(g, m);
  const double flx_x1 = (Q(MomentIndx::X1) - L(MomentIndx::X1) - 0.5*invs_dx*Tx - invs_dx*invs_dx*D*ax1)/Er;
  flux_x1_(g, m) = flx_x1;

  const double old_flx_y1 = flux_y1_(g, m);
  const double flx_y1 = (Q(MomentIndx::Y1) - L(MomentIndx::Y1) - 0.5*invs_dy*Ty - invs_dy*invs_dy*D*ay1)/Er;
  flux_y1_(g, m) = flx_y1;

  const double old_flx_z1 = flux_z1_(g, m);
  const double flx_z1 = (Q(MomentIndx::Z1) - L(MomentIndx::Z1) - 0.5*invs_dz*Tz - invs_dz*invs_dz*D*az1)/Er;
  flux_z1_(g, m) = flx_z1;

  // Calculate the second flux moments
  const double old_flx_x2 = flux_x2_(g, m);
  const double flx_x2 = (Q(MomentIndx::X2) - L(MomentIndx::X2) - 0.5*invs_dx*Lx - 3.*invs_dx*invs_dx*D*ax2)/Er;
  flux_x2_(g, m) = flx_x2;

  const double old_flx_y2 = flux_y2_(g, m);
  const double flx_y2 = (Q(MomentIndx::Y2) - L(MomentIndx::Y2) - 0.5*invs_dy*Ly - 3.*invs_dy*invs_dy*D*ay2)/Er;
  flux_y2_(g, m) = flx_y2;

  const double old_flx_z2 = flux_z2_(g, m);
  const double flx_z2 = (Q(MomentIndx::Z2) - L(MomentIndx::Z2) - 0.5*invs_dz*Lz - 3.*invs_dz*invs_dz*D*az2)/Er;
  flux_z2_(g, m) = flx_z2;

  //----------------------------------------------------------------------------
  // UPDATE INCOMING CURRENTS IN NEIGHBORING NODES / B.C.
  update_Jin_from_Jout(g, m);

  //----------------------------------------------------------------------------
  // RETURN MAXIMUM RELATIVE ERROR IN FLUX
  const double rd_flx_avg = std::abs((flx_avg - old_flx_avg) / flx_avg);
  const double rd_flx_x1 = std::abs((flx_x1 - old_flx_x1) / flx_x1);
  const double rd_flx_x2 = std::abs((flx_x2 - old_flx_x2) / flx_x2);
  const double rd_flx_y1 = std::abs((flx_y1 - old_flx_y1) / flx_y1);
  const double rd_flx_y2 = std::abs((flx_y2 - old_flx_y2) / flx_y2);
  const double rd_flx_z1 = std::abs((flx_z1 - old_flx_z1) / flx_z1);
  const double rd_flx_z2 = std::abs((flx_z2 - old_flx_z2) / flx_z2);
  return var_abs_max(7, rd_flx_avg, rd_flx_x1, rd_flx_x2, rd_flx_y1, rd_flx_y2, rd_flx_z1, rd_flx_z2);
}

void NEMDiffusionDriver::inner_iteration(xt::xtensor<double, 2>& errors) {
  // Iterate through all nodes
  for (std::size_t m = 0; m < NM_; m++) {
    const auto geom_indx = geom_->geom_indx(m);
    const double dx = geom_->dx(geom_indx[0]);
    const double dy = geom_->dy(geom_indx[1]);
    const double dz = geom_->dz(geom_indx[2]);
    const double invs_dx = 1. / dx;
    const double invs_dy = 1. / dy;
    const double invs_dz = 1. / dz;
    const auto& xs = *geom_->mat(m);

    for (std::size_t g = 0; g < NG_; g++) {
      errors(g, m) = calc_node(g, m, geom_indx, invs_dx, invs_dy, invs_dz, xs);
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
  xt::xtensor<double, 2> errors = flux_avg_;

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
    
    // Perform inner iterations until convergence
    double max_flx_diff = 100.;
    for(std::size_t ii = 0; ii < 10; ii++) {
      inner_iteration(errors);
      max_flx_diff = xt::amax(errors)();
      std::exit(1);
    }

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
}

}