#include <utils/bickley.hpp>
#include <utils/chebyshev.hpp>
#include <utils/constants.hpp>
#include <utils/gauss_kronrod.hpp>
#include <utils/scarabee_exception.hpp>

#include <array>
#include <cmath>
#include <sstream>

double Ki3(double x) {
  if (x < 1.) {
    constexpr double a{0x0p+0};  // a = 0.000000
    constexpr double b{0x1p+0};  // b = 1.000000
    constexpr std::array<double, 25> c{
        0x1.df97acf54c7a8p-1,   -0x1.12e17623b053bp-2,  0x1.5d0f8db96959ep-5,
        -0x1.56fc22a584252p-8,  0x1.4d302f53da333p-11,  -0x1.9807d3b547ae1p-14,
        0x1.64748942c70a4p-16,  -0x1.a482183151eb8p-18, 0x1.2ff395d3051ecp-19,
        -0x1.f977af7c66667p-21, 0x1.d22129487ae15p-22,  -0x1.d1df7559eb852p-23,
        0x1.f0b80eap-24,        -0x1.17447955c28f6p-24, 0x1.484d22acccccdp-25,
        -0x1.90b8f81d70a3ep-26, 0x1.f8f1aa5eb851fp-27,  -0x1.46c02b5c28f5cp-27,
        0x1.b021f3999999ap-28,  -0x1.224f30eb851ecp-28, 0x1.89318023d70a4p-29,
        -0x1.0916f46b851ecp-29, 0x1.5bb61bp-30,         -0x1.a4f449851eb85p-31,
        0x1.8d23ae3d70a3ep-32};

    return chebyshev_eval(x, a, b, c);
  } else if (x < 3.) {
    constexpr double a{0x1p+0};    // a = 1.000000
    constexpr double b{0x1.8p+1};  // b = 3.000000
    constexpr std::array<double, 25> c{
        0x1.a825725f0d20dp-3,   -0x1.9c9d68433ec12p-4,
        0x1.c08f993e2357ep-6,   -0x1.586950cdd77a4p-8,
        0x1.a514233905785p-11,  -0x1.bdd377225cb85p-14,
        0x1.b7dfcbaf6c29p-17,   -0x1.b136cbe1b5c29p-20,
        0x1.c27ed51eb851fp-23,  -0x1.fdbdabf851eb9p-26,
        0x1.3a55cf68f5c29p-28,  -0x1.a0b56a147ae15p-31,
        0x1.2429be147ae15p-33,  -0x1.abeaf33333333p-36,
        0x1.4471c28f5c28fp-38,  -0x1.fa2bd70a3d70ap-41,
        0x1.942a3d70a3d71p-43,  -0x1.47b3333333333p-45,
        0x1.121eb851eb852p-47,  -0x1.919999999999ap-50,
        0x1.43d70a3d70a3ep-52,  -0x1.cp-53,
        -0x1.599999999999ap-54, -0x1.e51eb851eb852p-53,
        0x1.a8a3d70a3d70ap-53};

    return chebyshev_eval(x, a, b, c);
  } else if (x < 5.) {
    constexpr double a{0x1.8p+1};  // a = 3.000000
    constexpr double b{0x1.4p+2};  // b = 5.000000
    constexpr std::array<double, 25> c{
        0x1.76cd1a526853ep-6,   -0x1.62810320acc1fp-7,  0x1.700d794405fd5p-9,
        -0x1.08020c91ab4cap-11, 0x1.23b7f93ab059ap-14,  -0x1.091ddb3d4870ap-17,
        0x1.9f503628c51ecp-21,  -0x1.23413ee766667p-24, 0x1.7aa1c260f5c29p-28,
        -0x1.d7b0308p-32,       0x1.2295df5c28f5cp-35,  -0x1.6bb5a147ae148p-39,
        0x1.d6db333333333p-43,  -0x1.3d9999999999ap-46, 0x1.bf9999999999ap-50,
        -0x1.447ae147ae148p-53, 0x1.0a3d70a3d70a4p-58,  0x1.ab851eb851eb8p-56,
        0x1.999999999999ap-60,  0x1.851eb851eb852p-56,  -0x1.a3d70a3d70a3ep-58,
        -0x1.dc28f5c28f5c3p-57, -0x1.8a3d70a3d70a4p-57, -0x1.770a3d70a3d71p-56,
        0x1.68f5c28f5c28fp-56};

    return chebyshev_eval(x, a, b, c);
  } else if (x < 7.) {
    constexpr double a{0x1.4p+2};  // a = 5.000000
    constexpr double b{0x1.cp+2};  // b = 7.000000
    constexpr std::array<double, 25> c{
        0x1.60f36b1c82ea4p-9,   -0x1.494fcce6897afp-10, 0x1.4f0b72616720ap-12,
        -0x1.d3e57c3d0e9b3p-15, 0x1.f2eea20abac52p-18,  -0x1.b05d2c885a51fp-21,
        0x1.3d7e56fb8p-24,      -0x1.979812661eb85p-28, 0x1.d518ddec28f5cp-32,
        -0x1.eedefd8f5c29p-36,  0x1.e8a95c28f5c29p-40,  -0x1.cd0a666666667p-44,
        0x1.a81d70a3d70a4p-48,  -0x1.83b851eb851ecp-52, 0x1.723d70a3d70a4p-56,
        -0x1.2e147ae147ae1p-60, -0x1.51eb851eb851fp-60, 0x1.9eb851eb851ecp-59,
        0x1.1eb851eb851ecp-62,  0x1.63d70a3d70a3ep-59,  -0x1.947ae147ae148p-61,
        -0x1.947ae147ae148p-60, -0x1.970a3d70a3d71p-60, -0x1.370a3d70a3d71p-59,
        0x1.3051eb851eb85p-59};

    return chebyshev_eval(x, a, b, c);
  } else if (x < 9.) {
    constexpr double a{0x1.cp+2};  // a = 7.000000
    constexpr double b{0x1.2p+3};  // b = 9.000000
    constexpr std::array<double, 25> c{
        0x1.572a4e0ec6729p-12,  -0x1.3d8ac8e362ac5p-13, 0x1.3f51b61c69bc2p-15,
        -0x1.b76a6af9ab338p-18, 0x1.cc0910b57aae1p-21,  -0x1.8593de9bdb99ap-24,
        0x1.15d69395a5c29p-27,  -0x1.578f0743d70a4p-31, 0x1.789dd6ecccccdp-35,
        -0x1.74c91ed70a3d7p-39, 0x1.528d58f5c28f6p-43,  -0x1.1e47851eb851fp-47,
        0x1.c950a3d70a3d7p-52,  -0x1.5e147ae147ae1p-56, 0x1.18f5c28f5c28fp-60,
        0x1.1eb851eb851ecp-66,  -0x1.a3d70a3d70a3ep-63, 0x1.999999999999ap-62,
        0x1.eb851eb851eb8p-66,  0x1.599999999999ap-62,  -0x1.8a3d70a3d70a4p-64,
        -0x1.7851eb851eb85p-63, -0x1.b5c28f5c28f5cp-63, -0x1.28f5c28f5c28fp-62,
        0x1.4e147ae147ae1p-62};

    return chebyshev_eval(x, a, b, c);
  } else if (x < 11.) {
    constexpr double a{0x1.2p+3};  // a = 9.000000
    constexpr double b{0x1.6p+3};  // b = 11.000000
    constexpr std::array<double, 25> c{
        0x1.5448e0162d652p-15,  -0x1.391de02664243p-16,
        0x1.38776b8f10e54p-18,  -0x1.a9fc6267118fbp-21,
        0x1.b8fdea0ac7ca4p-24,  -0x1.706d3fe32dd71p-27,
        0x1.0277b124e5c29p-30,  -0x1.3939d4ffd70a4p-34,
        0x1.4ee53b6e147aep-38,  -0x1.414a6fc28f5c3p-42,
        0x1.18801f5c28f5cp-46,  -0x1.c32fd70a3d70ap-51,
        0x1.5211eb851eb85p-55,  -0x1.dcp-60,
        0x1.55c28f5c28f5cp-64,  0x1.999999999999ap-69,
        -0x1.6147ae147ae15p-66, 0x1.70a3d70a3d70ap-65,
        0x1.851eb851eb852p-68,  0x1.4f5c28f5c28f6p-65,
        -0x1.6b851eb851eb8p-67, -0x1.a666666666667p-66,
        -0x1.c28f5c28f5c29p-66, -0x1.f0a3d70a3d70ap-66,
        0x1.1333333333333p-65};

    return chebyshev_eval(x, a, b, c);
  } else if (x < 13.) {
    constexpr double a{0x1.6p+3};  // a = 11.000000
    constexpr double b{0x1.ap+3};  // b = 13.000000
    constexpr std::array<double, 25> c{
        0x1.5600391eea6ap-18,   -0x1.396b7beb50297p-19, 0x1.371354985c579p-21,
        -0x1.a55a529f158dcp-24, 0x1.b0e58ef75a429p-27,  -0x1.66755f5081333p-30,
        0x1.f1b31f60d999ap-34,  -0x1.29d68fd628f5cp-37, 0x1.39bb4031eb852p-41,
        -0x1.279e1e6666667p-45, 0x1.f8fadae147ae1p-50,  -0x1.8b56f5c28f5c3p-54,
        0x1.1e8a3d70a3d71p-58,  -0x1.83ae147ae147bp-63, 0x1.1c28f5c28f5c3p-67,
        0x1.47ae147ae147bp-73,  -0x1.428f5c28f5c29p-69, 0x1.6147ae147ae15p-68,
        0x1.0a3d70a3d70a4p-71,  0x1.6666666666667p-68,  -0x1.851eb851eb852p-70,
        -0x1.851eb851eb852p-69, -0x1.170a3d70a3d71p-68, -0x1.ce147ae147ae1p-69,
        0x1.20f5c28f5c28fp-68};

    return chebyshev_eval(x, a, b, c);
  } else if (x < 15.) {
    constexpr double a{0x1.ap+3};  // a = 13.000000
    constexpr double b{0x1.ep+3};  // b = 15.000000
    constexpr std::array<double, 25> c{
        0x1.5b1e66cc700f6p-21,  -0x1.3d212941edbc3p-22, 0x1.397bb0dd51b33p-24,
        -0x1.a699754bca68p-27,  0x1.afce97fb87733p-30,  -0x1.635045877799ap-33,
        0x1.e9c45a973eb85p-37,  -0x1.22a2dacb851ecp-40, 0x1.2f27b38a3d70ap-44,
        -0x1.1a602e3333333p-48, 0x1.dbd0a7ae147aep-53,  -0x1.6e7b5c28f5c29p-57,
        0x1.0475c28f5c28fp-61,  -0x1.57d70a3d70a3ep-66, 0x1.eb851eb851eb8p-71,
        0x1.47ae147ae147bp-78,  -0x1.28f5c28f5c28fp-72, 0x1.6b851eb851eb8p-71,
        0x1.1eb851eb851ecp-74,  0x1.63d70a3d70a3ep-71,  -0x1.c7ae147ae147bp-73,
        -0x1.68f5c28f5c28fp-72, -0x1.07ae147ae147bp-71, -0x1.f5c28f5c28f5cp-72,
        0x1.3333333333333p-71};

    return chebyshev_eval(x, a, b, c);
  } else {
    return 0.;
  }
}

double Ki3_quad(double x) {
  auto integrand = [x](double theta) {
    const double cos_theta = std::cos(theta);
    const double cos_theta_sqrd = cos_theta * cos_theta;

    return cos_theta_sqrd * std::exp(-x / cos_theta);
  };

  GaussKronrodQuadrature<61> gk;
  auto integral = gk.integrate(integrand, 0., PI_2, 1.E-9, 100000);

  return integral.first;
}
