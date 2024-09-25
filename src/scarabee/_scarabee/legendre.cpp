#include <moc/quadrature/legendre.hpp>

#include <cmath>

namespace scarabee {

// N = 2
template <>
const std::array<double, 1> Legendre<2>::sin_ = {
    std::sqrt(1. - std::pow(5.77350269189625764509148780501957456e-01, 2.))};
template <>
const std::array<double, 1> Legendre<2>::invs_sin_ = {1. / sin_[0]};
template <>
const std::array<double, 1> Legendre<2>::wgt_ = {
    1.00000000000000000000000000000000000e+00};
template <>
const std::array<double, 1> Legendre<2>::wsin_ = {wgt_[0] * sin_[0]};
template <>
const std::array<double, 1> Legendre<2>::polar_angle_ = {std::asin(sin_[0])};

// N = 4
template <>
const std::array<double, 2> Legendre<4>::sin_ = {
    std::sqrt(1. - std::pow(3.39981043584856264802665759103244687e-01, 2.)),
    std::sqrt(1. - std::pow(8.61136311594052575223946488892809505e-01, 2.))};
template <>
const std::array<double, 2> Legendre<4>::invs_sin_ = {1. / sin_[0],
                                                      1. / sin_[1]};
template <>
const std::array<double, 2> Legendre<4>::wgt_ = {
    6.52145154862546142626936050778000593e-01,
    3.47854845137453857373063949221999407e-01};
template <>
const std::array<double, 2> Legendre<4>::wsin_ = {wgt_[0] * sin_[0],
                                                  wgt_[1] * sin_[1]};
template <>
const std::array<double, 2> Legendre<4>::polar_angle_ = {std::asin(sin_[0]), std::asin(sin_[1])};

// N = 6
template <>
const std::array<double, 3> Legendre<6>::sin_ = {
    std::sqrt(1. - std::pow(2.38619186083196908630501721680711935e-01, 2.)),
    std::sqrt(1. - std::pow(6.61209386466264513661399595019905347e-01, 2.)),
    std::sqrt(1. - std::pow(9.32469514203152027812301554493994609e-01, 2.))};
template <>
const std::array<double, 3> Legendre<6>::invs_sin_ = {
    1. / sin_[0], 1. / sin_[1], 1. / sin_[2]};
template <>
const std::array<double, 3> Legendre<6>::wgt_ = {
    4.67913934572691047389870343989550995e-01,
    3.60761573048138607569833513837716112e-01,
    1.71324492379170345040296142172732894e-01};
template <>
const std::array<double, 3> Legendre<6>::wsin_ = {
    wgt_[0] * sin_[0], wgt_[1] * sin_[1], wgt_[2] * sin_[2]};
template <>
const std::array<double, 3> Legendre<6>::polar_angle_ = {std::asin(sin_[0]), std::asin(sin_[1]), std::asin(sin_[2])};

// N = 8
template <>
const std::array<double, 4> Legendre<8>::sin_ = {
    std::sqrt(1. - std::pow(1.83434642495649804939476142360183981e-01, 2.)),
    std::sqrt(1. - std::pow(5.25532409916328985817739049189246349e-01, 2.)),
    std::sqrt(1. - std::pow(7.96666477413626739591553936475830437e-01, 2.)),
    std::sqrt(1. - std::pow(9.60289856497536231683560868569472990e-01, 2.))};
template <>
const std::array<double, 4> Legendre<8>::invs_sin_ = {
    1. / sin_[0], 1. / sin_[1], 1. / sin_[2], 1. / sin_[3]};
template <>
const std::array<double, 4> Legendre<8>::wgt_ = {
    3.62683783378361982965150449277195612e-01,
    3.13706645877887287337962201986601313e-01,
    2.22381034453374470544355994426240884e-01,
    1.01228536290376259152531354309962190e-01};
template <>
const std::array<double, 4> Legendre<8>::wsin_ = {
    wgt_[0] * sin_[0], wgt_[1] * sin_[1], wgt_[2] * sin_[2], wgt_[3] * sin_[3]};
template <>
const std::array<double, 4> Legendre<8>::polar_angle_ = {std::asin(sin_[0]), std::asin(sin_[1]), std::asin(sin_[2]), std::asin(sin_[3])};

// N = 10
template <>
const std::array<double, 5> Legendre<10>::sin_ = {
    std::sqrt(1. - std::pow(1.48874338981631210884826001129719985e-01, 2.)),
    std::sqrt(1. - std::pow(4.33395394129247190799265943165784162e-01, 2.)),
    std::sqrt(1. - std::pow(6.79409568299024406234327365114873576e-01, 2.)),
    std::sqrt(1. - std::pow(8.65063366688984510732096688423493049e-01, 2.)),
    std::sqrt(1. - std::pow(9.73906528517171720077964012084452053e-01, 2.))};
template <>
const std::array<double, 5> Legendre<10>::invs_sin_ = {
    1. / sin_[0], 1. / sin_[1], 1. / sin_[2], 1. / sin_[3], 1. / sin_[4]};
template <>
const std::array<double, 5> Legendre<10>::wgt_ = {
    2.95524224714752870173892994651338329e-01,
    2.69266719309996355091226921569469353e-01,
    2.19086362515982043995534934228163192e-01,
    1.49451349150580593145776339657697332e-01,
    6.66713443086881375935688098933317929e-02};
template <>
const std::array<double, 5> Legendre<10>::wsin_ = {
    wgt_[0] * sin_[0], wgt_[1] * sin_[1], wgt_[2] * sin_[2], wgt_[3] * sin_[3],
    wgt_[4] * sin_[4]};
template <>
const std::array<double, 5> Legendre<10>::polar_angle_ = {std::asin(sin_[0]), std::asin(sin_[1]), std::asin(sin_[2]), std::asin(sin_[3]), std::asin(sin_[4])};

// N = 12
template <>
const std::array<double, 6> Legendre<12>::sin_ = {
    std::sqrt(1. - std::pow(1.25233408511468915472441369463853130e-01, 2.)),
    std::sqrt(1. - std::pow(3.67831498998180193752691536643717561e-01, 2.)),
    std::sqrt(1. - std::pow(5.87317954286617447296702418940534280e-01, 2.)),
    std::sqrt(1. - std::pow(7.69902674194304687036893833212818076e-01, 2.)),
    std::sqrt(1. - std::pow(9.04117256370474856678465866119096193e-01, 2.)),
    std::sqrt(1. - std::pow(9.81560634246719250690549090149280823e-01, 2.))};
template <>
const std::array<double, 6> Legendre<12>::invs_sin_ = {
    1. / sin_[0], 1. / sin_[1], 1. / sin_[2],
    1. / sin_[3], 1. / sin_[4], 1. / sin_[5]};
template <>
const std::array<double, 6> Legendre<12>::wgt_ = {
    2.49147045813402785000562436042951211e-01,
    2.33492536538354808760849898924878056e-01,
    2.03167426723065921749064455809798377e-01,
    1.60078328543346226334652529543359072e-01,
    1.06939325995318430960254718193996224e-01,
    4.71753363865118271946159614850170603e-02};
template <>
const std::array<double, 6> Legendre<12>::wsin_ = {
    wgt_[0] * sin_[0], wgt_[1] * sin_[1], wgt_[2] * sin_[2],
    wgt_[3] * sin_[3], wgt_[4] * sin_[4], wgt_[5] * sin_[5]};
template <>
const std::array<double, 6> Legendre<12>::polar_angle_ = {std::asin(sin_[0]), std::asin(sin_[1]), std::asin(sin_[2]), std::asin(sin_[3]), std::asin(sin_[4]), std::asin(sin_[5])};

}  // namespace scarabee
