#include <moc/flat_source_region.hpp>
#include <utils/scarabee_exception.hpp>

#include <xtensor/xbuilder.hpp>

void FlatSourceRegion::initialize() {
  if (xs_ == nullptr) {
    throw ScarabeeException(
        "Cannot initialize FlatSourceRegion with no cross sections.");
  }

  flux_ = xt::zeros<double>({xs_->ngroups()});
  source_ = xt::zeros<double>({xs_->ngroups()});
}