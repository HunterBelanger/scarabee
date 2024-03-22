#include <moc/flat_source_region.hpp>
#include <utils/scarabee_exception.hpp>

void FlatSourceRegion::initialize() {
  if (xs_ == nullptr) {
    throw ScarabeeException("Cannot initialize FlatSourceRegion with no cross sections.");
  }

  flux_.resize(xs_->ngroups(), 0.);
  source_.resize(xs_->ngroups(), 0.);
}