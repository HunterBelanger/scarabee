#include <data/depletion_chain.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <algorithm>

namespace scarabee {

BranchingTargets::BranchingTargets(const std::vector<Branch>& branches)
    : branches_(branches) {
  if (branches_.size() < 2) {
    const auto mssg =
        "Must have at least 2 targets for a BranchingTargets instance.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Get sum of branch ratios
  double ratios_sum = 0.;
  for (const auto& branch : branches_) {
    if (branch.branch_ratio <= 0.) {
      const auto mssg =
          "Branching ratio for target \"" + branch.target + "\" is <= 0.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }

    ratios_sum += branch.branch_ratio;
  }

  for (auto& branch : branches_) {
    branch.branch_ratio /= ratios_sum;
  }
}

FissionYields::FissionYields(const std::vector<std::string>& targets,
                             const std::vector<double>& incident_energies,
                             const xt::xtensor<double, 2>& yields)
    : targets_(targets),
      incident_energies_(incident_energies),
      yields_(yields) {
  // Make sure we have targets and energies
  if (targets_.size() == 0) {
    const auto mssg = "No fission yield targets provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (incident_energies_.size() == 0) {
    const auto mssg = "No incident energies provided.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure energies are sorted
  if (std::is_sorted(incident_energies_.begin(), incident_energies_.end()) ==
      false) {
    const auto mssg = "Incident energies must be sorted.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (incident_energies_.front() <= 0.) {
    const auto mssg = "Incident energies must be > 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure shape of yields is correct
  if (yields_.shape()[0] != incident_energies_.size()) {
    const auto mssg =
        "First dimension of yields must have same length as number of incident "
        "energies.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (yields_.shape()[1] != targets_.size()) {
    const auto mssg =
        "Second dimension of yields must have same length as number of "
        "targets.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Make sure all yields are >= 0
  const double yield_min = xt::amin(yields_)();
  if (yield_min < 0.) {
    const auto mssg = "Fission yields must be >= 0.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }
}

double FissionYields::yield(std::size_t t, double E) const {
  if (t >= this->size()) {
    const auto mssg = "Target index is out of range.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (incident_energies_.size() == 1) {
    return yields_(0, t);
  }

  std::size_t iE = 0;
  double fE = 0.;

  if (E >= incident_energies_.back()) {
    iE = incident_energies_.size()-2;
    fE = 1.;
  } else if (E > incident_energies_.front()) {
    for (iE = 0; iE < incident_energies_.size()-1; iE++) {
      const double El = incident_energies_[iE];
      const double Eh = incident_energies_[iE+1];

      if (El <= E && E <= Eh) {
        fE = (E - El) / (Eh - El);
        break;
      }
    }
  }

  return yields_(iE, t)*(1. - fE) + yields_(iE+1, t)*fE;
}

}  // namespace scarabee