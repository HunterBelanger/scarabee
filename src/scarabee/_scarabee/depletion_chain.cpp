#include <data/depletion_chain.hpp>
#include <utils/logging.hpp>
#include <utils/nuclide_names.hpp>
#include <utils/scarabee_exception.hpp>

#include <cereal/archives/portable_binary.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>

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
    iE = incident_energies_.size() - 2;
    fE = 1.;
  } else if (E > incident_energies_.front()) {
    for (iE = 0; iE < incident_energies_.size() - 1; iE++) {
      const double El = incident_energies_[iE];
      const double Eh = incident_energies_[iE + 1];

      if (El <= E && E <= Eh) {
        fE = (E - El) / (Eh - El);
        break;
      }
    }
  }

  return yields_(iE, t) * (1. - fE) + yields_(iE + 1, t) * fE;
}

bool DepletionChain::holds_nuclide_data(const std::string& nuclide) const {
  if (data_.find(nuclide) == data_.end()) {
    return false;
  }

  return true;
}

const ChainEntry& DepletionChain::nuclide_data(
    const std::string& nuclide) const {
  const auto it = data_.find(nuclide);
  if (it == data_.end()) {
    const auto mssg =
        "Nuclide \"" + nuclide + "\" is not present in the depletion chain.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  return it->second;
}

void DepletionChain::insert_entry(const std::string& nuclide,
                                  const ChainEntry& entry) {
  if (this->holds_nuclide_data(nuclide)) {
    const auto mssg = "Nuclide \"" + nuclide +
                      "\" is already present in the depletion chain.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  data_[nuclide] = entry;
}

// Helper functions for extracting new targets from various types of targets
void extract_targets(const NoTarget& /*t*/,
                     std::set<std::string>& /*found_nuclides*/,
                     std::set<std::string>& /*next_nuclides*/) {}

void extract_targets(const SingleTarget& t,
                     std::set<std::string>& found_nuclides,
                     std::set<std::string>& next_nuclides) {
  if (found_nuclides.contains(t.target()))
    return;
  else if (next_nuclides.contains(t.target()))
    return;
  next_nuclides.insert(t.target());
}

void extract_targets(const BranchingTargets& bt,
                     std::set<std::string>& found_nuclides,
                     std::set<std::string>& next_nuclides) {
  for (const auto& branch : bt.branches()) {
    if (found_nuclides.contains(branch.target))
      continue;
    else if (next_nuclides.contains(branch.target))
      continue;
    next_nuclides.insert(branch.target);
  }
}

void extract_targets(const FissionYields& fy,
                     std::set<std::string>& found_nuclides,
                     std::set<std::string>& next_nuclides) {
  for (const auto& target : fy.targets()) {
    if (found_nuclides.contains(target))
      continue;
    else if (next_nuclides.contains(target))
      continue;
    next_nuclides.insert(target);
  }
}

std::vector<std::string> DepletionChain::descend_chains(
    std::set<std::string> nuclides, bool decay_only) const {
  std::set<std::string> found_nuclides = nuclides;
  std::set<std::string> next_nuclides;

  // While we have a new set of nuclides to follow
  while (nuclides.empty() == false) {
    // Check all products of each nuclide
    for (const auto& nuc_name : nuclides) {
      // Save this nuclide as having been "found"
      found_nuclides.insert(nuc_name);

      // Get the chain entry for the nuclides
      if (this->holds_nuclide_data(nuc_name) == false) {
        const auto mssg =
            "The nuclide \"" + nuc_name + "\" has no depletion chain entry.";
        spdlog::error(mssg);
        throw ScarabeeException(mssg);
      }
      const auto& entry = this->nuclide_data(nuc_name);

      // For each decay or transmutation reaction (including fission), save all
      // the newly encountered target nuclides.
      auto target_extractor = [&found_nuclides, &next_nuclides](const auto& t) {
        extract_targets(t, found_nuclides, next_nuclides);
      };

      if (entry.decay_targets())
        std::visit(target_extractor, *entry.decay_targets());
      if (decay_only) continue;

      if (entry.n_gamma()) std::visit(target_extractor, *entry.n_gamma());
      if (entry.n_2n()) std::visit(target_extractor, *entry.n_2n());
      if (entry.n_3n()) std::visit(target_extractor, *entry.n_3n());
      if (entry.n_p()) std::visit(target_extractor, *entry.n_p());
      if (entry.n_alpha()) std::visit(target_extractor, *entry.n_alpha());
      if (entry.n_fission())
        extract_targets(*entry.n_fission(), found_nuclides, next_nuclides);
    }

    nuclides.swap(next_nuclides);
    next_nuclides.clear();
  }

  // Make a vector with the found nuclides
  std::vector<std::string> found_vec(found_nuclides.begin(),
                                     found_nuclides.end());

  // Sort nuclides by Z then A
  std::sort(found_vec.begin(), found_vec.end(),
            [](const std::string& n1, const std::string& n2) {
              return nuclide_name_to_za(n1) < nuclide_name_to_za(n2);
            });

  return found_vec;
}

void DepletionChain::save(const std::string& fname) const {
  if (std::filesystem::exists(fname)) {
    std::filesystem::remove(fname);
  }

  std::ofstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryOutputArchive arc(file);

  arc(*this);
}

std::shared_ptr<DepletionChain> DepletionChain::load(const std::string& fname) {
  if (std::filesystem::exists(fname) == false) {
    std::stringstream mssg;
    mssg << "The file \"" << fname << "\" does not exist.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  std::shared_ptr<DepletionChain> out = std::make_shared<DepletionChain>();

  std::ifstream file(fname, std::ios_base::binary);

  cereal::PortableBinaryInputArchive arc(file);

  arc(*out);

  return out;
}

}  // namespace scarabee