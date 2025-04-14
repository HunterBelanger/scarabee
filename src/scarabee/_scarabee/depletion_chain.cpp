#include <data/depletion_chain.hpp>
#include <utils/logging.hpp>
#include <utils/nuclide_names.hpp>
#include <utils/scarabee_exception.hpp>

#include <cereal/archives/portable_binary.hpp>

#include <xtensor/views/xview.hpp>

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

  // Make sure sum isn't greater than 1 + a small tolerance
  if (ratios_sum > 1. + 1.E-5) {
    const auto mssg = "Sum of branching ratios exceeds unity.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // We don't normalize the branching ratios and don't care if the sum < 1.
  // This is due to the fact that if we remove short lived nuclides, and they
  // don't have a decay target, we then the ratios wouldn't sum to unity
  // anymore ! Normalizing after would artificially increase production of the
  // other possible targets in the branch which we don't want either.
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

void FissionYields::remove_nuclide(const std::string& nuclide) {
  auto has_nuclide = [this, &nuclide]() {
    for (const auto& target : this->targets_) {
      if (target == nuclide) return true;
    }
    return false;
  };

  while (has_nuclide()) {
    // Get the index of the nuclide
    std::size_t i = 0;
    auto it = targets_.begin();
    for (i = 0; i < targets_.size(); i++) {
      if (nuclide == targets_[i]) break;
      it++;
    }
    targets_.erase(it);  // Remove target name

    // We have the index, now we remove it from the array.
    // To do this, we first make a view
    auto view = xt::view(yields_, xt::all(), xt::drop(i));
    xt::xtensor<double, 2> temp_yields = view;
    yields_ = temp_yields;
  }
}

void FissionYields::replace_nuclide(const std::string& nuclide,
                                    const std::string& new_nuclide) {
  for (auto& target : targets_) {
    if (target == nuclide) target = new_nuclide;
  }
}

void FissionYields::replace_nuclide(const std::string& nuclide,
                                    const BranchingTargets& targets) {
  // Get the total yield of the nuclide we want to remove
  std::vector<double> nuclide_yield_sums(incident_energies_.size(), 0.);
  for (std::size_t i = 0; i < targets_.size(); i++) {
    if (nuclide == targets_[i]) {
      for (std::size_t e = 0; e < incident_energies_.size(); e++) {
        nuclide_yield_sums[e] += yields_(e, i);
      }
    }
  }

  // Remove all instances of the nuclide
  this->remove_nuclide(nuclide);

  // Make copy of yields and then zero yields
  xt::xtensor<double, 2> old_yields = yields_;
  yields_ =
      xt::zeros<double>({old_yields.shape()[0],
                         old_yields.shape()[1] + targets.branches().size()});

  // Re-fill yields
  for (std::size_t e = 0; e < yields_.shape()[0]; e++) {
    std::size_t ti = 0;
    for (std::size_t i = 0; i < yields_.shape()[1]; i++) {
      if (i < old_yields.shape()[1]) {
        // Use old data if in that region
        yields_(e, i) = old_yields(e, i);
      } else {
        // If in region of new target, get info from branch ratios
        yields_(e, i) =
            targets.branches()[ti].branch_ratio * nuclide_yield_sums[e];
        ti++;
      }
    }
  }

  // Add the new targets to target list
  for (const auto& branch : targets.branches()) {
    targets_.push_back(branch.target);
  }
}

void ChainEntry::remove_nuclide(const std::string& nuclide,
                                const Target& new_target) {
  auto target_eliminator = [this, &nuclide](const auto& t) {
    this->remove_nuclide(nuclide, t);
  };
  std::visit(target_eliminator, new_target);
}

void ChainEntry::remove_nuclide(const std::string& nuclide,
                                const NoTarget& new_target) {
  if (decay_targets_)
    remove_nuclide(nuclide, new_target, decay_targets_.value());
  if (n_gamma_) remove_nuclide(nuclide, new_target, n_gamma_.value());
  if (n_2n_) remove_nuclide(nuclide, new_target, n_2n_.value());
  if (n_3n_) remove_nuclide(nuclide, new_target, n_3n_.value());
  if (n_p_) remove_nuclide(nuclide, new_target, n_p_.value());
  if (n_alpha_) remove_nuclide(nuclide, new_target, n_alpha_.value());
  if (n_fission_) remove_nuclide(nuclide, new_target, n_fission_.value());
}

void ChainEntry::remove_nuclide(const std::string& nuclide,
                                const SingleTarget& new_target) {
  if (decay_targets_)
    remove_nuclide(nuclide, new_target, decay_targets_.value());
  if (n_gamma_) remove_nuclide(nuclide, new_target, n_gamma_.value());
  if (n_2n_) remove_nuclide(nuclide, new_target, n_2n_.value());
  if (n_3n_) remove_nuclide(nuclide, new_target, n_3n_.value());
  if (n_p_) remove_nuclide(nuclide, new_target, n_p_.value());
  if (n_alpha_) remove_nuclide(nuclide, new_target, n_alpha_.value());
  if (n_fission_) remove_nuclide(nuclide, new_target, n_fission_.value());
}

void ChainEntry::remove_nuclide(const std::string& nuclide,
                                const BranchingTargets& new_targets) {
  if (decay_targets_)
    remove_nuclide(nuclide, new_targets, decay_targets_.value());
  if (n_gamma_) remove_nuclide(nuclide, new_targets, n_gamma_.value());
  if (n_2n_) remove_nuclide(nuclide, new_targets, n_2n_.value());
  if (n_3n_) remove_nuclide(nuclide, new_targets, n_3n_.value());
  if (n_p_) remove_nuclide(nuclide, new_targets, n_p_.value());
  if (n_alpha_) remove_nuclide(nuclide, new_targets, n_alpha_.value());
  if (n_fission_) remove_nuclide(nuclide, new_targets, n_fission_.value());
}

void ChainEntry::remove_nuclide(const std::string& nuclide,
                                const NoTarget& new_target, Target& target) {
  if (std::holds_alternative<NoTarget>(target)) {
    // Transfering NoTarget to NoTarget. Nothing to do.
  } else if (std::holds_alternative<SingleTarget>(target)) {
    const auto target_name = std::get<SingleTarget>(target).target();
    if (target_name == nuclide) {
      target = new_target;
    }
  } else {
    // BranchingTargets
    auto& branch_targets = std::get<BranchingTargets>(target);
    branch_targets.remove_nuclide(nuclide);
  }
}

void ChainEntry::remove_nuclide(const std::string& nuclide,
                                const SingleTarget& new_target,
                                Target& target) {
  if (std::holds_alternative<NoTarget>(target)) {
    // Nothing to do. Target is empty so don't need to replace anything.
  } else if (std::holds_alternative<SingleTarget>(target)) {
    const auto& target_name = std::get<SingleTarget>(target).target();
    if (target_name == nuclide) {
      target = new_target;
    }
  } else {
    // BranchingTargets
    auto& branch_targets = std::get<BranchingTargets>(target);
    branch_targets.replace_nuclide(nuclide, new_target.target());
  }
}

void ChainEntry::remove_nuclide(const std::string& nuclide,
                                const BranchingTargets& new_targets,
                                Target& target) {
  if (std::holds_alternative<NoTarget>(target)) {
    // Nothing to do. Target is empty so don't need to replace anything.
  } else if (std::holds_alternative<SingleTarget>(target)) {
    const auto& target_name = std::get<SingleTarget>(target).target();
    if (target_name == nuclide) {
      target = new_targets;
    }
  } else {
    // BranchingTargets
    auto& branch_targets = std::get<BranchingTargets>(target);
    branch_targets.replace_nuclide(nuclide, new_targets);

    if (branch_targets.branches().size() == 0) {
      target = NoTarget();
    }
  }
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

std::set<std::string> DepletionChain::nuclides() const {
  std::set<std::string> nucs;
  for (const auto& nuc : data_) nucs.insert(nuc.first);
  return nucs;
}

void DepletionChain::remove_nuclide(const std::string& nuclide) {
  // If nuclide isn't in the chain, we can't remove it
  if (this->holds_nuclide_data(nuclide) == false) {
    const auto mssg = "Cannot remove nuclide \"" + nuclide +
                      "\", as not present in depletion chain.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Get the ChainEntry for the nuclide
  auto& nuc = this->data_[nuclide];

  // Does the nuclide have decay data ? If not, we can't exactly replace it.
  if (nuc.decay_targets().has_value() == false) {
    const auto mssg = "Nuclide \"" + nuclide + "\" has no decay targets.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  auto& decay_targets = nuc.decay_targets().value();

  // We need to make sure the decay target doesn't contain itself !!
  if (std::holds_alternative<SingleTarget>(decay_targets)) {
    if (std::get<SingleTarget>(decay_targets).target() == nuclide) {
      // If we only decay to ourself, set the decay to no target.
      decay_targets = NoTarget();
    }
  } else if (std::holds_alternative<BranchingTargets>(decay_targets)) {
    // If we have branches, remove ourself
    auto& dt = std::get<BranchingTargets>(decay_targets);
    dt.remove_nuclide(nuclide);

    if (dt.branches().size() == 0) {
      // We once we removed ourself there are no branches,
      // just set to no target.
      decay_targets = NoTarget();
    }
  }

  // Go through all other chain entries (but not ourself!) and replace
  // instances of this nuclide in targets with its decay targets.
  for (auto& entry : data_) {
    if (&entry.second == &nuc) continue;  // Don't do ourself !

    entry.second.remove_nuclide(nuclide, decay_targets);
  }

  // Delete nuclide from the chain
  data_.erase(nuclide);
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

      // Get the chain entry for the nuclides. If it doesn't have an entry,
      // go to the next nuclide.
      if (this->holds_nuclide_data(nuc_name) == false) continue;
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