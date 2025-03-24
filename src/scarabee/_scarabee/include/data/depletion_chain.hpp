#ifndef SCARABEE_DEPLETION_CHAIN_H
#define SCARABEE_DEPLETION_CHAIN_H

#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <xtensor/xtensor.hpp>

namespace scarabee {

class NoTarget {};

class SimpleTarget {
 public:
  SimpleTarget(const std::string& target) : target_(target) {}

  const std::string& target() const { return target_; }

 private:
  std::string target_;
};

class BranchingTargets {
 public:
  struct Branch {
    std::string target;
    double branch_ratio;
  };

  const std::vector<Branch>& branches() const { return branches_; }

 private:
  std::vector<Branch> branches_;
};

using Target = std::variant<NoTarget, SimpleTarget, BranchingTargets>;

class FissionYields {
 public:

  const std::vector<std::string>& targets() const { return targets_; }
  const std::vector<double>& incident_energies() const { return incident_energies_; }

  double branch_ratio(std::size_t t, double E) const;

 private:
  std::vector<std::string> targets_;
  std::vector<double> incident_energies_;
  // First index incident energy, second is target
  xt::xtensor<double, 2> branching_ratios_;
};

class ChainEntry {
 public:
  double half_life() const { return half_life_; }
  const std::optional<BranchingTargets>& decay_targets() const {
    return decay_targets_;
  }

  const std::optional<Target>& n_gamma() const { return n_gamma_; }
  const std::optional<Target>& n_2n() const { return n_2n_; }
  const std::optional<Target>& n_3n() const { return n_3n_; }
  const std::optional<Target>& n_p() const { return n_p_; }
  const std::optional<Target>& n_alpha() const { return n_alpha_; }
  const std::optional<FissionYields>& n_fission() const { return n_fission_; }

 private:
  // Radioactive Decay
  double half_life_;
  std::optional<BranchingTargets> decay_targets_;

  // Transmutation Reactions
  std::optional<Target> n_gamma_;
  std::optional<Target> n_2n_;
  std::optional<Target> n_3n_;
  std::optional<Target> n_p_;
  std::optional<Target> n_alpha_;
  std::optional<FissionYields> n_fission_;
};

class DepletionChain {
 public:
  bool holds_nuclide_data(const std::string& nuclide) const;
  const ChainEntry& nuclide_data(const std::string& nuclide) const;

  void gather_all_targets(const std::string& base,
                          std::set<std::string> targets) const;
  void recursive_gather_all_targets(const std::string& base,
                                    std::set<std::string> targets) const;

 private:
  std::map<std::string, ChainEntry> data_;
};

}  // namespace scarabee

#endif