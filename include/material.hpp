#ifndef SCARABEE_MATERIAL_H
#define SCARABEE_MATERIAL_H

#include <cstdint>
#include <map>
#include <vector>

#include <yaml-cpp/yaml.h>

class Material {
  public:
    Material(const YAML::Node& node);

    std::uint32_t id() const { return id_; }

    std::uint32_t ngroups() const { return ngroups_; }
    
    std::uint32_t ndelayed_families() const { return ndelayed_families_; }

    float Et(std::uint32_t g) const { return Et_[g]; }
    
    float Es(std::uint32_t gin, std::uint32_t gout) const { return Es_[gin][gout]; }

    float D(std::uint32_t g) const { return D_[g]; }

    float Ea(std::uint32_t g) const { return Ea_[g]; }

    float Ef(std::uint32_t g) const { return Ef_[g]; }

    float nu_prompt(std::uint32_t g) const {
      return nu_prompt_[g];
    }

    float nu_delayed(std::uint32_t g) const {
      if (nu_delayed_.empty() == false) {
        return nu_delayed_[g];
      }

      return 0.;
    }

    float chi(std::uint32_t gin, std::uint32_t gout) const {
      return chi_[gin][gout];
    }

    float nu(std::uint32_t g) const {
      return nu_prompt(g) + nu_delayed(g);
    }

    float decay_constant(std::uint32_t df) const { return decay_constants_[df]; }

    float delayed_family_probabliity(std::uint32_t df) const { return P_delayed_group_[df]; }

    float velocity(std::uint32_t g) const { return velocity_[g]; }

  private:
    std::uint32_t id_, ngroups_, ndelayed_families_;
    std::vector<float> Et_;
    std::vector<std::vector<float>> Es_;
    std::vector<float> D_;
    std::vector<float> Ea_;
    std::vector<float> Ef_;
    std::vector<std::vector<float>> chi_;
    std::vector<float> nu_prompt_;
    std::vector<float> nu_delayed_;
    std::vector<float> P_delayed_group_;
    std::vector<float> decay_constants_;
    std::vector<float> velocity_;
};


#endif
