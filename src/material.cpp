#include <material.hpp>

#include <exception>
#include <vector>

Material::Material(const YAML::Node& node) {
  // Get material id
  if (node["id"] && node["id"].IsScalar()) {
    id_ = node["id"].as<std::uint32_t>();
  } else {
    throw std::runtime_error("Material is missing a valid id.");
  }

  //===========================================================================
  // Get Total XS
  if (!node["total"] || !node["total"].IsSequence()) {
    std::stringstream mssg;
    mssg << "Invalid total xs entry in material " << id_ << ".";
    throw std::runtime_error(mssg.str());
  }
  Et_ = node["total"].as<std::vector<float>>();
  ngroups_ = static_cast<std::uint32_t>(Et_.size());

  //===========================================================================
  // Get Diffusion Coefficient
  if (!node["diffusion"] || !node["diffusion"].IsSequence()) {
    std::stringstream mssg;
    mssg << "Invalid diffusion coefficient entry in material " << id_ << ".";
    throw std::runtime_error(mssg.str());
  }
  D_ = node["diffuison"].as<std::vector<float>>();

  //===========================================================================
  // Get Absorption XS
  std::vector<double> Ea;
  if (!node["absorption"] || !node["absorption"].IsSequence() ||
      node["absorption"].size() != ngroups_) {
    std::stringstream mssg;
    mssg << "Invalid absorption xs entry in material " << id_ << ".";
    throw std::runtime_error(mssg.str());
  }
  Ea_ = node["absorption"].as<std::vector<float>>();

  //===========================================================================
  // Get the scatter matrix
  Es_ = std::vector<std::vector<float>>(ngroups_, std::vector<float>(ngroups_, 0.));
  if (!node["scatter"] || !node["scatter"].IsSequence() ||
      node["scatter"].size() != ngroups_) {
    std::stringstream mssg;
    mssg << "Invalid scatter matrix entry in material " << id_ << ".";
    throw std::runtime_error(mssg.str());
  }

  // Go through all rows (incoming energy groups)
  for (std::size_t ei = 0; ei < ngroups_; ei++) {
    // Make sure the row is the right size
    if (!node["scatter"][ei].IsSequence() ||
        node["scatter"][ei].size() != ngroups_) {
      std::stringstream mssg;
      mssg << "Row " << ei << " of the scatter matrix for material " << id_;
      mssg << " is invalid.";
      throw std::runtime_error(mssg.str());
    }

    // Go through all columns (outgoing energy groups)
    for (std::size_t eo = 0; eo < ngroups_; eo++) {
      Es_[ei][eo] = node["scatter"][ei][eo].as<float>();

      if (Es_[ei][eo] < 0.) {
        std::stringstream mssg;
        mssg << "Negative scattering component at row " << ei;
        mssg << ", column " << eo << " of the scattering matrix";
        mssg << " for material " << id_ << ".";
        throw std::runtime_error(mssg.str());
      }
    }
  } 

  //===========================================================================
  // Get Fission XS if Present
  Ef_ = std::vector<float>(ngroups_, 0.);
  bool fissile = false;
  if (node["fission"]) {
    if (!node["fission"].IsSequence() || node["fission"].size() != ngroups_) {
      std::stringstream mssg;
      mssg << "Invalid fission xs entry in material " << id_ << ".";
      throw std::runtime_error(mssg.str());
    }
    Ef_ = node["fission"].as<std::vector<float>>();

    for (const auto& xs_f : Ef_) {
      if (xs_f < 0.) {
        std::stringstream mssg;
        mssg << "Negative fission xs in material " << id_ << ".";
        throw std::runtime_error(mssg.str());
      }

      if (xs_f > 0.) fissile = true;
    }
  }

  //===========================================================================
  // Get Nu if we have fission
  std::vector<float> nu_prmpt(ngroups_, 0.);
  std::vector<float> nu_dlyd(ngroups_, 0.);
  if (fissile) {
    if (node["nu"]) {
      // Treat nu as nu_prmpt
      if (!node["nu"].IsSequence() || node["nu"].size() != ngroups_) {
        std::stringstream mssg;
        mssg << "Invalid nu entry in material " << id_ << ".";
        throw std::runtime_error(mssg.str());
      }
      nu_prmpt = node["nu"].as<std::vector<float>>();
    } else if (node["nu_prompt"] && node["nu_delayed"]) {
      // Prompt
      if (!node["nu_prompt"].IsSequence() ||
          node["nu_prompt"].size() != ngroups_) {
        std::stringstream mssg;
        mssg << "Invalid nu_prompt entry in material " << id_ << ".";
        throw std::runtime_error(mssg.str());
      }
      nu_prmpt = node["nu_prompt"].as<std::vector<float>>();

      // Delayed
      if (!node["nu_delayed"].IsSequence() ||
          node["nu_delayed"].size() != ngroups_) {
        std::stringstream mssg;
        mssg << "Invalid nu_delayed entry in material " << id_ << ".";
        throw std::runtime_error(mssg.str());
      }
      nu_dlyd = node["nu_delayed"].as<std::vector<float>>();
    } else {
      if (node["nu_prompt"]) {
        // No nu_delayed data is given. This is bad
        std::stringstream mssg;
        mssg << "No nu_delayed data is provided in material " << id_ << ".";
        throw std::runtime_error(mssg.str());
      } else if (node["nu_delayed"]) {
        // No nu_prompt data is given. This is bad
        std::stringstream mssg;
        mssg << "No nu_prompt data is provided in material " << id_ << ".";
        std::runtime_error(mssg.str());
      } else {
        // No nu data is given at all. This is really bad
        std::stringstream mssg;
        mssg << "No nu data is provided in material " << id_ << ".";
        throw std::runtime_error(mssg.str());
      }
    }
  }

  //===========================================================================
  // Get Chi Matrix if we have fission
  chi_ = std::vector<std::vector<float>>(ngroups_, std::vector<float>(ngroups_, 0.));
  if (fissile) {
    if (!node["chi"] || !node["chi"].IsSequence()) {
      std::stringstream mssg;
      mssg << "Invalid chi matrix entry in material " << id_ << ".";
      throw std::runtime_error(mssg.str());
    }

    if (node["chi"].size() != 1 && node["chi"].size() != ngroups_) {
      std::stringstream mssg;
      mssg << "Invalid chi matrix entry in material " << id_ << ".";
      throw std::runtime_error(mssg.str());
    }

    const bool chi_matrix = (node["chi"].size() == ngroups_);
    
    if (chi_matrix) {
      // Make sure all rows have the right length
      for (std::size_t ei = 0; ei < ngroups_; ei++) {
        if (!node["chi"][ei].IsSequence() ||
            node["chi"][ei].size() != ngroups_) {
          std::stringstream mssg;
          mssg << "Invalid length for chi[" << ei << "] in material " << id_
               << ".";
          throw std::runtime_error(mssg.str());
        }

        for (std::size_t eo = 0; eo < ngroups_; eo++) {
          chi_[ei][eo] = node["chi"][ei][eo].as<float>();

          if (chi_[ei][eo] < 0.) {
            std::stringstream mssg;
            mssg << "chi[" << ei << "][" << eo << "] is negative in material "
                 << id_ << ".";
            throw std::runtime_error(mssg.str());
          }
        }
      }
    } else {
      if (!node["chi"][0].IsSequence() ||
          node["chi"][0].size() != ngroups_) {
        std::stringstream mssg;
        mssg << "Invalid length for chi[0] in material " << id_ << ".";
        throw std::runtime_error(mssg.str());
      }

      for (std::size_t ei = 0; ei < ngroups_; ei++) {
        for (std::size_t eo = 0; eo < ngroups_; eo++) {
          chi_[ei][eo] = node["chi"][0][eo].as<float>();

          if (chi_[ei][eo] < 0.) {
            std::stringstream mssg;
            mssg << "chi[" << ei << "][" << eo << "] is negative in material "
                 << id_ << ".";
            throw std::runtime_error(mssg.str());
          }
        }
      }
    }
  }

  //===========================================================================
  // Get the delayed group info
  std::vector<float> P_delayed_grp;
  std::vector<float> delayed_constants;
  if (node["delayed_groups"] && node["delayed_groups"].IsMap()) {
    // Get the probability of each group
    if (!node["delayed_groups"]["probabilities"] ||
        !node["delayed_groups"]["probabilities"].IsSequence()) {
      std::stringstream mssg;
      mssg << "No probabilities entry in delayed_groups for material " << id_
           << ".";
      throw std::runtime_error(mssg.str());
    }
    P_delayed_grp =
        node["delayed_groups"]["probabilities"].as<std::vector<float>>();

    // Get the decay constant of each group
    if (!node["delayed_groups"]["constants"] ||
        !node["delayed_groups"]["constants"].IsSequence()) {
      std::stringstream mssg;
      mssg << "No constants entry in delayed_groups for material " << id_ << ".";
      throw std::runtime_error(mssg.str());
    }
    delayed_constants =
        node["delayed_groups"]["constants"].as<std::vector<float>>();

    // Make sure both have the same size
    if (P_delayed_grp.size() != delayed_constants.size()) {
      std::stringstream mssg;
      mssg << "In delayed_groups entry for material " << id_ << ", ";
      mssg << "probabilities and constants entries have different sizes.";
      throw std::runtime_error(mssg.str());
    }
  } else if (node["delayed_groups"]) {
    std::stringstream mssg;
    mssg << "Invalid delayed_groups entry in material " << id_ << ".";
    throw std::runtime_error(mssg.str());
  }

  //===========================================================================
  // Get the group speeds
  velocity_ = std::vector<float>(ngroups_, 1.);
  if (node["group-speeds"]) {
    if (!node["group-speeds"].IsSequence() ||
        node["group-speeds"].size() != ngroups_) {
      std::stringstream mssg;
      mssg << "Invalid group-speeds entry in material " << id_ << ".";
      throw std::runtime_error(mssg.str());
    }
    velocity_ = node["group-speeds"].as<std::vector<float>>();
  } else {
    std::stringstream mssg;
    mssg << "Missing group-speeds entry in material " << id_ << ".";
    throw std::runtime_error(mssg.str());
  }
}