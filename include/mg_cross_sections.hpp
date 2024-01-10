#ifndef SCARABEE_MG_CROSS_SECTIONS_H
#define SCARABEE_MG_CROSS_SECTIONS_H

#include <utils/scarabee_exception.hpp>

#include <cstdint>
#include <memory>
#include <optional>
#include <sstream>
#include <vector>

#include <ndarray.hpp>

class MGCrossSections {
 public:
  MGCrossSections();

  std::uint32_t ngroups() const {
    return static_cast<std::uint32_t>(Et_.size());
  }

  std::uint32_t ndelayed_families() const {
    if (fissile_) return static_cast<std::uint32_t>(nu_.shape()[0] - 1);
    return 0;
  }

  std::uint32_t max_legendre_order() const {
    return static_cast<std::uint32_t>(Es_.shape()[0] - 1);
  }

  bool fissile() const { return fissile_; }

  double Et(std::uint32_t g) const { return Et_[g]; }

  double Etr(std::uint32_t g) const { return Etr_[g]; }

  double Ea(std::uint32_t g) const { return Ea_[g]; }

  double Ef(std::uint32_t g) const {
    if (fissile_) return Ef_[g];
    return 0.;
  }

  double Er(std::uint32_t g) const { return Ea(g) + Es(g); }

  double Er_tr(std::uint32_t g) const { return Ea(g) + Es_tr(g); }

  double nu_prompt(std::uint32_t g) const {
    if (fissile_) return nu_(0, g);
    return 0.;
  }

  double nu_delayed(std::uint32_t g, std::uint32_t i) const {
    if (fissile_) {
      if (i < 1 || i > this->ndelayed_families()) {
        std::stringstream mssg;
        mssg << "Invalid delayed group " << i << ".";
        throw ScarabeeException(mssg.str());
      }
      return nu_(i, g);
    }
    return 0.;
  }

  double nu_delayed(std::uint32_t g) const {
    if (fissile_) {
      double nu_d = 0.;
      for (std::uint32_t i = 1; i <= this->ndelayed_families(); i++) {
        nu_d += nu_(i, g);
      }
      return nu_d;
    }
    return 0.;
  }

  double nu(std::uint32_t g) const { return nu_prompt(g) + nu_delayed(g); }

  double chi_prompt(std::uint32_t g) const {
    if (fissile_) return chi_(0, g);
    return 0.;
  }

  double chi_delayed(std::uint32_t g, std::uint32_t i) const {
    if (fissile_) {
      if (i < 1 || i > this->ndelayed_families()) {
        std::stringstream mssg;
        mssg << "Invalid delayed group " << i << ".";
        throw ScarabeeException(mssg.str());
      }
      return chi_(i, g);
    }
    return 0.;
  }

  double Es(std::uint32_t gin, std::uint32_t gout) const {
    return Es_(0, gin, gout);
  }

  double Es(std::uint32_t gin, std::uint32_t gout, std::uint32_t l) const {
    return Es_(l, gin, gout);
  }

  double Es(std::uint32_t gin) const {
    double Es_ret = 0.;
    for (std::uint32_t gout = 0; gout < this->ngroups(); gout++) {
      Es_ret += Es_(0, gin, gout);
    }
    return Es_ret;
  }

  double Es_tr(std::uint32_t gin, std::uint32_t gout) const {
    return Es_tr_(gin, gout);
  }

  double Es_tr(std::uint32_t gin) const {
    double Es_ret = 0.;
    for (std::uint32_t gout = 0; gout < this->ngroups(); gout++) {
      Es_ret += Es_tr_(gin, gout);
    }
    return Es_ret;
  }

  // Operators for constructing compound cross sections
  MGCrossSections operator+(const MGCrossSections& R) const;
  MGCrossSections operator*(double N) const;
  MGCrossSections& operator+=(const MGCrossSections& R);
  MGCrossSections& operator*=(double N);

 //private:
  NDArray<double> nu_;       // Fission yields
  NDArray<double> chi_;      // Fission spectrum
  NDArray<double> Es_;       // Scattering matrix
  NDArray<double> Es_tr_;    // Transport Corrected Scattering matrix
  std::vector<double> Et_;   // Total xs
  std::vector<double> Etr_;  // Transport xs
  std::vector<double> Ea_;   // Absorption xs
  std::vector<double> Ef_;   // Fission xs
  bool fissile_;             // Fissile indicator

  void generate_transport_corrected_data();
};

#endif
