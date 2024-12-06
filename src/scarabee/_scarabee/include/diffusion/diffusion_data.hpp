#ifndef SCARABEE_DIFFUSION_DATA_H
#define SCARABEE_DIFFUSION_DATA_H

#include <data/diffusion_cross_section.hpp>

#include <xtensor/xtensor.hpp>

#include <memory>
#include <string>

namespace scarabee {

class DiffusionData {
 public:
  enum class ADF : std::size_t { YP = 0, XP = 1, YN = 2, XN = 3 };
  enum class CDF : std::size_t { I = 0, II = 1, III = 2, IV = 3 };

  DiffusionData(std::shared_ptr<DiffusionCrossSection> xs);

  DiffusionData(std::shared_ptr<DiffusionCrossSection> xs,
                const xt::xtensor<double, 2>& form_factors);

  DiffusionData(std::shared_ptr<DiffusionCrossSection> xs,
                const xt::xtensor<double, 2>& form_factors,
                const xt::xtensor<double, 2>& adf,
                const xt::xtensor<double, 2>& cdf);

  std::size_t ngroups() const { return xs_->ngroups(); }

  bool fissile() const { return xs_->fissile(); }

  double D(std::size_t g) const { return xs_->D(g); }

  double Ea(std::size_t g) const { return xs_->Ea(g); }

  double Ef(std::size_t g) const { return xs_->Ef(g); }

  double vEf(std::size_t g) const { return xs_->vEf(g); }

  double nu(std::size_t g) const { return xs_->nu(g); }

  double Er(std::size_t g) const { return xs_->Er(g); }

  double chi(std::size_t g) const { return xs_->chi(g); }

  double Es(std::size_t gin, std::size_t gout) const {
    return xs_->Es(gin, gout);
  }

  double Es(std::size_t gin) const { return xs_->Es(gin); }

  double adf_xp(std::size_t g) const {
    return adf_.size() > 0 ? adf_(g, ADF::XP) : 1.;
  }

  double adf_xn(std::size_t g) const {
    return adf_.size() > 0 ? adf_(g, ADF::XN) : 1.;
  }

  double adf_yp(std::size_t g) const {
    return adf_.size() > 0 ? adf_(g, ADF::YP) : 1.;
  }

  double adf_yn(std::size_t g) const {
    return adf_.size() > 0 ? adf_(g, ADF::YN) : 1.;
  }

  double cdf_I(std::size_t g) const {
    if (cdf_.size() > 0) {
      return cdf_(g, CDF::I);
    } else if (adf_.size() > 0) {
      return 0.5 * (adf_xp(g) + adf_yp(g));
    } else {
      return 1.;
    }
  }

  double cdf_II(std::size_t g) const {
    if (cdf_.size() > 0) {
      return cdf_(g, CDF::II);
    } else if (adf_.size() > 0) {
      return 0.5 * (adf_xn(g) + adf_yp(g));
    } else {
      return 1.;
    }
  }

  double cdf_III(std::size_t g) const {
    if (cdf_.size() > 0) {
      return cdf_(g, CDF::III);
    } else if (adf_.size() > 0) {
      return 0.5 * (adf_xn(g) + adf_yn(g));
    } else {
      return 1.;
    }
  }

  double cdf_IV(std::size_t g) const {
    if (cdf_.size() > 0) {
      return cdf_(g, CDF::IV);
    } else if (adf_.size() > 0) {
      return 0.5 * (adf_xp(g) + adf_yn(g));
    } else {
      return 1.;
    }
  }

  void rotate_clockwise();
  void rotate_counterclockwise();

  const std::string& xs_name() const { return xs_->name(); }

  std::shared_ptr<DiffusionCrossSection> xs() const { return xs_; }

  const std::string& name() const { return name_; }
  void set_name(const std::string& new_name) { name_ = new_name; }

  const xt::xtensor<double, 2>& form_factors() const { return form_factors_; }
  void set_form_factors(const xt::xtensor<double, 2>& ff);

  const xt::xtensor<double, 2>& adf() const { return adf_; }
  void set_adf(const xt::xtensor<double, 2>& adf);

  const xt::xtensor<double, 2>& cdf() const { return cdf_; }
  void set_cdf(const xt::xtensor<double, 2>& cdf);

  void save(const std::string& fname) const;
  static std::shared_ptr<DiffusionData> load(const std::string& fname);

 private:
  std::shared_ptr<DiffusionCrossSection> xs_;
  xt::xtensor<double, 2> form_factors_;
  xt::xtensor<double, 2> adf_;  // group then ADF direction
  xt::xtensor<double, 2> cdf_;  // group then CDF direction
  std::string name_;
};

}  // namespace scarabee

#endif
