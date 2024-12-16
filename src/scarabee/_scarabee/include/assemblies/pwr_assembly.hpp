#ifndef PWR_ASSEMBLY_H
#define PWR_ASSEMBLY_H

#include <data/nd_library.hpp>
#include <data/material.hpp>
#include <data/cross_section.hpp>
#include <cylindrical_cell.hpp>
#include <cylindrical_flux_solver.hpp>
#include <data/diffusion_cross_section.hpp>
#include <diffusion/diffusion_data.hpp>
#include <assemblies/pins/fuel_pin.hpp>
#include <assemblies/pins/guide_tube.hpp>
#include <assemblies/pins/burnable_poison_pin.hpp>
#include <moc/quadrature/polar_quadrature.hpp>
#include <moc/quadrature/yamamoto_tabuchi.hpp>
#include <moc/cartesian_2d.hpp>
#include <moc/moc_driver.hpp>
#include <utils/criticality_spectrum.hpp>
#include <utils/serialization.hpp>

#include <xtensor/xtensor.hpp>

#include <cereal/cereal.hpp>
#include <cereal/types/variant.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/optional.hpp>
#include <cereal/types/string.hpp>

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace scarabee {

using Pin = std::variant<std::shared_ptr<FuelPin>, std::shared_ptr<GuideTube>,
                         std::shared_ptr<BurnablePoisonPin>>;

class PWRAssembly {
 public:
  PWRAssembly(double pitch, std::shared_ptr<Material> moderator,
              std::pair<std::size_t, std::size_t> shape,
              std::shared_ptr<NDLibrary> ndl);

  std::shared_ptr<NDLibrary> ndl() const { return ndl_; }
  void set_ndl(std::shared_ptr<NDLibrary> ndl) { ndl_ = ndl; }

  const std::optional<std::string> criticality_spectrum_method() const {
    return criticality_spectrum_method_;
  }
  void set_criticality_spectrum_method(const std::optional<std::string>& csm);

  const std::vector<std::pair<std::size_t, std::size_t>>& condensation_scheme()
      const {
    return condensation_scheme_;
  }
  void set_condensation_scheme(
      const std::vector<std::pair<std::size_t, std::size_t>>& cs) {
    condensation_scheme_ = cs;
  }

  const std::vector<std::pair<std::size_t, std::size_t>>&
  few_group_condensation_scheme() const {
    return few_group_condensation_scheme_;
  }
  void set_few_group_condensation_scheme(
      const std::vector<std::pair<std::size_t, std::size_t>>& cs) {
    few_group_condensation_scheme_ = cs;
  }

  const std::vector<Pin>& pins() const { return pins_; }
  void set_pins(const std::vector<Pin>&);

  std::pair<std::size_t, std::size_t> shape() const { return shape_; }

  std::shared_ptr<Material> moderator() const { return moderator_; }
  void set_moderator(std::shared_ptr<Material> mod);

  std::shared_ptr<CrossSection> moderator_xs() const { return moderator_xs_; }

  std::shared_ptr<CrossSection> average_fuel_pin() const { return avg_fp_; }

  double pitch() const { return pitch_; }

  std::uint32_t num_azimuthal_angles() const { return num_azimuthal_angles_; }
  void set_num_azimuthal_angles(std::uint32_t n);

  double track_spacing() const { return track_spacing_; }
  void set_track_spacing(double t);

  bool anisotropic() const { return anisotropic_; }
  void set_anisotropic(bool a) { anisotropic_ = a; }

  PolarQuadrature polar_quadrature() const { return polar_quadrature_; }
  void set_polar_quadrature(PolarQuadrature pq) { polar_quadrature_ = pq; }

  BoundaryCondition boundary_conditions() const { return boundary_conditions_; }
  void set_boundary_conditions(BoundaryCondition bc) {
    boundary_conditions_ = bc;
  }

  std::uint32_t dancoff_num_azimuthal_angles() const {
    return dancoff_num_azimuthal_angles_;
  }
  void set_dancoff_num_azimuthal_angles(std::uint32_t n);

  double dancoff_track_spacing() const { return dancoff_track_spacing_; }
  void set_dancoff_track_spacing(double t);

  PolarQuadrature dancoff_polar_quadrature() const {
    return dancoff_polar_quadrature_;
  }
  void set_dancoff_polar_quadrature(PolarQuadrature pq) {
    dancoff_polar_quadrature_ = pq;
  }

  bool plot_assembly() const { return plot_assembly_; }
  void set_plot_assembly(bool pa) { plot_assembly_ = pa; }

  double flux_tolerance() const { return flux_tolerance_; }
  void set_flux_tolerance(double ftol);

  double keff_tolerance() const { return keff_tolerance_; }
  void set_keff_tolerance(double ftol);

  const std::vector<double>& fuel_dancoff_corrections() const {
    return fuel_dancoff_corrections_;
  };
  const std::vector<double>& clad_dancoff_corrections() const {
    return clad_dancoff_corrections_;
  };

  const xt::xtensor<double, 2>& form_factors() const { return form_factors_; }
  const xt::xtensor<double, 2>& adf() const { return adf_; }
  const xt::xtensor<double, 2>& cdf() const { return cdf_; }
  std::shared_ptr<DiffusionCrossSection> diffusion_xs() const {
    return diffusion_xs_;
  }
  std::shared_ptr<MOCDriver> moc() const { return moc_; }
  std::shared_ptr<Cartesian2D> moc_geom() const { return moc_geom_; }

  void solve();
  void save_diffusion_data(const std::string& fname) const;

  void save(const std::string& fname) const;
  static std::shared_ptr<PWRAssembly> load(const std::string& fname);

 private:
  double pitch_;
  std::pair<std::size_t, std::size_t> shape_;
  std::shared_ptr<NDLibrary> ndl_;
  std::shared_ptr<Material> moderator_{nullptr};
  std::shared_ptr<CrossSection> moderator_xs_{nullptr};

  std::vector<std::pair<std::size_t, std::size_t>> condensation_scheme_;
  std::vector<std::pair<std::size_t, std::size_t>>
      few_group_condensation_scheme_;

  std::vector<Pin> pins_;

  // MOC parameters for computing dancoff corrections
  std::uint32_t dancoff_num_azimuthal_angles_{64};
  double dancoff_track_spacing_{0.05};
  PolarQuadrature dancoff_polar_quadrature_{YamamotoTabuchi<6>()};

  // MOC parameters for assembly calculation
  std::uint32_t num_azimuthal_angles_{64};
  double track_spacing_{0.05};
  double keff_tolerance_{1.0e-5};
  double flux_tolerance_{1.0e-5};
  PolarQuadrature polar_quadrature_{YamamotoTabuchi<6>()};
  BoundaryCondition boundary_conditions_{BoundaryCondition::Reflective};
  bool anisotropic_{false};

  bool plot_assembly_{false};
  std::shared_ptr<Cartesian2D> moc_geom_{nullptr};
  std::shared_ptr<MOCDriver> moc_{nullptr};
  std::shared_ptr<DiffusionCrossSection> diffusion_xs_{nullptr};
  xt::xtensor<double, 2> form_factors_;
  xt::xtensor<double, 2> adf_;
  xt::xtensor<double, 2> cdf_;
  std::vector<double> fuel_dancoff_corrections_;
  std::vector<double> clad_dancoff_corrections_;
  std::optional<std::string> criticality_spectrum_method_{"P1"};
  std::shared_ptr<CriticalitySpectrum> criticality_spectrum_{nullptr};

  std::vector<std::shared_ptr<CylindricalCell>> pin_1d_cells;
  std::vector<std::shared_ptr<CylindricalFluxSolver>> pin_1d_fluxes;
  std::shared_ptr<CrossSection> avg_fp_{nullptr};

  void get_fuel_dancoff_corrections();
  void get_clad_dancoff_corrections();
  void pin_cell_calc();
  void condense_xs();
  void moc_calc();
  void criticality_spectrum();
  void compute_form_factors();
  void few_group_xs();
  void compute_adf_cdf();

  enum class DancoffMaterial { Fuel, Clad };
  double isolated_fuel_pin_flux(DancoffMaterial dm) const;
  double isolated_guide_tube_flux() const;
  double isolated_burnable_poison_tube_flux() const;

  std::vector<double> compute_avg_surface_flx(
      const std::vector<std::pair<std::size_t, double>>& segments) const;
  std::vector<double> compute_avg_flx(const Vector& r,
                                      const Direction& u) const;

  friend class cereal::access;
  PWRAssembly() {}
  template <class Archive>
  void serialize(Archive& arc) {
    arc(CEREAL_NVP(pitch_), CEREAL_NVP(shape_), CEREAL_NVP(moderator_),
        CEREAL_NVP(moderator_xs_), CEREAL_NVP(condensation_scheme_),
        CEREAL_NVP(few_group_condensation_scheme_), CEREAL_NVP(pins_),
        CEREAL_NVP(dancoff_num_azimuthal_angles_),
        CEREAL_NVP(dancoff_track_spacing_),
        CEREAL_NVP(dancoff_polar_quadrature_),
        CEREAL_NVP(num_azimuthal_angles_), CEREAL_NVP(track_spacing_),
        CEREAL_NVP(keff_tolerance_), CEREAL_NVP(flux_tolerance_),
        CEREAL_NVP(polar_quadrature_), CEREAL_NVP(boundary_conditions_),
        CEREAL_NVP(anisotropic_), CEREAL_NVP(plot_assembly_),
        CEREAL_NVP(moc_geom_), CEREAL_NVP(moc_), CEREAL_NVP(diffusion_xs_),
        CEREAL_NVP(form_factors_), CEREAL_NVP(adf_), CEREAL_NVP(cdf_),
        CEREAL_NVP(fuel_dancoff_corrections_),
        CEREAL_NVP(clad_dancoff_corrections_),
        CEREAL_NVP(criticality_spectrum_method_),
        CEREAL_NVP(criticality_spectrum_), CEREAL_NVP(pin_1d_cells),
        CEREAL_NVP(pin_1d_fluxes), CEREAL_NVP(avg_fp_));
  }
};

};  // namespace scarabee

#endif