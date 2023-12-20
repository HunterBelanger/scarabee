#ifndef SCARABEE_MACRO_FLUX_SPECTRUM_H
#define SCARABEE_MACRO_FLUX_SPECTRUM_H

#include <self_shielding/maxwellian.hpp>
#include <self_shielding/watt.hpp>
#include <self_shielding/one_over_e.hpp>

#include <exception>
#include <variant>

class MacroFluxSpectrum {
  public:
    using Component = std::variant<Maxwellian,Watt,OneOverE>;

    MacroFluxSpectrum():
      thermal_{Maxwellian(0.0256)},
      epithermal_{OneOverE()},
      fission_{Watt(0.965E6, 2.29E-6)},
      thermal_point_{0.1},
      fission_point_{820.3E3} {
      set_constants();
    }

    MacroFluxSpectrum(Component thermal, Component epithermal, Component fission,
             double thermal_point, double fission_point):
             thermal_{thermal},
             epithermal_{epithermal},
             fission_{fission},
             thermal_point_{thermal_point},
             fission_point_{fission_point} {
      if (thermal_point_ < 0.) {
        throw std::runtime_error("Thermal - epithermal transition energy must be > 0.");
      }

      if (fission_point_ < 0.) {
        throw std::runtime_error("Epithermal - fission transition energy must be > 0.");
      }

      if (fission_point_ <= thermal_point_) {
        throw std::runtime_error("Epithermal-fission transition energy must be > thermal-epithermal transition energy.");
      }

      set_constants();
    }

    double operator()(const double& E) const {
      if (E < thermal_point_) {
        return std::visit(EvalVisitor{E}, thermal_);
      } else if (E < fission_point_) {
        return std::visit(EvalVisitor{E}, epithermal_);
      } else {
        return std::visit(EvalVisitor{E}, fission_);
      }
    }

  private:
    Component thermal_;
    Component epithermal_;
    Component fission_;
    double thermal_point_, fission_point_;

    struct EvalVisitor {
      EvalVisitor(const double& E): E{E} {}

      double operator()(const Maxwellian& comp) const {
        return comp(E);
      }

      double operator()(const OneOverE& comp) const {
        return comp(E);
      }

      double operator()(const Watt& comp) const {
        return comp(E);
      }

      double E;
    };

    struct SetConstVisitor {
      SetConstVisitor(const double& C): C{C} {}

      void operator()(Maxwellian& comp) const {
        comp.set_constant(C);
      }

      void operator()(OneOverE& comp) const {
        comp.set_constant(C);
      }

      void operator()(Watt& comp) const {
        comp.set_constant(C);
      }

      double C;
    };

    void set_constants() {
      // Set constants so that things align
      const double thermal_thermal_point_val = std::visit(EvalVisitor{thermal_point_}, thermal_);
      const double epithermal_thermal_point_val = std::visit(EvalVisitor{thermal_point_}, epithermal_);
      const double epi_const = thermal_thermal_point_val / epithermal_thermal_point_val;
      std::visit(SetConstVisitor{epi_const}, epithermal_);

      const double epithermal_fission_point_val = std::visit(EvalVisitor{fission_point_}, epithermal_);
      const double fission_fission_point_val = std::visit(EvalVisitor{fission_point_}, fission_);
      const double fis_const = epithermal_fission_point_val / fission_fission_point_val;
      std::visit(SetConstVisitor{fis_const}, fission_);
    }
};

#endif