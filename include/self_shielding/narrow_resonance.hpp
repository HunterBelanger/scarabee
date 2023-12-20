#ifndef SCARABEE_NARROW_RESONANCE_H
#define SCARABEE_NARROW_RESONANCE_H

#include <PapillonNDL/tabulated_1d.hpp>

#include <self_shielding/macro_flux_spectrum.hpp>

class NarrowResonance {
  public:
    NarrowResonance(const MacroFluxSpectrum& spectrum,
                    const pndl::Tabulated1D& total_xs,
                    const double& background_xs):
                    spectrum_(spectrum),
                    total_xs_(total_xs),
                    background_xs_(background_xs) {}

    double operator()(const double& E) const {
      return spectrum_(E) / (total_xs_(E) + background_xs_);
    }

    double background_xs() const { return background_xs_; }

    void set_background_xs(const double& bkgrnd_xs)  { background_xs_ = bkgrnd_xs; }

    const pndl::Tabulated1D& total_xs() const { return total_xs_; }

    const MacroFluxSpectrum& spectrum() const { return spectrum_; }

  private:
    MacroFluxSpectrum spectrum_;
    pndl::Tabulated1D total_xs_;
    double background_xs_;
};

#endif