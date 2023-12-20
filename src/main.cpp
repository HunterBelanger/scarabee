#include <material.hpp>
#include <geometry.hpp>
#include <solver.hpp>

#include <exception>
#include <iostream>
#include <map>
#include <string>

#include <docopt.h>

#include <yaml-cpp/yaml.h>

static const std::string help =
"Scarabée\n\n"

"Usage:\n"
"   scarabée <input>\n"
"   scarabée (-h | --help)\n"
"   scarabée --version\n\n"

" Options:\n"
"   <input>       Name of the input file.\n"
"   --version     Show the version information.\n"
"   -h --help     Show this help message.\n";

/*
int main(int argc, char** argv) {
  // Initialize docopt
  std::map<std::string, docopt::value> args = docopt::docopt(help, {argv + 1, argv + argc}, true, "Scarabée 0.1.0");

  // Read intput file into YAML node
  YAML::Node input = YAML::LoadFile(args["<input>"].asString());

  // Read all materials
  const auto& mats = input["materials"];
  if (mats.IsSequence() == false) {
    throw std::runtime_error("materials must be a sequence of material entries.");
  }

  for (std::size_t m = 0; m < mats.size(); m++) {
    Material mat(mats[m]);
    materials.insert({mat.id(), mat});
  }
  std::cout << "Made materials...\n";

  // Read the geometry
  const auto& geom_node = input["geometry"];
  Geom1D geom(geom_node);
  std::cout << "Made geometry...\n";

  // Make problem solver
  StaticSolver1D solver(std::move(geom));
  std::cout << "Made solver...\n";

  // Solve problem
  solver.solve();

  return 0;
}*/


#include <self_shielding/macro_flux_spectrum.hpp>
#include <self_shielding/narrow_resonance.hpp>
#include <tools/gauss_kronrod.hpp>

#include <ImApp/imapp.hpp>

#include <PapillonNDL/st_neutron.hpp>

#include <cmath>
#include <vector>
#include <memory>

const std::vector<double> xmas_172 {
    1.00001e-05, 3.00000e-03, 5.00000e-03, 6.90000e-03, 1.00000e-02,
    1.50000e-02, 2.00000e-02, 2.50000e-02, 3.00000e-02, 3.50000e-02,
    4.20000e-02, 5.00000e-02, 5.80000e-02, 6.70000e-02, 7.70000e-02,
    8.00000e-02, 9.50000e-02, 1.00001e-01, 1.15000e-01, 1.34000e-01,
    1.40000e-01, 1.60000e-01, 1.80000e-01, 1.89000e-01, 2.20000e-01,
    2.48000e-01, 2.80000e-01, 3.00000e-01, 3.14500e-01, 3.20000e-01,
    3.50000e-01, 3.91000e-01, 4.00000e-01, 4.33000e-01, 4.85000e-01,
    5.00000e-01, 5.40000e-01, 6.25000e-01, 7.05000e-01, 7.80000e-01,
    7.90000e-01, 8.50000e-01, 8.60000e-01, 9.10000e-01, 9.30000e-01,
    9.50000e-01, 9.72000e-01, 9.86000e-01, 9.96000e-01, 1.02000e+00,
    1.03500e+00, 1.04500e+00, 1.07100e+00, 1.09700e+00, 1.11000e+00,
    1.12535e+00, 1.15000e+00, 1.17000e+00, 1.23500e+00, 1.30000e+00,
    1.33750e+00, 1.37000e+00, 1.44498e+00, 1.47500e+00, 1.50000e+00,
    1.59000e+00, 1.67000e+00, 1.75500e+00, 1.84000e+00, 1.93000e+00,
    2.02000e+00, 2.10000e+00, 2.13000e+00, 2.36000e+00, 2.55000e+00,
    2.60000e+00, 2.72000e+00, 2.76792e+00, 3.30000e+00, 3.38075e+00,
    4.00000e+00, 4.12925e+00, 5.04348e+00, 5.34643e+00, 6.16012e+00,
    7.52398e+00, 8.31529e+00, 9.18981e+00, 9.90555e+00, 1.12245e+01,
    1.37096e+01, 1.59283e+01, 1.94548e+01, 2.26033e+01, 2.49805e+01,
    2.76077e+01, 3.05113e+01, 3.37201e+01, 3.72665e+01, 4.01690e+01,
    4.55174e+01, 4.82516e+01, 5.15780e+01, 5.55951e+01, 6.79041e+01,
    7.56736e+01, 9.16609e+01, 1.36742e+02, 1.48625e+02, 2.03995e+02,
    3.04325e+02, 3.71703e+02, 4.53999e+02, 6.77287e+02, 7.48518e+02,
    9.14242e+02, 1.01039e+03, 1.23410e+03, 1.43382e+03, 1.50733e+03,
    2.03468e+03, 2.24867e+03, 3.35463e+03, 3.52662e+03, 5.00451e+03,
    5.53084e+03, 7.46586e+03, 9.11882e+03, 1.11378e+04, 1.50344e+04,
    1.66156e+04, 2.47875e+04, 2.73944e+04, 2.92830e+04, 3.69786e+04,
    4.08677e+04, 5.51656e+04, 6.73795e+04, 8.22975e+04, 1.11090e+05,
    1.22773e+05, 1.83156e+05, 2.47235e+05, 2.73237e+05, 3.01974e+05,
    4.07622e+05, 4.50492e+05, 4.97871e+05, 5.50232e+05, 6.08101e+05,
    8.20850e+05, 9.07180e+05, 1.00259e+06, 1.10803e+06, 1.22456e+06,
    1.35335e+06, 1.65299e+06, 2.01897e+06, 2.23130e+06, 2.46597e+06,
    3.01194e+06, 3.67879e+06, 4.49329e+06, 5.48812e+06, 6.06531e+06,
    6.70320e+06, 8.18731e+06, 1.00000e+07, 1.16183e+07, 1.38403e+07,
    1.49182e+07, 1.73325e+07, 1.96403e+07};

std::pair<std::vector<double>, std::vector<double>> group_xs(const std::vector<double>& group_bounds,
                                                             const std::function<double(double)>& flux,
                                                             const std::vector<double>& energy,
                                                             const std::vector<double>& tot_xs) {
  std::vector<double> group_bounds_;
  std::vector<double> group_xs;
  
  std::vector<double> flux_vals(energy.size());
  std::vector<double> flux_xs_vals(energy.size());
  for (std::size_t i = 0; i < energy.size(); i++) {
    flux_vals[i] = flux(energy[i]);
    flux_xs_vals[i] = flux(energy[i]) * tot_xs[i];
  }

  // Create functions which will be responsible for the integrations
  pndl::Tabulated1D flux_func(pndl::Interpolation::LinLin, energy, flux_vals);
  pndl::Tabulated1D flux_xs_func(pndl::Interpolation::LinLin, energy, flux_xs_vals);

  // Generate each group constant
  for (std::size_t i = 0; i < group_bounds.size()-1; i++) {
    double Elow = group_bounds[i];
    double Ehi = group_bounds[i+1];

    double grp_xs = flux_xs_func.integrate(Elow, Ehi) / flux_func.integrate(Elow, Ehi);

    group_bounds_.push_back(Elow);
    group_xs.push_back(grp_xs);

    group_bounds_.push_back(Ehi);
    group_xs.push_back(grp_xs);
  }

  return {group_bounds_, group_xs};
}

struct CosineAngleCM {
  double operator()(const double& E, const double& Eout) const {
    const double R = A_ * std::sqrt(1. + (((A_+1.)*Q_)/(A_*E)));

    const double num = (Eout * (1. + A_) * (1. + A_)) - (E * (1. + (R * R)));
    const double denom = 2. * R * E;

    return num / denom;
  }

  double A_, Q_;
};

double R(const double& E, const double& A, const double& Q) {
  return A * std::sqrt(1. + (((A+1.)*Q)/(A*E)));
}

double lab_cosine(const double& R, const double& omega) {
  return (1. + R*omega) / std::sqrt(1. + R*R + 2.*R*omega);
}

double cm_cosine(const double& E, const double& Eout, const double& A, const double& R) {
  const double num = (Eout * (1. + A) * (1. + A)) - (E * (1. + (R * R)));
  const double denom = 2. * R * E;
  return num / denom;
}

struct ScatteringFeedFunc {

  double operator()(const double& E) const {
    const double R_ = R(E, A, Q);

    // Get bounds of integration
    double omega_low = cm_cosine(E, Eout_low, A, R_);
    double omega_hi = cm_cosine(E, Eout_hi, A, R_);

    // If scattering angles are outside of [-1,1], we need to truncate the
    // domain of integration.
    if (omega_low < -1. && omega_hi < -1.) { return 0.; }
    else if (omega_low > 1. && omega_hi > 1.) { return 0.; }
    else if (omega_low < -1.) { omega_low = -1.; }
    if (omega_hi > 1.) { omega_hi = 1.; }

    // Create integrand function
    auto func = [this, omega_low, omega_hi, R_, E](const double& omega) {
      const double mu = lab_cosine(R_, omega);
      return this->angle_dist_.pdf(E, omega) * std::legendre(this->l, mu);
    };

    // Integrate function
    GaussKronrodQuadrature<21> gk;
    auto I = gk.integrate(func, omega_low, omega_hi, 0.001, 100);

    // TODO check error on integral

    return I.first;
  }

  pndl::AngleDistribution angle_dist_; // Scattering distribution in CM frame for reaction
  double A, Q;
  double Eout_low, Eout_hi; // Outgoing energy group bounds
  unsigned int l; // Legendre order
};

class PlotLayer : public ImApp::Layer {
  public:
    PlotLayer() {
      pndl::ACE ace("/home/hunter/projects/scarabée/U238.293.6.ace");
      pndl::STNeutron nuc(ace);

      energy = nuc.total_xs().energy();
      for (auto& engy : energy) { engy *= 1.E6; }

      MacroFluxSpectrum spec;
      pndl::Tabulated1D tot_xs(pndl::Interpolation::LinLin, energy, nuc.total_xs().xs());
      NarrowResonance nr(spec, tot_xs, 1.0);

      std::cout << " START \n";
      grouped_xss.push_back(group_xs(xmas_172, spec, energy, tot_xs.y()));
      grouped_xss.push_back(group_xs(xmas_172, nr, energy, tot_xs.y()));
      std::cout << " STOP \n";
    }

    void render() override final {
      ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

      ImGui::SetNextWindowSize({500, 500}, ImGuiCond_Once);
      ImGui::Begin("XS Plot");
      ImPlot::BeginPlot("Grouped XS", ImVec2(-1,-1));
      ImPlot::SetupAxis(ImAxis_X1, "Energy [eV]");
      ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
      ImPlot::SetupAxis(ImAxis_Y1, "Cross Section [barns]");
      ImPlot::SetupAxisScale(ImAxis_Y1, ImPlotScale_Log10);

      ImPlot::PlotLine("Smooth Spectrum", grouped_xss[0].first.data(), grouped_xss[0].second.data(), static_cast<int>(grouped_xss[0].first.size()));
      ImPlot::PlotLine("NR Spectrum", grouped_xss[1].first.data(), grouped_xss[1].second.data(), static_cast<int>(grouped_xss[1].first.size()));

      ImPlot::EndPlot();
      ImGui::End();
    }

  private:
    std::vector<double> energy;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> grouped_xss;
};

int main() {
  std::unique_ptr<ImApp::Layer> layer = std::make_unique<PlotLayer>();
  
  ImApp::App app(1920, 1080, "Scarabee");
  app.push_layer(std::move(layer));
  app.enable_docking();
  app.run();

  return 0;
}