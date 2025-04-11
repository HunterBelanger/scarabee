// Disable parallelization of Eigen in this translation unit, because we will
// be running depletion of different materials in parallel already !
#define EIGEN_DONT_PARALLELIZE

#include <data/depletion_matrix.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>
#include <utils/nuclide_names.hpp>
#include <utils/scarabee_exception.hpp>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <algorithm>
#include <set>

namespace scarabee {

DepletionMatrix::DepletionMatrix(const std::vector<std::string>& nuclides)
    : nuclides_(nuclides),
      matrix_(static_cast<int>(size()), static_cast<int>(size())) {
  // Make sure nuclides are sorted
  std::sort(nuclides_.begin(), nuclides_.end(),
            [](const std::string& n1, const std::string& n2) {
              return nuclide_name_to_za(n1) < nuclide_name_to_za(n2);
            });
}

bool DepletionMatrix::has_nuclide(const std::string& nuclide) const {
  for (const auto& nuc : nuclides_) {
    if (nuc == nuclide) return true;
  }
  return false;
}

std::size_t DepletionMatrix::get_nuclide_index(
    const std::string& nuclide) const {
  for (std::size_t i = 0; i < nuclides_.size(); i++) {
    if (nuclides_[i] == nuclide) return i;
  }

  const auto mssg = "Depletion matrix does not contain \"" + nuclide + "\".";
  spdlog::error(mssg);
  throw ScarabeeException(mssg);

  // NEVER GETS HERE
  return 0;
}

void DepletionMatrix::exponential_product(std::span<double> N,
                                          bool cram48) const {
  if (N.size() != this->size()) {
    const auto mssg =
        "The number of provided number densities does not agree with the size "
        "of the matrix.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  std::span<const std::complex<double>> theta;
  std::span<const std::complex<double>> alpha;
  double alpha0;

  if (cram48) {
    // 48th order CRAM
    theta = std::span<const std::complex<double>>(cram48_theta_.begin(),
                                                  cram48_theta_.end());
    alpha = std::span<const std::complex<double>>(cram48_alpha_.begin(),
                                                  cram48_alpha_.end());
    alpha0 = cram48_alpha0_;
  } else {
    // 16th order CRAM
    theta = std::span<const std::complex<double>>(cram16_theta_.begin(),
                                                  cram16_theta_.end());
    alpha = std::span<const std::complex<double>>(cram16_alpha_.begin(),
                                                  cram16_alpha_.end());
    alpha0 = cram16_alpha0_;
  }

  detail_exp_product(N, theta, alpha, alpha0);
}

void DepletionMatrix::detail_exp_product(
    std::span<double> N, std::span<const std::complex<double>> theta,
    std::span<const std::complex<double>> alpha, double alpha0) const {
  // Initialize the complex matrix
  Eigen::SparseMatrix<std::complex<double>> Acmplx(this->size(), this->size());

  // Complex A is initialized as the current real matrix.
  // Effectively completes Acmplx = matrix_;
  for (int k = 0; k < matrix_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(matrix_, k); it; ++it) {
      Acmplx.coeffRef(it.row(), it.col()) = {it.value(), 0.};
    }
  }
  Acmplx.makeCompressed();

  // Create the solver object
  Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;

  // Initialize the complex vector for N
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> Ncmplx(this->size());
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> Ncmplx_out(
      this->size());
  for (std::size_t n = 0; n < this->size(); n++) {
    Ncmplx(n).real(N[n]);
    Ncmplx(n).imag(0.);
  }

  // Do all iterations for CRAM
  for (std::size_t i = 0; i < theta.size(); i++) {
    // Subtract theta[i] from diagonal
    for (std::size_t j = 0; j < this->size(); j++) {
      Acmplx.coeffRef(j, j) = matrix_.coeff(j, j) - theta[i];
    }

    // Solve the system for Acmplx @ x = N
    solver.analyzePattern(Acmplx);
    solver.factorize(Acmplx);
    if (solver.info() != Eigen::ComputationInfo::Success) {
      const auto mssg =
          "Could not factorize complex depletion matrix in CRAM iterations.";
      spdlog::error(mssg);
      throw ScarabeeException(mssg);
    }
    Ncmplx_out = solver.solve(Ncmplx);

    // Multiply by complex alpha[i]
    Ncmplx_out *= 2. * alpha[i];
    for (std::size_t j = 0; j < this->size(); j++) {
      const double val = Ncmplx(j).real() + Ncmplx_out(j).real();
      Ncmplx(j).real(val);
    }
  }

  // Reassign number densities for the input span
  for (std::size_t n = 0; n < this->size(); n++) {
    N[n] = Ncmplx(n).real() * alpha0;

    if (N[n] < 0.) N[n] = 0.;
  }
}

double compute_loss_term(const ChainEntry& nucinfo,
                         const DepletionReactionRates& nucrr) {
  double loss = 0.;

  // Get losses due to radioactive decay.
  // Must convert half life to decay constant !
  if (nucinfo.half_life()) loss += LN_2 / nucinfo.half_life().value();

  // Accumulate losses due to transmutation.
  if (nucinfo.n_gamma()) loss += nucrr.n_gamma;
  if (nucinfo.n_2n()) loss += nucrr.n_2n;
  if (nucinfo.n_3n()) loss += nucrr.n_3n;
  if (nucinfo.n_p()) loss += nucrr.n_p;
  if (nucinfo.n_alpha()) loss += nucrr.n_alpha;
  if (nucinfo.n_fission()) loss += nucrr.n_fission;

  return loss;
}

void fill_target_gains(DepletionMatrix& /*matrix*/, const NoTarget& /*target*/,
                       const double /*rate*/,
                       const std::size_t /*parent_index*/) {
  // Nothing to do when there is no target
}

void fill_target_gains(DepletionMatrix& matrix, const SingleTarget& target,
                       const double rate, const std::size_t parent_index) {
  const std::size_t target_index = matrix.get_nuclide_index(target.target());
  matrix.ref(target_index, parent_index) += rate;
}

void fill_target_gains(DepletionMatrix& matrix, const BranchingTargets& targets,
                       const double rate, const std::size_t parent_index) {
  for (const auto& branch : targets.branches()) {
    const std::size_t target_index = matrix.get_nuclide_index(branch.target);
    matrix.ref(target_index, parent_index) += rate * branch.branch_ratio;
  }
}

void fill_target_gains(DepletionMatrix& matrix, const FissionYields& fy,
                       const double rate, const double E,
                       const std::size_t parent_index) {
  for (std::size_t fyi = 0; fyi < fy.size(); fyi++) {
    const std::size_t target_index =
        matrix.get_nuclide_index(fy.targets()[fyi]);
    matrix.ref(target_index, parent_index) += rate * fy.yield(fyi, E);
  }
}

void fill_target_gains(DepletionMatrix& matrix, const Target& target,
                       const double rate, const std::size_t parent_index) {
  if (std::holds_alternative<NoTarget>(target)) {
    fill_target_gains(matrix, std::get<NoTarget>(target), rate, parent_index);
  } else if (std::holds_alternative<SingleTarget>(target)) {
    fill_target_gains(matrix, std::get<SingleTarget>(target), rate,
                      parent_index);
  } else {
    // Must hold BranchingTargets
    fill_target_gains(matrix, std::get<BranchingTargets>(target), rate,
                      parent_index);
  }
}

void fill_gain_terms(DepletionMatrix& matrix, const ChainEntry& nucinfo,
                     const DepletionReactionRates& nucrr, const std::size_t i) {
  if (nucinfo.decay_targets())
    fill_target_gains(matrix, nucinfo.decay_targets().value(),
                      LN_2 / nucinfo.half_life().value(), i);

  if (nucinfo.n_2n() && nucrr.n_2n > 0.)
    fill_target_gains(matrix, nucinfo.n_2n().value(), nucrr.n_2n, i);

  if (nucinfo.n_3n() && nucrr.n_3n > 0.)
    fill_target_gains(matrix, nucinfo.n_3n().value(), nucrr.n_3n, i);

  if (nucinfo.n_p() && nucrr.n_p > 0.)
    fill_target_gains(matrix, nucinfo.n_p().value(), nucrr.n_p, i);

  if (nucinfo.n_alpha() && nucrr.n_alpha > 0.)
    fill_target_gains(matrix, nucinfo.n_alpha().value(), nucrr.n_alpha, i);

  if (nucinfo.n_fission() && nucrr.n_fission > 0.)
    fill_target_gains(matrix, nucinfo.n_fission().value(), nucrr.n_fission,
                      nucrr.average_fission_energy, i);
}

std::shared_ptr<DepletionMatrix> build_depletion_matrix(
    std::shared_ptr<DepletionChain> chain, std::shared_ptr<Material> mat,
    std::span<const double> flux, std::shared_ptr<NDLibrary> ndl) {
  // Start by making a set of initial nuclides.
  // These are only the depletable nuclides !!
  std::set<std::string> initial_dep_nuclides;
  for (const auto& nuc : mat->composition().nuclides) {
    initial_dep_nuclides.insert(nuclide_name_to_simple_name(nuc.name));
  }

  // Get the sorted list of all possible targets
  std::vector<std::string> all_targets =
      chain->descend_chains(initial_dep_nuclides);

  // Get the vector of all reaction rate objects for the nuclides in the
  // material
  std::vector<DepletionReactionRates> nuc_rrs =
      mat->compute_depletion_reaction_rates(flux, ndl);

  // With this, we can build the matrix object
  std::shared_ptr<DepletionMatrix> matrix_ptr =
      std::make_shared<DepletionMatrix>(all_targets);
  DepletionMatrix& matrix = *matrix_ptr;

  // Now, go through each initial nuclide that is in the material and has
  // reaction rates. We fill the matrix elements for targets of this nuclide.
  for (const auto& nucrr : nuc_rrs) {
    // Remove the nuclide from the all_targets list if it has reaction rates.
    // This way, after the reaction rate terms, we can go through all the
    // remaining nuclides in all_targets which don't have a reaction rate, but
    // which almost certainly have radioactive decay !
    std::remove(all_targets.begin(), all_targets.end(), nucrr.nuclide);

    // Skip the nuclide if it isn't in the depletion chain, even if it might
    // have depletion cross section nuclear data and reaction rates.
    if (chain->holds_nuclide_data(nucrr.nuclide) == false) continue;

    // Get the depletion chain entry for the nuclide
    const auto& nucinfo = chain->nuclide_data(nucrr.nuclide);

    // Get the index of the nuclide
    const std::size_t i = matrix.get_nuclide_index(nucrr.nuclide);

    // First, account for all losses due to decay and transmutation.
    matrix.ref(i, i) -= compute_loss_term(nucinfo, nucrr);

    // Now we fill the gain terms for all targets
    fill_gain_terms(matrix, nucinfo, nucrr, i);
  }

  // Now we need to go over all the nuclides that remain in all_targets, which
  // weren't in the material and therefore don't have reaction rates.
  DepletionReactionRates zero_rr;
  for (const auto& nuclide : all_targets) {
    // Skip the nuclide if it isn't in the depletion chain, even if it might
    // have depletion cross section nuclear data and reaction rates.
    if (chain->holds_nuclide_data(nuclide) == false) continue;

    // Get the depletion chain entry for the nuclide
    const auto& nucinfo = chain->nuclide_data(nuclide);

    // Get the index of the nuclide
    const std::size_t i = matrix.get_nuclide_index(nuclide);

    // First, account for all losses due to decay and transmutation.
    matrix.ref(i, i) -= compute_loss_term(nucinfo, zero_rr);

    // Now we fill the gain terms for all targets
    fill_gain_terms(matrix, nucinfo, zero_rr, i);
  }

  matrix.compress();
  return matrix_ptr;
}

// Constants for taking matrix exponential with CRAM.
// Tables can be found in reference [1].

const std::array<std::complex<double>, 8> DepletionMatrix::cram16_alpha_{
    std::complex<double>{5.464930576870210e+3, -3.797983575308356e+4},
    std::complex<double>{9.045112476907548e+1, -1.115537522430261e+3},
    std::complex<double>{2.344818070467641e+2, -4.228020157070496e+2},
    std::complex<double>{9.453304067358312e+1, -2.951294291446048e+2},
    std::complex<double>{7.283792954673409e+2, -1.205646080220011e+5},
    std::complex<double>{3.648229059594851e+1, -1.155509621409682e+2},
    std::complex<double>{2.547321630156819e+1, -2.639500283021502e+1},
    std::complex<double>{2.394538338734709e+1, -5.650522971778156e+0}};

const std::array<std::complex<double>, 8> DepletionMatrix::cram16_theta_{
    std::complex<double>{+3.509103608414918, 8.436198985884374},
    std::complex<double>{+5.948152268951177, 3.587457362018322},
    std::complex<double>{-5.264971343442647, 16.22022147316793},
    std::complex<double>{+1.419375897185666, 10.92536348449672},
    std::complex<double>{+6.416177699099435, 1.194122393370139},
    std::complex<double>{+4.993174737717997, 5.996881713603942},
    std::complex<double>{-1.413928462488886, 13.49772569889275},
    std::complex<double>{-10.84391707869699, 19.27744616718165}};

const double DepletionMatrix::cram16_alpha0_{2.124853710495224e-16};

const std::array<std::complex<double>, 24> DepletionMatrix::cram48_alpha_{
    std::complex<double>{6.387380733878774e+2, -6.743912502859256e+2},
    std::complex<double>{1.909896179065730e+2, -3.973203432721332e+2},
    std::complex<double>{4.236195226571914e+2, -2.041233768918671e+3},
    std::complex<double>{4.645770595258726e+2, -1.652917287299683e+3},
    std::complex<double>{7.765163276752433e+2, -1.783617639907328e+4},
    std::complex<double>{1.907115136768522e+3, -5.887068595142284e+4},
    std::complex<double>{2.909892685603256e+3, -9.953255345514560e+3},
    std::complex<double>{1.944772206620450e+2, -1.427131226068449e+3},
    std::complex<double>{1.382799786972332e+5, -3.256885197214938e+6},
    std::complex<double>{5.628442079602433e+3, -2.924284515884309e+4},
    std::complex<double>{2.151681283794220e+2, -1.121774011188224e+3},
    std::complex<double>{1.324720240514420e+3, -6.370088443140973e+4},
    std::complex<double>{1.617548476343347e+4, -1.008798413156542e+6},
    std::complex<double>{1.112729040439685e+2, -8.837109731680418e+1},
    std::complex<double>{1.074624783191125e+2, -1.457246116408180e+2},
    std::complex<double>{8.835727765158191e+1, -6.388286188419360e+1},
    std::complex<double>{9.354078136054179e+1, -2.195424319460237e+2},
    std::complex<double>{9.418142823531573e+1, -6.719055740098035e+2},
    std::complex<double>{1.040012390717851e+2, -1.693747595553868e+2},
    std::complex<double>{6.861882624343235e+1, -1.177598523430493e+1},
    std::complex<double>{8.766654491283722e+1, -4.596464999363902e+3},
    std::complex<double>{1.056007619389650e+2, -1.738294585524067e+3},
    std::complex<double>{7.738987569039419e+1, -4.311715386228984e+1},
    std::complex<double>{1.041366366475571e+2, -2.777743732451969e+2}};

const std::array<std::complex<double>, 24> DepletionMatrix::cram48_theta_{
    std::complex<double>{-4.465731934165702e+1, 6.233225190695437e+1},
    std::complex<double>{-5.284616241568964e+0, 4.057499381311059e+1},
    std::complex<double>{-8.867715667624458e+0, 4.325515754166724e+1},
    std::complex<double>{+3.493013124279215e+0, 3.281615453173585e+1},
    std::complex<double>{+1.564102508858634e+1, 1.558061616372237e+1},
    std::complex<double>{+1.742097597385893e+1, 1.076629305714420e+1},
    std::complex<double>{-2.834466755180654e+1, 5.492841024648724e+1},
    std::complex<double>{+1.661569367939544e+1, 1.316994930024688e+1},
    std::complex<double>{+8.011836167974721e+0, 2.780232111309410e+1},
    std::complex<double>{-2.056267541998229e+0, 3.794824788914354e+1},
    std::complex<double>{+1.449208170441839e+1, 1.799988210051809e+1},
    std::complex<double>{+1.853807176907916e+1, 5.974332563100539e+0},
    std::complex<double>{+9.932562704505182e+0, 2.532823409972962e+1},
    std::complex<double>{-2.244223871767187e+1, 5.179633600312162e+1},
    std::complex<double>{+8.590014121680897e-1, 3.536456194294350e+1},
    std::complex<double>{-1.286192925744479e+1, 4.600304902833652e+1},
    std::complex<double>{+1.164596909542055e+1, 2.287153304140217e+1},
    std::complex<double>{+1.806076684783089e+1, 8.368200580099821e+0},
    std::complex<double>{+5.870672154659249e+0, 3.029700159040121e+1},
    std::complex<double>{-3.542938819659747e+1, 5.834381701800013e+1},
    std::complex<double>{+1.901323489060250e+1, 1.194282058271408e+0},
    std::complex<double>{+1.885508331552577e+1, 3.583428564427879e+0},
    std::complex<double>{-1.734689708174982e+1, 4.883941101108207e+1},
    std::complex<double>{+1.316284237125190e+1, 2.042951874827759e+1}};

const double DepletionMatrix::cram48_alpha0_{2.258038182743983e-47};

// [1] P. Maria, “Higher-Order Chebyshev Rational Approximation Method and
//     Application to Burnup Equations,” Nucl. Sci. Eng., vol. 182, no. 3,
//     pp. 297–318, 2016, doi: 10.13182/nse15-26.

}  // namespace scarabee