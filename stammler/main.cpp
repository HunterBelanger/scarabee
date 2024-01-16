#include <array>
#include <cmath>
#include <iostream>
#include <vector>

constexpr float PI {3.1415927};

float Ki3(float x) {
  constexpr std::array<float, 45> A {0.7853961, -0.9990226, 0.7266088,
    0.7852024, -0.9912340, 0.6466375, 0.7845986, -0.9791293, 0.5856605,
    0.7834577, -0.9638914, 0.5346648, 0.7817094, -0.9463843, 0.4907827,
    0.7793031, -0.9271152, 0.4521752, 0.7762107, -0.9064822, 0.4177388,
    0.7724519, -0.8849865, 0.3869945, 0.7679903, -0.8626685, 0.3590753,
    0.7628988, -0.8400133, 0.3338676, 0.7540982, -0.8054172, 0.2998569,
    0.7401279, -0.7587821, 0.2609154, 0.7239594, -0.7125290, 0.2278226,
    0.7058777, -0.6672761, 0.1994999, 0.6861762, -0.6234536, 0.1751248};

  constexpr std::array<float, 51> B {0.7247294,   -0.7538355,   0.3203223,
    -5.337485E-2, 0.6663720,   -0.6279752,   0.2295280,   -3.146833E-2, 
    0.5956163,   -0.5094124,   0.1631667,   -1.906198E-2,
    0.5191031,   -0.4046007,   0.1152418,   -1.174752E-2,
    0.4425954,   -0.3159648,   8.097913E-2, -7.328415E-3,
    0.3703178,   -0.2434341,   5.669960E-2, -4.617254E-3,
    0.1684022,   -7.158569E-2, 7.923547E-3,
    0.1278307,   -5.016344E-2, 5.095111E-3,
    9.611422E-2, -3.501524E-2, 3.286040E-3,
    7.170491E-2, -2.437465E-2, 2.126242E-3,
    4.616317E-2, -1.425519E-2, 1.123687E-3,
    2.475115E-2, -6.810124E-3, 4.762937E-4,
    1.302864E-2, -3.232035E-3, 2.031843E-4,
    6.749972E-3, -1.524126E-3, 8.701440E-5,
    3.454768E-3, -7.157367E-4, 3.742673E-5};

  constexpr std::array<int, 20> INDEXA {3, 6, 9, 12, 15, 18, 21, 24, 27, 30,
    33, 33, 36, 36, 39, 39, 42, 42, 45, 45};

  constexpr std::array<int, 20> INDEXB {4, 8, 12, 16, 20, 24, 27, 30, 33, 36,
    39, 39, 42, 42, 45, 45, 48, 48, 51, 51};

  x = std::abs(x);
  int I;
  if (x < 1.) {
    I = static_cast<int>(20. * x);
    I = INDEXA[I]-1;
    return x*(x*A[I]+A[I-1])+A[I-2];
  }

  I = static_cast<int>(2.5*(x-1.));
  if (x < 3.4) {
    I = INDEXB[I]-1;
    return x*(x*(x*B[I]+B[I-1])+B[I-2])+B[I-3];
  }

  if (x < 9.) {
    I = INDEXB[I]-1;
    return x*(x*B[I]+B[I-1])+B[I-2];
  }

  return 0.;
}

void copran(const std::vector<float>& R,
            const std::vector<float>& Et,
            std::vector<std::vector<float>>& P,
            std::vector<float>& GAM,
            std::size_t IG) {
  
  // Constants for calculation
  constexpr std::array<std::array<float, 6>, 6> GJC {
    std::array<float,6>({0.00000000,.55555556,.87393877,.95491150,.98046718,.99029084}),
    std::array<float,6>({2.00000000,.00000000,.28606124,.65127016,.82660307,.90725799}),
    std::array<float,6>({0.72783448,1.27216552,.0000000,.16932809,.47704397,.68412769}),
    std::array<float,6>({0.27930792,.91696442,.80372766,.00000000,.11094751,.35681753}),
    std::array<float,6>({0.12472388,.51939018,.81385828,.54202764,.00000000,.07803490}),
    std::array<float,6>({0.06299166,.29563548,.58554794,.66869856,.38712636,.00000000})
  };

  const std::size_t N = R.size(); // Number of regions
  std::vector<float> R2(R.size(), 0.);
  std::vector<float> SIGV(R.size(), 0.);
  std::vector<float> TAU(R.size(), 0.);

  // Prepare
  float Y = 0.;
  for (std::size_t i = 0; i < N; i++) {
    R2.at(i) = R.at(i) * R.at(i);
    SIGV[i] = PI * Et.at(i) * (R2.at(i) - Y);
    Y = R2.at(i);

    for (std::size_t j = 0; j < N; j++) {
      P.at(i).at(j) = 0.;
    }
  }

  // Fill P[i][j] with the integrals of (Ki3(T+) - Ki3(T-))
  for (std::size_t i = 0; i < N; i++) {
    float rstart = 0.;
    if (i > 0) rstart = R[i-1];
    float dr = R[i] - rstart;

    // Loop over Gauss-Jacobi points
    for (std::size_t k = 0; k < IG; k++) {
      Y = rstart + dr*GJC[k][IG];
      float fac = dr * GJC[IG][k];
      if (i > 0) TAU[i-1] = 0.;
      float TPLUS = 0.;
      float Y2 = Y*Y;
      for (std::size_t j = i; j < N; j++) {
        float TMINUS = TPLUS;
        TPLUS = std::sqrt(R2[j] - Y2);

        if (j > 0) TAU[j] = TAU[j-1] + Et[j]*(TPLUS - TMINUS);
        else TAU[j] = Et[j]*(TPLUS - TMINUS);

        for (std::size_t l = i; l <= j; l++) {
          P[l][j] += fac*(Ki3(TAU[j]+TAU[l]) - Ki3(TAU[j]-TAU[l]));
        }
      }
    }
  }

  // Compose P[i][j] for black boundaries
  for (std::size_t i = 0; i < N; i++) {
    std::size_t j = N - i - 1;

    for (std::size_t k = j; k < N; k++) {
      std::size_t l = j + N - k - 1;

      if (l == 0) goto l30;
      
      if (l == j) P[j][l-1] = P[l-1][j];
      P[j][l] -= P[j][l-1];

      if (j == 0) goto l30;

      P[j][l] = P[j][l] - P[j-1][l] + P[j-1][l-1];

      l30: P[l][j] = P[j][l];
    }

    P[j][j] += SIGV[j];
  }

  // First flight blackness
  const float S = 2. / (PI * R.back());
  for (std::size_t i = 0; i < N; i++) {
    float sum = SIGV[i];
    for (std::size_t j = 0; j < N; j++) {
      sum -= P[i][j];
    }
    GAM[i] = S * sum;
  }
}

void gauscp(const std::vector<float>& SEC, // Secondaries Es / E
            const std::vector<float>& ST1, // Inverse sigma total
            const std::vector<float>& VSI, // Vol * Sigma
            std::vector<std::vector<float>>& P, // Symmetric col probs, then matrix
            std::vector<float>& Y, // At first, partial blackness, then solution
            std::vector<std::vector<float>>& X) { // At first RHS, then solution

  const std::size_t N = SEC.size();

  // Set up all RHSs
  for (std::size_t i = 0; i < N; i++) {
    for (std::size_t j = 0; j < N; j++) {
      X[i][j] = ST1[j] * P[j][i];
      P[j][i] = -SEC[i] * P[j][i];
    }
    P[i][i] += VSI[i];
  }

  for (std::size_t i = 0; i < P.size(); i++) {
    for (std::size_t j = 0; j < P.size(); j++) {
      std::cout << P[i][j];
      if (j != P.size()-1) std::cout << ", ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  for (std::size_t i = 0; i < Y.size(); i++) {
    std::cout << Y[i] << "\n";
  }
  std::cout << "\n";

  // Sweep lower diagonal
  for (std::size_t i = 0; i < N; i++) {
    float F = 1. / P[i][i];
    Y[i] *= F;

    for (std::size_t j = 0; j < N; j++) {
      X[i][j] *= F;
    }

    if (i == N-1) goto l70;

    for (std::size_t j = i+1; j < N; j++) {
      P[i][j] *= F;
    }

    // Do actual sweeping
    for (std::size_t j = i+1; j < N; j++) {
      Y[j] -= Y[j]*P[j][i];

      for (std::size_t k = 0; k < N; k++) {
        X[j][k] -= X[i][k]*P[j][i];
      }

      for (std::size_t k = i+1; k < N; k++) {
        P[j][k] -= P[i][k] * P[j][i];
      }
    }
  }

  // Solve. The last RHS is already the last component
  l70: std::size_t IMAX = N -1;
  if (IMAX == 0) return;

  for (std::size_t i = 0; i < IMAX; i++) {
    std::size_t j = N - i - 1;
    std::size_t j1 = j + 1;

    for (std::size_t k = j1; k < N; k++) {
      Y[j] -= Y[k]*P[j][k];

      for (std::size_t kk = 0; kk < N; kk++) {
        X[j][kk] -= X[k][kk]*P[j][k];
      }
    }
  }
}

int main() {
  /*
  std::vector<float> R {0.5, 0.61, 0.8, 4.55};
  std::vector<float> Et {2., 0.5, 1.5, 1.5};
  std::vector<std::vector<float>> P(4, std::vector<float>(4, 0.));
  std::vector<float> Gam(4, 0.);
  */

  std::size_t N = 4;
  std::vector<float> R {0.27, 0.54, 0.625439, 0.710879};
  std::vector<float> Et {0.177949, 0.177949, 0.159206, 0.159206};
  std::vector<float> Es {1.27537E-01, 1.27537E-01, 4.44777E-02, 4.44777E-02};
  std::vector<std::vector<float>> P(N, std::vector<float>(N, 0.));
  std::vector<float> gam(N, 0.);

  std::vector<float> ST1 = Et; // Inverse Sigma Total
  for (auto& v : ST1) v = 1. / v;

  std::vector<float> VSI(R.size(), 0.); // Vol * Sigma
  for (std::size_t i = 0; i < N; i++) {
    if (i == 0) {
      VSI[i] = Et[i] * PI * R[i] * R[i];
    } else {
      VSI[i] = Et[i] * PI * (R[i]*R[i] - R[i-1]*R[i-1]);
    }
  }

  std::vector<float> SEC = Es;
  for (std::size_t i = 0; i < N; i++) SEC[i] /= Et[i];

  std::vector<std::vector<float>> X(N, std::vector<float>(N, 0.));

  copran(R, Et, P, gam, 5);

  for (std::size_t i = 0; i < P.size(); i++) {
    for (std::size_t j = 0; j < P.size(); j++) {
      std::cout << P[i][j];
      if (j != P.size()-1) std::cout << ", ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  gauscp(SEC, ST1, VSI, P, gam, X);

  /*
  for (std::size_t i = 0; i < P.size(); i++) {
    for (std::size_t j = 0; j < P.size(); j++) {
      std::cout << P[i][j];
      if (j != P.size()-1) std::cout << ", ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  */

  /*
  // Print X
  for (std::size_t i = 0; i < X.size(); i++) {
    for (std::size_t j = 0; j < X.size(); j++) {
      std::cout << X[i][j];
      if (j != X.size()-1) std::cout << ", ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  // Print Y
  for (std::size_t j = 0; j < gam.size(); j++) {
    std::cout << gam[j];
    if (j != gam.size()-1) std::cout << ", ";
  }
  std::cout << "\n";
  */

  return 0;
}