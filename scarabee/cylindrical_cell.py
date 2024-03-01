import numpy as np
from numba import njit
from scipy.integrate import quad
from scipy.linalg import solve as mat_solve
from typing import Iterable

from scarabee.transport_cross_sections import TransportCrossSections
from scarabee.bickley import Ki3

class CylindricalCell:
  def __init__(self, radii: Iterable[float], mats: Iterable[TransportCrossSections]):
    self._radii = np.array(radii)
    self._mats = np.array(mats)

    # Do logic checks on input
    if len(self._radii) < 2:
      raise RuntimeError("Number of radii must be >= 2.")
    
    if np.all(self._radii[:-1] <= self._radii[1:]) is False:
      raise RuntimeError("Radii must be sorted.")

    if self._radii[0] <= 0.:
      raise RuntimeError("All radii must be > 0.")

    if self._radii.shape != self._mats.shape:
      raise RuntimeError("Must have the same number of radii and materials.")

    ngrps = self._mats[0].ngroups
    for i in range(1, self.nregions):
      if self._mats[i].ngroups != ngrps:
        raise RuntimeError("Materials must have the same number of groups.")

    # Calculate the volume of each rin
    self._vols = np.zeros(self._radii.shape)
    for i in range(self.nregions):
      self._vols[i] = np.pi * self._radii[i] * self._radii[i]
      if i != 0:
        self._vols[i] -= np.pi * self._radii[i-1] * self._radii[i-1]

    # Calculate all the collision probabilities
    self._p = np.zeros((self.ngroups, self.nregions, self.nregions))
    self._calc_collision_probabilities()

    # Solve system    
    self._X = np.zeros((self.ngroups, self.nregions, self.nregions))
    self._Y = np.zeros((self.ngroups, self.nregions))
    self._Gamma = np.zeros(self.ngroups)
    self._solve_system()

  @property
  def nregions(self) -> int:
    return len(self._radii)

  @property
  def ngroups(self) -> int:
    return self._mats[0].ngroups

  @property
  def Sb(self) -> float:
    return 2. * np.pi * self._radii[-1]

  def V(self, i: int) -> float:
    return self._vols[i]

  def R(self, i: int) -> float:
    return self._radii[i]

  def mat(self, i: int) -> TransportCrossSections:
    return self._mats[i]

  def x(self, g: int, i: int) -> float:
    return 0.25 * self.Sb() * self.V(i) * self._Y[g, i]

  def Y(self, a: float, g: int, i: int) -> float:
    Yi = self._Y[g, i]
    Gamma = self._Gamma[g]
    return Yi / (1. - a * (1. - Gamma))

  def X(self, a: float, g: int, i: int, k: int) -> float:
    xk = self.x(g, k)
    Xik = self._X[g, i, k]
    return Xik + a * xk * self.Y(a, g, i)

  def Gamma(self, g: int) -> float:
    return self._Gamma[g]

  def p(self, g: int, i: int, j: int) -> float:
    return self._p[g, i, j]

  def _calc_collision_probabilities(self):
    S = np.zeros((self.nregions, self.nregions))

    # Calculate the matrix for each energy group
    for g in range(self.ngroups):
      # Load S
      for i in range(self.nregions):
        for j in range(i, self.nregions):
          S[i, j] = self._calc_S(i, j, g)

          # These are symmetric, so we can assign the diagonal term too
          if j != i:
            S[j, i] = S[i, j]

      # S has now been filled. We can now load p
      for i in range(self.nregions):
        for j in range(i, self.nregions):
          self._p[g, i, j] = 0.
          if i == j:
            self._p[g, i, j] += self._vols[i] * self._mats[i].Et(g)

          self._p[g, i, j] += 2. * S[i, j]
          if i > 0 and j > 0:
            self._p[g, i, j] += 2. * S[i - 1, j - 1]
          if i > 0:
            self._p[g, i, j] -= 2. * S[i - 1, j]
          if j > 0:
            self._p[g, i, j] -= 2. * S[i, j - 1]

          if i != j:
            self._p[g, j, i] = self._p[g, i, j]

      S.fill(0.)

  def _calc_S(self, i: int, j: int, g: int):
    if i > j:
      tmp = j
      j = i
      i = tmp
    
    # We must solve S_{i,j} = int_{0}^{R_i} (Ki3(Tmax(y)) - Ki3(Tmin(y))) dy.
    # This integral is performed by doing the integral out to each annular ring
    # and adding it to a sum.

    S_ij = 0.

    for k in range(i+1):
      # We do the integral from radii_[k-1] to radii_[k]. If k == 0, then we
      # start the integral from the center of the cylinder.
      # Get min and max radii for determining y points
      Rmin = 0 if k == 0 else self._radii[k-1]
      Rmax = self._radii[k]

      # Initialize an array for the x coordinates
      x = np.zeros((j+1-k))

      # Initialize an array to hold Et for each material in group g
      Ets = np.zeros(self.nregions)
      for i in range(self.nregions):
        Ets[i] = self._mats[i].Et(g)

      # We now integrate from Rmin to Rmax
      integral = quad(_cyl_col_prob_integrand, Rmin, Rmax, args=(x, k, i, j, self._radii, Ets))

      # TODO check integral

      S_ij += integral[0]

    return S_ij

  def _solve_system(self):
    # The Y and X terms all depend on the energy group, and are described by
    # the same matrix. We therefore load the matrix once, and solve it for
    # all equations.
    for g in range(self.ngroups):
      # Initialize and fill matrix for group g
      M = np.zeros((self.nregions, self.nregions))
      for i in range(self.nregions):
        for j in range(self.nregions):
          Et = self._mats[j].Et(g)
          Es_gg = self._mats[j].Es(g, g)
          c_j = Es_gg / Et

          M[i, j] = -c_j * self._p[g, j, i]

          if j == i:
            M[i, j] += Et * self._vols[i]

      # There are nregions systems for solve for X, and 1 for Y
      b = np.zeros((self.nregions, self.nregions+1))

      # Here, we fill b for all nregion X systems
      for k in range(self.nregions):
        # Load the b vector
        Et_k = self._mats[k].Et(g)
        for i in range(self.nregions):
          b[i, k] = self._p[g, k, i] / Et_k

      # Here we fill the last index of b for the Y system
      coeff = 4. / self.Sb
      for i in range(self.nregions):
        Et_i = self._mats[i].Et(g)
        Vol_i = self._vols[i]
        sum_p = 0.
        for j in range(self.nregions):
          sum_p += self._p[g, i, j]
        b[i, self.nregions] = coeff * (Et_i * Vol_i - sum_p)

      # Solve for this set of X and Y
      mat_sol = mat_solve(M, b)
      self._X[g,:,:] = mat_sol[:,:self.nregions]
      self._Y[g,:] = mat_sol[:,-1]

      # Calculate multicollision blackness for this group
      self._Gamma[g] = 0.
      for i in range(self.nregions):
        self._Gamma[g] += self._mats[i].Er(g) * self._vols[i] * self._Y[g, i]

def _cyl_col_prob_integrand(y, x, k, i, j, radii, Ets):
  y2 = y * y

  # Get the x coordinates for the given y
  for n in range(len(x)):
    x[n] = np.sqrt(radii[k + n] * radii[k + n] - y2)

  # Calculate tau_pls and tau_min by iterating through all segments
  tau_pls = 0.
  tau_min = 0.
  for s in range(len(x)):
    # Get the length of the segment
    t = x[s] if s == 0 else x[s] - x[s - 1]

    # Get optical depth contribution
    dtau = t * Ets[s+k]

    # Add to the tau variables
    if s + k <= i:
      tau_pls += 2. * dtau
    else:
      tau_pls += dtau
      tau_min += dtau

  if i == j:
    tau_min = 0.0

  return Ki3(tau_pls) - Ki3(tau_min)