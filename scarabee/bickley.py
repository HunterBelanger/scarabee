import numpy as np
from scipy.integrate import quad
from numba import njit

@njit
def _Ki3_integrand(theta: float, x: float) -> float:
  cos = np.cos(theta)
  cos_sqrd = cos * cos
  return cos_sqrd * np.exp(-x / cos)

def Ki3(x: float) -> float:
  return quad(_Ki3_integrand, 0., 0.5*np.pi, args=(x), epsabs=0)[0]