import numpy as np
from typing import Optional

class TransportCrossSections:
  def __init__(self):
    self._Et  = np.array([])
    self._Ea  = np.array([])
    self._Ef  = np.array([])
    self._Es  = np.array([])
    self._nu  = np.array([])
    self._chi = np.array([])
    self._fissile = False

  @property
  def ngroups(self):
    return len(self._Et)

  @property
  def fissile(self):
    return self._fissile

  def Et(self, g: int):
    return self._Et[g]

  def Ea(self, g: int):
    return self._Ea[g]

  def Ef(self, g: int):
    return self._Ef[g]

  def Er(self, g: int):
    return self.Ea(g) + self.Es(g) - self.Es(g,g)

  def Es(self, gin: int, gout: Optional[int] = None):
    if gout is not None:
      return self._Es[gin][gout]
    else:
      return np.sum(self._Es[gin,:])