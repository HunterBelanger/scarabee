from enum import Enum

class CriticalLeakage(Enum):
  """
  Defines a critical leakage model which can be applied to an assembly.
  """

  P1 = 1
  """
  Homogeneous P1 approximation.
  """

  B1 = 2
  """
  Homogeneous B1 approximation.
  """

  NoLeakage = 3
  """
  No leakage (assembly calculation is not modified).
  """
