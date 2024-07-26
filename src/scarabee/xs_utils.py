import numpy as np
from _scarabee import DiffusionCrossSection

# Utility functions for saving/loading diffusion xs

def save_diffusion_xs(fname, xs):
    """
    Saves a set of diffusion cross sections to a numpy file.

    Parameters
    ----------
    fname : str
        Name of file in which to save data.
    xs : DiffusionCrossSection
        The diffusion cross sections to save.
    """
    # Make giant numpy array
    NG = xs.ngroups

    data = np.zeros((5+NG, NG))
    
    # Save diffusion coefficients
    for g in range(NG):
        data[0, g] = xs.D(g)

    # Save Ea
    for g in range(NG):
        data[1, g] = xs.Ea(g)

    # Save Ef
    for g in range(NG):
        data[2, g] = xs.Ef(g)

    # Save vEf
    for g in range(NG):
        data[3, g] = xs.vEf(g)

    # Save chi
    for g in range(NG):
        data[4, g] = xs.chi(g)

    # Save Es
    for gin in range(NG):
        for gout in range(NG):
            data[5+gin, gout] = xs.Es(gin, gout)

    np.save(fname, data)

def load_diffusion_xs(fname):
    """
    Loads a set of diffusion cross sections from a numpy file.

    Parameters
    ----------
    fname : str
        Name of file containing the saved data.

    Returns
    -------
    DiffusionCrossSection
        Diffusion cross sections from the file.
    """

    data = np.load(fname)

    D   = data[0,:]
    Ea  = data[1,:]
    Ef  = data[2,:]
    vEf = data[3,:]
    chi = data[4,:]
    Es = data[5:,:]

    return DiffusionCrossSection(D, Ea, Es, Ef, vEf, chi)
