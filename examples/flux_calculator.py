import numpy as np
import matplotlib.pyplot as plt
import pyPapillonNDL as pndl
from scarabee import FluxCalculator

def main():
    # Load the continuous energy data for U238
    #ace = pndl.ACE("C:\\Users\\hunte\\Documents\\nuclear_data\\ace\\endf8.0\\abeille_binary\\data\\neutron_dir\\U\\U238\\U238.293.6.ace.bin", pndl.ACEType.BINARY)
    ace = pndl.ACE("/mnt/c/Users/hunte/Documents/nuclear_data/ace/endf8.0/abeille_binary/data/neutron_dir/U/U238/U238.293.6.ace.bin", pndl.ACEType.BINARY)
    Nuc = pndl.STNeutron(ace)

    # Make an array for the ultra-fine slowing down calculation
    energy_bounds = np.geomspace(1.E-5, 20.E6, 500000)
    sig_t = np.zeros(energy_bounds.size-1)
    sig_s = np.zeros(energy_bounds.size-1)
    for g in range(sig_t.size):
        E = 0.5*(energy_bounds[g+1] + energy_bounds[g]) * 1.E-6
        sig_t[g] = Nuc.total_xs()(E)
        sig_s[g] = Nuc.elastic_xs()(E)
    
    # Atomic Weight Ratio of H1
    AWR_H1 = 0.9991673

    # Perform slowing down calculation for U238 with a 50 barn background xs from H1
    calc = FluxCalculator(energy_bounds, sig_t, sig_s, Nuc.awr())
    calc.add_background_nuclide(50., AWR_H1)
    calc.solve()

    plt.plot(calc.avg_energy, calc.flux)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [eV]')
    plt.ylabel('Flux [Arb. Units]')
    plt.tight_layout()
    plt.show()
    
if __name__ == '__main__':
    main()
