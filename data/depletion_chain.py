import xml.etree.ElementTree as ET
from scarabee._scarabee import (
    NoTarget,
    SingleTarget,
    Branch,
    BranchingTargets,
    FissionYields,
    ChainEntry,
    DepletionChain,
)
import numpy as np

def strip_name(name):
    if '_' in name:
        return name.replace('_', '')
    return name

def read_fission_yields(node):
    if 'parent' in node.attrib:
        #raise RuntimeError("Fission yields refers to parent nuclide {:}.".format(node.attrib['parent']))
        return strip_name(node.attrib['parent'])

    energies = node[0].text.split()
    energies = [float(E) for E in energies]
    NE = len(energies) # Number of energies

    products = None
    yields = []

    for i in range(1, NE+1):
        fiss_yields = node[i]

        prods = fiss_yields[0].text.split()
        prods = [strip_name(n) for n in prods]
        if products is None:
            products = prods
        else:
            # Make sure products are in the same order !
            if len(products) != len(prods):
                raise RuntimeError("Product lists at different energies are of different lengths.")
            for i in range(len(products)):
                if products[i] != prods[i]:
                    raise RuntimeError("Product lists at different energies do not agree.")

        i_yields = fiss_yields[1].text.split()
        yields.append([float(y) for y in i_yields])

    yields = np.array(yields)

    return FissionYields(products, energies, yields)

def main():
    root = ET.parse("chain_casl_pwr.xml").getroot()

    dc = DepletionChain()

    # Parse all nuclides in the chain
    for nuc in root:
        nuclide_name = strip_name(nuc.attrib["name"])
        print("Running {:}".format(nuclide_name))

        half_life = None
        decay_branches_list = []

        n_gamma = None
        n_2n = None
        n_3n = None
        n_p = None
        n_alpha = None
        n_fission = None

        if "half_life" in nuc.attrib:
            half_life = float(nuc.attrib["half_life"])

        for entry in nuc:
            if entry.tag == "decay":
                branch = Branch()
                branch.target = strip_name(entry.attrib["target"])
                branch.branch_ratio = float(entry.attrib["branching_ratio"])
                if branch.branch_ratio > 0.:
                    decay_branches_list.append(branch)

            elif entry.tag == "reaction":
                reaction = entry.attrib["type"]
                transmut_target = NoTarget()
                if "target" in entry.attrib:
                    transmut_target = SingleTarget(strip_name(entry.attrib["target"]))
                
                
                if reaction == '(n,gamma)':
                    n_gamma = transmut_target
                elif reaction == '(n,2n)':
                    n_2n = transmut_target
                elif reaction == '(n,3n)':
                    n_3n = transmut_target
                elif reaction == '(n,p)':
                    n_p = transmut_target
                elif reaction == '(n,a)':
                    n_alpha = transmut_target
                else:
                    # OpenMC considers others like (n,4n), but Scarabée doesn't
                    # track these other transmutation reactions.
                    pass
            
            elif entry.tag == 'neutron_fission_yields':
                n_fission = read_fission_yields(entry)

        # Create decay targets (if there are any)
        decay_targets = None
        if len(decay_branches_list) == 1:
            decay_targets = SingleTarget(decay_branches_list[0].target)
        elif len(decay_branches_list) > 1:
            decay_targets = BranchingTargets(decay_branches_list)

        # Build the chain entry for the nuclide
        ce = ChainEntry()
        ce.half_life = half_life
        ce.decay_targets = decay_targets
        ce.n_gamma = n_gamma
        ce.n_2n = n_2n
        ce.n_3n = n_3n
        ce.n_p = n_p
        ce.n_alpha = n_alpha
        ce.n_fission = n_fission

        dc.insert_entry(nuclide_name, ce)

    dc.save_xml("chain.xml")
    dc.save_json("chain.json")
    dc.save_bin("chain.bin")

if __name__ == "__main__":
    main()