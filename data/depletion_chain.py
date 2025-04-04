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
import sys


def strip_name(name):
    if "_" in name:
        return name.replace("_", "")
    return name

def read_fission_yields(node):
    if "parent" in node.attrib:
        # raise RuntimeError("Fission yields refers to parent nuclide {:}.".format(node.attrib['parent']))
        return strip_name(node.attrib["parent"])

    energies = node[0].text.split()
    energies = [float(E) for E in energies]
    NE = len(energies)  # Number of energies

    products = None
    yields = []

    for i in range(1, NE + 1):
        fiss_yields = node[i]

        prods = fiss_yields[0].text.split()
        prods = [strip_name(n) for n in prods]
        if products is None:
            products = prods
        else:
            # Make sure products are in the same order !
            if len(products) != len(prods):
                raise RuntimeError(
                    "Product lists at different energies are of different lengths."
                )
            for i in range(len(products)):
                if products[i] != prods[i]:
                    raise RuntimeError(
                        "Product lists at different energies do not agree."
                    )

        i_yields = fiss_yields[1].text.split()
        yields.append([float(y) for y in i_yields])

    yields = np.array(yields)

    return FissionYields(products, energies, yields)

def main():
    if len(sys.argv) != 2:
        raise RuntimeError("Script takes exactly 1 argument, the file name of an xml depletion chain.")

    root = ET.parse(sys.argv[1]).getroot()

    dc = DepletionChain()

    # First, go through all nuclides and obtain all fission yields which don't
    # just refer to a different nuclide.
    fiss_yields = {}
    for nuc in root:
        if nuc[-1].tag == "neutron_fission_yields" and "parent" not in nuc[-1].attrib:
            nuclide_name = strip_name(nuc.attrib["name"])
            fiss_yields[nuclide_name] = read_fission_yields(nuc[-1])

    # Parse all nuclides in the chain
    for nuc in root:
        nuclide_name = strip_name(nuc.attrib["name"])

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
                if branch.branch_ratio > 0.0:
                    decay_branches_list.append(branch)

            elif entry.tag == "reaction":
                reaction = entry.attrib["type"]
                transmut_target = None
                branch_ratio = None
                if "target" in entry.attrib:
                    transmut_target = strip_name(entry.attrib["target"])
                if "branching_ratio" in entry.attrib:
                    branch_ratio = float(entry.attrib["branching_ratio"])

                if reaction == "(n,gamma)":
                    if branch_ratio is None and transmut_target is None:
                        n_gamma = NoTarget()
                    elif branch_ratio is None:
                        n_gamma = SingleTarget(transmut_target)
                    else:
                        if n_gamma is None:
                            n_gamma = []
                        branch = Branch()
                        branch.target = transmut_target
                        branch.branch_ratio = branch_ratio
                        n_gamma.append(branch)

                elif reaction == "(n,2n)":
                    if branch_ratio is None and transmut_target is None:
                        n_2n = NoTarget()
                    elif branch_ratio is None:
                        n_2n = SingleTarget(transmut_target)
                    else:
                        if n_2n is None:
                            n_2n = []
                        branch = Branch()
                        branch.target = transmut_target
                        branch.branch_ratio = branch_ratio
                        n_2n.append(branch)

                elif reaction == "(n,3n)":
                    if branch_ratio is None and transmut_target is None:
                        n_3n = NoTarget()
                    elif branch_ratio is None:
                        n_3n = SingleTarget(transmut_target)
                    else:
                        if n_3n is None:
                            n_3n = []
                        branch = Branch()
                        branch.target = transmut_target
                        branch.branch_ratio = branch_ratio
                        n_3n.append(branch)

                elif reaction == "(n,p)":
                    if branch_ratio is None and transmut_target is None:
                        n_p = NoTarget()
                    elif branch_ratio is None:
                        n_p = SingleTarget(transmut_target)
                    else:
                        if n_p is None:
                            n_p = []
                        branch = Branch()
                        branch.target = transmut_target
                        branch.branch_ratio = branch_ratio
                        n_p.append(branch)                   

                elif reaction == "(n,a)":
                    if branch_ratio is None and transmut_target is None:
                        n_alpha = NoTarget()
                    elif branch_ratio is None:
                        n_alpha = SingleTarget(transmut_target)
                    else:
                        if n_alpha is None:
                            n_alpha = []
                        branch = Branch()
                        branch.target = transmut_target
                        branch.branch_ratio = branch_ratio
                        n_alpha.append(branch)                   
                else:
                    # OpenMC considers others like (n,4n), but Scarabée doesn't
                    # track these other transmutation reactions.
                    pass

            elif entry.tag == "neutron_fission_yields":
                parent = nuclide_name
                if 'parent' in entry.attrib:
                    parent = strip_name(entry.attrib['parent'])
                n_fission = fiss_yields[parent]

        # Create decay targets (if there are any)
        decay_targets = None
        if len(decay_branches_list) == 1:
            decay_targets = SingleTarget(decay_branches_list[0].target)
        elif len(decay_branches_list) > 1:
            decay_targets = BranchingTargets(decay_branches_list)

        if isinstance(n_gamma, list):
            n_gamma = BranchingTargets(n_gamma)
        if isinstance(n_2n, list):
            n_2n = BranchingTargets(n_2n)
        if isinstance(n_3n, list):
            n_3n = BranchingTargets(n_3n)
        if isinstance(n_p, list):
            n_p = BranchingTargets(n_p)
        if isinstance(n_alpha, list):
            n_alpha = BranchingTargets(n_alpha)

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

    # Expunge all short-lives nuclides (i.e. half life less than 24 hrs) that
    # aren't in important chains.
    nuclide_names = dc.nuclides
    for nuc_name in nuclide_names:
        if nuc_name in ["I135", "Xe135"]:
            continue

        nuc = dc.nuclide_data(nuc_name)

        if nuc.half_life is not None and nuc.half_life < 60.*60.*24.:
            print("Removing {:} from chain".format(nuc_name))
            dc.remove_nuclide(nuc_name)
        elif nuc_name in ["Cd115", "Rh102", "Rh102m1", "Sb127"]: # No Evals in ENDF-8.0 for these. RIP
            print("Removing {:} from chain".format(nuc_name))
            dc.remove_nuclide(nuc_name)

    dc.save("chain.bin")

if __name__ == "__main__":
    main()
