import mcsmrff_opt
import itertools

LJ = [[2.0, -1.0, 1.0], [3.98892, 4.15625, 4.0], [0.1, 0.1, 0.1]]

Pb = 907
Cl = 344
Cs = 352
elems = ["Pb", "Cl", "Cs"]
tersoff_atoms = [Pb, Cl, Cs]

three_body_interactions = list(itertools.product(elems, repeat=3))

mcsmrff_opt.run_mcsmrff_optimizer(
    n_sets=5,
    sim_name="test",
    lennard_jones=LJ,
    atom_list=three_body_interactions,
    tersoff=None,
    lj_coul=None,
    constant_charge=True,
    training_set_pickle_path=(
        "/fs/home/hch54/"
        "Grad-MCSMRFF/PbCl3Cs/set2/"
        "set2.pickle"),
    tersoff_atoms=tersoff_atoms
)
