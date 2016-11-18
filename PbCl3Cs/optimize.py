import mcsmrff_opt

LJ = [[2.1, -1.05, 1.0], [3.98892, 4.15625], [0.1, 0.1]]
three_body_interactions = ['Pb', 'Pb', 'Pb',
                           'Pb', 'Pb', 'Cl',
                           'Pb', 'Pb', 'Cs',
                           'Pb', 'Cl', 'Pb',
                           'Pb', 'Cl', 'Cl',
                           'Pb', 'Cl', 'Cs',
                           'Pb', 'Cs', 'Pb',
                           'Pb', 'Cs', 'Cl',
                           'Pb', 'Cs', 'Cs',
                           'Cl', 'Pb', 'Pb',
                           'Cl', 'Pb', 'Cl',
                           'Cl', 'Pb', 'Cs',
                           'Cl', 'Cl', 'Pb',
                           'Cl', 'Cl', 'Cl',
                           'Cl', 'Cl', 'Cs',
                           'Cl', 'Cs', 'Pb',
                           'Cl', 'Cs', 'Cl',
                           'Cl', 'Cs', 'Cs',
                           'Cs', 'Pb', 'Pb',
                           'Cs', 'Pb', 'Cl',
                           'Cs', 'Pb', 'Cs',
                           'Cs', 'Cl', 'Pb',
                           'Cs', 'Cl', 'Cl',
                           'Cs', 'Cl', 'Cs',
                           'Cs', 'Cs', 'Pb',
                           'Cs', 'Cs', 'Cl',
                           'Cs', 'Cs', 'Cs']

Pb_, Pb = 111, 907
Cl_, Cl = 21, 344
Cs_, Cs = 72, 352
tersoff_atoms = [Pb, Cl, Cs]

mcsmrff_opt.run_mcsmrff_optimizer(
    n_sets=2,
    sim_name="test",
    lennard_jones=LJ,
    atom_list=three_body_interactions,
    training_set_pickle_path=(
        "/fs/home/hch54/"
        "Grad-MCSMRFF/PbCl3Cs/set2/"
        "set2.pickle")
)
