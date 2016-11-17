import mcsmrff_opt

LJ = [[2.1, -1.05], [3.98892, 4.15625], [0.1, 0.1]]
three_body_interactions = ['Pb', 'Pb', 'Pb',
                           'Pb', 'Pb', 'Cl',
                           'Pb', 'Pb', 'H3',
                           'Pb', 'Cl', 'Pb',
                           'Pb', 'Cl', 'Cl',
                           'Pb', 'Cl', 'H3',
                           'Pb', 'H3', 'Pb',
                           'Pb', 'H3', 'Cl',
                           'Pb', 'H3', 'H3',
                           'Cl', 'Pb', 'Pb',
                           'Cl', 'Pb', 'Cl',
                           'Cl', 'Pb', 'H3',
                           'Cl', 'Cl', 'Pb',
                           'Cl', 'Cl', 'Cl',
                           'Cl', 'Cl', 'H3',
                           'Cl', 'H3', 'Pb',
                           'Cl', 'H3', 'Cl',
                           'Cl', 'H3', 'H3',
                           'H3', 'Pb', 'Pb',
                           'H3', 'Pb', 'Cl',
                           'H3', 'Pb', 'H3',
                           'H3', 'Cl', 'Pb',
                           'H3', 'Cl', 'Cl',
                           'H3', 'Cl', 'H3',
                           'H3', 'H3', 'Pb',
                           'H3', 'H3', 'Cl',
                           'H3', 'H3', 'H3']

mcsmrff_opt.run_mcsmrff_optimizer(
    n_sets=5,
    sim_name="test",
    lennard_jones=LJ,
    atom_list=three_body_interactions,
    training_set_pickle_path=(
        "/fs/home/hch54/"
        "Grad-MCSMRFF/PbCl3Cs/set2/"
        "set2.pickle")
)
