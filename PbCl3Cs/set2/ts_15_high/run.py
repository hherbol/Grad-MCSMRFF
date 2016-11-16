import mcsmrff_train

mcsmrff_train.create_training_set("set2",
                                  N_perterbations_per_seed=5,
                                  perterbation_dx=0.3,
                                  perterbation_dr=10,
                                  expansion_perterbation_dx=0.1,
                                  expansion_perterbation_dr=0.5,
                                  expansion_step=0.5,
                                  N_expansion=5,
                                  N_perterbations_per_expansion=5,
                                  min_seeds=False,
                                  training_sets_folder="training_set",
                                  queue=None,
                                  procs=2,
                                  persist=False)
