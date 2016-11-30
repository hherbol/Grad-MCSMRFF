import structures

# This holds indices to OPLS parameters for corresponding
# Underscored variables depict the bonding indice while
# no underscore indicates the atom type index
I_, I = 66, 838
Cl_, Cl = 21, 344
Pb_, Pb = 111, 907
Cs_, Cs = 72, 352
Si_, Si = 112, 908
O_, O = 113, 909
H_ = 54
N_ = 53
HN = 233

tersoff_atoms = [Pb, Cl, Cs]

extra_Pb = {
    Pb: structures.Struct(index=Pb,
                          index2=Pb_,
                          element_name='Pb',
                          element=82,
                          mass=207.2,
                          charge=0.4,
                          vdw_e=10.1,
                          vdw_r=3.0),
}

extra_Si = {
    Si: structures.Struct(index=Si,
                          index2=Si_,
                          element_name='Si',
                          element=14,
                          mass=28.0855,
                          charge=4,
                          vdw_e=10.1,
                          vdw_r=3.3),
}

extra_O = {
    O: structures.Struct(index=O,
                         index2=O_,
                         element_name='O',
                         element=8,
                         mass=15.999,
                         charge=-2,
                         vdw_e=10.1,
                         vdw_r=1.52),
}

opls_path = '/fs/europa/g_pc/Forcefields/OPLS/oplsaa.prm'
lammps_mcsmrff = '/fs/home/hch54/lammps/lammps-7Dec15/src/lmp_serial'
