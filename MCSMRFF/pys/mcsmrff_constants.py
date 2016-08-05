import utils
# This holds indices to OPLS parameters for corresponding
# Underscored variables depict the bonding indice while no underscore indicates the atom type index
I_, I = 66, 838
Cl_, Cl = 21, 344
Pb_, Pb = 111, 907
H_ = 54
N_ = 53
HN = 233

extra_Pb = {
		Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.4, vdw_e=10.1, vdw_r=3.0),
	}

opls_path = '/fs/europa/g_pc/Forcefields/OPLS/oplsaa.prm'
lammps_mcsmrff = '/fs/home/hch54/lammps/lammps-7Dec15/src/lmp_serial'