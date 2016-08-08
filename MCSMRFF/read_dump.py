import os, sys
import numpy as np

def read_dump(fptr, data=["element","x","y","z"]):
	ids = None
	atoms_flag = False

	frames, frame = [], []
	for line in open(fptr,'r'):
		if line.startswith("ITEM: ATOMS"):
			atoms_flag = True
			frame = []
			if ids is None:
				ids = {name:index for index,name in enumerate(line.split()[2:])}
				continue
		elif line.startswith("ITEM:"):
			atoms_flag = False
			if frame != []: frames.append(frame)
			frame = []
			continue
		elif atoms_flag:
			# Read in the data
			atom = []
			values = line.split()
			for d in data:
				try:
					atom.append( np.float64(values[ids[d]]) )
				except:
					atom.append( values[ids[d]] )
					print values
			frame.append(atom)

	return frames

def frames_to_xyz(filename, frames):
	fptr = open(filename,'w')
	for frame in frames:
		fptr.write("%d\nAtoms\n" % len(frame))
		for atom in frame:
			fptr.write("%s\t%lg\t%lg\t%lg\n" % tuple(atom))
	fptr.close()

fptr = sys.argv[1]
frames = read_dump(fptr, data=["type","xu","yu","zu"])

swap = {1:"Pb",2:"Cl",3:"N",4:"C",5:"H",6:"H"}

for i,frame in enumerate(frames):
	for j,atom in enumerate(frame):
		frames[i][j][0] = swap[atom[0]]

frames_to_xyz("out.xyz", frames)
