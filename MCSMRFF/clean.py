import sys, os
for directory in os.listdir("training_sets"):
	for file in os.listdir("training_sets/%s"%directory):
		if ".out" in file or ".cml" in file or ".engrad" in file: continue
		os.system("rm training_sets/%s/%s" % (directory,file))

