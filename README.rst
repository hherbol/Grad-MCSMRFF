A LAMMPS fix for parameterizing Tersoff potentials. Optional Coulombic and Lennard-Jones contributions are available for the range outside the Tersoff cutoffs. 

# Getting started

1.	Download [LAMMPS](http://lammps.sandia.gov/download.html) ([older versions here](http://lammps.sandia.gov/tars/)). Some versions will not work - the 7 Dec 2015 is the one we've tested the most, so use that one. 
	
2.	Open a terminal in the src directory in your lammps folder

3.	In your LAMMPS src directory run `make yes-manybody`
	
4.	If you don't have an SSH key, generate one like this:

		ssh-keygen -t rsa -b 4096 -C "your_email@example.com"  #Creates an ssh key, using your GitHub e-mail as a label
		
	When prompted to "Enter a file in which to save the key," press Enter
	If asked to Overwrite, enter 'y'
	At the prompt, type a secure passphrase
	Retype your secure passphrase

5.	Adding a new SSH key to your GitHub account

		gedit ~/.ssh/id_rsa.pub
	
	In the top right corner of any GitHub page in your browser, click on your profile photo, the click 'Settings'
	In the user settings sidebar, click 'SSH keys'
	Click 'New SSH key'
	In the "Title" field, add a descriptive label for the new key
	Copy and Paste the contents from the 'id_rsa.pub' file into the "Key" field
	Click 'Add SSH key'

6.	Load your keys into your SSH agent
	
		eval "$(ssh-agent -s)"
		
		ssh-add
		
	Enter passphrase
	
7.	Test your SSH connection

		ssh -T git@github.com
		
	You should see the message "Hi 'username'! You've successfully authenicated, but GitHub does not provide shell access."

8.	Clone repository wherever
	
		git clone git@github.com:hherbol/Grad-MCSMRFF.git

9.  Copy over contents of Grad-MCSMRFF/LAMMPS to your lammps src folder.  Note, this overwrites the min.h and pair_tersoff.h files.

10. Make using the new Makefile.mcsmrff

		make mcsmrff -j 4

11. You're done! Now you can use lammps with the lmp_mcsmrff file in your lammps/src directory.

# To parameterize a force field

1.	Build example structures for your system in Avogadro (including any OPLS bonds)

2.	Optimize with orca (e.g. Opt B97-D3 def2-TZVP)

3.	Analyze at higher level (e.g. SP RI-B2PLYP D3BJ def2-TZVP)

4.	Great a training set directory that houses directories for each orca result.  Each directory must include the .out, .engrad and system.cml file for that system

5.	Run `python gradient_tersoff.py RUN_NAME`. The output parameters will start to appear in lammps/RUN_NAME_best.tersoff.

6.  ...