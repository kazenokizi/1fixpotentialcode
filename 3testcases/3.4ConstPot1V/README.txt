This folder contains:

1	Input file: Ion1.in, Ion2.in
2	Data file:  Ion0.restart
These two files could be put together to run the 
simulation with 1V potential criteria on atoms of type
5 & 6. Further soft coding on atom types is required. 

3	Bash file:  job.sh
This file is used to call loop runs for obtaining voltages
on atoms. It requires files from folder Input. 

When the simulation is done, files as follows will be
generated:
1	charges.txt
It contains the atom # and its charge. It includes both the
system atoms and the ghost atoms.

2	list1.txt
This file is processed during the simulation due to fix_potential
development. It contains the list and the atom #, and they 
will be used for mapping during the fix.

3	list2.txt
This file is processed during the simulation due to fix_potential
development. It contains the atom# and its corresponding local
atom #, and they will be used for back-mapping during the fix.

4	mapresult.txt
This file is processed during the simulation due to fix_potential
development. It contains the atom# and its corresponding local
atom #, and they are used to check the success of mapping. 

5	charge.dump
This file includes the atom positions and their charges during 
each timestep. 

6	Ion1.data/.restart
These files are generated after the fix_potential test.

7	Ion2.data/.restart
These files are generated after the rerun. 

8	log.lammps
Self-explanatory. 