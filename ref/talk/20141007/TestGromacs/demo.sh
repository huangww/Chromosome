#!/usr/bin/env bash

# step 1: Dowload pdb file 
# http://www.rcsb.org/pdb/home/home.do

# step 2: generate some intermediated files
gmx pdb2gmx -f 1AKI.pdb -o 1AKI_processed.gro -water spce
# select the force field: 15

# step 3: define box and solvent
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

# step 4: setting MD simulation by writting a .mdp file(ions.mdp) and assemble .tpr file
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -nn 8
# choose group 13 "SOL" for embedding ions

# step 5: Assemble the binary input using grompp using this input parameter file: (minim.mdp)
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

# step 6: We are now ready to invoke mdrun to carry out the EM:
gmx mdrun -v -deffnm em

# step 7: analysis the results
gmx energy -f em.edr -o potential.xvg
# type "10 0" to select Potential (10); zero (0) terminates input. 
