# VASP2QE
A simple code to convert a VASP POSCAR to a Quantum ESPRESSO (QE) input file.
The code was used in a project testing the time-reversal symmtery feature of LOBSTER as reported in this publication:
https://s3-eu-west-1.amazonaws.com/itempdf74155353254prod/12024756/LOBSTER__Local_Orbital_Projections__Atomic_Charges__and_Chemical_Bonding_Analysis_from_Projector-Augmented-Wave-Based_DF_v1.pdf

The code needs the boost library and can be compiled by the following command:
g++ -o vasp2qe.x vasp2qe.cpp -lboost_filesystem  -lboost_system

Once the binary is generated, one can just simplly use the following command if POSCAR/CONTCAR is already available in the working dir:

./vasp2qe.x

or if POSCAR/CONTCAR in a different directory:

./vasp2qe.x <path-to-poscar-or-contcar>/POSCAR

The result will be output to the stdout. Once the QE input file is produced, you still need to set the values for the variables ecutwfc and nbnd, the pseudo-potential filenames, and the number of KPOINTS in order for the input file to be useable for a QE calculation.
