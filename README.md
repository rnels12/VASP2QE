# VASP2QE
A simple code to covert a VASP POSCAR to a Quantum ESPRESSO input file.
The code was used for a project that test the time-reversal symmtery feature of LOBSTER as reported in this publication:
https://s3-eu-west-1.amazonaws.com/itempdf74155353254prod/12024756/LOBSTER__Local_Orbital_Projections__Atomic_Charges__and_Chemical_Bonding_Analysis_from_Projector-Augmented-Wave-Based_DF_v1.pdf

The code needs the boost library and can be compiled by the following command:
g++ -o vasp2qe.x vasp2qe.cpp -lboost_filesystem  -lboost_system
