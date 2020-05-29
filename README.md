# autofit
## automated PIP fitting code

## Authors
Ahren W. Jasper

## Functionality

## Contact
Ahren Jasper [ajasper@anl.gov]

## Install

## Input Files
# fit.in
Contains a series of records, each a single line of input parameters (see example fit.in files)
Below are listed each record and associated parameters, with data type in parantheses

Record 1: cut0, cut1, cut2, cut3
cut0 (dp)	Values of energy ranges in kcal/mol with cut0 > cut1 > cut2 > cut3

Record 2: epsilon
epsilon (dp)	Range parameter in kcal/mol for the weight function (see )

Record 3: sx(3)
sx(3) (dp)	Sets the "zeroes" for the weight function for reactant, product, and saddle point (R,P,SP) batches

Record 4: sc(3)
sc(3) (dp)	Sets the relative weights for each batch (R,P,SP) -- leave as 1,1,1 for default weighting

Record 5: natom
natom (int)	Number of atoms in system

Record 6: iagroup(natom)
iagroup(natom) (int)	Atom group of each atom

Record 7: symb(natom)
symb(natom) (char * 2)	Atom symbol for each atom

Record 8: ipow, ipowt
ipow (int)	Maximum allowed order for each factor
ipowt (int)	Total order of the term

Record 9: lreadbasis, lreaddisc
lreadbasis (logical)	= TRUE, read basis from basis.dat
			= FALSE, compute basis and write to basis.dat
lreaddisc (logical)	= TRUE, determine disconnected atom groups to remove from basis
			= FALSE, do not consider disconnected groups

Record 10: dsymb1,dsymb2 (if lreaddisc=TRUE)
dsymb1 (char * 2)	Atom symbol of bond type to monitor for disconnected terms 

# ai.dat, ai2.dat, ai3.dat

