# autofit
## automated PIP fitting code

## Authors
Ahren W. Jasper

## Functionality

## Contact
Ahren Jasper [ajasper@anl.gov]

# Install
./comp.x in autofit/src

For pyfit wrapper input generator:\
   Mako Templates module for Python is required\
   see [https://www.makotemplates.org/]

# Input Files
## Fortran input file
The Fortran input file contains a series of records, each a single line of input parameters (see example fit.in files)
Alternatively, input can be formatted in more user-friendly file and converted with accompanying Python wrapper (located in /pyfit)
Below are listed each record and associated parameters, with data type in parantheses

**Record 1**: EnergyRanges\
cut0 cut1 cut2 cut3 (dp)\
-Values of energy ranges in kcal/mol with cut0 > cut1 > cut2 > cut3

**Record 2**: Epsilon\
epsilon (dp)\
-Range parameter in kcal/mol for the weight function (see )

**Record 3**: NumBatches\
batches (int)\
-Number of user defined batches of training set data in ai.all

**Record 4**: BatchZeroes\
z1 z2... (dp)\
-Sets the "zeroes" for the weight function for each batch

**Record 5**: BatchWeights\
w1 w2... (dp)\
-Sets the relative weights for each batch -- leave as 1 1 1 for default weighting

**Record 6**: NumAtoms\
natom (int)\
-Number of atoms in system

**Record 7**: AtomGroups\
iagroup(natom) (int)\
-Atom group of each atom

**Record 8**: Symbols\
symb(natom) (char * 2)\
-Atom symbol for each atom

**Record 9**: FactorOrder, TotalOrder\
ipow ipowt (int)\
-Maximum allowed order for each factor and total order of the term

**Record 10**: ReadBasis\
lreadbasis (logical)\
    = TRUE, read basis from basis.dat\
    = FALSE, compute basis and write to basis.dat

**Record 11**: FindDisconnected\
lreaddisc (logical)\
    = TRUE, determine disconnected atom groups to remove from basis\
    = FALSE, do not consider disconnected groups

if lreaddisc=TRUE:\
**Record 12**: DisconnectedGroups\
dgroup(natom) (int)\
-Disconnected molecular groups in the product system. Currently only supports a single product channel.

## ai.all, ai.test

