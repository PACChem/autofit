"""
Sets up Autofit
"""

import os
import sys
import fitparser

# Set paths
DRIVE_PATH = os.getcwd()

# Read the input file into a string
with open(os.path.join(DRIVE_PATH, 'input.dat'), 'r') as infile:
    INPUT_STRING = infile.read()

# Check training data parameters
fitparser.check_training_data_keywords(INPUT_STRING)

# Read the training data
ENERGY_RANGES = fitparser.read_energy_ranges(INPUT_STRING)
ENERGY_WEIGHTS = fitparser.read_energy_weights(INPUT_STRING)
EPSILON = fitparser.read_epsilon(INPUT_STRING)
DATA_BATCHES = fitparser.read_data_batches(INPUT_STRING)
BATCH_ZEROES = fitparser.read_batch_zeroes(INPUT_STRING)
BATCH_WEIGHTS = fitparser.read_batch_weights(INPUT_STRING)
ENERGY_UNITS = fitparser.read_energy_units(INPUT_STRING)

# Check functional form parameters
fitparser.check_functional_form_keywords(INPUT_STRING)

# Read the functional form
NATOMS = fitparser.read_num_atoms(INPUT_STRING)
SYMB_LST = fitparser.read_symbols(INPUT_STRING,NATOMS)
GROUPS = fitparser.read_groups(INPUT_STRING,NATOMS)
TOTAL_ORDER = fitparser.read_total_order(INPUT_STRING)
FACTOR_ORDER = fitparser.read_factor_order(INPUT_STRING)
READ_BASIS_FLAG = fitparser.read_read_basis(INPUT_STRING)
FIND_DISCONNECTED_FLAG = fitparser.read_find_disconnected(INPUT_STRING)
DISCONNECTED_GROUPS = fitparser.read_disconnected_groups(INPUT_STRING,NATOMS)
MONITOR_GROUPS = fitparser.read_monitor_groups(INPUT_STRING)

# Write paramaters to input file
print('Writing input file to fit.in ...')
fitparser.inp_setup.write_fit_input(
    job_dir_path=DRIVE_PATH, energy_ranges=ENERGY_RANGES, 
    energy_weights=ENERGY_WEIGHTS, epsilon=EPSILON,
    data_batches=DATA_BATCHES, batch_zeroes=BATCH_ZEROES, 
    batch_weights=BATCH_WEIGHTS, energy_units=ENERGY_UNITS, 
    num_atoms=NATOMS, symbols=SYMB_LST, atom_groups=GROUPS, 
    total_order=TOTAL_ORDER, factor_order=FACTOR_ORDER,
    read_basis=READ_BASIS_FLAG, 
    find_disconnected=FIND_DISCONNECTED_FLAG,
    disconnected_groups=DISCONNECTED_GROUPS,
    monitor_groups=MONITOR_GROUPS)
