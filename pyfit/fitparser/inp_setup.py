"""
input file writing routines for parallel OneDMin runs
"""

import os
import subprocess
import fitparser


def write_fit_input_mako(job_dir_path,
            energy_ranges, energy_weights, epsilon, data_batches, 
            batch_zeroes, batch_weights, energy_units, num_atoms, symbols,  
            atom_groups, total_order, factor_order, read_basis, 
            find_disconnected, disconnected_groups, 
            monitor_groups):
    """ write autofit input file
    """

    inp_str = fitparser.writer.autofit_input(
        energy_ranges, energy_weights, epsilon,
        data_batches, batch_zeroes, batch_weights, energy_units,
        num_atoms, symbols, atom_groups, total_order, factor_order,
        read_basis, find_disconnected, disconnected_groups, 
        monitor_groups)

    job_file_path = os.path.join(job_dir_path, 'fit.in')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)

def write_fit_input(job_dir_path,
            energy_ranges, energy_weights, epsilon, data_batches, 
            batch_zeroes, batch_weights, energy_units, num_atoms, symbols,  
            atom_groups, total_order, factor_order, read_basis, 
            find_disconnected, disconnected_groups, 
            monitor_groups):
    """ write autofit input file
    """

    inp_str = fitparser.writer.autofit_input(
        energy_ranges, energy_weights, epsilon,
        data_batches, batch_zeroes, batch_weights, energy_units,
        num_atoms, symbols, atom_groups, total_order, factor_order,
        read_basis, find_disconnected, disconnected_groups, 
        monitor_groups)

    job_file_path = os.path.join(job_dir_path, 'fit.in')
    with open(job_file_path, 'w') as input_file:
        input_file.write(inp_str)

