"""
Executes the automation part of 1DMin
"""

import os
from mako.template import Template


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def autofit_input(energy_ranges, energy_weights, epsilon,
                  data_batches, batch_zeroes, batch_weights, energy_units,
                  num_atoms, symbols, groups, total_order, factor_order,
                  gen_basis, gen_disconnected, check_disconnected):
    """ writes the autofit input file for each instance
    """

    # Set the dictionary for the autofit input file
    fill_vals = {
        "EnergyRanges": energy_ranges,
        "EnergyWeights": energy_weights,
        "Epsilon": epsilon,
        "DataBatches": data_batches,
        "BatchZeroes": batch_zeroes,
        "BatchWeights": batch_weights,
        "EnergyUnits": energy_units,
        "NumAtoms": num_atoms,
        "Symbols": symbols,
        "Groups": groups,
        "TotalOrder": total_order,
        "FactorOrder": factor_order,
        "GenerateBasis": gen_basis,
        "GenerateDisconnected": gen_disconnected,
        "CheckDisconnected": check_disconnected
    }

    # Set template name and path for the autofit input file
    template_file_name = 'autofit_inp.mako'
    template_file_path = os.path.join(TEMPLATE_PATH, template_file_name)

    # Build the 1dmin input string
    input_str = Template(filename=template_file_path).render(**fill_vals)

    return input_str
