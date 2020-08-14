"""
parses the input file for keywords
"""

from .find import first_capture
from .find import all_captures
from .pattern import capturing
from .pattern import zero_or_more
from .pattern import one_or_more
from .pattern import series
from .pattern import escape
from .pattern import NONSPACE
from .pattern import SPACE
from .pattern import WILDCARD
from .pattern import INTEGER
from .pattern import FLOAT
from .pattern import LOGICAL
from .pattern import LINE_FILL
from .pattern import NONNEWLINE
from .pattern import NEWLINE


INPUT_SUPPORTED_SECTIONS = [
    'training_data',
    'functional_form'
]
INPUT_REQUIRED_SECTIONS = [
    'training_data',
    'functional_form'
]

TD_SUPPORTED_KEYWORDS = [
    'EnergyRanges',
    'EnergyWeights',
    'Epsilon',
    'DataBatches',
    'BatchZeroes',
    'BatchWeights',
    'EnergyUnits',
]
FF_SUPPORTED_KEYWORDS = [
    'NumAtoms',
    'Symbols', 
    'AtomGroups',
    'TotalOrder',
    'FactorOrder',
    'ReadBasis',
    'FindDisconnected',
    'DisconnectedGroups',
    'MonitorGroups',
]

TD_REQUIRED_KEYWORDS = [
    'EnergyRanges',
    'Epsilon',
    'DataBatches',
    'BatchZeroes',
    'EnergyUnits',
]
FF_REQUIRED_KEYWORDS = [
    'NumAtoms',
    'Symbols',
    'AtomGroups',
    'TotalOrder',
    'FactorOrder',
    'ReadBasis',
    'FindDisconnected',
]

# Read the targets and baths sections and species


def read_energy_ranges(input_string):
    """ obtain 
    """

    pattern = ('EnergyRanges' +
               one_or_more(SPACE) + capturing(FLOAT) + # cut0
               one_or_more(SPACE) + capturing(FLOAT) + # cut1
               one_or_more(SPACE) + capturing(FLOAT) + # cut2
               one_or_more(SPACE) + capturing(FLOAT))  # cut3
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    out='  '.join(keyword)
#    print(keyword)
#    print(out)

    return out


def read_energy_weights(input_string):
    """ 
    """

    pattern = ('EnergyWeights' +
               one_or_more(SPACE) + capturing(FLOAT) + # weight0
               one_or_more(SPACE) + capturing(FLOAT) + # weight1
               one_or_more(SPACE) + capturing(FLOAT) + # weight2
               one_or_more(SPACE) + capturing(FLOAT))  # weight3
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    out='  '.join(keyword)

    return out


def read_epsilon(input_string):
    """ 
    """

    pattern = ('Epsilon' +
               one_or_more(SPACE) + 
               capturing(FLOAT))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
#    keyword = int(keyword)

    return keyword


def read_data_batches(input_string):
    """ 
    """

    pattern = ('DataBatches' +
               one_or_more(SPACE) + capturing(one_or_more(NONSPACE)) + 
               one_or_more(SPACE) + capturing(one_or_more(NONSPACE)) + 
               one_or_more(SPACE) + capturing(one_or_more(NONSPACE)))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    out='  '.join(keyword)

    return out


def read_batch_zeroes(input_string):
    """ 
    """

    pattern = ('BatchZeroes' +
               one_or_more(SPACE) + capturing(FLOAT) +
               one_or_more(SPACE) + capturing(FLOAT) +
               one_or_more(SPACE) + capturing(FLOAT))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    out='  '.join(keyword)

    return out


def read_batch_weights(input_string):
    """ 
    """

    pattern = ('BatchWeights' +
               one_or_more(SPACE) + capturing(FLOAT) +
               one_or_more(SPACE) + capturing(FLOAT) +
               one_or_more(SPACE) + capturing(FLOAT))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

#    assert keyword is not None
    if keyword is None:
        print("No BatchWeights found, setting to 1 1 1")
        keyword = ('1.','1.','1.')

    out='  '.join(keyword)

    return out
    

def read_energy_units(input_string):
    """ 
    """

    pattern = ('EnergyUnits' +
               one_or_more(SPACE) + capturing(one_or_more(NONSPACE)))
    block = _get_training_data_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def _get_training_data_section(input_string):
    """ grabs the section of text containing all of the job keywords
        for training data
    """
    pattern = (escape('$training_data') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section

#FUNCTIONAL FORM FUNCTIONS
def read_num_atoms(input_string):
    """ 
    """

    pattern = ('NumAtoms' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_symbols(input_string,natoms):
    """ 
    """

    symb_line = _get_symbols_line(input_string,natoms)

    out=' '.join(symb_line)

    assert symb_line is not None
    out=' '.join(symb_line)

    return out

def _get_symbols_line(input_string,natoms):
    """ grabs the line of text containing atom symbols
    """
    tmp = []
    for _ in range(int(natoms)):
        tmp.append(one_or_more(SPACE))
        tmp.append(capturing(NONSPACE))

    pattern = ('Symbols' + ''.join(tmp)
        )
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


def read_groups(input_string,natoms):
    """ 
    """

    groups_line = _get_groups_line(input_string,natoms)

    assert groups_line is not None
    out=' '.join(groups_line)

    return out

def _get_groups_line(input_string,natoms):
    """ grabs the line of text containing atom symbols
    """
    tmp = []
    for _ in range(int(natoms)):
        tmp.append(one_or_more(SPACE))
        tmp.append(capturing(INTEGER))

    pattern = ('AtomGroups' + ''.join(tmp))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


def read_total_order(input_string):
    """ 
    """

    pattern = ('TotalOrder' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_factor_order(input_string):
    """ 
    """

    pattern = ('FactorOrder' +
               one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_read_basis(input_string):
    """ 
    """

    pattern = ('ReadBasis' +
               one_or_more(SPACE) + capturing(LOGICAL))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_find_disconnected(input_string):
    """ 
    """

    pattern = ('FindDisconnected' +
               one_or_more(SPACE) + capturing(LOGICAL))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None

    return keyword


def read_disconnected_groups(input_string,natoms):
    """ 
    """

    groups_line = _get_disc_groups_line(input_string,natoms)

    assert groups_line is not None
    out=' '.join(groups_line)

    return out

def _get_disc_groups_line(input_string,natoms):
    """ grabs the line of text containing atom symbols
    """
    tmp = []
    for _ in range(int(natoms)):
        tmp.append(one_or_more(SPACE))
        tmp.append(capturing(INTEGER))

    pattern = ('DisconnectedGroups' + ''.join(tmp))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section


def read_monitor_groups(input_string):
    """ 
    """

    pattern = ('MonitorGroups'
            + one_or_more(SPACE) + capturing(INTEGER)
            + one_or_more(SPACE) + capturing(INTEGER))
    block = _get_functional_form_section(input_string)

    keyword = first_capture(pattern, block)

    assert keyword is not None
    out=' '.join(keyword)

    return out


def _get_functional_form_section(input_string):
    """ grabs the section of text containing all of the job keywords
        for functional form of PIPs
    """
    pattern = (escape('$functional_form') + LINE_FILL + NEWLINE +
               capturing(one_or_more(WILDCARD, greedy=False)) +
               escape('$end'))
    section = first_capture(pattern, input_string)

    assert section is not None

    return section

# Functions to check for errors in the input file

def check_training_data_keywords(input_string):
    """ obtains the keywords defined in the input by the user
    """
    section_string = _get_training_data_section(input_string)
    defined_keywords = _get_defined_keywords(section_string)

    # Check if keywords are supported
    if not all(keyword in TD_SUPPORTED_KEYWORDS
               for keyword in defined_keywords):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(keyword in defined_keywords
               for keyword in TD_REQUIRED_KEYWORDS):
        raise NotImplementedError

    print("Training Data Input:")
    print(section_string)


def check_functional_form_keywords(input_string):
    """ obtains the keywords defined in the input by the user
    """
    section_string = _get_functional_form_section(input_string)
    defined_keywords = _get_defined_keywords(section_string)

    # Check if keywords are supported
    if not all(keyword in FF_SUPPORTED_KEYWORDS
               for keyword in defined_keywords):
        raise NotImplementedError

    # Check if elements of keywords
    if not all(keyword in defined_keywords
               for keyword in FF_REQUIRED_KEYWORDS):
        raise NotImplementedError

    print("Functional Form Input:")
    print(section_string)


def _get_defined_keywords(section_string):
    """ gets a list of all the keywords defined in a section
    """

    defined_keys = []
    for line in section_string.splitlines():
        tmp = line.strip().split(' ')[0]
        defined_keys.append(tmp.strip())

    return defined_keys
