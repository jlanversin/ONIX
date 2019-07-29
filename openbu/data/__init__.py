# This file contain diverse dictionary that are called by the functions and class of Salameche

import os
import math as m
from .read_lib_functions import *
from . import list_and_dict


### Generate the default libraries

# Default libraries path

default_atm_mass_lib_path = os.path.join(os.path.dirname(__file__), 'default_libs/mass.mas12')
default_decay_b_lib_path = os.path.join(os.path.dirname(__file__), 'default_libs/decay_lib')
default_xs_lib_path = os.path.join(os.path.dirname(__file__), 'default_libs/xs_lib')
default_fy_lib_path = os.path.join(os.path.dirname(__file__), 'default_libs/fy_lib')

# Default libraries data

default_atm_mass_lib = read_mass_lib(default_atm_mass_lib_path)
default_decay_lib_b = read_decay_lib(default_decay_b_lib_path)
default_decay_lib_a = conv_decay_b_a(default_decay_lib_b)
default_xs_lib = read_xs_lib(default_xs_lib_path)
default_fy_lib = read_fy_lib(default_fy_lib_path)


### Generate the default matrices and nuclide lists

default_C = default_decay_mat_from_Ctxt()
default_B = default_xs_mat_from_Btxt()
default_nucl_list = default_nucl_list_from_txt()