import openbu.utils as utils

fy_lib = '/home/julien/Open-Burnup.dev/openbu/data/default_libs/fy_lib'
decay_lib = '/home/julien/Open-Burnup.dev/openbu/data/default_libs/decay_lib'

utils.plot_nuclide_chart_color_per_nuclear_data(decay_lib, fy_lib)