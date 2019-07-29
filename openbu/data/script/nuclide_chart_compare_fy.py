import openbu.utils as utils

path1 = '/home/julien/Open-Burnup.dev/openbu/data/default_libs/fy_lib'
path2 = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/jeff33/jeff33_fy_lib'

utils.plot_nuclide_chart_compare_fy(path1, path2, '922350')