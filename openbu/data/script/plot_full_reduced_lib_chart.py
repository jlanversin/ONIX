import openbu.utils as utils

decay_path1 = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/decay_lib'
decay_path2 = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/decay_lib_reduced'
fy_path1 = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/fy_lib'
fy_path2 = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/fy_lib_reduced'

utils.plot_compare_two_nuclear_data_on_nuclide_chart(decay_path1,fy_path1,decay_path2, fy_path2)