import openbu.utils as utils

path1 = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/fy_lib'
path2 = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/jeff33/fy_lib'


parent_list = ['902320',
'922330',
'922340',
'922350',
'922360',
'922380',
'932370',
'932380',
'942380',
'942390',
'942400',
'942410',
'942420',
'952410',
'952430',
'962430',
'962440',
'962450']

#utils.plot_compare_libs_sum_over_parents(path1, path2, parent_list)
utils.plot_compare_libs(path1, path2, '922350')