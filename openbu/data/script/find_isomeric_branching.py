import openbu.data as d
import openbu.utils as utils

path = '/home/julien/Open-Burnup.dev/openbu/data/default_libs/xs_lib'

xs_dic = d.read_xs_lib(path)

selected_zamid = []
selected_name = []
for nucl in xs_dic:
	ngamma = xs_dic[nucl]['(n,gamma)'][0]
	ngammaX = xs_dic[nucl]['(n,gamma)X'][0]
	total = ngamma+ngammaX
	name = utils.zamid_to_name(nucl)
	if ngamma == 0:
		continue
	# if ngamma ==0 and ngammaX != 0: # There are apparently no ngammaX with ngamma  = 0
	# 	print (nucl)
	# 	print (name)
	# 	print('ngamma {}, ngammaX {}\n\n\n\n'.format(ngamma, ngammaX))
	# 	selected_0ngamma_zamid.append(nucl)
	# 	#selected_0ngamma_name.append(name)
	elif ngammaX/total > 0.1:
		# print (nucl)
		# print (name)
		# print('ngamma {}, ngammaX {}'.format(ngamma, ngammaX))
		# print ('ratio:{}\n\n\n\n'.format(ngammaX/total))
		selected_zamid.append(nucl)
		#selected_name.append(name)

order_selected_zamid = utils.order_nuclide_per_z(selected_zamid)
selected_name = utils.zamid_list_to_name_list(order_selected_zamid)


print ('\n\n\n\n\n')
print (len(selected_name))
for nucl in selected_name:
	print (nucl)
