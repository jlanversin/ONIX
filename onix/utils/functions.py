# Small tricks and calculations to make the life of an nuclear engineer easier
from math import log
import os
import shutil
from onix.data import time_dic
import onix.data as d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.colors as colors
import xml.etree.ElementTree as ET

NA = 6.02214086e+23

def decay_to_halflife(decay_constant, unit):

	half_life_s = log(2)/decay_constant
	half_life = half_life_s/time_dic[unit]

	return half_life


def halflife_to_decay(half_life, unit):

	half_life_s = half_life*time_dic[unit]
	decay_constant = log(2)/half_life_s

	return decay_constant

def halflife_to_second(half_life, unit):

	half_life_s = half_life*time_dic[unit]

	return half_life_s


def is_int(s):
		try:
			int(s)
			return True
		except ValueError:
			return False

def is_zamid(string):

	statement = True

	for number in string:
		if not is_number(number):
			statement = False

	return statement

def is_name(string):

	if '-' in string:
		return True
	else:
		return False

def is_list_redundant(l):

	result = False
	l_set = set(l)
	if len(l) > len(l_set):
		result = True

def get_list_redundant_elt(l):

	count_elt = []
	for elt in l:
		if elt in count_elt:
			redundant_elt.append(elt)
		else:
			count_elt.append(elt)

	return redundant_elt

def zamid_list_to_name_list(zamid_list):

	name_list = []
	for zamid in zamid_list:
		name = zamid_to_name(zamid)
		name_list.append(name)

	return name_list

def name_list_to_zamid_list(name_list):

	zamid_list = []
	for name in name_list:
		zamid = name_to_zamid(name)
		zamid_list.append(zamid)

	return zamid_list

def get_zamid_z(zamid):

	z = int(zamid[:-4])

	return z

def get_name_z(name):

	zamid = name_to_zamid(name)
	z = int(zamid[:-4])

	return z

def get_zamid_a(zamid):

	a = int(zamid[-4:-1])

	return a

def get_zamid_n(zamid):

	a = int(zamid[-4:-1])
	z = int(zamid[:-4])

	return a - z

def get_zamid_s(zamid):

	s = int(zamid[-1])

	return s

def zamid_to_name(zamid):
	""""Finds and returns the name of the nuclide"""
	dic = d.nuc_name_dic

	if len(zamid) == 5:
		nz = int(zamid[0:1])
		na = int(zamid[1:4])
		state = int(zamid[4])
	if len(zamid) == 6:
		nz = int(zamid[0:2])
		na = int(zamid[2:5])
		state = int(zamid[5])
	if len(zamid) == 7:
		nz = int(zamid[0:3])
		na = int(zamid[3:6])
		state = int(zamid[6])

	if state == 0:
		nuc_name = '{}-{}'.format(d.nuc_zz_dic[nz], na) 
	else:
		nuc_name = '{}-{}*'.format(d.nuc_zz_dic[nz], na)

	return nuc_name

def name_to_zamid(name):
	"""Finds and returns the zzaaam id of the nuclide"""
	dic = d.nuc_name_dic

	elt_name = name.split('-')[0]
	na = int(name.split('-')[1].replace('*',''))
	if '*' in name:
		state = 1
	else:
		state = 0
	zzaaam = 10000*d.nuc_name_dic[elt_name] + na*10 + state
	zamid = str(zzaaam)

	return zamid


def get_hm(passlist, hm_vol):

	hmmd = 0 # Heavy Metal Mass Density
	for i in passlist.passport_list:
		if int(i.zamid[:-4]) > 90:
			hmmd += i.current_dens*1E+24*i.mass/NA

	hm = hmmd*hm_vol

	return hm

def convert_mass_to_atom(mass, nuclide):

	if is_name(nuclide):
		zamid = name_to_zamid(nuclide)
		zaid =zamid[:-1]
	else:
		zamid = nuclide
		zaid =zamid[:-1]

	molar_mass = d.default_atm_mass_lib[zaid]

	atom = mass*NA/molar_mass

	return atom

def convert_atom_to_mass(atom, nuclide):

	if is_name(nuclide):
		zamid = name_to_zamid(nuclide)
		zaid =zamid[:-1]
	else:
		zamid = nuclide
		zaid =zamid[:-1]

	molar_mass = d.default_atm_mass_lib[zaid]

	mass = atom*molar_mass/NA

	return mass

def get_bu_sec_conv_factor(vol, ihm):

	bu_sec_conv_factor = vol*1e-3/(ihm*24*3600) # Unit in L/g

	return bu_sec_conv_factor

def get_keylist_from_dict(dict):

	keylist = list(dict.keys())

	return keylist

def get_decay_nucl(decay_a_lib):

	decay_nucl = []
	for zamid in decay_a_lib:
		decay_nucl.append(zamid)

	return decay_nucl

def get_xs_nucl(xs_lib):

	xs_nucl = []
	for zamid in xs_lib:
		xs_nucl.append(zamid)

	return xs_nucl

def get_fy_nucl(fy_lib):

	fy_nucl = []
	for zamid in fy_lib:
		fy_nucl.append(zamid)

	return fy_nucl

def get_all_nucl(list_of_dict):

	unfiltered_list = []
	for dictionary in list_of_dict:
		unfiltered_list += get_keylist_from_dict(dictionary)

	all_nucl = list(set(unfiltered_list))

	return all_nucl

def is_lista_in_listb(lista, listb):

	result = all(elem in listb for elem in lista)

	return result

def get_fy_parent_nucl(fy_lib):

	fy_nucl = get_fy_nucl(fy_lib)
	fy_parent = []
	sample_zamid = fy_nucl[0]
	sample = fy_lib[sample_zamid]

	for fission_parent in sample:
		fy_parent.append(fission_parent)

	return fy_parent


def get_cell_folder_path(file_name, *dir_path):

	folder_name = '{}_cell'.format(file_name)

	if not dir_path:
		dir_path = os.getcwd()

	folder_path = dir_path + '/' + folder_name

	return folder_path

def gen_cell_folder(name, *dir_path):

	if dir_path:
		folder_path = get_cell_folder_path(name, dir_path)
	elif not dir_path:
		folder_path = get_cell_folder_path(name)

	if os.path.exists(folder_path):
		shutil.rmtree(folder_path)

	os.makedirs(folder_path)

def get_folder_path(folder_name, *dir_path):

	if not dir_path:
		dir_path = os.getcwd()

	folder_path = dir_path + '/' + folder_name

	return folder_path

def gen_folder(folder_name, *dir_path):

	if dir_path:
		folder_path = get_folder_path(folder_name, dir_path)
	elif not dir_path:
		folder_path = get_folder_path(folder_name)

	if os.path.exists(folder_path):
		shutil.rmtree(folder_path)

	os.makedirs(folder_path)

# Convert a cell dictionary into a cell list.
# Dict keys must be cell IDs
# List is ordered according to cells ID
def cell_dict_to_cell_list(cell_dict):

	cell_list = []
	for ID in cell_dict:
		cell_list.append(cell_dict[ID])

	changed = True
	while changed:
		changed = False
		for i in range(len(cell_list) - 1):
			cell_id0 = cell_list[i].id
			cell_id1 = cell_list[i+1].id
			if cell_id0 > cell_id1:
				cell_list[i], cell_list[i+1] = cell_list[i+1], cell_list[i]
				changed = True


	return cell_list

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

def openmc_name_to_onix_name(name):

	i = 0
	while not is_number(name[i]):
		i += 1

	am = name[i:]
	am = am.replace('_m1', '*')
	am = am.replace('m', '*') # used in jeff3.3
	am = am.replace('n', '*') # used in jeff3.3
	onix_name = name[:i] + '-' + am



	return onix_name

def onix_name_to_openmc_name(name):

	openmc_name = name.replace('-', '')
	openmc_name = openmc_name.replace('*','_m1')

	return openmc_name

def bu_namelist_to_mc_namelist(name_list):

	new_name_list = []
	for name in name_list:
		new_name = onix_name_to_openmc_name(name)
		new_name_list.append(new_name)

	return new_name_list

def mc_namelist_to_bu_namelist(name_list):

	new_name_list = []
	for name in name_list:
		new_name = openmc_name_to_onix_name(name)
		new_name_list.append(new_name)

	return new_name_list

def order_nuclide_per_z(nucl_list):

		changed = True
		while changed:
			changed = False
			for i in range(len(nucl_list) - 1):
				zamid0 = int(nucl_list[i])
				zamid1 = int(nucl_list[i+1])

				if zamid0 > zamid1:
					nucl_list[i], nucl_list[i+1] = nucl_list[i+1], nucl_list[i]
					changed = True

		return nucl_list

def order_nuclide_name_per_z(nucl_name_list):

	nucl_name_list_old_format = mc_namelist_to_bu_namelist(nucl_name_list)
	zamid_list = name_list_to_zamid_list(nucl_name_list_old_format)
	ordered_zamid_list = order_nuclide_per_z(zamid_list)
	ordered_nucl_name_list_old_format = zamid_list_to_name_list(ordered_zamid_list)
	ordered_nucl_name_list = bu_namelist_to_mc_namelist(ordered_nucl_name_list_old_format)

	return ordered_nucl_name_list

def order_nuclide_per_a(nucl_list):

		changed = True
		while changed:
			changed = False
			for i in range(len(nucl_list) - 1):
				zamid0 = nucl_list[i]
				zamid1 = nucl_list[i+1]


				if zamid0 > zamid1:
					nucl_list[i], nucl_list[i+1] = nucl_list[i+1], nucl_list[i]
					changed = True

		return nucl_list


def plot_compare_libs(lib1_path, lib2_path, fissile_parent):

	fy_dict1 = d.read_fy_lib(lib1_path)
	fy_dict2 = d.read_fy_lib(lib2_path)

	unordered_keys1 = get_keylist_from_dict(fy_dict1)
	unordered_keys2 = get_keylist_from_dict(fy_dict2)

	ordered_keys1 = order_nuclide_per_z(unordered_keys1)
	ordered_keys2 = order_nuclide_per_z(unordered_keys2)

	rel_diff = []
	abs_diff = []
	nucl_list = []
	fy_list1 = []
	fy_list2 = []

	# sum over fy in list 1
	for nuclide in ordered_keys1:
		fy_1 = fy_dict1[nuclide][fissile_parent][0]
		fy_list1.append(fy_1)

	# sum over fy in list 2
	for nuclide in ordered_keys2:
		fy_2 = fy_dict2[nuclide][fissile_parent][0]
		fy_list2.append(fy_2)


	for nuclide in ordered_keys1:
		if nuclide in ordered_keys2:
			fy_1 = fy_dict1[nuclide][fissile_parent][0]
			fy_2 = fy_dict2[nuclide][fissile_parent][0]

			if fy_1 == 0 and fy_2 == 0:
				rel_diff_val = 0
			elif fy_1 == 0 and fy_2 != 0: # if fy1 = 0 but fy2 != 1, set label to 100%. The inverse situation yield - 100%
				rel_diff_val = 100
			else:
				rel_diff_val = (fy_2 - fy_1)*100/fy_1

			abs_diff_val = fy_2 -fy_1

			rel_diff.append(rel_diff_val)
			abs_diff.append(abs_diff_val)


			nucl_list.append(nuclide)

	x_vect = [i for i in range(len(rel_diff))]

	z_vect = []
	x_z_vect = []
	count = 0
	for i in range(len(nucl_list)):
		if count == 0:
			z_vect.append(nucl_list[i][:-4])
			x_z_vect.append(i)
		count += 1

		if count == 49:
			count = 0

	plt.figure(1)
	plt.title('Relative percent difference between fission yield data')
	plt.xticks(x_z_vect, z_vect)
	plt.plot(x_vect, rel_diff, 'ro')
	plt.plot(x_vect, [100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')
	plt.ylabel('Relative percent difference')
	plt.xlabel('Atomic number')
	#plt.yscale('log')

	plt.figure(2)
	plt.xticks(x_vect, nucl_list)
	plt.plot(x_vect, rel_diff, 'ro')
	plt.plot(x_vect, [100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')
	#plt.yscale('log')

	plt.figure(3)
	plt.xticks(x_vect, nucl_list)
	plt.plot(x_vect, abs_diff, 'ro')
	# plt.plot(x_vect, [100]*len(rel_diff), 'k')
	# plt.plot(x_vect, [-100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')



	plt.show()


def plot_compare_libs_sum_over_parents(lib1_path, lib2_path, parent_list):

	fy_dict1 = d.read_fy_lib(lib1_path)
	fy_dict2 = d.read_fy_lib(lib2_path)

	unordered_keys1 = get_keylist_from_dict(fy_dict1)
	unordered_keys2 = get_keylist_from_dict(fy_dict2)

	ordered_keys1 = order_nuclide_per_z(unordered_keys1)
	ordered_keys2 = order_nuclide_per_z(unordered_keys2)

	rel_diff = []
	abs_diff = []
	nucl_list = []
	# fy_list1 = []
	# fy_list2 = []

	# # sum over fy in list 1
	# for nuclide in ordered_keys1:
	# 	sum_fy = 0
	# 	for parent in parent_list:
	# 		sum_fy += fy_dict1[nuclide][parent][0]
	# 	fy_list1.append(sum_fy)

	# # sum over fy in list 2
	# for nuclide in ordered_keys2:
	# 	sum_fy = 0
	# 	for parent in parent_list:
	# 		sum_fy += fy_dict2[nuclide][parent][0]
	# 	fy_list2.append(sum_fy)


	for nuclide in ordered_keys1:
		if nuclide in ordered_keys2:
			sum_fy1 = 0
			for parent in parent_list:
				sum_fy1 += fy_dict1[nuclide][parent][0]
			sum_fy2 = 0
			for parent in parent_list:
				sum_fy2 += fy_dict2[nuclide][parent][0]


			if sum_fy1 == 0 and sum_fy2 == 0:
				rel_diff_val = 0
			elif sum_fy1 == 0 and sum_fy2!= 0: # if fy1 = 0 but fy2 != 1, set label to 100%. The inverse situation yield - 100%
				rel_diff_val = 100
			else:
				rel_diff_val = (sum_fy2 - sum_fy1)*100/sum_fy1

			abs_diff_val = sum_fy2 -sum_fy1

			rel_diff.append(rel_diff_val)
			abs_diff.append(abs_diff_val)


			nucl_list.append(nuclide)

	x_vect = [i for i in range(len(rel_diff))]

	z_vect = []
	x_z_vect = []
	count = 0
	for i in range(len(nucl_list)):
		if count == 0:
			z_vect.append(nucl_list[i][:-4])
			x_z_vect.append(i)
		count += 1

		if count == 49:
			count = 0

	plt.figure(1)
	plt.title('Relative percent difference between fission yield data')
	plt.xticks(x_z_vect, z_vect)
	plt.plot(x_vect, rel_diff, 'ro')
	plt.plot(x_vect, [100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')
	plt.ylabel('Relative percent difference')
	plt.xlabel('Atomic number')
	#plt.yscale('log')

	plt.figure(2)
	plt.xticks(x_vect, nucl_list)
	plt.plot(x_vect, rel_diff, 'ro')
	plt.plot(x_vect, [100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')
	#plt.yscale('log')

	plt.figure(3)
	plt.xticks(x_vect, nucl_list)
	plt.plot(x_vect, abs_diff, 'ro')
	# plt.plot(x_vect, [100]*len(rel_diff), 'k')
	# plt.plot(x_vect, [-100]*len(rel_diff), 'k')
	plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')

	plt.show()

# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))


def plot_nuclide_chart_compare_fy(lib1_path, lib2_path, fissile_parent):

	fy_dict1 = d.read_fy_lib(lib1_path)
	fy_dict2 = d.read_fy_lib(lib2_path)

	unordered_keys1 = get_keylist_from_dict(fy_dict1)
	unordered_keys2 = get_keylist_from_dict(fy_dict2)

	ordered_keys1 = order_nuclide_per_z(unordered_keys1)
	ordered_keys2 = order_nuclide_per_z(unordered_keys2)

	rel_diff = []
	label_list = []
	nucl_list = []
	fy_list1 = []
	fy_list2 = []

	# sum over fy in list 1
	for nuclide in ordered_keys1:
		fy_1 = fy_dict1[nuclide][fissile_parent][0]
		fy_list1.append(fy_1)

	# sum over fy in list 2
	for nuclide in ordered_keys2:
		fy_2 = fy_dict2[nuclide][fissile_parent][0]
		fy_list2.append(fy_2)


	for nuclide in ordered_keys1:
		if nuclide in ordered_keys2:
			fy_1 = fy_dict1[nuclide][fissile_parent][0]
			fy_2 = fy_dict2[nuclide][fissile_parent][0]

			if fy_1 == 0 and fy_2 == 0:
				label = 'both = 0'
				rel_diff_val = 0
			elif fy_1 == 0 and fy_2 != 0: # if fy1 = 0 but fy2 != 1, set label to 100%. The inverse situation yield - 100%
				label = 'ori = 0'
				rel_diff_val = 0
			elif fy_1 != 0 and fy_2 == 0: # if fy1 = 0 but fy2 != 1, set label to 100%. The inverse situation yield - 100%
				label = 'jef33 = 0'
				rel_diff_val = 0
			else:
				rel_diff_val = (fy_2 - fy_1)*100/fy_1

				if abs(rel_diff_val) > 500:
					rel_diff_val = 500*abs(rel_diff_val)/rel_diff_val # set value to 500 or -500

				label = '{:4.2E}'.format((fy_2 - fy_1)*100/fy_1)

			# if rel_diff_val == 0:
			# 	rel_diff_log.append(abs(rel_diff_val))
			# 	rel_diff.append(rel_diff_val)
			# else:
			rel_diff.append(rel_diff_val)
			label_list.append(label)

			nucl_list.append(nuclide)



	not_in_jeff33_zamid = []
	for zamid in ordered_keys1:
		if zamid not in ordered_keys2:
			not_in_jeff33_zamid.append(zamid)

	not_in_endf4_zamid = []
	for zamid in ordered_keys2:
		if zamid not in ordered_keys1:
			not_in_endf4_zamid.append(zamid)

	# test if not_in list overlap
	for zamid in nuclide:
		if zamid in not_in_endf4_zamid:
			print (zamid)
		if zamid in not_in_jeff33_zamid:
			print (zamid)


	# loop over the z
	# z list for lib 1
	z_list1 = []
	for zamid in ordered_keys1:
		if zamid in ordered_keys2:
			z = get_zamid_z(zamid)
			z_list1.append(z)


	# loop over the z for the nuclide in endf4 but not jeff33
	not_in_jeff33_z_list = []
	for zamid in not_in_jeff33_zamid:
		z = get_zamid_z(zamid)
		not_in_jeff33_z_list.append(z)

	# loop over the z for the nuclide in jeff33 but not endf4
	not_in_endf4_z_list = []
	for zamid in not_in_endf4_zamid:
		z = get_zamid_z(zamid)
		not_in_endf4_z_list.append(z)

	# # z list for lib 2
	# z_list2 = []
	# for zamid in ordered_keys2:
	# 	z = get_zamid_z(zamid)
	# 	z_list2.append(z)

	# loop over the a
	# a list for lib 1
	a_list1 = []
	for zamid in ordered_keys1:
		if zamid in ordered_keys2:
			a = get_zamid_a(zamid)
			a_list1.append(a)

	# loop over the a for the nuclide in endf4 but not jeff33
	not_in_jeff33_a_list = []
	for zamid in not_in_jeff33_zamid:
		a = get_zamid_a(zamid)
		not_in_jeff33_a_list.append(a)

	# loop over the z for the nuclide in jeff33 but not endf4
	not_in_endf4_a_list = []
	for zamid in not_in_endf4_zamid:
		a = get_zamid_a(zamid)
		not_in_endf4_a_list.append(a)

	# loop over the n
	# a list for lib 1
	n_list1 = []
	for zamid in ordered_keys1:
		if zamid in ordered_keys2:
			n = get_zamid_n(zamid)
			n_list1.append(n)

	# loop over the n for the nuclide in endf4 but not jeff33
	not_in_jeff33_n_list = []
	for zamid in not_in_jeff33_zamid:
		n = get_zamid_n(zamid)
		not_in_jeff33_n_list.append(n)

	# loop over the n for the nuclide in jeff33 but not endf4
	not_in_endf4_n_list = []
	for zamid in not_in_endf4_zamid:
		n = get_zamid_n(zamid)
		not_in_endf4_n_list.append(n)


	# check if list overlap

	# # z list for lib 2
	# a_list2 = []
	# for zamid in ordered_keys2:
	# 	a = get_zamid_a(zamid)
	# 	a_list2.append(a)


	N = len(rel_diff) # Number of labels

	# setup the plot
	plt.figure(1)

	tag = rel_diff# Tag each point with a corresponding label    

	# define the colormap
	cmap = plt.cm.Spectral
	# extract all colors from the .jet map
	# cmaplist = [cmap(i) for i in range(cmap.N)]
	# # create the new map
	# cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)	 

	# define the bins and normalize
	bounds = np.linspace(0,N,N+1)
	#norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

	# make the scatter
	scat_a_z = plt.scatter(n_list1,z_list1,c=tag,cmap=cmap, marker='s', s = [20]*len(z_list1), norm=MidpointNormalize(midpoint=0,vmin=-100, vmax=max(rel_diff)))
	#scat_a_z = plt.scatter(a_list1,z_list1,c=tag,cmap=cmap, marker='s', s = [20]*len(z_list1), norm=MidpointNormalize(midpoint=0,vmin=-100, vmax=max(rel_diff)))
	
	scat2 = plt.scatter(not_in_jeff33_a_list,not_in_jeff33_z_list, marker = 'o', color='k', s = [20]*len(not_in_jeff33_z_list), label='Not in JEF33')
	scat3 = plt.scatter(not_in_endf4_a_list,not_in_endf4_z_list, marker = '*', color='k', s = [20]*len(not_in_endf4_z_list), label='Not in ENDF4')
	
	#cb = plt.colorbar(scat, spacing='proportional',ticks=bounds)
	#plt.imshow(ras, cmap=cmap, clim=(elev_min, elev_max), norm=MidpointNormalize(midpoint=mid_val,vmin=elev_min, vmax=elev_max))

	cb = plt.colorbar(scat_a_z, spacing='proportional')
	cb.set_label('JEF33 < ENDF4              JEF33 > ENDF4')
	# create the colorbar
	plt.title('Relative Difference in Percent between ENDF4 and JEF33 Fission Yields')
	plt.xlabel('A')
	plt.ylabel('Z')
	plt.grid('on')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.legend()

	# for i, txt in enumerate(label_list):
	# 	plt.annotate(txt, (a_list1[i], z_list1[i]))

	plt.show()

def get_openmc_xs_nucl_list():

	#path_to_xs_xml = os.environ['OPENMC_CROSS_SECTIONS']
	path_to_xs_xml = '/home/julien/Open-Burnup.dev/ENDFVIII_cross_sections/cross_sections.xml'
	MC_XS_nucl_list = []

	tree = ET.parse(path_to_xs_xml)
	root = tree.getroot()

	for child in root:
		if child.attrib['type'] == 'neutron':
			MC_XS_nucl_list.append(child.attrib['materials'])

	# Remove trouble makers		
	# For some reason, OpenMC can't find these nuclides in jeff lib at 800K
	MC_XS_nucl_list.remove('Cu63')
	MC_XS_nucl_list.remove('Cu65')
	MC_XS_nucl_list.remove('Mn55')

	try:
		MC_XS_nucl_list.remove('C0')
	except ValueError:
		pass
	try:
		MC_XS_nucl_list.remove('V0')
	except ValueError:
		pass
	try:
		MC_XS_nucl_list.remove('Zn0')
	except ValueError:
		pass

	return MC_XS_nucl_list


def plot_nuclide_chart_color_per_nuclear_data(decay_path, fy_path):

	fy_dict = d.read_fy_lib(fy_path)
	fy_unordered_keys = get_keylist_from_dict(fy_dict)
	fy_nucl_list = order_nuclide_per_z(fy_unordered_keys)

	decay_dict = d.read_decay_lib(decay_path)
	decay_unordered_keys = get_keylist_from_dict(decay_dict)
	decay_nucl_list = order_nuclide_per_z(decay_unordered_keys)

	MC_xs_nucl_name_list = get_openmc_xs_nucl_list()
	xs_nucl_list_name_list = mc_namelist_to_bu_namelist(MC_xs_nucl_name_list)
	xs_nucl_list = name_list_to_zamid_list(xs_nucl_list_name_list)
	xs_nucl_list = order_nuclide_per_z(xs_nucl_list)

	NAX_nucl_list = d.NAX_nucl_list

	orphan_NAX_z = []
	orphan_NAX_n = []

	only_decay_z = []
	only_decay_n = []
	only_decay_NAX_z = []
	only_decay_NAX_n = []	

	only_fy_z = []
	only_fy_n = []
	only_fy_NAX_z = []
	only_fy_NAX_n = []

	only_xs_z = []
	only_xs_n = []
	only_xs_NAX_z = []
	only_xs_NAX_n = []
	only_xs_zamid = []

	decay_xs_z = []
	decay_xs_n = []
	decay_xs_NAX_z = []
	decay_xs_NAX_n = []

	decay_fy_z = []
	decay_fy_n = []
	decay_fy_NAX_z = []
	decay_fy_NAX_n = []

	xs_fy_z = []
	xs_fy_n = []
	xs_fy_NAX_z = []
	xs_fy_NAX_n = []

	decay_xs_fy_z = []
	decay_xs_fy_n = []
	decay_xs_fy_NAX_z = []
	decay_xs_fy_NAX_n = []

	for decay_zamid in decay_nucl_list:

		if get_zamid_s(decay_zamid) == 1:
			continue

		in_fy = 'no'
		for fy_zamid in fy_nucl_list:
			if fy_zamid == decay_zamid:
				in_fy = 'yes'
				break

		in_xs = 'no'
		for xs_zamid in xs_nucl_list:
			if xs_zamid == decay_zamid:
				in_xs = 'yes'
				break

		in_NAX = 'no'
		for NAX_zamid in NAX_nucl_list:
			if NAX_zamid == decay_zamid:
				in_NAX = 'yes'
				break

		if in_fy =='no' and in_xs =='no' and in_NAX =='no':
			only_decay_z.append(get_zamid_z(decay_zamid))
			only_decay_n.append(get_zamid_n(decay_zamid))

		if in_fy =='no' and in_xs =='no' and in_NAX =='yes':
			only_decay_NAX_z.append(get_zamid_z(decay_zamid))
			only_decay_NAX_n.append(get_zamid_n(decay_zamid))

		if in_fy =='yes' and in_xs =='no' and in_NAX =='no':
			decay_fy_z.append(get_zamid_z(decay_zamid))
			decay_fy_n.append(get_zamid_n(decay_zamid))

		if in_fy =='yes' and in_xs =='no' and in_NAX =='yes':
			decay_fy_NAX_z.append(get_zamid_z(decay_zamid))
			decay_fy_NAX_n.append(get_zamid_n(decay_zamid))

		if in_fy =='no' and in_xs =='yes' and in_NAX =='no':
			decay_xs_z.append(get_zamid_z(decay_zamid))
			decay_xs_n.append(get_zamid_n(decay_zamid))

		if in_fy =='no' and in_xs =='yes' and in_NAX =='yes':
			decay_xs_NAX_z.append(get_zamid_z(decay_zamid))
			decay_xs_NAX_n.append(get_zamid_n(decay_zamid))

		if in_fy =='yes' and in_xs =='yes' and in_NAX =='no':
			decay_xs_fy_z.append(get_zamid_z(decay_zamid))
			decay_xs_fy_n.append(get_zamid_n(decay_zamid))

		if in_fy =='yes' and in_xs =='yes' and in_NAX =='yes':
			decay_xs_fy_NAX_z.append(get_zamid_z(decay_zamid))
			decay_xs_fy_NAX_n.append(get_zamid_n(decay_zamid))

	for fy_zamid in fy_nucl_list:

		if get_zamid_s(fy_zamid) == '1':
			continue

		in_xs = 'no'
		for xs_zamid in xs_nucl_list:
			if fy_zamid == xs_zamid:
				in_xs = 'yes'
				break

		in_decay = 'no'
		for decay_zamid in decay_nucl_list:
			if fy_zamid == decay_zamid:
				in_decay = 'yes'
				break

		in_NAX = 'no'
		for NAX_zamid in NAX_nucl_list:
			if fy_zamid == NAX_zamid:
				in_NAX = 'yes'
				break

		if in_xs =='no' and in_decay =='no' and in_NAX =='no':
			only_fy_z.append(get_zamid_z(fy_zamid))
			only_fy_n.append(get_zamid_n(fy_zamid))

		if in_xs =='no' and in_decay =='no' and in_NAX =='yes':
			only_fy_NAX_z.append(get_zamid_z(fy_zamid))
			only_fy_NAX_n.append(get_zamid_n(fy_zamid))

		if in_xs =='yes' and in_decay =='no' and in_NAX =='no':
			xs_fy_z.append(get_zamid_z(fy_zamid))
			xs_fy_n.append(get_zamid_n(fy_zamid))

		if in_xs =='yes' and in_decay =='no' and in_NAX =='yes':
			xs_fy_NAX_z.append(get_zamid_z(fy_zamid))
			xs_fy_NAX_n.append(get_zamid_n(fy_zamid))

	for xs_zamid in xs_nucl_list:

		if get_zamid_s(xs_zamid) == '1':
			continue

		in_decay = 'no'
		for decay_zamid in decay_nucl_list:
			if xs_zamid == decay_zamid:
				in_decay = 'yes'
				break

		in_fy = 'no'
		for fy_zamid in fy_nucl_list:
			if xs_zamid == fy_zamid:
				in_fy = 'yes'
				break

		in_NAX = 'no'
		for NAX_zamid in NAX_nucl_list:
			if xs_zamid == NAX_zamid:
				in_NAX = 'yes'
				break

		if in_fy =='no' and in_decay =='no' and in_NAX =='no':
			only_xs_z.append(get_zamid_z(xs_zamid))
			only_xs_n.append(get_zamid_n(xs_zamid))

		if in_fy =='no' and in_decay =='no' and in_NAX =='yes':
			only_xs_NAX_z.append(get_zamid_z(xs_zamid))
			only_xs_NAX_n.append(get_zamid_n(xs_zamid))

	for NAX_zamid in NAX_nucl_list:

		if get_zamid_s(NAX_zamid) == '1':
			continue

		in_decay = 'no'
		for decay_zamid in decay_nucl_list:
			if NAX_zamid == decay_zamid:
				in_decay = 'yes'
				break

		in_fy = 'no'
		for fy_zamid in fy_nucl_list:
			if NAX_zamid == fy_zamid:
				in_fy = 'yes'
				break

		in_xs = 'no'
		for xs_zamid in xs_nucl_list:
			if NAX_zamid == xs_zamid:
				in_xs = 'yes'
				break

		if in_fy =='no' and in_decay =='no' and in_xs =='no':
			orphan_NAX_z.append(get_zamid_z(NAX_zamid))
			orphan_NAX_n.append(get_zamid_n(NAX_zamid))




	# for i in range(len(decay_fy_n)):
	# 	if decay_fy_n[i] == 52 and decay_fy_z[i] ==41:
	# 		print (i)
	# 		print ('410920 in decay_fy')

	# for i in range(len(decay_xs_fy_NAX_n)):
	# 	if decay_xs_fy_NAX_n[i] == 52 and decay_xs_fy_NAX_z[i] ==41:
	# 		print (i)
	# 		print ('410920 in decay_xs_fy_NAX')


	size = 13

	# setup the plot
	plt.figure(1)

	plt.scatter(only_decay_n, only_decay_z, marker = 's', color='k', label = 'decay only', s = size)
	plt.scatter(only_decay_NAX_n, only_decay_NAX_z, marker = '.', color='k', label = 'decay only NAX', s = size)
	
	plt.scatter(only_xs_n, only_xs_z, marker = 's', color='gold', label = 'xs only', s = size)
	plt.scatter(only_xs_NAX_n, only_xs_NAX_z, marker = '.', color='gold', label = 'xs only NAX', s = size)
	
	plt.scatter(only_fy_n, only_fy_z, marker = 's', color='pink', label = 'fy only', s = size)
	plt.scatter(only_fy_NAX_n, only_fy_NAX_z, marker = '.', color='pink', label = 'fy only NAX', s = size)


	plt.scatter(decay_fy_n, decay_fy_z, marker = 's', color='b', label = 'decay & fy', s = size)
	plt.scatter(decay_fy_NAX_n, decay_fy_NAX_z, marker = '.', color='b', label = 'decay & fy NAX', s = size)


	plt.scatter(decay_xs_n, decay_xs_z, marker = 's', color='lime', label = 'decay & xs', s = size)
	plt.scatter(decay_xs_NAX_n, decay_xs_NAX_z, marker = '.', color='lime', label = 'decay & xs NAX', s = size)


	plt.scatter(xs_fy_n, xs_fy_z, marker = 's', color='fuchsia', label = 'xs & fy', s = size)
	plt.scatter(xs_fy_NAX_n, xs_fy_NAX_z, marker = '.', color='fuchsia', label = 'xs & fy NAX', s = size)


	plt.scatter(decay_xs_fy_n, decay_xs_fy_z, marker = 's', color='r', label = 'decay xs fy', s = size)
	plt.scatter(decay_xs_fy_NAX_n, decay_xs_fy_NAX_z, marker = '.', color='r', label = 'decay xs fy NAX', s = size)

	plt.scatter(orphan_NAX_n, orphan_NAX_z, marker = '.', color='grey', label = 'orphan NAX', s = size)


	plt.title('Nuclide colored per data availability\nNAX nuclides included')
	plt.xlabel('N')
	plt.ylabel('Z')
	plt.grid(b=True, which='major', color='k', linestyle='-')
	plt.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2, markevery=50)
	plt.minorticks_on()
	# plt.grid('on')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.legend()
	#plt.set_minor_frequency(2)

	# for i, txt in enumerate(label_list):
	# 	plt.annotate(txt, (a_list1[i], z_list1[i]))

	plt.show()

def plot_compare_two_nuclear_data_on_nuclide_chart(decay_path1, fy_path1,decay_path2, fy_path2):

	fy_dict1 = d.read_fy_lib(fy_path1)
	fy_unordered_keys1 = get_keylist_from_dict(fy_dict1)
	fy_nucl_list1 = order_nuclide_per_z(fy_unordered_keys1)

	decay_dict1 = d.read_decay_lib(decay_path1)
	decay_unordered_keys1 = get_keylist_from_dict(decay_dict1)
	decay_nucl_list1 = order_nuclide_per_z(decay_unordered_keys1)

	fy_dict2 = d.read_fy_lib(fy_path2)
	fy_unordered_keys2 = get_keylist_from_dict(fy_dict2)
	fy_nucl_list2 = order_nuclide_per_z(fy_unordered_keys2)

	decay_dict2 = d.read_decay_lib(decay_path2)
	decay_unordered_keys2 = get_keylist_from_dict(decay_dict2)
	decay_nucl_list2 = order_nuclide_per_z(decay_unordered_keys2)

	MC_xs_nucl_name_list = get_openmc_xs_nucl_list()
	xs_nucl_list_name_list = mc_namelist_to_bu_namelist(MC_xs_nucl_name_list)
	xs_nucl_list = name_list_to_zamid_list(xs_nucl_list_name_list)
	xs_nucl_list = order_nuclide_per_z(xs_nucl_list)

	NAX_nucl_list = d.NAX_nucl_list

	decay_z1 = []
	decay_n1 = []
	fy_z1 = []
	fy_n1 = []
	decay_z2 = []
	decay_n2 = []
	fy_z2 = []
	fy_n2 = []
	xs_z2 = []
	xs_n2 = []
	NAX_z2 = []
	NAX_n2 = []

	for decay_zamid1 in decay_nucl_list1:

		if get_zamid_s(decay_zamid1) == 1:
			continue

		decay_z1.append(get_zamid_z(decay_zamid1))
		decay_n1.append(get_zamid_n(decay_zamid1))


	for fy_zamid1 in fy_nucl_list1:

		if get_zamid_s(fy_zamid1) == '1':
			continue

		fy_z1.append(get_zamid_z(fy_zamid1))
		fy_n1.append(get_zamid_n(fy_zamid1))

	for decay_zamid2 in decay_nucl_list2:

		if get_zamid_s(decay_zamid2) == 1:
			continue

		decay_z2.append(get_zamid_z(decay_zamid2))
		decay_n2.append(get_zamid_n(decay_zamid2))


	for fy_zamid2 in fy_nucl_list2:

		if get_zamid_s(fy_zamid2) == '1':
			continue

		fy_z2.append(get_zamid_z(fy_zamid2))
		fy_n2.append(get_zamid_n(fy_zamid2))

	for xs_zamid2 in xs_nucl_list:

		if get_zamid_s(xs_zamid2) == '1':
			continue

		xs_z2.append(get_zamid_z(xs_zamid2))
		xs_n2.append(get_zamid_n(xs_zamid2))

	for NAX_zamid2 in NAX_nucl_list:

		if get_zamid_s(NAX_zamid2) == '1':
			continue

		NAX_z2.append(get_zamid_z(NAX_zamid2))
		NAX_n2.append(get_zamid_n(NAX_zamid2))



	size = 13

	# setup the plot
	plt.figure(1)

	# Plot nuclear data 1
	plt.scatter(decay_n1, decay_z1, marker = 's', color='k', label = 'Full', s = size)
	plt.scatter(fy_n1, fy_z1, marker = 's', color='k', s = size)
	plt.scatter(xs_n2, xs_z2, marker = 's', color='k',s = size)
	#plt.scatter(NAX_n2, NAX_z2, marker = 's', color='k', s = size)

	# Plot nuclear data 2	
	plt.scatter(decay_n2, decay_z2, marker = 's', color='r', label = 'Reduced', s = size)
	plt.scatter(fy_n2, fy_z2, marker = 's', color='r', s = size)
	plt.scatter(xs_n2, xs_z2, marker = 's', color='r',s = size)
	#plt.scatter(NAX_n2, NAX_z2, marker = 's', color='b', s = size, label = 'Archaeology')


	plt.xlabel('Number of neutrons', fontsize=16)
	plt.ylabel('Number of protons', fontsize=16)
	plt.grid(b=True, which='major', color='k', linestyle='-')
	plt.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2, markevery=50)
	plt.minorticks_on()
	plt.tick_params(labelsize=15)
	# plt.grid('on')
	plt.gca().set_aspect('equal', adjustable='box')
	plt.legend(prop={'size': 15})
	#plt.set_minor_frequency(2)

	# for i, txt in enumerate(label_list):
	# 	plt.annotate(txt, (a_list1[i], z_list1[i]))

	plt.show()


def convert_spectrum_to_janis_weighting_format(path_to_simulation, bucell, BU):

	path = path_to_simulation +'/output_summary/{}_flux_spectrum'.format(bucell)
	spectrum_file = open(path)

	lines = spectrum_file.readlines()

	#Energy bins
	energy_bins = lines[0].split()[1:]

	#Energy mid points
	energy_mid_points = lines[2].split()[1:]

	# The data starts at the 7th line
	for line in lines[6:]:
		line = line.split()
		if float(line[1]) == BU:
			spectrum_lethargy = line[3:]

	spectrum = [float(x)/float(y) for x,y in zip(spectrum_lethargy, energy_mid_points)]

	energy_bin_file = open('energy_bin.gst', 'w')

	txt_bin = 'neutron group structure......anl 299 group\n'
	for i in range(len(energy_bins)-1):
		Emin = energy_bins[i]
		Emax = energy_bins[i+1]
		txt_bin += '{} {} {}\n'.format(i+1, Emin, Emax)

	energy_bin_file.write(txt_bin)
	energy_bin_file.close()

	spectrum_file = open('spectrum.txt', 'w')
	txt_spectrum = ''
	for i in range(len(spectrum)):
		E = energy_mid_points[i]
		flux = spectrum[i]
		txt_spectrum += '{} {}\n'.format(E, flux)

	spectrum_file.write(txt_spectrum)
	spectrum_file.close()


def get_zamid_natural_abundance(zamid):

	name_old_format = zamid_to_name(zamid)
	name_new_format = onix_name_to_openmc_name(name_old_format)
	nat_abun_dict = d.NATURAL_ABUNDANCE
	if name_new_format in nat_abun_dict:
		nat_abun = d.NATURAL_ABUNDANCE[name_new_format]
	else:
		nat_abun = 0.0

	return nat_abun

def get_name_natural_abundance(name):

	name_new_format = onix_name_to_openmc_name(name)
	nat_abun_dict = d.NATURAL_ABUNDANCE
	if name_new_format in nat_abun_dict:
		nat_abun = d.NATURAL_ABUNDANCE[name_new_format]
	else:
		nat_abun = 0.0

	return nat_abun

def find_zamid_precursor(zamid, reaction):

	xs_prod_fromS_toS = d.xs_prod_fromS_toS
	zamid_shift = xs_prod_fromS_toS[reaction]
	precursor_zamid = int(zamid) - 10000*zamid_shift[0] - 10*zamid_shift[1] - zamid_shift[2]

	return str(precursor_zamid)


def smooth_triangle(data, degree, dropVals=False):
    triangle=np.array(list(range(degree)) + [degree] + list(range(degree)[::-1])) + 1
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(sum(point)/sum(triangle))
    if dropVals:
        return smoothed
    #smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    # while len(smoothed) < len(data):
    #     smoothed.append(smoothed[-1])

        
    return smoothed

def moving_average(data, window):

	weights = np.repeat(1.0, window)/window
	smas = np.convolve(data, weights, 'valid')

	print (smas)
	print (data[:int(window/2)])

	return list(data[:int(window/2)]) + list(smas) + list(data[-int(window/2):-1])

def read_BUCell_vol(path, cell):

	path_to_parameters = path +'/system_parameters'

	parameter_file = open(path_to_parameters)

	lines = parameter_file.readlines()

	search = 'BuCell'

	for line in lines:
		print (line)
		if line == 'BuCell {}\n'.format(cell):
			search = 'volume'

		if search == 'volume':
			if line.split()[0] == 'Volume':
				vol = line.split()[3]
				break

	return float(vol)


class Empty_argument(Exception):
	"""Raise when the user calls decay_halflife_conv without entering any argument """
	pass