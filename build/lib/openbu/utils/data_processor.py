from .functions import *
import ast
import re
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
import openbu.data as d
#import pylab

def read_nuclide_reac_rank(nuclide, step, path):

	file = open(path, 'r')
	lines = file.readlines()

	zamid = name_to_zamid(nuclide)

	search = 'nuclide'
	for i in range(0,len(lines)):
		l = lines[i]
		if l == ' ==={}({}) ===\n'.format(nuclide, zamid):
			search = 'step'
		elif l == 'STEP {}\n'.format(step) and search == 'step':
			data = lines[i+2]
			break

	data = ast.literal_eval(data)


	destruction = {}

	dest_total = 0
	for tuples in data:
		if '-' not in tuples[0]:
			if tuples[1] == 0.0:
				continue
			dest_total += tuples[1]


	for tuples in data:
		if '-' not in tuples[0]:
			if tuples[1] == 0.0:
				continue
			destruction[tuples[0]] = [tuples[1], tuples[1]/dest_total] # val and percent


	production = {}

	prod_total = 0
	for tuples in data:
		if '-' in tuples[0]:
			reaction_val = tuples[1]
			if reaction_val == 0.0:
				continue	
			prod_total += reaction_val

	for tuples in data:
		if '-' in tuples[0]:
			parent = tuples[0].split()[0]
			reaction = tuples[0].split()[1]
			reaction_val = tuples[1]
			if reaction_val == 0.0:
				continue	
			production[parent] = [reaction, reaction_val, reaction_val/prod_total]

	return [destruction, production]

def plot_bucell_nuclide_network(nuclide, step, path, threshold):

	file = open(path, 'r')
	lines = file.readlines()

	zamid = name_to_zamid(nuclide)

	search = 'nuclide'
	for i in range(0,len(lines)):
		l = lines[i]
		if l == ' ==={}({}) ===\n'.format(nuclide, zamid):
			search = 'step'
		elif l == 'STEP {}\n'.format(step) and search == 'step':
			data = lines[i+2]
			break

	data = ast.literal_eval(data)


	destruction = {}

	dest_total = 0
	for tuples in data:
		if '-' not in tuples[0]:
			# if tuples[1] == 0.0:
			# 	continue
			if tuples[1] < threshold:
				continue
			dest_total += tuples[1]


	for tuples in data:
		if '-' not in tuples[0]:
			# if tuples[1] == 0.0:
			# 	continue
			if tuples[1] < threshold:
				continue
			destruction[tuples[0]] = [tuples[1], tuples[1]/dest_total] # val and percent


	production = {}

	prod_total = 0
	for tuples in data:
		if '-' in tuples[0]:
			reaction_val = tuples[1]
			# if reaction_val == 0.0:
			# 	continue	
			if reaction_val < threshold:
				continue	
			prod_total += reaction_val

	for tuples in data:
		if '-' in tuples[0]:
			parent = tuples[0].split()[0]
			reaction = tuples[0].split()[1]
			reaction_val = tuples[1]
			# if reaction_val == 0.0:
			# 	continue	
			if reaction_val < threshold:
				continue
			production[parent] = [reaction, reaction_val, reaction_val/prod_total]


	G = nx.MultiDiGraph()

	for parent in production:
		label = '{}\n{:.2E}[{:.2%}]'.format(production[parent][0], production[parent][1], production[parent][2])
		G.add_edge(parent, nuclide, label = label, length = 10)

	for edge in destruction:
		G.add_edge(nuclide, edge, label = '{:.2E}[{:.2%}]'.format(destruction[edge][0],destruction[edge][1] ), length = 10)

	# Get target nuclide index
	index = 0
	for node in G.nodes():
		if node == nuclide:
			break
		index += 1


	edges = G.edges()
	edge_labels = []

	edge_labels=dict([((u,v,),d['label'])
	                 for u,v,d in G.edges(data=True)])
	node_color = []
	for node in G.nodes():
		if node == nuclide:
			node_color.append('mediumaquamarine')
		if node in production:
			node_color.append('lightskyblue')
		if node in destruction:
			node_color.append('darksalmon')

	edges = G.edges()
	# edge_weights = [G[u][v]['weight'] for u,v in edges]
	# red_edges = [('C','D'),('D','A')]
	# edge_colors = ['black' if not edge in red_edges else 'red' for edge in G.edges()]
	pos=nx.circular_layout(G, scale = 2)
	pos[nuclide] = np.array([0, 0])
	nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)
	nx.draw_networkx_labels(G,pos)
	nx.draw_networkx_edges(G,pos, edges = edges)
	nx.draw(G,pos, node_size=4000, node_color = node_color, fontsize = 6)

	plt.show()

def plot_nuclide_dens(bucell, nuclide):

	path = os.getcwd() +'/{}_subdens'.format(bucell)

	time_seq = read_time_seq(path)
	dens_seq = read_dens(path)

	plt.figure(1)
	plt.plot(time_seq, dens_seq, color = 'orange', marker = 'o')
	plt.xlabel('Time [day]')
	plt.ylabel('Density [atm/cm3]')
	plt.grid()
	plt.title('{} {} density evolution'.format(bucell, nuclide))

	plt.show()

def plot_nuclide_dens_from_path(bucell, nuclide, path):

	time_seq = read_time_seq(path)
	dens_seq = read_dens(path)

	plt.figure(1)
	plt.plot(time_seq, dens_seq, color = 'orange', marker = 'o')
	plt.xlabel('Time [day]')
	plt.ylabel('Density [atm/cm3]')
	plt.grid()
	plt.title('{} {} density evolution'.format(bucell, nuclide))

	plt.show()


def plot_xs_time_evolution(bucell, nuclide, xs_name):

	path = os.getcwd() +'/{}_xs_lib'.format(bucell)

	time_seq = read_time_seq(path)
	xs_seq = read_xs(nuclide, xs_name, path)

	marker_list = ['x', '+', 'o', '*', '^', 's']
	color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

	plt.figure(1)
	plt.plot(time_seq, xs_seq, color = 'teal', marker = 'o')

	plt.xlabel('Time in Days')
	plt.ylabel('Eff. XS in barn')
	plt.legend()
	plt.grid()
	plt.title('{} {} effective cross sections evolution'.format(bucell, nuclide))

	plt.show()

def plot_xs_bu_evolution(bucell_list, nuclide, xs_name):

	index = 0
	for bucell in bucell_list:
		path = os.getcwd() +'/{}_xs_lib'.format(bucell)

		xs_seq = read_xs_seq(nuclide, xs_name, path)

		marker_list = ['x', '+', 'o', '*', '^', 's']
		color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

		plt.figure(1)
		plt.plot(bu_seq, xs_seq, color = color_list[index], marker = marker_list[index], label = bucell)

		index += 1

	bu_seq = read_bu_seq(path)

	plt.xlabel('BU in MWd/kg')
	plt.ylabel('Eff. XS in barn')
	plt.legend()
	plt.grid()
	plt.title('{} {} effective cross sections evolution'.format(nuclide, xs_name))

	plt.show()


def plot_xs_time_evolution_from_path(bucell, nuclide, xs_name, path):

	time_seq = read_time_seq()
	xs_seq = read_xs_seq(nuclide, xs_name, path)

	marker_list = ['x', '+', 'o', '*', '^', 's']
	color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

	plt.figure(1)
	plt.plot(time_seq, xs_seq, color = 'teal', marker = 'o')

	plt.xlabel('Time in Days')
	plt.ylabel('Eff. XS in barn')
	plt.legend()
	plt.grid()
	plt.title('{} {} effective cross sections evolution'.format(bucell, nuclide))

	plt.show()

def plot_xs_bu_evolution_from_path(bucell_list, nuclide, xs_name, path):

	index = 0
	for bucell in bucell_list:
		path_xs = path +'/{}_xs_lib'.format(bucell)

		xs_seq = read_xs_seq(nuclide, xs_name,path_xs)

		marker_list = ['x', '+', 'o', '*', '^', 's']
		color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

		plt.figure(1)
		bu_seq = read_bu_seq(path_xs)
		plt.plot(bu_seq, xs_seq, color = color_list[index], marker = marker_list[index], label = bucell)

		index += 1


	plt.xlabel('BU in MWd/kg')
	plt.ylabel('Eff. XS in barn')
	plt.legend()
	plt.grid()
	plt.title('{} {} effective cross sections evolution'.format(nuclide, xs_name))

	plt.show()

# Compare xs evolution for the same nuclide, for various xs for various different runs
def compare_xs_bu_evolution_from_path(bucell, nuclide, xs_name_list, path_list, name_list):

	bu_seq_list = []
	xs_seq_list = []
	for path in path_list:

		bu_seq_list.append(read_bu_seq(path))
		xs_seq_list.append([])

		for xs in xs_name_list:
			xs_seq_list[-1].append(read_xs_seq(nuclide, xs, path))


	marker_list = ['x', '+', 'o', '*', '^', 's']
	color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']
	
	# plt.figure(1)

	# for i in range(len(path_list)):

	# 	for j in range(len(xs_name_list)):

	# 		plt.plot(bu_seq_list[i], xs_seq_list[i][j], label = '{} {}'.format(name_list[i], xs_name_list[j]))

	# Three subplots sharing both x/y axes
	f, ax_tuple = plt.subplots(len(xs_name_list), sharex=True)
	ax_tuple[0].set_title('{} {} {} effective cross sections evolution'.format(bucell, nuclide, xs_name_list))
	for i in range(len(xs_name_list)):
		for j in range(len(path_list)):
			ax_tuple[i].plot(bu_seq_list[j], xs_seq_list[j][i], label = name_list[j])
		ax_tuple[i].set_ylabel('{} Eff. XS [barn]'.format(xs_name_list[i]))
		ax_tuple[i].grid()
		ax_tuple[i].set_xlabel('BU [MWd/kg]')
	# Fine-tune figure; make subplots close to each other and hide x ticks for
	# all but bottom plot.
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
	plt.legend()
	plt.show()


def plot_flux(bucell):

	path = os.getcwd() +'/{}_xs_lib'.format(bucell)

	time_seq = read_time_seq(path)
	flux_seq = read_flux(path)

	plt.figure(1)
	plt.step(time_seq, flux_seq, where = 'pre', color = 'blue')
	plt.xlabel('Time [day]')
	plt.ylabel('Flux [neutron/cm3.s-1]')
	plt.grid()
	plt.title('{} neutron flux evolution'.format(bucell))

	plt.show()

def plot_flux_from_path(bucell, path):

	time_seq = read_time_seq(path)
	flux_seq = read_flux(path)

	plt.figure(1)
	plt.step(time_seq, flux_seq, where = 'pre', color = 'blue')
	plt.xlabel('Time [day]')
	plt.ylabel('Flux [neutron/cm3.s-1]')
	plt.grid()
	plt.title('{} neutron flux evolution'.format(bucell))

	plt.show()

def plot_flux_spectrum_bu_evolution_from_path(bucell_list, steps_list, path):

	bucell_index = 0
	for bucell in bucell_list:
		path_flux_spectrum = path +'/{}_flux_spectrum'.format(bucell)

		flux_spectrum_list = read_flux_spectrum(path_flux_spectrum, steps_list)

		energy_mid_points = read_energy_mid_points(path_flux_spectrum)

		marker_list = ['x', '+', 'o', '*', '^', 's']
		color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

		plt.figure(bucell_index)
		index = 0
		for flux_spectrum in flux_spectrum_list:
			plt.plot(energy_mid_points, flux_spectrum, color = color_list[index], marker = marker_list[index], label=steps_list[index])
			index += 1


		plt.xlabel('BU in MWd/kg')
		plt.ylabel('Eff. XS in barn')
		plt.legend()
		plt.grid()
		plt.title('neutron spectrum in cell {} for steps {} '.format(bucell, steps_list))
		
	plt.show()

def plot_lethargy_spectrum_bu_evolution_from_path(bucell_list, steps_list, path):

	bucell_index = 0
	for bucell in bucell_list:
		path_flux_spectrum = path +'/{}_flux_spectrum'.format(bucell)

		flux_spectrum_list = read_flux_spectrum(path_flux_spectrum, steps_list)

		energy_mid_points = read_energy_mid_points(path_flux_spectrum)
		energy_bin_length = read_energy_bin_length(path_flux_spectrum)

		marker_list = ['x', '+', 'o', '*', '^', 's']
		color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

		plt.figure(bucell_index)
		index = 0
		for flux_spectrum in flux_spectrum_list:
			lethargy_spectrum = [x*y/z for x,y,z in zip(flux_spectrum, energy_mid_points, energy_bin_length)]
			plt.plot(energy_mid_points, lethargy_spectrum, color = color_list[index], marker = marker_list[index], label=steps_list[index])
			index += 1


		plt.xlabel('BU in MWd/kg')
		plt.ylabel('Eff. XS in barn')
		plt.legend()
		plt.xscale('log')
		plt.yscale('log')
		plt.grid()
		plt.title('neutron spectrum in cell {} for steps {} '.format(bucell, steps_list))
		
	plt.show()

def plot_xs_dens_flux(bucell, xs_nuclide, xs_name, dens_nuclide, xs_path, dens_path):

	time_subseq = read_time_seq(dens_path)
	time_seq = read_time_seq(xs_path)
	xs_seq = read_xs_seq(xs_nuclide, xs_name, xs_path)
	dens_subseq = read_dens(dens_nuclide, dens_path)
	flux_subseq = read_flux(dens_path)

	# Three subplots sharing both x/y axes
	f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
	ax1.plot(time_seq, xs_seq, color = 'teal', marker = 'o')
	ax1.set_ylabel('{}\n{} Eff. XS [barn]'.format(xs_nuclide, xs_name))
	ax1.set_title('{} cell'.format(bucell))
	ax1.grid()
	ax2.plot(time_subseq, dens_subseq, color = 'orange', marker = 'o')
	ax2.set_ylabel('{}\nDensity [atm/barn-cm]'.format(dens_nuclide))
	ax2.grid()
	ax3.step(time_subseq, flux_subseq, color = 'blue')
	ax3.set_ylabel('Neutron flux [neutron/cm2-s]')
	ax3.set_xlabel('Time [day]')
	ax3.grid()
	# Fine-tune figure; make subplots close to each other and hide x ticks for
	# all but bottom plot.
	f.subplots_adjust(hspace=0)
	plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
	plt.show()


def read_time_seq(path):

	time_file = open(path, 'r')

	lines = time_file.readlines()

	# Find and store time
	for line in lines:
		if line != '\n':
			if line.split()[0] == 'TIME':
				time_seq = [float(x) for x in line.split()[1:]]
				break

	return time_seq

def get_step_time_length_seq(path):

	time_seq = read_time_seq(path)

	if time_seq[0] == 0:
		step_time_length = [x-y for x,y in zip(time_seq[1:], time_seq[:-1])]
	else:
		step_time_length = [time_seq[0]] +  [x-y for x,y in zip(time_seq[1:], time_seq[:-1])]
	
	return step_time_length


def read_bu_seq(path):

	bu_file = open(path, 'r')

	lines = bu_file.readlines()

	# Find and store time
	for line in lines:
		if line != '\n':
			if line.split()[0] == 'SYSTEM-BU':
				bu_seq = [float(x) for x in line.split()[1:]]
				break

	return bu_seq


def read_flux(path):

	flux_file = open(path, 'r')

	lines = flux_file.readlines()

	for line in lines:
		if line != '\n':
			if line.split()[0] == 'FLUX':
				flux_seq = [float(x) for x in line.split()[1:]]
				break
	# Add an initial 0 for the flux to have same len for flux_seq and time_seq
	flux_seq = [0] + flux_seq

	return flux_seq

def read_flux_subseq(path):

	flux_file = open(path, 'r')

	lines = flux_file.readlines()

	for line in lines:
		if line != '\n':
			if line.split()[0] == 'FLUX':
				flux_seq = [float(x) for x in line.split()[1:]]
				break

	return flux_seq

def get_fluence_seq(path, cell):

	dens_file = path +'/{}_dens'.format(cell)

	flux_seq = read_flux(dens_file)
	time_seq = read_time_seq(dens_file)

	#Reminder: flux starts with 0

	fluence_seq = [0]
	pre_time = 0
	for i in range(1, len(time_seq)):
		time = time_seq[i]
		flux = flux_seq[i]
		time_int = (time-pre_time)*24*3600
		fluence_seq.append(time_int*flux + fluence_seq[i-1])
		pre_time = time

	return fluence_seq

def get_fluence_subseq(path, cell):

	subdens_file = path +'/{}_subdens'.format(cell)

	flux_subseq = read_flux(subdens_file)
	time_subseq = read_time_seq(subdens_file)

	#Reminder: flux starts with 0

	fluence_subseq = [0]
	pre_time = 0
	for i in range(1, len(time_subseq)):
		time = time_subseq[i]
		flux = flux_subseq[i]
		time_int = (time-pre_time)*24*3600
		fluence_subseq.append(time_int*flux + fluence_subseq[i-1])
		pre_time = time

	return fluence_subseq

# This method build the fluence sequence until a final time
# Warning this function is only valid for non-fissile material as the flux
# is the same within each step
def get_fluence_seq_until_time(path, cell, final_time):

	final_step = find_step_from_time(path, cell, final_time)
	extra_fluence = get_extra_fluence_from_time(path, cell, final_time)
	fluence_seq = get_fluence_seq(path, cell)

	fluence_seq_until_time = fluence_seq[:final_step+1] + [fluence_seq[final_step] + extra_fluence]

	return fluence_seq_until_time

def get_fluence_subseq_until_time(path, cell, final_time):

	final_substep =find_substep_from_time(path, cell, final_time)
	extra_fluence = get_extra_subfluence_from_time(path, cell, final_time)
	fluence_subseq = get_fluence_subseq(path, cell)

	fluence_subseq_until_time = fluence_subseq[:final_substep+1] + [fluence_subseq[final_substep] + extra_fluence]

	return fluence_subseq_until_time

# This function calculate the additional fluence from the previous time point
# to where the final time is set
def get_extra_fluence_from_time(path, cell, time):

	step = find_step_from_time(path, cell, time)
	dens_file = path +'/{}_dens'.format(cell)
	# flux seq has an added zero at the beginning of the arry
	flux = read_flux(dens_file)[step+1]
	previous_time_point = read_time_seq(dens_file)[step]

	time_int = time-previous_time_point

	return flux*time_int*24*3600

# This function calculate the additional fluence from the previous time point
# to where the final time is set
# From subdens file
def get_extra_subfluence_from_time(path, cell, time):

	substep = find_substep_from_time(path, cell, time)
	subdens_file = path +'/{}_subdens'.format(cell)
	# flux seq has an added zero at the beginning of the arry

	flux_subseq = read_flux_subseq(subdens_file)
	flux= flux_subseq[substep]
	previous_time_point = read_time_seq(subdens_file)[substep]

	time_int = time-previous_time_point

	return flux*time_int*24*3600

def find_step_from_time(path, cell, time):

	dens_file = path +'/{}_dens'.format(cell)
	time_seq = read_time_seq(dens_file)

	step = 0
	for t in time_seq[1:]:
		if time <= t:
			break
		step += 1

	return step

def find_substep_from_time(path, cell, time):

	subdens_file = path +'/{}_subdens'.format(cell)
	time_seq = read_time_seq(subdens_file)

	substep = 0
	for t in time_seq[1:]:
		if time <= t:
			break
		substep += 1

	return substep


def get_step_fluence_length(path, cell):

	fluence_seq = get_fluence_seq(path, cell)

	step_fluence_length = [x-y for x,y in zip(fluence_seq[1:], fluence_seq[:-1])]
	
	return step_fluence_length

def read_flux_spectrum(path, steps_list):

	flux_spectrum_file = open(path, 'r')

	lines = flux_spectrum_file.readlines()

	flux_spectrum_list = []

	step_count = 0
	for line in lines[6:]: # The flux spectrum data always start at the 6th line
		if step_count in steps_list:
			flux_spectrum_list.append([float(x) for x in line.split()[3:]])
		step_count += 1

	return flux_spectrum_list

def read_energy_mid_points(path):

	flux_spectrum_file = open(path, 'r')

	lines = flux_spectrum_file.readlines()

	energy_mid_points = [float(x) for x in lines[2].split()[1:]]

	return energy_mid_points

def read_energy_bin_length(path):

	flux_spectrum_file = open(path, 'r')

	lines = flux_spectrum_file.readlines()

	energy_bin_length = [float(x) for x in lines[1].split()[1:]]

	return energy_bin_length

def read_dens(nuclide, path):

	zamid = name_to_zamid(nuclide)

	dens_seq = []

	dens_file = open(path, 'r')

	lines = dens_file.readlines()

	for line in lines:
		if line != '\n':
			if line.split()[0] == zamid:
				dens_seq = [float(x) for x in line.split()[1:]]
				break

	return dens_seq

# cumulative dens
def get_cum_dens(nuclide, path):

	dens_seq = read_dens(nuclide, path)

	cum_dens_seq = []
	cum_dens_seq.append(dens_seq[0])
	for i in range(1, len(dens_seq)):
		cum_dens_seq.append(dens_seq[i] + cum_dens_seq[i-1])

	return cum_dens_seq

def convert_dens_seq_to_cum_dens_seq(dens_seq):

	cum_dens_seq = []
	cum_dens_seq.append(dens_seq[0])
	for i in range(1, len(dens_seq)):
		cum_dens_seq.append(dens_seq[i] + cum_dens_seq[i-1])

	return cum_dens_seq

def get_nucl_atomic_mass(nucl):

	zamid = name_to_zamid(nucl)
	zaid = zamid[:-1]
	if zaid in d.default_atm_mass_lib:
		M = d.default_atm_mass_lib[zaid]
	else:
		M = int(get_zamid_a(zamid))

	return M

# calculate total mass density at certain step
def get_total_mass_density(path, cell, step):

	nucl_name_list = read_dens_nucl(path, cell)
	dens_path = path+'/{}_dens'.format(cell)
	dens_file = open(dens_path )

	NA = d.NA

	total_mass_density = 0

	for nucl in nucl_name_list:
		dens = read_dens(nucl, dens_path)[step]
		M = get_nucl_atomic_mass(nucl)
		mass_density = dens*(M/NA)
		total_mass_density += mass_density

	return total_mass_density


def get_pu_subseq_mat(path, cell, EFPD):

	final_substep =find_substep_from_time(path, cell, EFPD)

	path = path +'/{}_subdens'.format(cell)

	name_list = d.Pu_isotopes_name
	
	time_subseq = read_time_seq(path)

	t_before = time_subseq[final_substep]
	t_after = time_subseq[final_substep+1]

	pu_subseq_mat = []
	for name in name_list:
		dens_subseq = read_dens(name, path)
		dens_subseq_until_time = dens_subseq[:final_substep+1]
		dens_before = dens_subseq[final_substep]
		dens_after = dens_subseq[final_substep+1]
		pair1 = [t_before, dens_before]
		pair2 = [t_after, dens_after]
		interpolated_dens = interpolation_between_two_points(pair1, pair2, EFPD)
		dens_subseq = dens_subseq_until_time + [interpolated_dens]

		pu_subseq_mat.append(dens_subseq)

	return pu_subseq_mat

# cumulative plutonium production
def get_cum_pu_subseq_mat(path, cell, EFPD):

	path = path +'/{}_subdens'.format(cell)

	name_list = ['Pu-238', 'Pu-239', 'Pu-240', 'Pu-241', 'Pu-242', 'Pu-243']

	final_substep =find_substep_from_time(path, cell, EFPD)
	
	time_subseq = read_time_seq(path)

	t_before = time_subseq[final_substep]
	t_after = time_subseq[final_substep+1]

	cum_pu_subseq_mat = []
	for name in name_list:
		dens_subseq = read_dens(name, path)
		dens_subseq_until_time = dens_subseq[:final_substep+1]
		dens_before = dens_subseq[final_substep]
		dens_after = dens_subseq[final_substep+1]
		pair1 = [t_before, dens_before]
		pair2 = [t_after, dens_after]
		interpolated_dens = interpolation_between_two_points(pair1, pair2, EFPD)
		dens_subseq = dens_subseq_until_time + [interpolated_dens]
		cum_dens_subseq = convert_dens_seq_to_cum_dens_seq(dens_sbuseq)

		cum_pu_subseq_mat.append(cum_dens_subseq)

	return cum_pu_subseq_mat


# linear interpolation between two points
def interpolation_between_two_points(pair1, pair2, x):

	a = (pair2[1] - pair1[1])/(pair2[0] - pair1[0])
	b = (pair1[1]*pair2[0] - pair1[0]*pair2[1])/(pair2[0] - pair1[0])

	y = a*x+b

	return y


def read_xs_seq(nuclide, xs_name, path, cell):

	path = path + '/{}_xs_lib'.format(cell)

	zamid = name_to_zamid(nuclide)
	xs_name_found = 'no'

	xs_file = open(path, 'r')

	lines = xs_file.readlines()

	# Search for the line
	line_index = 0
	for line in lines:
		if line != '\n':
			if line.split()[0] == nuclide:
				break
		line_index += 1

	# First line needs to be treated differently
	if lines[line_index].split()[2] == xs_name:
		xs_name_found = 'yes'
		xs_seq= [float(x) for x in lines[line_index].split()[3:]]

	xs_loop = line_index+1

	while lines[xs_loop].split()[0] == zamid:
		if lines[xs_loop].split()[1] == xs_name:
			xs_name_found = 'yes'
			xs_seq = [float(x) for x in lines[xs_loop].split()[2:]]
		xs_loop += 1

	if xs_name_found == 'no':
		raise xs_name_not_found("nuclide {} has no data for cross section {}".format(nuclide, xs_name))
	else:
		return xs_seq

def get_time_averaged_xs(nuclide, xs_name, path, cell):

	xs_seq = read_xs_seq(nuclide, xs_name, path, cell)
	xs_lib_path = path + '/{}_dens'.format(cell)
	time_seq = read_time_seq(xs_lib_path)
	tot_time = time_seq[-1]

	av_xs = 0

	for i in range(len(xs_seq)):

		xs = xs_seq[i]
		time_bos = time_seq[i]
		time_eos = time_seq[i+1]
		time_coeff = (time_eos - time_bos)/tot_time
		av_xs += xs*time_coeff

	return av_xs

def get_time_averaged_flux(path, cell):

	xs_lib_path = path + '/{}_dens'.format(cell)
	flux_seq = read_flux(xs_lib_path)
	time_seq = read_time_seq(xs_lib_path)
	tot_time = time_seq[-1]

	av_flux = 0

	for i in range(len(flux_seq)-1):

		flux = flux_seq[i+1]
		time_bos = time_seq[i]
		time_eos = time_seq[i+1]
		time_coeff = (time_eos - time_bos)/tot_time
		av_flux += flux*time_coeff

	return av_flux


def get_tot_xs(nuclide, path, cell):

	xs_name_list = ['fission','(n,gamma)','(n,2n)','(n,3n)','(n,p)','(n,a)','(n,gamma)X']

	tot_xs_seq = []
	i = 0
	for xs_name in xs_name_list:
		try:
			xs_seq = read_xs_seq(nuclide, xs_name, path, cell)
		except xs_name_not_found:
			continue
		if i == 0:
			tot_xs_seq = xs_seq
		else:
			tot_xs_seq = [x+y for x,y in zip(tot_xs_seq,xs_seq)]

		i += 1

	return tot_xs_seq




# This method list all nuclides that are present in xs lib
def read_xs_nucl(xs_lib_file_path):

	xs_lib_file = open(xs_lib_file_path)

	lines = xs_lib_file.readlines()

	nucl_name_list = []
	for line in lines:
		if line == '\n':
			continue
		line = line.split()
		if line[0].split('-')[0] in d.nuc_name_dic:
			nucl_name_list.append(line[0])

	xs_lib_file.close()


	return nucl_name_list

# make a list of all nuclide present in the dens_file
def read_dens_nucl(path, cell):

	dens_path = path + '/{}_dens'.format(cell)
	dens_file = open(dens_path, 'r')

	lines = dens_file.readlines()

	nucl_zamid_list = []
	for line in lines[7:]:
		line = line.split()
		nucl_zamid_list.append(line[0])

	dens_file.close()

	nucl_name_list = zamid_list_to_name_list(nucl_zamid_list)

	return nucl_name_list


def rank_nuclide_per_dens(bucell, step_list, path):

	dens_path = path +'/output_summary' + '/{}_dens'.format(bucell)
	dens_file = open(dens_path, 'r')

	lines = dens_file.readlines()

	dens_list_per_step = []

	for step in step_list:

		dens_dict = {}
		# Data starts at 8th line
		for line in lines[7:]:
			line = line.split()
			dens_dict[line[0]] = float(line[step+1])

		# for key, value in sorted(dens_dict.iteritems(), key=lambda (k,v): (v,k)):
		#     print ("%s: %s" % (key, value))

		sorted_dens_dict = sorted(dens_dict.items(), key=lambda kv: kv[1], reverse=True)
		dens_list_per_step.append(sorted_dens_dict)

	dens_file.close()

	cwd = os.getcwd()
	sorted_dens_file = open('ranked dens', 'w')

	txt = ''
	for step in step_list:
		txt += '{:<20}'.format(step)
	txt += '\n\n'

	for i in range(len(dens_list_per_step[0])):
		for step in step_list:
			dens_list = dens_list_per_step[step]
			txt += '{:<8}{:<12.2E}'.format(dens_list[i][0], dens_list[i][1])
		txt += '\n'	

	sorted_dens_file.write(txt)
	sorted_dens_file.close()


def rank_nuclide_per_reac_rate(bucell, step_list, path, file_name):

	dens_path = path +'/output_summary' + '/{}_dens'.format(bucell)
	dens_file = open(dens_path, 'r')
	xs_path = path +'/output_summary' + '/{}_xs_lib'.format(bucell)
	xs_file = open(xs_path, 'r')

	# Read densities
	lines = dens_file.readlines()
	dens_dict_per_step = []
	for step in step_list:
		dens_dict = {}
		# Data starts at 8th line
		for line in lines[7:]:
			line = line.split()
			name = zamid_to_name(line[0])
			dens_dict[name] = float(line[step+1])

		dens_dict_per_step.append(dens_dict)

	# Read xs
	lines = xs_file.readlines()
	xs_dict_per_step = []
	for step in step_list:
		xs_dict = {}
		count = 0
		# Data starts at 8th line
		for line in lines[7:]:
			if line == '\n':
				continue
			line = line.split()

			if line[0].split('-')[0] in d.nuc_name_dic:
				# For first data, you need to read data first
				if count == 0:
					nucl_name = line[0]
					abs_xs = float(line[step+3])

				# Reached new nuclide, need to store data
				else:
					xs_dict[nucl_name] = abs_xs
					nucl_name = line[0]
					abs_xs = float(line[step+3])

			else:

				abs_xs += float(line[step+2])
				count += 1

		xs_dict_per_step.append(xs_dict)

	flux_per_step = []
	flux_seq = read_flux(dens_path)
	for step in step_list:
		flux_per_step.append(flux_seq[step+1]) # flux starts with 0


	# Now we will go over each nuclide in the dict per step and multiply xs with density
	#print (xs_dict_per_step)

	sorted_reac_dict_per_step = []
	total_abs_per_step = []
	for step in range(len(step_list)):
		dens_dict = dens_dict_per_step[step]
		xs_dict = xs_dict_per_step[step]
		reac_dict = {}
		for nucl in xs_dict:
			if nucl in dens_dict:
				nucl_dens = dens_dict[nucl]
				nucl_xs = xs_dict[nucl]
				reac_rate = nucl_dens*nucl_xs
				reac_dict[nucl] = reac_rate

		# Create a sorted list of tuples
		sorted_reac_tuple = sorted(reac_dict.items(), key=lambda kv: kv[1], reverse=True)
		sorted_reac_dict_per_step.append(sorted_reac_tuple)
		total_abs = 0
		for i in sorted_reac_tuple:
			total_abs += i[1]
		total_abs_per_step.append(total_abs)

	dens_file.close()
	xs_file.close()


	cwd = os.getcwd()
	sorted_reac_file = open('{} ranked react'.format(file_name), 'w')

	txt = ''
	for step in step_list:
		txt += '{:<20}'.format(step)
	txt += '\n\n'

	for step in range(len(step_list)):
		txt += 'flux={:<10.5E}'.format(flux_per_step[step])

	txt += '\n'

	for step in range(len(step_list)):
		txt += 'tot-abs={:<10.5E}'.format(total_abs_per_step[step])

	txt += '\n\n'

	for i in range(len(sorted_reac_dict_per_step[0])):
		for step in range(len(step_list)):
			reac_tuple_list = sorted_reac_dict_per_step[step]
			txt += '{:<8}{:<12.2E}'.format(reac_tuple_list[i][0], reac_tuple_list[i][1])
		txt += '\n'	

	sorted_reac_file.write(txt)
	sorted_reac_file.close()

def plot_matrix_from_compressed_matrix(path, step, cell):

	path_to_xs = path +'/step_{}'.format(step) +'/{}_cell'.format(cell) +'/matrix/xs_mat'
	path_to_decay = path +'/step_{}'.format(step) +'/{}_cell'.format(cell) +'/matrix/decay_mat'
	file_xs = open(path_to_xs, 'r')
	file_decay = open(path_to_decay, 'r')
	lines_xs = file_xs.readlines()
	lines_decay = file_decay.readlines()

	plt.figure(1)

	count = 0
	# x_vect = [i for i in range(len(lines))]
	for i in range(len(lines_xs)):
		line_xs = lines_xs[i]
		line_decay = lines_decay[i]
		zamid = line_xs.split('|')[0]
		line_elt_xs = line_xs.split(':')[1]
		elts_xs = line_elt_xs.split(',')[:-1] # Last element empty because of last coma in each line
		current_line_xs = []
		current_x_vect_xs = []
		elt_val_xs = len(lines_xs) - count  # The value assigned is the index of the nuclide starting from the last
		for elt_xs in elts_xs:
			elt_index = int(elt_xs.split()[0])
			current_x_vect_xs .append(elt_index)
			current_line_xs.append(elt_val_xs)
		if i == len(lines_xs)-1:
			plt.scatter(current_x_vect_xs , current_line_xs , marker='s', color = 'k', s = 4, label = 'cross section')
		else:
			plt.scatter(current_x_vect_xs , current_line_xs , marker='s', color = 'k', s = 4)
		line_elt_decay = line_decay.split(':')[1]
		elts_decay = line_elt_decay.split(',')[:-1] # Last element empty because of last coma in each line
		current_line_decay = []
		current_x_vect_decay = []
		elt_val_decay = len(lines_decay) - count  # The value assigned is the index of the nuclide starting from the last
		for elt_decay in elts_decay:
			elt_index = int(elt_decay.split()[0])
			current_x_vect_decay.append(elt_index)
			current_line_decay.append(elt_val_decay)
		if i == len(lines_xs)-1:	
			plt.scatter(current_x_vect_decay , current_line_decay , marker='+', color = 'r', s = 4, label = 'decay')
		else:
			plt.scatter(current_x_vect_decay , current_line_decay , marker='+', color = 'r', s = 4)
		#mat.append(current_line)
		count += 1

	file_xs.close()
	file_decay.close()
	plt.legend()
	plt.show()

def plot_matrix_bysign_from_compressed_matrix(path, step, cell):

	path_to_xs = path +'/step_{}'.format(step) +'/{}_cell'.format(cell) +'/matrix/xs_mat'
	path_to_decay = path +'/step_{}'.format(step) +'/{}_cell'.format(cell) +'/matrix/decay_mat'
	file_xs = open(path_to_xs, 'r')
	file_decay = open(path_to_decay, 'r')
	lines_xs = file_xs.readlines()
	lines_decay = file_decay.readlines()

	size = 0.8

	plt.figure(1)

	count = 0
	# x_vect = [i for i in range(len(lines))]
	for i in range(len(lines_xs)):
		line_xs = lines_xs[i]
		line_decay = lines_decay[i]
		zamid = line_xs.split('|')[0]
		line_elt_xs = line_xs.split(':')[1]
		elts_xs = line_elt_xs.split(',')[:-1] # Last element empty because of last coma in each line
		current_line_xs = []
		current_x_vect_xs = []
		#elt_val_xs = len(lines_xs) - count  # The value assigned is the index of the nuclide starting from the last
		elt_val_xs = count
		for elt_xs in elts_xs[1:]:
			elt_index = int(elt_xs.split()[0])
			current_x_vect_xs .append(elt_index)
			current_line_xs.append(elt_val_xs)
		if i == len(lines_xs)-1:
			plt.scatter(current_x_vect_xs , current_line_xs , marker='s', color = 'r', s = size, label = 'Production Terms')
		else:
			plt.scatter(current_x_vect_xs , current_line_xs , marker='s', color = 'r', s = size)
		line_elt_decay = line_decay.split(':')[1]
		elts_decay = line_elt_decay.split(',')[:-1] # Last element empty because of last coma in each line
		current_line_decay = []
		current_x_vect_decay = []
		#elt_val_decay = len(lines_decay) - count  # The value assigned is the index of the nuclide starting from the last
		elt_val_decay = count
		for elt_decay in elts_decay:
			elt_index = int(elt_decay.split()[0])
			current_x_vect_decay.append(elt_index)
			current_line_decay.append(elt_val_decay)
		if i == len(lines_xs)-1:	
			plt.scatter(current_x_vect_decay , current_line_decay , marker='s', color = 'r', s = size)
		else:
			plt.scatter(current_x_vect_decay , current_line_decay , marker='s', color = 'r', s = size)
		#mat.append(current_line)
		count += 1

	# cover diag elt
	for i in range(len(lines_xs)):
		if i == 0:
			plt.scatter([i], [i],  marker='s', color = 'k', s = size, label = 'Destruction Terms')
		else:
			plt.scatter([i], [i],  marker='s', color = 'k', s = size)
	file_xs.close()
	file_decay.close()
	plt.gca().invert_yaxis()
	plt.ylabel('Row index')
	plt.xlabel('Column index')
	plt.legend()
	plt.grid()
	plt.show()



class xs_name_not_found(Exception):
    """Raise when the user tries to access fission XS for a nuclide which fission XS have not been set yet """
    pass





# def bucells_average_density(bucell_list, step_list):

# 	bucell_vol_list = []

# 	system_parameters_file = open(os.getcwd() + '/output_summary/system_parameters')
	







 