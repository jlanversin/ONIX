import onix.utils as utils
from . import salameche
import os

class System(object):

	def __init__(self, system_id):

		self._id = system_id
		self._bucell_dict = None
		self._bounding_box = None
		self._sequence = None
		self._output_summary_path = None

	@property
	def id(self):
		return self._id
	
	@property
	def bucell_dict(self):

		return self._bucell_dict

	@bucell_dict.setter
	def bucell_dict(self, bucell_dict):

		self._bucell_dict = bucell_dict

	def add_bucell_dict(self, new_bucell_dict):

		# You need to make sure that there are no cells with same ids
		if self.bucell_dict == None:
			self.bucell_dict = new_bucell_dict
		else:
			old_bucell_dict = self.bucell_dict.copy()
			updated_bucell_dict = {**old_bucell_dict, **new_bucell_dict}
			self.bucell_dict = updated_bucell_dict

	def add_bucell(self, new_bucell):

		if self.bucell_dict == None:
			self.bucell_dict = {new_bucell.id : new_bucell}
		else:
			# You need to make sure this id is not already taken by other cell
			updated_bucell_dict = self.bucell_dict.copy()
			updated_bucell_dict[new_bucell.id] = new_bucell

			self.bucell_dict = updated_bucell_dict

	def get_bucell(self, name):

		result = None
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			if bucell.name == name:
				result = bucell
		
		if result ==  None:
			raise Cell_name_not_found('BUCell named {} has not been found in system no. {}'.format(name, self.id))

		return result

	def get_bucell_list(self):

		bucell_dict = self.bucell_dict
		bucell_list = utils.cell_dict_to_cell_list(bucell_dict)

		return bucell_list

	
	@property
	def bounding_box(self):

		return self._bounding_box

	@bounding_box.setter
	def bounding_box(self, bounding_box):

		self._bounding_box = bounding_box

	@property
	def total_vol(self):

		return self._total_vol

	@total_vol.setter
	def total_vol(self, total_vol):

		self._total_vol = total_vol

	def get_tot_hm(self):

		tot_hm = 0
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			tot_hm += bucell.get_hm()

		return tot_hm

	def get_tot_ihm(self):

		tot_ihm = 0
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			tot_ihm += bucell.ihm

		return tot_ihm
	
	def set_default_decay_for_all(self):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.set_default_decay_lib()

	def set_default_decay_for_all_no_add(self):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.set_default_decay_lib_no_add()

	def set_decay_for_all(self, decay_lib_path):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.set_decay_lib(decay_lib_path)

	def set_default_xs_for_all(self):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.set_default_xs_lib()

	def set_xs_for_all(self, xs_lib_path):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.set_xs_lib(xs_lib_path)

	def set_default_fy_for_all(self):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.set_default_fy_lib()

	def set_default_fy_for_all_no_add(self):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.set_default_fy_lib_no_add()

	def set_fy_for_all(self, fy_lib_path):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.set_fy_lib(fy_lib_path)

	@property
	def sequence(self):

		return self._sequence

	@sequence.setter
	def sequence(self, sequence):

		self._sequence = sequence

	# In addition to set sequence for each cells
	# This method also converts the average power_dens sequence to tot_power sequence if constant_power mode
	def set_sequence(self, sequence, mode = 'stand alone'):

		# Create tot_pow, av_pow_dens sequence
		# set time seq or bu seq depending on input

		sequence._initial_system_conversion(self)

		self.sequence = sequence
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			if mode == 'stand alone':
				bucell.set_sequence(sequence)
			elif mode == 'couple':
				bucell.set_sequence(sequence, mode = 'couple')

	@property
	def bu_sec_conv_factor(self):
		
		return self._bu_sec_conv_factor
	

	def _set_bu_sec_conv_factor(self):

		total_vol = self.total_vol
		# Since this get_tot_hm is called before burn, hm is ihm
		total_ihm = self.get_tot_hm()

		print('total_ihm', total_ihm)

		if total_ihm == 0:
			self._bu_sec_conv_factor = 0
		else:
			self._bu_sec_conv_factor = utils.get_bu_sec_conv_factor(total_vol, total_ihm)

	def zam_order_passlist(self):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell.passlist.zam_order_passport_list_2()

	# Burn the system according to the sequence set to system
	# def burn_loop(self):

	# 	sequence = self.sequence
	# 	steps_number = sequence.steps_number
	# 	for s in range(steps_number):

	def _print_current_allreacs_rank(self):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell._print_current_allreacs_rank()

	def _print_summary_allreacs_rank(self):

		summary_path = self._output_summary_path
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell._print_summary_allreacs_rank(summary_path)

	def _print_summary_dens(self):

		summary_path = self._output_summary_path
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell._print_summary_dens(summary_path)

	def _print_summary_subdens(self):

		summary_path = self._output_summary_path
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell._print_summary_subdens(summary_path)

	def _print_summary_xs(self):

		summary_path = self._output_summary_path
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell._print_summary_xs(summary_path)

	def _print_summary_isomeric_branching_ratio(self):

		summary_path = self._output_summary_path
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell._print_summary_isomeric_branching_ratio(summary_path)

	def _print_summary_flux_spectrum(self, mg_energy_bin):

		summary_path = self._output_summary_path
		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			bucell._print_summary_flux_spectrum(summary_path, mg_energy_bin)

	def _print_summary_kinf(self):

		summary_path = self._output_summary_path
		file_name = summary_path + '/kinf'
		write_file = open(file_name, 'w')
		sequence = self.sequence
		time_seq = sequence.time_seq
		system_bu_seq = sequence.system_bu_seq
		kinf_seq = sequence.kinf_seq
		steps_number = sequence.macrosteps_number

		txt = ''
		txt += '{:<12}'.format('TIME')
		txt += '{:<13.5E}'.format(time_seq[0]/(24*3600))# in days
		for s in range(steps_number):
			txt += '{:<13.5E}'.format(time_seq[s+1]/(24*3600))# in days

		txt += '\n'

		txt += '{:<12}'.format('SYSTEM-BU')
		txt += '{:<13.5E}'.format(system_bu_seq[0])
		for s in range(steps_number):
			txt += '{:<13.5E}'.format(system_bu_seq[s+1])

		txt += '\n\n'

		txt += '{:<12}'.format('K-INF')
		txt += '{:<13.5E}'.format(kinf_seq[0][0])
		for s in range(steps_number):
			txt += '{:<13.5E}'.format(kinf_seq[s+1][0])

		txt += '\n'

		txt += '{:<12}'.format('UNCERTAINTY')
		txt += '{:<13.5E}'.format(kinf_seq[0][1])
		for s in range(steps_number):
			txt += '{:<13.5E}'.format(kinf_seq[s+1][1])

		write_file.write(txt)
		write_file.close()


	def _print_summary_param(self):

		summary_path = self._output_summary_path
		file_name = summary_path + '/system_parameters'
		write_file = open(file_name, 'w')
		txt = 'System Volume [cm³] = {}\n'.format(self.total_vol)
		txt += 'System IHM [g] = {}\n\n'.format(self.get_tot_ihm())
		for bucell_id in self.bucell_dict:
			bucell = self.bucell_dict[bucell_id]
			txt += 'BuCell {}\n'.format(bucell.name)
			txt += 'Volume [cm³] = {}\n'.format(bucell.vol)
			txt += 'IHM [g] = {}\n\n'.format(bucell.ihm)
		write_file.write(txt)
		write_file.close()

	def copy_cell_folders_to_step_folder(self, s):

		bucell_dict = self.bucell_dict
		for bucell_id in bucell_dict:
			bucell = bucell_dict[bucell_id]
			print(bucell.name, bucell.id)
			bucell.copy_cell_folders_to_step_folder(s)

	def _gen_output_summary_folder(self):

		name = 'output_summary'
		self._output_summary_path = utils.get_folder_path(name)
		utils.gen_folder(name)

	def burn(self):

		sequence = self.sequence
		steps_number = sequence.steps_number
		for s in range(steps_number):

			print ('\n\n\n\n STEP {}\n\n\n\n'.format(s))

			sequence._gen_step_folder(s)
			salameche.burn_step(self, s)

		system._gen_output_summary_folder()
		system._print_summary_allreacs_rank()
		system._print_summary_dens()




	def print_bucell_nuclides(self, bucell, step, nuclide_list):

		bucell = self.get_bucell(bucell)
		file_path = os.getcwd() + '/{}_step{}_nucl_print'.format(bucell.name, step)
		file = open(file_path, 'w')

		passport_list = bucell.passlist.passport_list

		nucl_name_density = []
		txt = ''
		for nucl in passport_list:
			nucl_name = utils.onix_name_to_openmc_name(nucl.name)
			if nucl_name in nuclide_list:
				txt += '{:<10}'.format(nucl_name)
				txt += '{:<13.5E}'.format(nucl.dens_seq[step])
				txt += '\n'

		file.write(txt)
		file.close()


class Cell_name_not_found(Exception):

	pass
