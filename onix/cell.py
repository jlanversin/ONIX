from .sequence import Sequence as seq
from .passlist import Passlist as pl
from .passport import Passport as pp
from . import data
import onix.utils as utils
import copy
import os
import shutil
import operator
import time



class Cell(object):

	_NA = 6.02214086e+23
	zero_dens_1_atm = 1E-24

	def __init__(self, cell_id, name):

		self._id = cell_id
		self._name = name

		self._passlist = None
		self._initial_nucl = None
		self._decay_b_lib = None
		self._decay_a_lib = None
		self._xs_lib = None
		self._MC_XS_nucl_list = None
		self._fy_lib = None
		self._nucl_set = None
		# self._MC_flux = None
		# self._MC_flux_seq = None
		# self._flux = None
		self._FMF = None
		self._output_summary_path = None

	# Mainly designed for the code. Used when the user use text input rather than python module
	# probably obsolete
	# def set_from_input(self, input):

	# 	cell_id = self._cell_id
	# 	self._vol = input.vol(cell_id)
	# 	self._hm_vol = input.hm_vol(cell_id)
	# 	self._initial_nuc = input.initial_nuc(cell_id)

	# 	self._nuc_list = input.nucl_list
	# 	self._passlist = pl.gen_passlist_from_input(input, cell_id)
	# 	self._passdic = pl.passlist_to_passdic(self._passlist)
	# 	self._index_dic = pl.index_dic(self._passlist)

	# 	# Generate the list of fission child for each fissile nuclide
	# 	pl.gen_fission_child(self._passdic, input)

	# 	self.set_ihm(self._passlist, self._hm_vol)
	# 	self.set_bu_sec_conv_factor(self._vol, self._ihm)

	# 	self.link_sequence_from_input(input, cell_id, self._bu_sec_conv_factor)

	# 	self._set_tree()
	# 	self._set_leaves()
	# 	fy_parent = input.fy_parent
	#	self._set_fission_tree()
	#	self._set_fission_leaves()
	# 	self._fission_tree = self.build_fission_tree(self._tree, fy_parent)
	# 	self._fission_leaves = self.build_fission_leaves()

	@property
	def id(self):

		return self._id

	@property
	def name(self):

		return self._name



	# Set the initial densities to nuclides and create the initial nuc list
	def set_initial_dens(self, dens_dict):

		nucl_list = list(dens_dict.keys())

		if utils.is_name(nucl_list[0]):
			nucl_list = utils.name_list_to_zamid_list(nucl_list)

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)


		passlist = self.passlist
		passlist._set_initial_dens(dens_dict)

		self._set_initial_nucl(nucl_list)

	def _check_nucl_list_consistency(self):

		# Return lib_nucl and will also stop code if there are no decay/xs libs set
		lib_nucl = self.get_lib_nucl()
		# If no nucl_set set, nucl_set is set to default to lib_nucl
		# If nucl_set set, check if it is included in lib_nucl
		if self.nucl_set == None:
			self.nucl_set = lib_nucl
		elif not utils.is_lista_in_listb(self.nucl_set, lib_nucl):
			raise Nucl_set_not_in_Lib_nucl('Some of bucell {} nuclide-set"s nuclides are not in data libraries'.format(self.name))

		initial_nucl = self.initial_nucl

		nucl_set = self.nucl_set

		if initial_nucl == None:
			raise Initial_nucl_not_set('Cell {} has not been attributed initial nuclides'.format(self.name))
		elif not utils.is_lista_in_listb(initial_nucl, nucl_set):
			raise Initial_nucl_not_in_Nucl_set('Some of bucell {} initial nuclides {} are not in nuclide set'.format(self.name, initial_nucl))

	def _set_step_dens(self):

		passlist = self.passlist
		passport_list = passlist.passport_list

		for i in range(len(passport_list)):
			nuc_pass = passport_list[i]
			nuc_pass._set_step_dens()

	# This is effectively set_substep_dens
	def _update_dens(self, N, ss, ssn):

		# Technically I don't need to use N_dic since passport_list
		# is actually automatically ordered in the same way as N via mb
		# (the argument passport_list is a pointer to the object, not a copy of the object)
		# You need to be sure that N and passport_list are ordered in the same
		passlist = self.passlist
		passport_list = passlist.passport_list

		for i in range(len(passport_list)):
			nuc_pass = passport_list[i]
			nuc_pass._set_substep_dens(N[i], ss)

	@property
	def initial_nucl(self):

		return self._initial_nucl

	# def _get_initial_nucl(self, dens_dic):

	# 	initial_nuc = []
	# 	for zamid in dens_dic:
	# 		initial_nuc.append(zamid)

	# 	return initial_nuc

	def _set_initial_nucl(self, nucl_id_list: str):

		# Here a copy is made because the list nucl_id_list is latter
		# modified by other functions (set_default_decay and others)
		# This will accidentally modify initial_nucl if it is simply equalled to nucl_id_list
		self._initial_nucl = nucl_id_list.copy()

	@property
	def MC_XS_nucl_list(self):
		return self._MC_XS_nucl_list

	@MC_XS_nucl_list.setter
	def MC_XS_nucl_list(self, MC_XS_nucl_list):

		self._MC_XS_nucl_list = MC_XS_nucl_list
	

	@property
	def nucl_set(self):

		return self._nucl_set

	@nucl_set.setter
	def nucl_set(self, nucl_set):

		self._nucl_set = nucl_set

	# User choose among the prepared nucl set
	# Development for latter
	# def choose_nucl_set(self, nucl_set_family):

	# This property is just to know which nuclides
	# are initially present (dens non-0) in the cell
	@property
	def init_nucl(self):

		return self._init_nucl

	@init_nucl.setter
	def init_nucl(self, init_nucl):

		self._init_nucl = init_nucl
	

	def get_total_dens(self):

		total_dens = 0
		passlist = self.passlist
		passport_list = passlist.passport_list
		for nuc_pass in passport_list:
			total_dens += nuc_pass.current_dens

		return total_dens

	def get_subtotal_dens(self, nucl_list):

		subtotal_dens = 0
		passlist = self.passlist
		passport_list = passlist.passport_list
		for nuc_pass in passport_list:
			if nuc_pass.name in nucl_list:

				subtotal_dens += nuc_pass.current_dens

		return subtotal_dens

	# This version of the method counts zero dens atom as 1E-24. It is used
	# when computing the density of material for OpenMC
	def get_subtotal_dens_counting_zero_dens(self, nucl_list):

		subtotal_dens = 0
		passlist = self.passlist
		passport_list = passlist.passport_list
		for nuc_pass in passport_list:
			if nuc_pass.name in nucl_list:

				# following is a possible solution to the xs drop problem
				# nuclide that are zero in onix and supposed to be tallied
				# are actually set to 1E-24 in openmc
				# therefore, we count them as 1E-24 not zero

				if nuc_pass.current_dens == 0:
					dens = self.zero_dens_1_atm
				else:
					dens = nuc_pass.current_dens

				subtotal_dens += dens

				#subtotal_dens += nuc_pass.current_dens


		return subtotal_dens

	def get_nucl_dens_for_openmc(self, nucl_id):

		nucl = self.get_nuclide(nucl_id)
		nucl_dens = nucl.current_dens

		if nucl_dens == 0:
			return self.zero_dens_1_atm
		else:
			return nucl_dens

	def get_nucl_ao(self, nucl_id):

		nucl = self.get_nuclide(nucl_id)
		nucl_dens = nucl.current_dens
		total_dens = self.get_total_dens()

		return nucl_dens*100/total_dens

	def get_nucl_subao(self, nucl_id, nucl_list):

		nucl = self.get_nuclide(nucl_id)
		nucl_dens = nucl.current_dens

		# OpenMC behaves in a bizzare way when certain nuclides are set to 0.0
		# in material
		# Here we replaced any 0.0 dens nuclide with 1E-24
		if nucl_dens == 0.0:

			return self.zero_dens_1_atm

		else:

			subtotal_dens = self.get_subtotal_dens_counting_zero_dens(nucl_list)

			return nucl_dens/subtotal_dens

	@property
	def vol(self):

		return self._vol

	@vol.setter
	def vol(self, vol):

		self._vol = vol

	@property
	def hm_vol(self):

		return self._hm_vol

	@hm_vol.setter
	def hm_vol(self, hm_vol):

		self._hm_vol = hm_vol


	def set_passlist(self, nucl_list):

		if utils.is_list_redundant(nucl_list) == True:
			redundant_elt = utils.get_list_redundant_elt(nucl_list)
			raise Nuclide_list_redundant('Cell {} passlist object has been given a redundant nuclide list with following redundant element {}'.format(self.id, redundant_elt))

		self.passlist = self.get_passlist(nucl_list)

	def get_passlist(self, nucl_list):

		passlist = pl(nucl_list)

		return passlist

	def get_nucl_list(self):

		passlist = self.passlist

		return passlist.nucl_list

	def get_nuclide(self, nuclide_id):

		if utils.is_name(nuclide_id):
			nuclide_id = utils.name_to_zamid(nuclide_id)

		passport_dict = self.passlist._get_zamid_passport_dict()

		return passport_dict[nuclide_id]

	# def update_passlist(self, new_nucl_list):

	# 	current_passlist = self.passlist
	# 	current_passlist.add_nucl_list(new_extra_nuclides)
	# 	merged_nucl_list = current_passlist_nucl_list + list(set(new_nucl_list) - set(current_passlist_nucl_list))

	# 	self.set_passlist(merged_nucl_list)

	@property
	def passlist(self):

		return self._passlist

	@passlist.setter
	def passlist(self, passlist):

		self._passlist = passlist

	@property
	def index_dict(self):

		return self._index_dict

	@index_dict.setter
	def index_dict(self, index_dict):

		self._index_dict = index_dict

	@property
	def passdic(self):

		return self._passdic

	@passdic.setter
	def passdic(self, passdic):

		self._passdic = passdic



	# def _set_sequence_from_input(self, sequence_dict):


	# 	cell_id = self._cell_id
	# 	sequence = seq(cell_id)

	# 	self._set_ihm()
	# 	self._set_bu_sec_conv_factor()

	# 	passlist = self.passlist
	# 	bu_sec_conv_factor = self._bu_sec_conv_factor

	# 	sequence._set_from_input(sequence_dict, passlist, bu_sec_conv_factor)

	# 	self._sequence = sequence

	def set_sequence(self, sequence, mode = 'stand_alone'):

		# A copy of the sequence is done because when setting bu_sec_conv_factor, we don't want the original to 
		# be changed too because it might be used for other	 cells
		sequence_copy = copy.deepcopy(sequence)

		self._set_ihm()
		self._set_bu_sec_conv_factor()

		passlist = self.passlist
		bu_sec_conv_factor = self._bu_sec_conv_factor

		#The _cell_conversion method was actually just calling _set_initial_bucell_bu()
		#sequence_copy._cell_conversion(passlist, bu_sec_conv_factor, mode)
		sequence_copy._set_initial_bucell_bu()
		self._sequence = sequence_copy

	@property
	def sequence(self):

		return self._sequence




	# WARNING: ihm is the initial mass
	# As long as it is set before burn that's good
	# But if it is set during or qfter burn, the value
	# retrieved by @property will not be the initial mass anymore

	@property
	def ihm(self):
		"""Returns the absolute values of the decay constant of the nuclide"""
		ihm = self._ihm
		if ihm is None:
			pass  # define exception for undefined variable

		return ihm

	def set_ihm(self, passlist, hm_vol):

		ihm = utils.get_hm(passlist, hm_vol)

		self._ihm = self.ihm

	def get_hm(self):

		passlist = self.passlist
		vol = self.vol
		hm = utils.get_hm(passlist, vol)

		return hm


	def _set_ihm(self):

		passlist = self._passlist
	#	hm_vol = self._hm # when user had to define the hm_vol
		vol = self.vol
	#	ihm = utils.get_ihm(passlist, hm_vol) # when user had to define the hm_vol
		ihm = utils.get_hm(passlist, vol)

		self._ihm = ihm

	def check_act_presence(self):

		hm = utils.get_hm(self.passlist, self.vol)
		if hm == 0.0:
			return 'no'
		elif hm > 0.0:
			return 'yes'

	@property
	def bu_sec_conv_factor(self):
		"""Returns the absolute values of the decay constant of the nuclide"""
		bu_sec_conv_factor = self._bu_sec_conv_factor
		if bu_sec_conv_factor is None:
			pass  # define exception for undefined variable
		return bu_sec_conv_factor

	def set_bu_sec_conv_factor(self, vol, ihm):

		if ihm == 0:
			self._bu_sec_conv_factor = 0
		else:
			self._bu_sec_conv_factor = utils.get_bu_sec_conv_factor(vol, ihm)

	def _set_bu_sec_conv_factor(self):

		vol = self._vol
		ihm = self._ihm

		if ihm == 0:
			self._bu_sec_conv_factor = 0
		else:
			self._bu_sec_conv_factor = utils.get_bu_sec_conv_factor(vol, ihm)

	# Update power density from flux in each cell at each substep
	def _update_pow_dens(self, flux):

		passlist = self.passlist
		fission_energy_rate = 0
		conv_Mev_J = 1.60218e-13
		for i in passlist.passport_list:
			if i.current_xs != None:
				if i.fission_E != None and 'fission' in i.current_xs:
					fission_xs = i.current_xs['fission'][0] # cm2 10^-24
					fission_E = i.fission_E # MeV
					dens = i.current_dens # cm-1 barn-1 (cm-3 10+24)

					# print (i.name)
					# print (fission_E)
					# print (fission_xs)
					# print (dens)

					fission_energy_rate += fission_xs*fission_E*dens


		# print ('fission_energy_rate',fission_energy_rate)

		new_pow_dens = flux*(fission_energy_rate*conv_Mev_J) # W/cm3 = kW/l

		#print ('new_pow_dens',new_pow_dens)

		return new_pow_dens

	# Update flux from power density in each cell at each substep
	def _update_flux(self, pow_dens):

		passlist = self.passlist
		passport_list = passlist.passport_list
		print('update_flux called')
		fission_energy_rate = 0
		conv_Mev_J = 1.60218e-13
		for i in passport_list:
			if i.fission_E != None and i.current_xs != None:
				if 'fission' in i.current_xs:
					fission_xs = i.current_xs['fission'][0]
					fission_E = i.fission_E
					dens = i.current_dens

					fission_energy_rate += fission_xs*fission_E*dens

		new_flux = pow_dens/(fission_energy_rate*conv_Mev_J)

		return new_flux

	def _change_total_density(self, s):

		sequence = self.sequence
		density_change_dict = sequence.density_change_dict
		if density_change_dict != None:

			if self.name in density_change_dict:
				bucell_density_change_dict = density_change_dict[self.name]
				# When user sets a new density for step s+1, the Monte Carlo simulation of step s+1
				# needs to run with this new density
				# For conveniency, ONIX thus changes the density at the end of step s
				# This means that if we are now at step s, we need to look if the user has specified a new 
				# density for step s+1
				if s+1 in bucell_density_change_dict:

					new_tot_dens = bucell_density_change_dict[s+1]
					current_total_dens = self.get_total_dens()
					factor = new_tot_dens/current_total_dens
					passport_list = self.passlist.passport_list
					for nuc_pass in passport_list:
						current_dens = nuc_pass.current_dens
						nuc_pass.current_dens = current_dens*factor

	def _change_isotope_density(self,s):

		sequence = self.sequence
		isotopic_change_dict = sequence.isotopic_change_dict
		if isotopic_change_dict != None:
			if self.name in isotopic_change_dict:
				bucell_isotopic_change_dict = isotopic_change_dict[self.name]
				unit = bucell_isotopic_change_dict['unit']
				current_total_dens = self.get_total_dens()
				for nucl_name in bucell_isotopic_change_dict:
					if nucl_name == 'unit':
						continue
					nuclide_isotopic_change_dict = bucell_isotopic_change_dict[nucl_name]

					# When user sets a new density for step s+1, the Monte Carlo simulation of step s+1
					# needs to run with this new density
					# For conveniency, ONIX thus changes the density at the end of step s
					# This means that if we are now at step s, we need to look if the user has specified a new 
					# density for step s+1
					if s+1 in nuclide_isotopic_change_dict:
						nucl_passport = self.get_nuclide(nucl_name)
						# Only the current_dens is set to the user-defined density
						# the dens_seq and dens_subseq_mat last value are left unchanged
						# However both set_dens_cells (which updates material for new openmc simulation)
						# and mat_builder (which prepares the matrix for new depletion calculation) uses
						# current_dens. Therefore, only changing current_dens will update calculation without
						# changing stored value that will be printed

						# If unit is number density
						if unit == 'number density':
							new_dens = nuclide_isotopic_change_dict[s+1]
						# else if unit is atom fraction
						if unit == 'atom fraction':
							new_dens = nuclide_isotopic_change_dict[s+1]*current_total_dens

						nucl_passport.current_dens = new_dens


	# Get the a sub passport list of all actinides
	def get_act_passport_list(self):

		passport_list = self.passlist.passport_list
		act_passport_list = []
		for nucl in passport_list:
			if nucl.get_FAM() == 'ACT':
				act_passport_list.append(nucl)

		return act_passport_list

	# Get the a sub passport list of all fission products
	def get_fp_passport_list(self):

		passport_list = self.passlist.passport_list
		fp_passport_list = []
		for nucl in passport_list:
			if nucl.get_FAM() == 'FP':
				fp_passport_list.append(nucl)

		return fp_passport_list

	# Get the a sub passport list of all activation products
	def get_avt_passport_list(self):

		passport_list = self.passlist.passport_list
		avt_passport_list = []
		for nucl in passport_list:
			if nucl.get_FAM() == 'AVT':
				avt_passport_list.append(nucl)

		return avt_passport_list

	def _set_libs_from_input(self, lib):

		self._decay_a_lib = lib['decay_a']
		self._decay_b_lib = lib['decay_b']
		self._xs_lib = lib['xs']
		self._fy_lib = lib['fy']


	# Set decay dictionary from a decay text library which path is indicated by the user
	def set_decay_lib(self, decay_lib_path):

		decay_b = data.read_lib_functions.read_decay_lib(decay_lib_path)
		decay_a = data.read_lib_functions.conv_decay_b_a(decay_b)

		nucl_list = list(decay_a.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)
		#	self.update_passlist(nucl_list)

		self._decay_b_lib = decay_b
		self._decay_a_lib = decay_a

		self.passlist._set_decay(decay_b, decay_a)

	# Set decay dictionary from the default decay text library
	# Will add to passlist those nuclides present in the lib that are not present in passlist
	def set_default_decay_lib(self):

		# This command will find the absolute path of cell.py
		# Since cell.py is located in onix, the default library is just in __file__path + ./data/default_libs/decay_lib
		#__file__path = os.path.abspath(os.path.dirname(__file__))

		#default_decay_lib_path = __file__path+ '/data/default_libs/decay_lib'
		default_decay_lib_path = data.default_decay_b_lib_path
		#default_decay_lib_path = '/home/julien/Open-Burnup.dev/onix/data/default_libs/decay_lib'

		decay_b = data.read_lib_functions.read_decay_lib(default_decay_lib_path)
		decay_a = data.read_lib_functions.conv_decay_b_a(decay_b)

		nucl_list = list(decay_a.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)

		self._decay_b_lib = decay_b
		self._decay_a_lib = decay_a

		# If the passlist has already been defined, pass the decay to each nuclide in passlist
		self.passlist._set_decay(decay_b, decay_a)

	# Set decay dictionary from the default decay text library
	# Will NOT add to passlist those nuclides present in the lib that are not present in passlist
	def set_default_decay_lib_no_add(self):

		# This command will find the absolute path of cell.py
		# Since cell.py is located in onix, the default library is just in __file__path + ./data/default_libs/decay_lib
		__file__path = os.path.abspath(os.path.dirname(__file__))

		default_decay_lib_path = __file__path+ '/data/default_libs/decay_lib'

		#default_decay_lib_path = '/home/julien/Open-Burnup.dev/onix/data/default_libs/decay_lib'

		decay_b = data.read_lib_functions.read_decay_lib(default_decay_lib_path)
		decay_a = data.read_lib_functions.conv_decay_b_a(decay_b)

		passlist = self.passlist
		if passlist == None:
			raise Passlist_not_defined('No passlist for this cell {} has been defined'.format(cell.id))

		self._decay_b_lib = decay_b
		self._decay_a_lib = decay_a

		self.passlist._set_decay(decay_b, decay_a)

	# Set the decay dictionary from a decay object defined by the user
	def set_decay(self, decay_object):

		decay_b = decay_object.decay_b
		decay_a = decay_object.decay_a

		nucl_list = list(decay_a.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)	

		self._decay_b_lib = decay_b
		self._decay_a_lib = decay_a

		# If the passlist has already been defined, pass the decay to each nuclide in passlist
		self.passlist._set_decay(decay_b, decay_a)

	@property
	def decay_b_lib(self):

		return self._decay_b_lib

	@property
	def decay_a_lib(self):

		return self._decay_a_lib



	def set_xs_lib(self, xs_lib_path):

		xs_dic = data.read_lib_functions.read_xs_lib(xs_lib_path)

		nucl_list = list(xs_dic.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)

		self._xs_lib = xs_dic

		self.passlist._set_xs(xs_dic)

	def set_default_xs_lib(self):

		# This command will find the absolute path of cell.py
		# Since cell.py is located in onix, the default library is just in __file__path + ./data/default_libs/decay_lib
		__file__path = os.path.abspath(os.path.dirname(__file__))

		default_xs_lib_path = __file__path+ '/data/default_libs/xs_lib'

		# default_xs_lib_path = '/home/julien/Open-Burnup.dev/onix/data/default_libs/xs_lib'
		# pupu_xs_lib_path = '/home/julien/Open-Burnup.dev/onix/data/default_libs/xs_lib_pupu'

		xs_dic = data.read_lib_functions.read_xs_lib(default_xs_lib_path)

		nucl_list = list(xs_dic.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)

		self._xs_lib = xs_dic

		self.passlist._set_xs(xs_dic)

	def set_default_xs_lib_no_add(self):

		# This command will find the absolute path of cell.py
		# Since cell.py is located in onix, the default library is just in __file__path + ./data/default_libs/decay_lib
		__file__path = os.path.abspath(os.path.dirname(__file__))

		default_xs_lib_path = __file__path+ '/data/default_libs/xs_lib'

		# default_xs_lib_path = '/home/julien/Open-Burnup.dev/onix/data/default_libs/xs_lib'
		# pupu_xs_lib_path = '/home/julien/Open-Burnup.dev/onix/data/default_libs/xs_lib_pupu'

		xs_dic = data.read_lib_functions.read_xs_lib(default_xs_lib_path)

		xs_dic = data.read_lib_functions.read_xs_lib(default_xs_lib_path)

		nucl_list = list(xs_dic.keys())

		passlist = self.passlist
		if passlist == None:
			raise Passlist_not_defined('No passlist for this cell {} has been defined'.format(cell.id))

		self._xs_lib = xs_dic

		self.passlist._set_xs(xs_dic)

	# Set the xs dictionary from a xs object defined by the user
	def set_xs(self, xs_object):

		xs_dic = xs_object.xs

		nucl_list = list(xs_dic.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)

		self._xs_lib = xs_dic

		self.passlist._set_xs(xs_dic)

	def overwrite_xs(self, xs_object):

		xs_dic = xs_object.xs

		nucl_list = list(xs_dic.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)

		self._xs_lib = xs_dic

		self.passlist._overwrite_xs(xs_dic)

	@property
	def xs_lib(self):

		return self._xs_lib


	def set_fy_lib(self, fy_lib_path):

		fy_dic = data.read_lib_functions.read_fy_lib(fy_lib_path)

		nucl_list = list(fy_dic.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)

		self._fy_lib = fy_dic

		self.passlist._set_fy(fy_dic)

	def set_default_fy_lib(self):

		# This command will find the absolute path of cell.py
		# Since cell.py is located in onix, the default library is just in __file__path + ./data/default_libs/decay_lib
		__file__path = os.path.abspath(os.path.dirname(__file__))

		#default_fy_lib_path = __file__path+ '/data/default_libs/fy_lib'

		default_fy_lib_path = data.default_fy_lib_path

		#default_fy_lib_path = '/home/julien/Open-Burnup.dev/onix/data/default_libs/fy_lib'

		fy_dic = data.read_lib_functions.read_fy_lib(default_fy_lib_path)

		nucl_list = list(fy_dic.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)

		self._fy_lib = fy_dic

		self.passlist._set_fy(fy_dic)

	# Warning, for fy, only the FP in the dic are automatically
	# Added to passlist nucl list. The parents nuclide defined are not
	# User must make be sure parents nuclide are already in nucl list

	def set_default_fy_lib_no_add(self):

		# This command will find the absolute path of cell.py
		# Since cell.py is located in onix, the default library is just in __file__path + ./data/default_libs/decay_lib
		__file__path = os.path.abspath(os.path.dirname(__file__))

		default_fy_lib_path = __file__path+ '/data/default_libs/fy_lib'

		#default_fy_lib_path = '/home/julien/Open-Burnup.dev/onix/data/default_libs/fy_lib'

		fy_dic = data.read_lib_functions.read_fy_lib(default_fy_lib_path)

		passlist = self.passlist
		if passlist == None:
			raise Passlist_not_defined('No passlist for this cell {} has been defined'.format(cell.id))

		self._fy_lib = fy_dic

		self.passlist._set_fy(fy_dic)

	def set_fy(self, fy_object):

		fy_dic = fy_object.fy

		nucl_list = list(fy_dic.keys())

		passlist = self.passlist
		if passlist == None:
			self.set_passlist(nucl_list)
		else:
			passlist._add_nucl_list(nucl_list)

		self._fy_lib = fy_dic

		self.passlist._set_fy(fy_dic)

	@property
	def fy_lib(self):

		return self._fy_lib




	def _print_dens(self, time_point, bu_point, s):

		self._print_dens_1(s)
		self._print_dens_2(time_point, bu_point, s)

	def _print_initial_dens(self):


		sequence = self.sequence
		initial_time = sequence.time_seq_vect[0]
		initial_bu = sequence.bu_seq_vect[0]
		i = 0 # dummy i
		s = -1

		self._print_dens_1(s)
		self._print_dens_2(s, i)


	def _print_xs_lib(self):

		cell_id = self.id
		passlist = self.passlist
		xs_lib = self.xs_lib
		passport_list = passlist.passport_list
		cell_folder_path = self.folder_path

		file_name = cell_folder_path + '/xs_lib'

		write_file = open(file_name, 'w')
		txt = ''
		txt += '\n\n--- Actinides ---\n\n'
		txt += utils.printer.xs_lib_header

		#loop over actinides
		for nucl in self.get_act_passport_list():
			nucl_name = nucl.name
			nucl_zamid = nucl.zamid
			nucl_xs = nucl.current_xs
			# If this nuclide has not been assigned cross section, continue
			if nucl_xs == None:
				continue
			count = 0
			for xs_name in nucl_xs:
				# removal should not be printed on the lib
				if xs_name ==  'removal':
					continue
				xs_val = nucl_xs[xs_name][0]
				xs_unc = nucl_xs[xs_name][1]
				if count == 0:
					txt += '{:^19}'.format(nucl_name)
				else:
					txt += '{:^19}'.format('')
				txt += '{:^12}'.format(nucl_zamid)
				txt += '{:^17}'.format(xs_name)	
				txt += '{:^17.16E}'.format(xs_val)
				txt += '{:^17}'.format(xs_unc)
				txt += '\n'
				count += 1

		txt += '\n\n--- Fission Products ---\n\n'
		txt += utils.printer.xs_lib_header

		#loop over fission products
		for nucl in self.get_fp_passport_list():
			nucl_name = nucl.name
			nucl_zamid = nucl.zamid
			nucl_xs = nucl.current_xs
			# If this nuclide has not been assigned cross section, continue
			if nucl_xs == None:
				continue
			count = 0
			for xs_name in nucl_xs:
				# removal should not be printed on the lib
				if xs_name ==  'removal':
					continue
				xs_val = nucl_xs[xs_name][0]
				xs_unc = nucl_xs[xs_name][1]
				if count == 0:
					txt += '{:^19}'.format(nucl_name)
				else:
					txt += '{:^19}'.format('')
				txt += '{:^12}'.format(nucl_zamid)
				txt += '{:^17}'.format(xs_name)
				txt += '{:^17.16E}'.format(xs_val)
				txt += '{:^17}'.format(xs_unc)
				txt += '\n'
				count += 1

		txt += '\n\n--- Activation Products ---\n\n'
		txt += utils.printer.xs_lib_header

		#loop over activation products
		for nucl in self.get_avt_passport_list():
			nucl_name = nucl.name
			nucl_zamid = nucl.zamid
			nucl_xs = nucl.current_xs
			# If this nuclide has not been assigned cross section, continue
			if nucl_xs == None:
				continue
			count = 0
			for xs_name in nucl_xs:
				# removal should not be printed on the lib
				if xs_name ==  'removal':
					continue
				xs_val = nucl_xs[xs_name][0]
				xs_unc = nucl_xs[xs_name][1]
				if count == 0:
					txt += '{:^19}'.format(nucl_name)
				else:
					txt += '{:^19}'.format('')
				txt += '{:^12}'.format(nucl_zamid)
				txt += '{:^17}'.format(xs_name)
				txt += '{:^17.16E}'.format(xs_val)
				txt += '{:^17}'.format(xs_unc)
				txt += '\n'
				count += 1

		write_file.write(txt)
		write_file.close()

	def _print_general_dens_1(self, s):

		sequence = self.sequence
		time_point = sequence.time_point(s)
		bu_point = sequence.bu_point(s)
		cell_id = self.id
		passlist = self.passlist
		passport_list = passlist.passport_list
		cell_folder_path = self.folder_path
		total_dens = self.get_total_dens()

		file_name = cell_folder_path + '/cell_{}_dens'.format(cell_id)

	   # if os.stat(file_name).st_size == 0:
		if s == -1:
			write_file = open(file_name, 'w')
			txt = ''
			txt += '       {}\n'.format(time_point)
			txt += '       {}\n'.format(bu_point)
			for i in passport_list:
				zamid = i.zamid
				dens = i.current_dens
				txt += '{}  {}\n'.format(zamid, dens)
			txt += 'Total  {}'.format(total_dens)
			write_file.write(txt)
			write_file.close()

	# #    elif os.stat(file_name).st_size != 0:
		elif s >= -1:
			print('file already exists')
			read_file = open(file_name, 'r')
			lines = read_file.readlines()
			append_file = open(file_name, 'w')
			for i in range(len(lines)):
				line = lines[i]
				if i == 0:
					append_txt = ' {}\n'.format(time_point)
					append_file.write(line.strip('\n') + append_txt)
				elif i == 1:
					append_txt = ' {}\n'.format(bu_point)
					append_file.write(line.strip('\n') + append_txt)
				elif i == len(lines) - 1:
					append_txt = ' {}'.format(total_dens)
					append_file.write(line + append_txt)
				else:
					j = i-2
					nuc_pass = passport_list[j]
					dens = nuc_pass.current_dens
					append_txt = ' {}\n'.format(nuc_pass.current_dens)
					append_file.write(line.strip('\n') + append_txt)

			read_file.close()
			append_file.close()

	def _print_substep_dens(self, s):

		cell_id = self.id
		passlist = self.passlist
		passport_list = passlist.passport_list
		cell_folder_path = self.folder_path
		sequence = self.sequence

		time_subseq = sequence.time_subseq_mat[s] # sequence time subsequence starts with initial point
		
		# you can't calculate system_bu_subseq_mat here if the norma is flux and sequence in time
		# because time to bu requires total power but total power is only known once you have the subseq
		# evolution of all cells
		system_bu_subseq = sequence.system_bu_subseq_mat[s]
		
		bucell_bu_subseq = sequence.bucell_bu_subseq_mat[s]
		flux_subseq = sequence.flux_subseq_mat[s]
		pow_dens_subseq = sequence.pow_dens_subseq_mat[s]
		# substeps length is s-1
		substeps = sequence.microsteps_number(s-1)

		file_name = cell_folder_path + '/subdens'

		write_file = open(file_name, 'w')
		txt = ''

		txt += '{:<10}'.format('TIME')
		for ss in range(substeps):
			txt += '{:^13}'.format(time_subseq[ss])

		txt += '\n'

		txt += '{:<10}'.format('SYS BU')
		for ss in range(substeps):
			txt += '{:<13.5E}'.format(system_bu_subseq[ss])

		txt += '\n'

		txt += '{:<10}'.format('CELL BU')
		for ss in range(substeps):
			txt += '{:<13.5E}'.format(bucell_bu_subseq[ss])

		txt += '\n'

		txt += '{:<10}'.format('FLUX')
		for ss in range(substeps):
			txt += '{:<13.5E}'.format(flux_subseq[ss])

		txt += '\n'

		txt += '{:<10}'.format('POW DENS')
		for ss in range(substeps):
			txt += '{:<13.5E}'.format(pow_dens_subseq[ss])

		txt += '\n'

		txt += '{:<10}'.format('POW')
		for ss in range(substeps):
			txt += '{:<13.5E}'.format(pow_dens_subseq[ss]*self.vol*1E-3)

		txt += '\n\n'

		for nucl in passport_list:
			zamid = nucl.zamid
			dens = nucl.get_dens_subseq(s)
			txt += '{:<10}'.format(zamid)

			for ss in range(substeps):
				txt += '{:<13.5E}'.format(dens[ss])

			txt += '\n'

		write_file.write(txt)
		write_file.close()

	def _print_summary_subdens(self, summary_path):

		cell_name = self.name
		passlist = self.passlist
		passport_list = passlist.passport_list
		cell_folder_path = self.folder_path
		sequence = self.sequence
		time_subseq_mat = sequence.time_subseq_mat
		system_bu_subseq_mat = sequence.system_bu_subseq_mat
		bucell_bu_subseq_mat = sequence.bucell_bu_subseq_mat
		flux_subseq_mat = sequence.flux_subseq_mat
		pow_dens_subseq_mat = sequence.pow_dens_subseq_mat
		steps_number = sequence.macrosteps_number
		file_name = summary_path + '/{}_subdens'.format(cell_name)

		write_file = open(file_name, 'w')
		txt = ''

		txt += '{:<10}'.format('TIME')
		txt += '{:^13}'.format(time_subseq_mat[0][0]/(24*3600))# in days
		for s in range(steps_number):
			substeps_number = sequence.microsteps_number(s)
			for ss in range(substeps_number):
				txt += '{:<13.5E}'.format(time_subseq_mat[s+1][ss]/(24*3600))# in days


		# Having onix storing the system bu subseq is a little laborious
		# I let it aside for now

		txt += '\n'

		txt += '{:<10}'.format('SYS BU')
		txt += '{:^13}'.format(system_bu_subseq_mat[0][0])
		for s in range(steps_number):
			substeps_number = sequence.microsteps_number(s)
			for ss in range(substeps_number):
				txt += '{:<13.5E}'.format(system_bu_subseq_mat[s+1][ss])

		txt += '\n'

		txt += '{:<10}'.format('CELL BU')
		txt += '{:^13}'.format(bucell_bu_subseq_mat[0][0])
		for s in range(steps_number):
			substeps_number = sequence.microsteps_number(s)
			for ss in range(substeps_number):
				txt += '{:<13.5E}'.format(bucell_bu_subseq_mat[s+1][ss])

		txt += '\n'

		txt += '{:<10}'.format('FLUX')
		txt += '{:^13}'.format('')
		for s in range(steps_number):
			substeps_number = sequence.microsteps_number(s)
			for ss in range(substeps_number):
				txt += '{:<13.5E}'.format(flux_subseq_mat[s+1][ss])

		txt += '\n'

		txt += '{:<10}'.format('POW DENS')
		txt += '{:^13}'.format('')
		for s in range(steps_number):
			substeps_number = sequence.microsteps_number(s)
			for ss in range(substeps_number):
				txt += '{:<13.5E}'.format(pow_dens_subseq_mat[s+1][ss])

		txt += '\n'

		txt += '{:<10}'.format('POW')
		txt += '{:^13}'.format('')
		for s in range(steps_number):
			substeps_number = sequence.microsteps_number(s)
			for ss in range(substeps_number):
				txt += '{:<13.5E}'.format(pow_dens_subseq_mat[s+1][ss]*self.vol*1E-3)

		txt += '\n\n'

		for nucl in passport_list:
			zamid = nucl.zamid

			dens = nucl.dens_subseq_mat[s+1]
			txt += '{:<10}'.format(zamid)
			init_dens = nucl.dens_seq[0]
			txt += '{:<13.5E}'.format(init_dens)
			for s in range(steps_number):
				substeps_number = sequence.microsteps_number(s)
				for ss in range(substeps_number):
					dens = nucl.dens_subseq_mat[s+1]
					txt += '{:<13.5E}'.format(dens[ss])

			txt += '\n'

		write_file.write(txt)
		write_file.close()

	def _print_summary_dens(self, summary_path):

		cell_name = self.name
		passlist = self.passlist
		passport_list = passlist.passport_list
		cell_folder_path = self.folder_path
		sequence = self.sequence
		time_seq = sequence.time_seq
		system_bu_seq = sequence.system_bu_seq
		bucell_bu_seq = sequence.bucell_bu_seq
		flux_seq = sequence.flux_seq
		pow_dens_seq = sequence.pow_dens_seq
		steps_number = sequence.macrosteps_number
		file_name = summary_path + '/{}_dens'.format(cell_name)

		write_file = open(file_name, 'w')
		txt = ''

		txt += '{:<10}'.format('TIME')
		txt += '{:<13.5E}'.format(time_seq[0]/(24*3600))# in days
		for s in range(steps_number):
			txt += '{:<13.5E}'.format(time_seq[s+1]/(24*3600))# in days

		txt += '\n'

		txt += '{:<10}'.format('SYSTEM-BU')
		txt += '{:<13.5E}'.format(system_bu_seq[0])
		for s in range(steps_number):
			txt += '{:<13.5E}'.format(system_bu_seq[s+1])

		txt += '\n'

		txt += '{:<10}'.format('CELL-BU')
		txt += '{:<13.5E}'.format(bucell_bu_seq[0])
		for s in range(steps_number):
			txt += '{:<13.5E}'.format(bucell_bu_seq[s+1])

		txt += '\n'

		txt += '{:<10}'.format('FLUX')
		txt += '{:^13}'.format('')
		for s in range(steps_number):
				txt += '{:<13.5E}'.format(flux_seq[s])

		txt += '\n'

		txt += '{:<10}'.format('POW-DENS')
		txt += '{:^13}'.format('')
		for s in range(steps_number):
			txt += '{:<13.5E}'.format(pow_dens_seq[s])

		txt += '\n'

		txt += '{:<10}'.format('POW')
		txt += '{:^13}'.format('')
		for s in range(steps_number):
				txt += '{:<13.5E}'.format(pow_dens_seq[s]*self.vol*1E-3)

		txt += '\n\n'

		for nucl in passport_list:
			zamid = nucl.zamid

			txt += '{:<10}'.format(zamid)
			init_dens = nucl.dens_seq[0]
			txt += '{:<13.5E}'.format(init_dens)
			for s in range(steps_number):
				dens = nucl.dens_seq[s+1]
				txt += '{:<13.5E}'.format(dens)

			txt += '\n'

		write_file.write(txt)
		write_file.close()


	def _print_summary_flux_spectrum(self, summary_path, mg_energy_bin):

		bucell_id = self.id
		passlist = self.passlist
		passport_list = passlist.passport_list
		sequence = self.sequence
		time_seq = sequence.time_seq
		system_bu_seq = sequence.system_bu_seq
		bucell_bu_seq = sequence.bucell_bu_seq
		flux_spectrum_seq = sequence.flux_spectrum_seq

		bin_len = [x-y for x,y in zip(mg_energy_bin[1:],mg_energy_bin[:-1])]
		mid_points = [(x+y)/2 for x,y in zip(mg_energy_bin[1:],mg_energy_bin[:-1])]

		file_name = summary_path + '/{}_flux_spectrum'.format(self.name)
		write_file = open(file_name, 'w')

		txt = '{:<34}'.format('ENERGY-BIN')
		for val in mg_energy_bin:
			txt += '{:^18.10E}'.format(val)
		txt += '\n'
		txt += '{:<52}'.format('ENERGY-BIN-LEN')
		for val in bin_len:
			txt += '{:^18.10E}'.format(val)
		txt += '\n'
		txt += '{:<52}'.format('ENERGY-MID-POINTS')
		for val in mid_points:
			txt += '{:^18.10E}'.format(val)
		txt += '\n\n'

		flux_spectrum_seq_lethargy = []
		for i in range(len(flux_spectrum_seq)):
			flux_spectrum_seq_lethargy.append([x*y/z for x,y,z in zip(flux_spectrum_seq[i], mid_points, bin_len)])

		txt += '  TIME    SYSTEM-BU    BUCELL-BU  \n\n'

		for i in range(len(time_seq)-1):
			time = time_seq[i]/(24*3600) # in day
			system_bu = system_bu_seq[i]
			bucell_bu = bucell_bu_seq[i]

			txt += '{:^8.6}'.format(float(time))
			txt += '{:^13.6}'.format(float(system_bu))
			txt += '{:^13.6}'.format(float(bucell_bu))
			txt += '{:18}'.format('')

			for val in flux_spectrum_seq_lethargy[i]:
				txt += '{:^18.10E}'.format(val)

			txt += '\n'


		write_file.write(txt)
		write_file.close()


	def _print_summary_xs(self, summary_path):

		bucell_id = self.id
		passlist = self.passlist
		passport_list = passlist.passport_list
		sequence = self.sequence
		time_seq = sequence.time_seq
		system_bu_seq = sequence.system_bu_seq
		bucell_bu_seq = sequence.bucell_bu_seq

		file_name = summary_path + '/{}_xs_lib'.format(self.name)
		write_file = open(file_name, 'w')

		# Write time, system burnup and cell burnup
		txt = '{:<51}'.format('TIME')
		for time in time_seq[1:]: # initial time is discarded as no XS are defined for it
			time = time/(24*3600)
			txt += '{:^17}'.format(time)
		txt += '\n'

		txt += '{:<51}'.format('SYSTEM-BU')
		for bu in system_bu_seq[1:]: # initial time is discarded as no XS are defined for it
			txt += '{:^17.5E}'.format(bu)
		txt += '\n'

		txt += '{:<51}'.format('BUCELL-BU')
		for bu in bucell_bu_seq[1:]: # initial time is discarded as no XS are defined for it
			txt += '{:^17.5E}'.format(bu)
		txt += '\n'

		txt += '\n\n--- Actinides ---\n\n'
		#txt += utils.printer.xs_lib_header

		#loop over actinides
		for nucl in self.get_act_passport_list():
			nucl_name = nucl.name
			nucl_zamid = nucl.zamid
			nucl_xs_seq = nucl.xs_seq
			# If this nuclide has not been assigned cross section, continue
			if nucl_xs_seq == None:
				continue
			count = 0
			sample_dict = nucl_xs_seq[0] # the xs names should be the same at each step
			for xs_name in sample_dict:
				# removal should not be printed on the lib
				if xs_name ==  'removal':
					continue
				if count == 0:
					txt += '{:^17}'.format(nucl_name)
				else:
					txt += '{:^17}'.format('')
				txt += '{:^17}'.format(nucl_zamid)
				txt += '{:^17}'.format(xs_name)
				for xs_dict in nucl_xs_seq:
					xs_val = xs_dict[xs_name][0]
					txt += '{:^17.8E}'.format(xs_val)

				txt += '\n'
				count += 1

		txt += '\n\n--- Fission Products ---\n\n'
		txt += utils.printer.xs_lib_header

		#loop over fission products
		for nucl in self.get_fp_passport_list():
			nucl_name = nucl.name
			nucl_zamid = nucl.zamid
			nucl_xs_seq = nucl.xs_seq
			# If this nuclide has not been assigned cross section, continue
			if nucl_xs_seq == None:
				continue
			count = 0
			sample_dict = nucl_xs_seq[0] # the xs names should be the same at each step
			for xs_name in sample_dict:
				# removal should not be printed on the lib
				if xs_name ==  'removal':
					continue
				if count == 0:
					txt += '{:^17}'.format(nucl_name)
				else:
					txt += '{:^17}'.format('')
				txt += '{:^17}'.format(nucl_zamid)
				txt += '{:^17}'.format(xs_name)
				for xs_dict in nucl_xs_seq:
					xs_val = xs_dict[xs_name][0]
					txt += '{:^17.8E}'.format(xs_val)

				txt += '\n'
				count += 1

		txt += '\n\n--- Activation Products ---\n\n'
		txt += utils.printer.xs_lib_header

		#loop over activation products
		for nucl in self.get_avt_passport_list():
			nucl_name = nucl.name
			nucl_zamid = nucl.zamid
			nucl_xs_seq = nucl.xs_seq
			# If this nuclide has not been assigned cross section, continue
			if nucl_xs_seq == None:
				continue
			count = 0
			sample_dict = nucl_xs_seq[0] # the xs names should be the same at each step
			for xs_name in sample_dict:
				# removal should not be printed on the lib
				if xs_name ==  'removal':
					continue
				if count == 0:
					txt += '{:^17}'.format(nucl_name)
				else:
					txt += '{:^17}'.format('')
				txt += '{:^17}'.format(nucl_zamid)
				txt += '{:^17}'.format(xs_name)
				for xs_dict in nucl_xs_seq:
					xs_val = xs_dict[xs_name][0]
					txt += '{:^17.8E}'.format(xs_val)

				txt += '\n'
				count += 1

		write_file.write(txt)
		write_file.close()

	# Decay lib shoult not change. But I still implement an on the fly get

	def _print_summary_isomeric_branching_ratio(self, summary_path):

		bucell_id = self.id
		sequence = self.sequence
		time_seq = sequence.time_seq
		system_bu_seq = sequence.system_bu_seq
		bucell_bu_seq = sequence.bucell_bu_seq
		isomeric_branching_ratio_seq = sequence.isomeric_branching_ratio_seq


		file_name = summary_path + '/{}_isomeric_branching_ratio'.format(self.name)
		write_file = open(file_name, 'w')

		# Write time, system burnup and cell burnup
		txt = '{:<17}'.format('TIME')
		for time in time_seq[1:]: # initial time is discarded as no XS are defined for it
			time = time/(24*3600)
			txt += '{:^17}'.format(time)
		txt += '\n'

		txt += '{:<17}'.format('SYSTEM-BU')
		for bu in system_bu_seq[1:]: # initial time is discarded as no XS are defined for it
			txt += '{:^17.5E}'.format(bu)
		txt += '\n'

		txt += '{:<17}'.format('BUCELL-BU')
		for bu in bucell_bu_seq[1:]: # initial time is discarded as no XS are defined for it
			txt += '{:^17.5E}'.format(bu)
		txt += '\n\n'

		#loop over nucl in isomeric branching 
		unordered_nucl_list = []
		for nucl in isomeric_branching_ratio_seq[0]:
			unordered_nucl_list.append(nucl)

		ordered_nucl_list = utils.order_nuclide_name_per_z(unordered_nucl_list)

		for nucl in ordered_nucl_list:

			txt += '{:<17}'.format(nucl)

			# First loop over the ground state ratio
			for isomeric_branching_ratio in isomeric_branching_ratio_seq:

				nucl_data = isomeric_branching_ratio[nucl]
				ratio_to_gs = nucl_data['(n,gamma)']
				txt += '{:^17.4}'.format(ratio_to_gs)

			txt += '\n'

			# then loop over the metastable state ratio
			txt += '{:<17}'.format('')
			for isomeric_branching_ratio in isomeric_branching_ratio_seq:

				nucl_data = isomeric_branching_ratio[nucl]
				ratio_to_ms = nucl_data['(n,gamma)X']
				txt += '{:^17.4}'.format(ratio_to_ms)

			txt += '\n'

		write_file.write(txt)
		write_file.close()

	@property
	def get_decay_nucl(self):

		decay_lib = self._decay_a_lib
		decay_nucl = utils.get_decay_nucl(decay_lib)

		return decay_nucl



	# xs lib can change so I implement a get on the fly


	@property
	def get_xs_nucl(self):

		xs_lib = self._xs_lib
		xs_nucl = utils.get_xs_nucl(xs_lib)

		return xs_nucl


	# fy lib can change so I implement a get on the fly

	def get_fy_nucl(self):

		fy_lib = self._fy_lib
		fy_nucl = utils.get_fy_nucl(fy_lib)

		return fy_nucl




	def get_lib_nucl(self):

		decay_a_lib = self._decay_a_lib

		# In couple mode, the xs_lib is only set after MC simulation. Therefore, it would only be possible
		# to build lib_nucl after MC call. This is not practical as lib_nucl is needed before the loop begins
		# Instead, MC_XS_nucl_list is used. It is not a lib but a list of nuclides.
		# This change is however not compatible if user wants to use constant cross section as lib nucl will rely
		# only on MC_nucl_list for nuclides with xs data
		#xs_lib = self._xs_lib
		MC_XS_nucl_list = self.MC_XS_nucl_list

		fy_lib = self.fy_lib

		decay_fy_lib_list = [decay_a_lib, fy_lib]

		#if decay_a_lib == None and xs_lib == None and fy_lib == None:
		if decay_a_lib == None and fy_lib == None:
			raise No_nuclear_lib_set('Cell {} has not been attributed decay and fy data library'.format(self.id))

		decay_fy_lib_list = [dic for dic in decay_fy_lib_list if dic != None]

		decay_fy_lib_nucl = utils.get_all_nucl(decay_fy_lib_list)
		#print (decay_fy_lib_nucl)
		unfiltered_lib_nucl = decay_fy_lib_nucl + MC_XS_nucl_list
		lib_nucl = list(set(unfiltered_lib_nucl))


		return lib_nucl


	def get_fy_parent_nucl(self):

		fy_lib = self._fy_lib
		fy_parent = utils.get_fy_parent_nucl(fy_lib)

		return fy_parent


	# This set somehow uptset the on the fly nature of the 
	# various nuclide list. Will have to deal with that later

	# def _set_nucl_list(self):

	# 	decay_a_lib = self._decay_a_lib
	# 	xs_lib = self._xs_lib
	# 	fy_lib = self._fy_lib

	# 	self._decay_nucl = self.get_decay_nucl
	# 	self._xs_nucl = self.get_xs_nucl
	# 	self._fy_nucl = self.get_fy_nucl
	# 	self._fy_parent_nucl = self.get_fy_parent_nucl
	# 	self._all_nucl = self.get_all_nucl


	def _set_MC_tallies(self, mc_nuclide_densities, flux_tally, flux_spectrum_tally, rxn_rate_tally, sampled_isomeric_branching_data, sampled_ng_cross_section_data, xs_mode, s):

		MC_flux = flux_tally.mean[0][0][0]
		self.sequence._set_macrostep_MC_flux(MC_flux)

		flux_spectrum = [x[0][0] for x in flux_spectrum_tally.mean]
		self.sequence._set_macrostep_flux_spectrum(flux_spectrum)

		self._set_step_isomeric_branching_ratio(flux_spectrum, sampled_isomeric_branching_data, sampled_ng_cross_section_data)

		# rxn_rate_tally are distributed in a xs_dict object
		xs_lib = utils.xs_lib('{} rxn rate'.format(self.name))
		nuclides = rxn_rate_tally.nuclides
		scores = rxn_rate_tally.scores
		macro_xs = rxn_rate_tally / flux_tally
		passlist = self.passlist
		nucl_dict = passlist._get_name_passport_dict()
		isomeric_branching_ratio = self.sequence.current_isomeric_branching_ratio


		for nucl in nuclides:

			onix_nucl_name = utils.openmc_name_to_onix_name(nucl)
			nucl_passport = nucl_dict[onix_nucl_name]

			# If the nuclide is artificially added to openmc material with 1 atm
			# Its density set for onix is 0
			# But the macro xs need to be divided by 1E-24 and not 0

			# For
			if nucl_passport.current_dens == 0.0:

				# mc_nuclide_densities gives a tuples where the first element is the nuclide name
				# and the sedonc element is the density
				#nucl_dens = mc_nuclide_densities[nucl][1]

				nucl_dens = self.zero_dens_1_atm
			else:
				nucl_dens = nucl_passport.current_dens
			xs_dict = {}
			for score in scores:
				macro_xs_val = macro_xs.get_values(scores = ['({} / flux)'.format(score)], nuclides = ['({} / total)'.format(nucl)])[0][0][0]
				xs_val = macro_xs_val/nucl_dens
				xs_dict[score] = xs_val
			xs_lib.add_xs_dict(nucl_passport.zamid, xs_dict)

		# This function screens the the xs_lib and look for which nuclide ngamma reactions needs to 
		# be branched
		xs_lib.isomeric_branching_weighting(isomeric_branching_ratio)

		# Here, in constant lib mode, if it is the first step, you need to call 
		# overwrite xs (because the xs will have been already set by constant lib)
		if xs_mode == 'constant lib' and s == 1:
			self.overwrite_xs(xs_lib)
		else:
			self.set_xs(xs_lib)

	def _set_step_isomeric_branching_ratio(self, flux_spectrum, sampled_isomeric_branching_data, sampled_ng_cross_section_data):

		isomeric_branching_ratio = {}
		for nucl in sampled_isomeric_branching_data:
			if nucl in sampled_ng_cross_section_data:
				isomeric_branching_ratio[nucl] = {}
				branching_data = sampled_isomeric_branching_data[nucl]
				xs_data = sampled_ng_cross_section_data[nucl]

				numerator_ground = sum([x*y*z for x,y,z in zip(flux_spectrum,branching_data['0'],xs_data)])
				numerator_excited = sum([x*y*z for x,y,z in zip(flux_spectrum,branching_data['1'],xs_data)])
				denominator = sum([x*y for x,y in zip(flux_spectrum,xs_data)])
				isomeric_branching_ratio[nucl]['(n,gamma)'] = numerator_ground/denominator
				isomeric_branching_ratio[nucl]['(n,gamma)X'] = numerator_excited/denominator

		self.sequence._set_macrostep_isomeric_branching_ratio(isomeric_branching_ratio)
		

	def _set_allreacs_dic(self, s, ss, ssn):


		passlist = self.passlist
		passport_list = passlist.passport_list
		passport_dict = passlist._get_zamid_passport_dict()

		sequence = self.sequence

		#flux = sequence.flux_point(s)
		flux = sequence.current_flux
		N = len(passport_list)

		# EOS = 1 means that the sequence reached end of step
		EOS = 0
		if ss == ssn - 1:
			EOS = 1

		for row in range(N):
			nuc_pass = passport_list[row]
			nuc_zamid = nuc_pass.zamid
			nuc_name = nuc_pass.name
			nuc_dens = nuc_pass.current_dens
			creation_dic = {}
			destruction_dic = {}
			allreacs_dic = {}

			if nuc_pass.decay_a != None and nuc_pass.decay_a != 'stable':
				decay = nuc_pass.decay_a
				decay_rate = decay.copy()
				del decay_rate['total decay']
				del decay_rate['half-life']
				for i in decay_rate:
					decay_rate[i] = decay_rate[i]*nuc_dens
				destruction_dic = decay_rate.copy()



			# Multiply xs by flux
			if nuc_pass.current_xs != None:
				xs = nuc_pass.current_xs
				xs_rate = {}
				for i in xs:
					xs_rate[i] = xs[i][0]*flux*1e-24*nuc_dens
				destruction_dic.update(xs_rate)
				del destruction_dic['removal']

			# Multiply fy by xs fission of father and phi
			if nuc_pass.fy != None:
				fy = nuc_pass.fy
				fy_rate = {}
				for i in fy:
					if i in passport_dict:
						father_pass = passport_dict[i]
						father_dens = father_pass.current_dens
						father_name = father_pass.name
						if father_pass.current_xs == None:
							continue
						if 'fission' not in father_pass.current_xs:
							continue
						father_fission_xs = father_pass.current_xs['fission'][0]
						entry = '{} fission'.format(father_name)
						fy_rate[entry] = fy[i][0]*1e-2*father_fission_xs*flux*1e-24*father_dens
				creation_dic = fy_rate.copy()

			# Now we gather the creation terms (excluding fission as it has been already collected in fyxsphi)
			xs_parent = nuc_pass.xs_parent
			decay_parent = nuc_pass.decay_parent

			xs_parent_val = {}
			for i in xs_parent:
				father_zamid = xs_parent[i]
				if father_zamid in passport_dict:
					father_pass = passport_dict[father_zamid]
					father_dens = father_pass.current_dens
					father_state = father_pass.state
					father_name = father_pass.name
					# The reacname from xs_parent is different from parent_xs in case parent is in excited state
					if father_state == 1:
						reac_name = i[1:]
					elif father_state == 0:
						reac_name = i
					if father_pass.current_xs != None:
						if reac_name in father_pass.current_xs:
							father_xs_rate = father_pass.current_xs[reac_name][0]*flux*1e-24*father_dens

							entry = '{} {}'.format(father_name, reac_name)

							xs_parent_val[entry] = father_xs_rate
			creation_dic.update(xs_parent_val)

			decay_parent_val = {}
			for i in decay_parent:
				father_zamid = decay_parent[i]
				if father_zamid in passport_dict:
					father_pass = passport_dict[father_zamid]
					father_dens = father_pass.current_dens
					father_state = father_pass.state
					father_name = father_pass.name
					# The reacname from xs_parent is different from parent_xs in case parent is in excited state
					if father_state == 1:
						reac_name = i[1:]
					elif father_state == 0:
						reac_name = i
					if father_pass.decay_a != None and father_pass.decay_a != 'stable':
						if reac_name in father_pass.decay_a:
							father_decay = father_pass.decay_a[reac_name]*father_dens

							entry = '{} {}'.format(father_name, reac_name)

							decay_parent_val[entry] = father_decay
			creation_dic.update(decay_parent_val)

			allreacs_dic = destruction_dic.copy()
			allreacs_dic.update(creation_dic)

			nuc_pass.destruction_dic = destruction_dic
			nuc_pass.creation_dic = creation_dic
			nuc_pass.allreacs_dic = allreacs_dic
			nuc_pass.allreacs_dic_list_append(allreacs_dic)

			sorted_allreacs = sorted(list(allreacs_dic.items()), key=operator.itemgetter(1))

			nuc_pass.append_current_sorted_allreacs_tuple_list(sorted_allreacs, ss)

			if EOS == 1:
				nuc_pass.append_sorted_allreacs_tuple_mat()

	def _print_current_allreacs_rank(self):

		passport_list = self.passlist.passport_list
		cell_id = self.id
		cell_folder_path = self.folder_path

		file_name = cell_folder_path +'/cell_{}_reacs_rank'.format(cell_id)
		file = open(file_name, 'w')

		txt = ''
		for nucl in passport_list:
			reacs_rank = nucl.current_sorted_allreacs_tuple_list
			txt += '\n{}({})\n'.format(nucl.name, nucl.zamid)
			for ss in range(len(reacs_rank)):
				txt += 'substep = {} | dens = {}\n'.format(ss, nucl.get_current_dens_subseq()[ss])
				txt += '{}\n'.format(reacs_rank[ss])

		file.write(txt)
		file.close()

	def _print_summary_allreacs_rank(self, summary_path):

		passport_list = self.passlist.passport_list
		cell_name = self.name

		file_name = summary_path +'/cell_{}_reacs_rank'.format(cell_name)
		file = open(file_name, 'w')

		txt = ''
		for nucl in passport_list:
			reacs_rank = nucl.sorted_allreacs_tuple_mat
			txt += '\n ==={}({}) ===\n'.format(nucl.name, nucl.zamid)
			for s in range(len(reacs_rank)):
				txt += '\nSTEP {}\n'.format(s+1)
				for ss in range(len(reacs_rank[s])):
					# subseq has s+1 elements because of the initial element 
					txt += 'substep = {} | dens = {}\n'.format(ss, nucl.dens_subseq_mat[s+1][ss])
					txt += '{}\n'.format(reacs_rank[s][ss])

		file.write(txt)
		file.close()

	def _reduce_nucl_set(self):

		self._set_all_leaves()
		total_leaves = self._total_leaves
		NAX_nucl_list = data.NAX_nucl_list
		reduced_nucl_set = list(set(total_leaves+NAX_nucl_list))
		ordered_reduced_set = utils.order_nuclide_per_z(reduced_nucl_set)


	def _set_all_leaves(self):

		start_time = time.time()
		self._set_tree()
		self._set_leaves()
		self._set_fission_tree()
		self._set_fission_leaves()
		run_time = time.time() - start_time
		total_leaves = list(set(self._leaves+self._fission_leaves))
		self._total_leaves = total_leaves


	@property
	def get_tree(self):

		nuc_list = self.nucl_set
		initial_nuc_list = self._initial_nucl
		# print (nuc_list)
		# print (initial_nuc_list)
		trees_dic = {}
		passdic = self.passlist._get_zamid_passport_dict()

		# This list will dynamically store the remaining nuclide to allocate in the tree.
		# If a nuclide has already been put in the tree but there is a chain loop, because the nuclide will not appear anymore in this list
		# the infinite looping will be prevented
		# Here we remove the initial nuclides
		#reduced_nuc_list = [x for x in nuc_list if x not in initial_nuc_list]
		reduced_nuc_list = list(nuc_list)

		# Initialize the dictionary where each entry is the name of the initial nuclide and the first element of the list is their direct child
		for nuc in initial_nuc_list:
			nuc_pass = passdic[nuc]
			nuc_zamid = nuc_pass.zamid
			nuc_name = nuc_pass.name
			
			reduced_nuc_list.remove(nuc_zamid)

			nuc_non0_child = nuc_pass.get_all_non0_child()
			nuc_child = list(set(nuc_non0_child).intersection(nuc_list))
			trees_dic[nuc_zamid] = [nuc_child]

			reduced_nuc_list_copy = list(reduced_nuc_list)
			reduced_nuc_list = [x for x in reduced_nuc_list_copy if x not in nuc_child]

			# Build the tree for each initial nuclide
			tree = trees_dic[nuc_zamid]

			end_of_tree = 0
			gen = 0

			# As long as there are child that exists
			while end_of_tree == 0:
				gen += 1
				next_child_gen = []
				for parent in tree[gen-1]:
					parent_pass = passdic[parent]
					parent_child = parent_pass.get_all_non0_child()

					reduced_parent_child = list(set(parent_child).intersection(reduced_nuc_list))

					reduced_nuc_list_copy = list(reduced_nuc_list)
					reduced_nuc_list = [x for x in reduced_nuc_list_copy if x not in reduced_parent_child]

					next_child_gen += reduced_parent_child

				if next_child_gen == []:
					end_of_tree = 1
					#reduced_nuc_list = 
					break

				tree.append(next_child_gen)

			reduced_nuc_list = list(nuc_list)

		return trees_dic

	def _set_tree(self):

		self._tree = self.get_tree





	@property
	def leaves(self):

		return self._leaves

	def _gen_leaves(self):

		tree_dic = self._tree

		conca_tree = []
		for initial_nuc in tree_dic:
			conca_tree += [initial_nuc]
			for gen in tree_dic[initial_nuc]:
					conca_tree += gen

		# Set is a set of unique element in the list. List convert that back into a new list
		conca_tree_no_redundance = list(set(conca_tree))

		return conca_tree_no_redundance

	def _set_leaves(self):

		self._leaves = self._gen_leaves()





# This trees are the one built from the fissile nuclide through their fission products and their own trees
# This trees need the trees_dic from build tree to know which fissile nuclide are actually in the original tree (i.e. which 
# fissile material have a chance to be created)



	@property
	def fission_tree(self):

		return self._fission_tree

	def _gen_fission_tree(self):

		passlist = self.passlist
		passdic = self.passlist._get_zamid_passport_dict()
		initial_nuc_list = self._initial_nucl
		nuc_list = self.nucl_set
		leaves = self._leaves
		fy_parent = self.get_fy_parent_nucl()
		fy_nucl = self.get_fy_nucl()
		fission_trees_dic = {}

		passlist._set_fission_child(fy_nucl, fy_parent)

		reduced_nuc_list = list(nuc_list)

		# First we enter the initial nuc that are also fissile
		for parent in fy_parent:
			if parent in leaves:
				nuc_pass = passdic[parent]
				nuc_zamid = nuc_pass.zamid

				reduced_nuc_list.remove(nuc_zamid)

				fission_child = nuc_pass.fission_child

				reduced_fission_child = list(set(fission_child).intersection(reduced_nuc_list))

				fission_trees_dic[parent] = [reduced_fission_child]

				reduced_nuc_list_copy = list(reduced_nuc_list)
				reduced_nuc_list = [x for x in reduced_nuc_list_copy if x not in reduced_fission_child]

			# Then we build the tree from the fission products
		#	for initial_nuc in fission_trees_dic:
				tree = fission_trees_dic[parent]
				end_of_tree = 0

				gen = 0

				# As long as there are child that exists
				while end_of_tree == 0:
					gen += 1
					next_child_gen = []
					for parent in tree[gen-1]:
						parent_pass = passdic[parent]
						parent_child = parent_pass.get_all_non0_child()

						reduced_parent_child = list(set(parent_child).intersection(reduced_nuc_list))

						reduced_nuc_list_copy = list(reduced_nuc_list)
						reduced_nuc_list = [x for x in reduced_nuc_list_copy if x not in reduced_parent_child]

						next_child_gen += reduced_parent_child

					if next_child_gen == []:
						end_of_tree = 1
						#reduced_nuc_list = 
						break

					tree.append(next_child_gen)

				reduced_nuc_list = list(nuc_list)

		return fission_trees_dic

	def _set_fission_tree(self):

		self._fission_tree = self._gen_fission_tree()





	@property
	def fission_leaves(self):

		return self._fission_leaves


	def _gen_fission_leaves(self):

		fission_tree_dic = self._fission_tree

		conca_fission_tree = []
		for initial_nuc in fission_tree_dic:
			conca_fission_tree += [initial_nuc]
			for gen in fission_tree_dic[initial_nuc]:
					conca_fission_tree += gen

		# Set is a set of unique element in the list. List convert that back into a new list
		conca_fission_tree_no_redundance = list(set(conca_fission_tree))

		return conca_fission_tree_no_redundance

	def _set_fission_leaves(self):

		self._fission_leaves = self._gen_fission_leaves()


	def _print_tree(self):

		folder_path = self.folder_path + '/tree_&_leaves'
		if os.path.exists(folder_path):
			shutil.rmtree(folder_path)

		os.makedirs(folder_path)

		tree = self.get_tree

		for branch_zamid in tree:

			branch_zamid_pass = pp(branch_zamid)
			branch_name = branch_zamid_pass.name
			branch_path = folder_path + '/{}_branch'.format(branch_name)
			txt =''
			branch = tree[branch_zamid]
			count = 0

			for leaves in branch:
				txt += '{}\n{}\n\n'.format(count, leaves)
				count = count + 1

			f = open(branch_path, 'w')
			f.write(txt)
			f.close()

		fission_tree = self.fission_tree

		for branch_zamid in fission_tree:

			branch_zamid_pass = pp(branch_zamid)
			branch_name = branch_zamid_pass.name

			branch_path = folder_path + '/{}_fission_branch'.format(branch_name)
			txt =''
			branch = fission_tree[branch_zamid]
			count = 0

			for leaves in branch:
				txt += '{}\n{}\n\n'.format(count, leaves)
				count = count + 1

			f = open(branch_path, 'w')
			f.write(txt)
			f.close()




	def _gen_allreacs_ranking(self):

		allreacs_ranking = pl._gen_allreacs_dic()


	@property
	def folder_path(self):

		return self._folder_path

	# Create a folder for the cell in the directory indicated in the argument
	# If no directory is passed as argument, the default directory is the current directory

	def gen_folder(self):

		cell_name = self.name

		utils.gen_cell_folder(cell_name)

	def _set_folder(self):

		cell_name = self.name

		utils.gen_cell_folder(cell_name)
		self._folder_path = utils.get_cell_folder_path(cell_name)

	def copy_cell_folders_to_step_folder(self, s):

		cell_folder_path = self.folder_path
		shutil.copytree(cell_folder_path, os.getcwd() + '/step_{}/{}_cell'.format(s, self.name))
		shutil.rmtree(cell_folder_path)

	# def _gen_output_summary_folder(self):

	# 	name = 'output_summary'
	# 	self._output_summary_path = utils.get_folder_path(name)
	# 	utils.gen_folder(name)


class Initial_nucl_not_set(Exception):
	"""Raise when the user forgot to set the initial nuclide of the cell and tries to burn cell"""
	pass

class Nucl_set_not_in_Lib_nucl(Exception):
	"""Raise when the user forgot to set the initial nuclide of the cell and tries to burn cell"""
	pass

class Initial_nucl_not_in_Nucl_set(Exception):
	"""Raise when the user forgot to set the initial nuclide of the cell and tries to burn cell"""
	pass

class Nuclide_list_redundant(Exception):

	pass

class Passlist_not_defined(Exception):
	"""Raise when the user forgot to defined passlist for a cell"""
	pass