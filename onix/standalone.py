import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
import copy
import xml.etree.ElementTree as ET
import glob
import pdb
import time

import openmc
import openmc.mgxs as mgxs
from onix.cell import Cell
from onix.system import System
from onix import salameche

from onix import utils
from onix import data

class Stand_alone(object):
	"""This class is used to execute standalone simulations

	Through this class the user can set the nuclear data libraries used for the simulation,
	manually set the volumes of each BUCell, set the burnup/time sequence and finally
	launch a simulation with the method "burn" """

	def __init__(self):

		self._fy_lib_set = 'no'
		self._decay_lib_set = 'no'
		self._xs_lib_set = 'no'

		self.system = System(1)

	@property
	def system(self):
		"""Returns an instantiation of the class system"""
		return self._system

	@system.setter
	def system(self, system):
		"""Sets an instantiation of the class system"""
		self._system = system

	def add_bucell(self, bucell):
		"""Adds a BUCell to the system"""
		system = self.system
		system.add_bucell(bucell)

	@property
	def total_vol(self):
		"""Returns the total volume of the system"""
		return self._total_vol

	@total_vol.setter
	def total_vol(self, total_vol):
		"""Sets the total volume of the system"""
		self._total_vol = total_vol

	def set_sequence(self, sequence):
		"""Sets the burnup/time sequence to the system"""
		self.sequence = sequence
		system = self.system
		system.set_sequence(sequence, mode = 'stand alone')

	def set_decay_lib(self, decay_lib_path):
		"""Sets a decay library chosen by the user that will be used in the simulation

		The user needs to specify the path of the chosen library"""
		system = self.system
		self._decay_lib_set = 'yes'
		self._decay_lib_path = decay_lib_path
		system.set_decay_for_all(decay_lib_path)


	def set_default_decay_lib(self):
		"""Sets the decay library to the default decay library (ENDF/B-VIII)"""
		system = self.system
		self._decay_lib_set = 'yes'
		#system.set_default_decay_for_all_no_add()
		self._decay_lib_path = 'default'
		system.set_default_decay_for_all()

	def set_decay_from_object(self, bucell, object):
		"""Sets the decay library from a decay object created by the user"""
		system = self.system
		# This should not set yes since it is only for one bucell
		# Need to be fixed later
		self._decay_lib_set = 'yes'
		self._decay_lib_path = 'decay_lib object defined by the user'
		bucell = system.get_bucell(bucell)
		bucell.set_decay(object)

	def set_xs_lib(self, xs_lib_path):
		"""Sets a cross section library chosen by the user that will be used in the simulation

		The user needs to specify the path of the chosen library"""
		system = self.system
		self._xs_lib_set = 'yes'
		self._xs_lib_path = xs_lib_path
		system.set_xs_for_all(xs_lib_path)

	def set_default_xs_lib(self):
		"""Sets the cross section library to the default cross section library (ENDF/B-VIII)"""
		system = self.system
		self._xs_lib_set = 'yes'
		#system.set_default_decay_for_all_no_add()
		system.set_default_xs_for_all()

	def set_xs_from_object(self, bucell, object):
		"""Sets the cross section library from a cross section object created by the user"""
		system = self.system
		# This should not set yes since it is only for one bucell
		# Need to be fixed later
		self._xs_lib_set = 'yes'
		self._xs_lib_path = 'xs_lib object defined by the user'
		bucell = system.get_bucell(bucell)
		bucell.set_xs(object)

	def set_fy_lib(self, fy_lib_path):
		"""Sets a fission yield library chosen by the user that will be used in the simulation

		The user needs to specify the path of the chosen library"""
		system = self.system
		self._fy_lib_set = 'yes'
		self._fy_lib_path = fy_lib_path
		system.set_fy_for_all(fy_lib_path)

	def set_default_fy_lib(self):
		"""Sets the fission yields library to the default fission yield library (ENDF/B-VIII)"""
		system = self.system
		self._fy_lib_set = 'yes'
		self._fy_lib_path = 'default'
		#system.set_default_fy_for_all_no_add()
		system.set_default_fy_for_all()

	def set_fy_from_object(self, bucell, object):
		"""Sets the fission yield library from a fission yeild object created by the user"""
		system = self.system
		# This should not set yes since it is only for one bucell
		# Need to be fixed later
		self._fy_lib_set = 'yes'
		self._fy_lib_path = 'fy_lib object defined by the user'
		bucell = system.get_bucell(bucell)
		bucell.set_fy(object)

	def set_vol(self, vol_dict):
		"""Sets the volume of each BUCell by providing a dictionnary of BUCell volumes
		where each key is the name of the BUCell and the entries are the volume in cm^3"""
		system = self.system
		bucell_dict = system.bucell_dict

		# Need to loop over bucell_dict because there might be more cells than bucells
		for i in bucell_dict:
			bucell = bucell_dict[i]
			if bucell.name in vol_dict:
				bucell.vol = vol_dict[bucell.name]

		# We treat total volume separately
		system.total_vol = vol_dict['total volume']

		self._volume_set = 'yes'

	def _step_normalization(self, s):

		system = self.system
		bucell_list = system.get_bucell_list()
		master_sequence = system.sequence
		system_flux = master_sequence._av_flux_seq[s]
		for bucell in bucell_list:
			bucell_sequence = bucell.sequence
			pow_dens = bucell._update_pow_dens(system_flux)
			bucell_sequence._set_macrostep_flux(system_flux)
			bucell_sequence._set_macrostep_pow_dens(pow_dens)


	def burn(self):
		"""Burns the system

		The method burn will set default nuclear data to the system if the user has not
		set nuclear data already

		A macrostep folder named "step_s" (with s the index of the macrostep) will be created per macrostep and
		contains various information on the system for the corresponding macrostep

		For each macrostep, burn calls salameche.burnstep to deplete the system until the next macrostep

		At the end of the simulation, burn will print various information on the system in the output_summary folder
		 """
		start_time = time.time()

		# If no decay libs and fy libs have been set, set default libs
		if self._decay_lib_set == 'no':
			self.set_default_decay_lib()
			print ('\n\n\n----  Default decay constants library set for system  ----\n')
		else:
			print ('\n\n\n----  User defined path for decay library  ----\n\n')
			print ('----  {}  ----\n\n\n'.format(self._decay_lib_path))
		
		if self._fy_lib_set == 'no':
			self.set_default_fy_lib()
			print ('\n\n\n----  Default fission yields library set for system  ----\n\n\n')
		else:
			print ('\n\n\n----  User defined path for fission yields library ----\n\n')
			print ('----  {}  ----\n\n\n'.format(self._fy_lib_path))
		
		if self._xs_lib_set == 'no':
			self.set_default_xs_lib()
			print ('\n\n\n----Default cross section library set for system----\n\n\n')
		else:
			print ('\n\n\n----  Path for cross sections library ----\n\n')
			print ('----  {}  ----\n\n\n'.format(self._xs_lib_path))

		system = self.system
		sequence = system.sequence
		norma_mode = sequence.norma_unit

		# This should be somewhere else but for now it is done here
		system.zam_order_passlist()

		macrosteps_number = sequence.macrosteps_number
		# Shift loop from 1 in order to align loop s and step indexes
		for s in range(1, macrosteps_number+1):

			print ('\n\n\n\n====== STEP {}======\n\n\n\n'.format(s))
			sequence._gen_step_folder(s)
			self._step_normalization(s)
			print (('\n\n\n=== Salameche Burn {}===\n\n\n'.format(s)))
			salameche.burn_step(system, s, 'stand alone')

			# To develop. Basially update bu against time when doing constant flux
			#sequence.dynamic_system_time_bu_conversion(system)
		
		
		system._gen_output_summary_folder()
		system._print_summary_allreacs_rank()
		system._print_summary_subdens()
		system._print_summary_dens()
		system._print_summary_param()

		run_time = time.time() - start_time
		print ('\n\n\n >>>>>> ONIX burn took {} seconds <<<<<<< \n\n\n'.format(run_time))
	