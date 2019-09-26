""" This module defines multiple Python class that are designed to be used by the user when using
the Python environment to define and launch an onix calculation"""
import math as m
from onix.data import list_and_dict
import onix.data as d
from .functions import *

class decay_lib(object):

	def __init__(self, id_number):

		self._id = id_number

		self._dic = {}
		self._decay_a = {}
		self._decay_b = {}


	def add_data(self, zamid, **kwargs):

		if not kwargs:
			raise Empty_data("No data for {} has been entered".format(zamid))

		self._dic[zamid] = kwargs

		if 'half_life' in kwargs:
			if kwargs['half_life'] == 'stable':
				self._dic[zamid]['total decay'] = 0
			else:
				total = m.log(2)/kwargs['half_life']
				self._dic[zamid]['total decay'] = total

		elif 'total' in kwargs:
			if kwargs['total decay'] == 0:
				self._dic[zamid]['half_life'] = 'stable'
			else:
				half_life = m.log(2)/kwargs['total decay']
				self._dic[zamid]['half_life'] = half_life

		elif 'total' not in kwargs and 'half_life' not in kwargs:
			total = 0
			for i in kwargs:
				total += kwargs[i]
			self._dic[zamid]['total decay'] = total

			half_life = m.log(2)/total
			self._dic[zamid]['half_life'] = half_life

		self.create_decay_a(zamid, self._dic[zamid])
		self.create_decay_b(zamid, self._dic[zamid])



	@property
	def dic(self):

		return self._dic

	def create_decay_a(self, zamid, dic):

		self._decay_a[zamid] = dic.copy()
		self._decay_a[zamid]['half-life'] = self._decay_a[zamid]['half_life']
		del self._decay_a[zamid]['half_life']

	@property
	def decay_a(self):

		return self._decay_a

	def create_decay_b(self, zamid,  dic):

		decay_b = dic.copy()

		for entries in decay_b:
			if entries not in ['half_life', 'total decay']:
				decay_b[entries] = decay_b[entries]/decay_b['total decay']

		self._decay_b[zamid] = decay_b
		self._decay_b[zamid]['half-life'] = self._decay_b[zamid]['half_life']


		if self._decay_b[zamid]['half-life'] == 'stable':
			self._decay_b[zamid]['unit'] = 'n/a'
		else:
			self._decay_b[zamid]['unit'] = 's'

		del self._decay_b[zamid]['half_life']

	@property
	def decay_b(self):

		return self._decay_b

class xs_lib(object):

	def __init__(self, name):

		self._name = name

		self._dict = {}
		self._xs = {}

	def add_data(self, zamid, **kwargs):

		if not kwargs:
			raise Empty_data("No data for {} has been entered".format(zamid))


		self._dict[zamid] = kwargs

		# If removal is not defined, calculate removal from xs values
		if 'removal' not in kwargs:
			removal = 0
			for i in kwargs:
				removal += kwargs[i]
			self._dict[zamid]['removal'] = removal


		#xs_dict = self._dict[zamid].copy()

		xs_dict = {}
		for key in self._dict[zamid]:
			correct_key = d.xs_lib_object_name_dict[key]
			xs_dict[correct_key] = self._dict[zamid][key]
			
		for key in xs_dict:
			xs_dict[key] = [xs_dict[key], 0.0] # Needs to reserve a spot for uncertainties
		
		#self.isomeric_branching(zamid, xs_dict)

		self._xs[zamid] = xs_dict

	def add_xs_dict(self, zamid, xs_dict):

		# If removal is not defined, calculate removal from xs values
		if 'removal' not in xs_dict:
			removal = 0
			for key in xs_dict:
				removal += xs_dict[key]
			xs_dict['removal'] = removal

		for key in xs_dict:
			xs_dict[key] = [xs_dict[key], 0.0] # Needs to reserve a spot for uncertainties

		#self.isomeric_branching(zamid, xs_dict)

		self._xs[zamid] = xs_dict


	# def isomeric_branching(self,zamid, xs_dict):

	# 	# This need to be developed later on
	# 	# For now, it only set the branching ration for Am241 to Am242/242m
	# 	if zamid == '952410' and '(n,gamma)' in xs_dict:
	# 		xs_val = xs_dict['(n,gamma)'][0].copy()
	# 		xs_dict['(n,gamma)'][0] = xs_val*0.89
	# 		xs_dict['(n,gamma)X'] = [xs_val*0.11, 0.0]

	# 	# Pm147 to Pm148/148m (very important for Sm149)
	# 	if zamid == '611470' and '(n,gamma)' in xs_dict:
	# 		xs_val = xs_dict['(n,gamma)'][0].copy()
	# 		xs_dict['(n,gamma)'][0] = xs_val*0.6953
	# 		xs_dict['(n,gamma)X'] = [xs_val*0.3046, 0.0]

	# This function screens the the xs_lib and look for which nuclide ngamma reactions needs to 
	# be branched
	def isomeric_branching_weighting(self, isomeric_branching_ratio):


		for zamid in self._xs:
			xs_dict = self._xs[zamid]
			name = zamid_to_name(zamid)
			mc_name = onix_name_to_openmc_name(name)
			if mc_name in isomeric_branching_ratio:
				ngamma_ratio = isomeric_branching_ratio[mc_name]['(n,gamma)']
				ngammaX_ratio = isomeric_branching_ratio[mc_name]['(n,gamma)X']
				xs_val = xs_dict['(n,gamma)'][0].copy()
				xs_dict['(n,gamma)'][0] = xs_val*ngamma_ratio
				xs_dict['(n,gamma)X'] = [xs_val*ngammaX_ratio, 0.0]




	@property
	def xs(self):

		return self._xs


class fy_lib(object):

	def __init__(self, id_number):

		self._id = id_number
		self._dic = {}
		self._fy = {}

	def add_data(self, zamid, dic):

		self._dic[zamid] = dic

		fy_dic = self._dic[zamid].copy()
		for key in fy_dic:
			fy_dic[key] = [fy_dic[key], 0.0] # Needs to reserve a spot for uncertainties

		self._fy[zamid] = fy_dic

	@property
	def fy(self):

		return self._fy





class Empty_data(Exception):
	"""Raise when the user does not enter any data while add_data has been called for a nuclide"""
	pass







