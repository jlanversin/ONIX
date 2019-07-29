import math as m
import os
import ast
import copy
from . import data
from openbu.cell import Cell

class Input(object):
	"""input reads, stores and process the input data in the input file provided by the user


	"""
	_file = None
	_time = None
	_cells = None
	_mode = None
	_lib = None
	_new_decay_lib_path = None
	_new_xs_lib_path = None
	_new_fy_lib_path = None
	_decay_lib_b = None
	_decay_lib_a = None
	_xs_lib = None
	_fy_lib = None
	_decay_nucl = None
	_xs_nucl = None
	_nucl_list = None
	_cell_dict = {}


	def __init__(self, input_path):

		self.input_path = input_path

		self._set_file()
		self._set_cell_id()
		self._set_mode()
		if self.mode in ['couple-openmc', 'couple-MCNP', 'couple-serpent']:
			self._set_MC_input_path()

		cell_list = []
		cell_dict = {}

		# Directly instantiate the cells
		for i in self._cell_id_list:


			cell = Cell(i) # instantiate a dummy cell with id of one of the real cell
			cell.vol = self._get_vol(i)
			cell.hm_vol = self._get_hm_vol(i)
		   # cell.set_initial_nuc(self.read_initial_nuc(i))

			cell._set_libs_from_input(self._get_lib(i))
			cell._set_nucl_list()
			cell.set_passlist(cell.get_lib_nucl())
			cell.set_initial_dens(self._get_dens(i))


			cell._set_sequence_from_input(self._get_sequence(i))

			#Should be optional
			cell.set_tree
			cell._set_leaves()
			cell._set_fission_tree()
			cell._set_fission_leaves()

			cell_list.append(copy.deepcopy(cell))
			cell_dict[i] = copy.deepcopy(cell)

		self._cell_list = cell_list
		self._cell_dict = cell_dict

	
	@property
	def file(self):

		return self._file

	def _set_file(self):

		path_to_file = self.input_path
		file = open(path_to_file, 'r')
		line = file.readlines()
		new_line = []
		for i in line:
			new_line.append(i.split('#')[0])
		if '' in new_line:
			new_line.remove('') # remove the remaining '' from line that start with #

		self._file = new_line



	@property
	def cells(self):
		"""Returns the absolute values of the decay constant of the nuclide"""
		if self._cells is None:
			pass  # define exception for undefined variable
		return self._cells

	@property
	def cell_id_list(self):
		"""Returns the absolute values of the decay constant of the nuclide"""
		if self._cells is None:
			pass  # define exception for undefined variable
		return self._cell_id_list

	def _set_cell_id(self):

		line = self._file
		cell_id_list = []
		for i in range(0,len(line)):
			l = line[i]
			l_split = l.split()
			if 'Cell' in l_split:
					cell_id_list.append(l_split[2])


		self._cell_id_list = cell_id_list
		self._cells = len(cell_id_list)


	@property
	def cell_list(self):

		return self._cell_list

	@property
	def cell_dict(self):

		return self._cell_dict


	@property
	def mode(self):
		"""Returns the absolute values of the decay constant of the nuclide"""
		if self._mode is None:
			pass  # define exception for undefined variable
		return self._mode

	def _set_mode(self):

		line = self._file
		search = 'var'
		for i in range(0,len(line)):
			l = line[i]
			if l == '=== Global Variables ===\n':
				search = 'mode'
			elif l not in ['\n', '\r\n'] and search == 'mode':
				if l.split()[0] == 'mode':
					mode = l.split()[2].replace('\n', '')

		self._mode = mode

	@property
	def MC_input_path(self):

		return self._MC_input_path

	def _set_MC_input_path(self):

		line = self._file
		search = 'var'
		for i in range(0,len(line)):
			l = line[i]
			if l == '=== Global Variables ===\n':
				search = 'MC_input_path'
			elif l not in ['\n', '\r\n'] and search == 'MC_input_path':
				if l.split()[0] == 'MC_input_path':
					MC_input_path = l.split()[2].replace('\n', '')
				elif l.split()[0] == '===':
					MC_input_path = os.getcwd()

		self._MC_input_path = MC_input_path


	@property
	def lib(self):
		"""Returns the absolute values of the decay constant of the nuclide"""
		if self._lib is None:
			pass  # define exception for undefined variable
		return self._lib

	def _get_lib(self, cell_id):

		line = self.file

		# Searching for the libraries path
		lib_path_dict = {}
		search = 'cell'
		for i in range(0, len(line)):
			l = line[i]
			if l == '=== Cell {} ===\n'.format(cell_id):
				search = 'lib'
			elif search == 'lib' and l == '--- Libraries ---\n':
				search = 'path'
			elif search == 'path':
				if l not in ['\n', '\r\n']:
					if l.split()[0] in ['===', '---']:
						break
					else:
						lib_path_dict[l.split()[0]] = l.split()[2]


		# Setting the libraries data
		lib_dict = {}

		for lib in lib_path_dict:
			if lib == 'decay':
				lib_dict[lib] = {}
				if lib_path_dict[lib] == 'default':
					lib_dict['decay_b'] = data.default_decay_lib_b
					lib_dict['decay_a'] = data.default_decay_lib_a
				else:
					lib_dict['decay_b'] = data.read_decay_lib(lib_path_dict[lib])
					lib_dict['decay_a'] = data.conv_decay_b_a(lib_dict['decay_b'])

			if lib == 'xs':
				if lib_path_dict[lib] == 'default':
					lib_dict[lib] = data.default_xs_lib
				else:
					lib_dict[lib] = data.read_xs_lib(lib_path_dict[lib])

			if lib == 'fy':
				if lib_path_dict[lib] == 'default':
					lib_dict[lib] = data.default_fy_lib
				else:
					lib_dict[lib] = data.read_fy_lib(lib_path_dict[lib])

		return lib_dict

	# Probably obsolete. Not use in the code anymore
	# def read_decay_lib_path(self):

	#     line = self._file
	#     exist = 0
	#     search = 'var'
	#     for i in range(0,len(line)):
	#         l = line[i]
	#         if l == '=== Global Variables ===\n':
	#             search = 'decay_library'
	#         elif l not in ['\n', '\r\n'] and search == 'decay_library':
	#             if l.split()[0] == 'decay_library':
	#                 decay_lib_path = l.split()[2]
	#                 self._new_decay_lib_path = decay_lib_path #only update _decay_library if it actually exists in the input file
	#                 exist = 1

	#     return exist

	# @property
	# def decay_lib_b(self):

	#     return self._decay_lib_b

	# @property
	# def decay_lib_a(self):

	#     return self._decay_lib_a


	# Probably obsolete. Not use in the code anymore
	# def read_xs_lib_path(self):

	#     line = self._file
	#     exist = 0
	#     search = 'var'
	#     for i in range(0,len(line)):
	#         l = line[i]
	#         if l == '=== Global Variables ===\n':
	#             search = 'xs_library'
	#         elif l not in ['\n', '\r\n'] and search == 'xs_library':
	#             if l.split()[0] == 'xs_library':
	#                 xs_lib_path = l.split()[2]
	#                 self._new_xs_lib_path = xs_lib_path #only update _xs_library if it actually exists in the input file
	#                 exist = 1

	#     return exist

	# @property
	# def xs_lib(self):

	#     return self._xs_lib


	# Probably obsolete. Not use in the code anymore
	# def read_fy_lib_path(self):

	#     line = self._file
	#     exist = 0
	#     search = 'var'
	#     for i in range(0,len(line)):
	#         l = line[i]
	#         if l == '=== Global Variables ===\n':
	#             search = 'fy_library'
	#         elif l not in ['\n', '\r\n'] and search == 'fy_library':
	#             if l.split()[0] == 'fy_library':
	#                 fy_lib_path = l.split()[2]
	#                 self._new_fy_lib_path = fy_lib_path #only update _xs_library if it actually exists in the input file
	#                 exist = 1

	#     return exist

	# @property
	# def fy_lib(self):

	#     return self._fy_lib


	def _get_vol(self, cell):

		line = self._file
		dens_dict = {}
		search = 'cell'
		for i in range(0, len(line)):
			l = line[i]
			if l == '=== Cell {} ===\n'.format(cell):
				search = 'vol'
			elif search == 'vol' and l == '--- Volume ---\n':
				search = 'val'
			elif search == 'val' and l not in ['\n', '\r\n']:
				if l.split()[0] in ['===', '---']:
					break
				vol = float(l.split()[0])

		return vol


	def _get_hm_vol(self, cell):

		line = self._file
		dens_dict = {}
		search = 'cell'
		for i in range(0, len(line)):
			l = line[i]
			if l == '=== Cell {} ===\n'.format(cell):
				search = 'hm_vol'
			elif search == 'hm_vol' and l == '--- Heavy Metal Volume ---\n':
				search = 'val'
			elif search == 'val' and l not in ['\n', '\r\n']:
				if l.split()[0] in ['===', '---']:
					break
				hm_vol = float(l.split()[0])

		return hm_vol


	def _get_dens(self, cell):

		line = self._file
		dens_dict = {}
		search = 'cell'
		for i in range(0, len(line)):
			l = line[i]
			if l == '=== Cell {} ===\n'.format(cell):
				search = 'dens'
			elif search == 'dens' and l == '--- Nuclides Densities ---\n':
				search = 'val'
			elif search == 'val' and l not in ['\n', '\r\n']:
				if l.split()[0] in ['===', '---']:
					break
				zamid = l.split()[0]
				dens = float(l.split()[2])
				dens_dict[zamid] = dens

		return dens_dict

	# Probably obsolete. Not use in the code anymore
	# def read_initial_nuc(self, cell):

	#     line = self._file
	#     initial_nuc_list = []
	#     search = 'cell'
	#     for i in range(0, len(line)):
	#         l = line[i]
	#         if l == '=== Cell {} ===\n'.format(cell):
	#             search = 'dens'
	#         elif search == 'dens' and l == '--- Nuclides Densities ---\n':
	#             search = 'val'
	#         elif search == 'val' and l not in ['\n', '\r\n']:
	#             if l.split()[0] in ['===', '---']:
	#                 break
	#             initial_nuc_list.append(l.split()[0])

	#     return initial_nuc_list


# initial density

	def _get_sequence(self, cell):

		line = self._file
		seq_dict = {}
		search = 'cell'
		for i in range(0, len(line)):
			l = line[i]
			if l == '=== Cell {} ===\n'.format(cell):
				search = 'sequence'
			elif search == 'sequence' and l == '--- Sequence ---\n':
				search = 'val'
			elif search == 'val' and l not in ['\n', '\r\n']:
				if l.split()[0] in ['===', '---']:
					break
				if l.split()[0] in ['unit_vector', 'norma_vector', 'substeps_vector']:
					seq_dict[l.split()[0]] = ast.literal_eval(l.split('=')[1].replace(' ',''))
				elif l.split()[0] in ['unit', 'normalization', 'flux_approximation']:
					seq_dict[l.split()[0]] = l.split()[2]

		# record the number of burnup steps
		seq_dict['bu_steps'] = len(seq_dict['unit_vector'])

		if seq_dict['unit'] in ['MWd/kg', 'kWd/kg', 'GWd/kg']:
			seq_dict['unit_type'] = 'bu'
		elif seq_dict['unit'] in ['s', 'd','m']:
			seq_dict['unit_type'] = 'time'

		# conversion of unit_vector in seconds if input unit is in another unit
		if seq_dict['unit'] != 's':
			seq_dict['unit_vector'] = [x*data.list_and_dict.time_dic[seq_dict['unit']] for x in seq_dict['unit_vector']]

		return seq_dict


