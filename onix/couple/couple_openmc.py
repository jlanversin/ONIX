import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
import copy
import xml.etree.ElementTree as ET
import glob
import pdb
import time
from copy import deepcopy

import openmc
import openmc.mgxs as mgxs
from onix.cell import Cell
from onix.system import System
from onix import salameche
from .openmc_fix import *

from onix import utils
from onix import data


class Couple_openmc(object):
    """This class is used to execute coupled-mode simulations

    Through this class, the user can:

    - chose to parallelize the Monte Carlo simulations via the "set_MPI" method

    - set the nuclear data libraries to be used in the simulation

    - set the burnup/time sequence for the simulation

    - manually set the volume of each BUCells (note that by default, ONIX uses the OpenMC stochastic volume calculation
      to find the volume of each cell. However, this method might produce important errors on the volume and it is
      advised to set the volume manually)
    - select which of the Cells defined in the OpenMC input should be depleted (these cells will be known as BUCells)
    - select a list of nuclides in each BUCell which cross sections will be tallied. If the user does not specify any such list of nuclides, ONIX will by default tally the cross sections of all
      nuclides which data are found in the cross section library.

    The Couple_openmc class also takes care of the coupling between ONIX and OpenMC. At the beginning of a simulation,
    it imports certain initial parameters from the OpenMC input to ONIX such as the name of the BUCells and the initial
    nuclide densities of each BUCell. During the simulation, it makes sure that information such as the flux tallies,
    neutron spectrum tallies, the reaction rates and the updated nuclide densities are correctly transfered between the burnup
    module and the neutron transport module.

    The Couple_openmc class will also sample the isomeric branching ratios and (n,gamma) cross sections on the
    the same energy points to prepare for the calculations of one-group isomeric branching ratios.

    **IMPORTANT: the root cell in the OpenMC input should be named "root cell"**

    Parameters
    ----------
    MC_input_path : str
        Specifies where the OpenMC input files in .xml format are located. This parameter should not be specified as ONIX cannot yet work without Python API.

    xs_mode : str
        Choice between 'constant lib' and 'no constant lib'. Default to 'no constant lib'.
        'constant lib' indicates to ONIX that the user has provided a constant one-group cross section libraries to be used.
        'no constant lib' indicates to ONIX that it should use one-group cross section computed by OpenMC.

    MPI : str
        Choice between 'on' and 'off'.
        Indicates whether OpenMC will be parrallelized or not.
    
    """



    # One-group energy bin
    energy_bin = openmc.EnergyFilter([0., 20.0e6])

    # Multigroup energy bin
    minorder = -3
    maxorder = 7
    mg_energy = np.logspace(minorder, maxorder, (maxorder - minorder) * 30 + 1)
    mg_energy_mid_points = [(x+y)/2 for x,y in zip(mg_energy[1:],mg_energy[:-1])]
    #mg_energy = np.logspace(-3, 7, num=300, base=10.0)
    mg_energy_bin = openmc.EnergyFilter(mg_energy)

    zero_dens_1_atm = 1E-24

    def __init__(self, MC_input_path = None, xs_mode = 'no constant lib', MPI = None):

        # If no MC_input_path is input, it is set to cwd
        if MC_input_path == None:
            MC_input_path = os.getcwd()
        self._MC_input_path = MC_input_path

        # If no mode is input by the user, mode is set by default to 'const_lib'
        if xs_mode == None:
            xs_mode == 'constant lib'
        self._xs_mode = xs_mode

        # # If MPI is not set to on by the user, MPI is set to off
        # if MPI == None:
        #   MPI == 'off'
        # self._MPI = MPI
        self._MPI = None

        self._volume_set = 'no'

        # List of selected cells to deplete
        self.selected_bucells_name_list = None
        # Dict of selected cells to deplete with their user defined nucl list
        self.selected_bucells_nucl_list_dict = None

        self._default_fy_lib_set = 'no'
        self._user_fy_lib_set = 'no'
        self._complete_user_fy_lib = 'no'
        self._complete_fy_lib = 'no'
        self._decay_lib_set = 'no'
        self._xs_lib_set = 'no'

        self._sampled_isomeric_branching_data = None
        self._sampled_ng_cross_section_data = None

        # This is the path set in OpenMC for the hdf5 point-wise cross sections
        self._cross_sections_path = None

        self._openmc_bin_path = None

        # By default all available isomeric branching ratio are calculated
        self._few_isomeric = 'off'

        # By default, ONIX does not compute reactions rates ranking (it takes a lot of memory)
        self._reac_rank = 'off'

        # By default tolerance value is set to 500 K for treating intermediate temperatures.
        self._tolerance = 500.0

        # Old way of defaulting MC_input_path to cwd
        # if args:
        #   self._MC_input_path = arg[0]
        # else:
        #   self._MC_input_path = os.getcwd()

# This method is used within the code to create a new material that is 
#  uniquely associated to a particular cell

    @property
    def MC_input_path(self):
        """Returns the path to OpenMC input files."""

        return self._MC_input_path

    @property
    def xs_mode(self):
        """Returns the mode for cross sections library."""

        return self._xs_mode
    
    # @property
    # def nucl_list_dict(self):

    #   return self._nucl_list_dict

    @property
    def root_cell(self):
        """Returns the root cell as define in the OpenMC geometry."""
        return self._root_cell

    @property
    def MPI(self):
        """Returns wether MPI is activated or not ('on' if the user chose to parallelize OpenMC with MPI, 'off' if not)"""
        return self._MPI

    def set_MPI(self, execu, tasks):
        """Sets the settings for the MPI parallelization

        Parameters
        ----------
        execu : str
            The MPI execute command, for example 'srun' on cluster using slurm

        tasks : str
            The number of tasks to parallelize"""
            
        self._MPI = 'on'
        self._tasks = tasks
        self._exec = execu

    def few_isomeric_on(self):
        """Calling this function will force ONIX to only calculate isomeric branching ratios of Pm147 and Am241
        By default ONIX calculates all isomeric branching ratios for which data
        can be found in the EAF-2010 multiplicities library
        """
        self._few_isomeric = 'on'

    @property
    def reac_rank(self):

        return self._reac_rank
    
    def reac_rank_on(self):
        """Calling this function will tell ONIX to produce reaction rates ranking and print them for each BUCells
        By default ONIX does not produce reaction rates ranking as it takes a lot of memory.

        Important: In the current version of ONIX, this function needs to be called before couple.import_openmc() is called other
        wise the function will not work.
        """
        self._reac_rank = 'on'

    def select_bucells(self, bucell_list):
        """Selects the cells from the OpenMC input that should be depleted.

        The user can choose to specify the nuclides for which cross sections should be updated by OpenMC simulations. To do so, the user must enter a tuple where the first 
        element is the OpenMC cell object and the second element is a list of the name of nuclides 
        which cross sections should be updated.

        If the user does not specify for which nuclides cross sections should be updated,
        only the OpenMC cell objects should be entered. ONIX will by default calculate the cross sections
        for all nuclides available in the cross section library.

        Parameters
        ----------
        bucell_list: list of openmc.Cell or list of tuples
            Defines the list of cells that should be depleted by ONIX.

            For example, the user can set *bucell_list = [(fuel_cell, ['U-235', 'Pu-239']), (cladding_cell, ['Al-27'])]*
        """
        self.selected_bucells_name_list = []
        self.selected_bucells_nucl_list_dict = {}

        for arg in bucell_list:
            # If this element is a tuple (bucell , nucl list)
            if isinstance(arg,tuple):   
                bucell_name = arg[0].name
                self.selected_bucells_name_list.append(bucell_name)
                self.selected_bucells_nucl_list_dict[bucell_name] = arg[1]
            else:
                self.selected_bucells_name_list.append(arg.name)

    def _get_nucl_to_be_tallied(self, bucell):

        # If the user has provided a list of nuclide to be tallied
        if bucell.name in self.selected_bucells_nucl_list_dict:
            nucl_list_input = self.selected_bucells_nucl_list_dict[bucell.name]
            if nucl_list_input == 'initial nuclides':
                nucl_list = bucell.init_nucl
            elif nucl_list_input == 'NAX':
                NAX_nucl_list_name = utils.zamid_list_to_name_list(data.NAX_nucl_list)
                NAX_nucl_list_name_new_format = utils.bu_namelist_to_mc_namelist(NAX_nucl_list_name)
                # I add init_nucl because I believe it is important for init nuclides to be tallied
                # Their density is usually high enough that their change can influence spectrum etc...
                nucl_list = [x for x in NAX_nucl_list_name_new_format if x in self.MC_XS_nucl_list] + bucell.init_nucl
                # Here we remove the potential duplicates
                nucl_list = list(dict.fromkeys(nucl_list))
            else:
                nucl_list = nucl_list_input
        else:
            nucl_list = self.MC_XS_nucl_list

        return nucl_list

    @property
    def system(self):
        """Returns the system object of the simulation.
        """
        return self._system

    @system.setter
    def system(self, system):
        """Sets the system object of the simulation
        """
        self._system = system
    

    def import_openmc(self, root_cell):
        """This method import certain parameters from OpenMC to ONIX

        This method launches the _pre_run script which runs a first Monte Carlo
        simulation in order for ONIX to get a handle on the various OpenMC objects
        (this feature is essential when the user decides to defines input parameters
        via a text input file rather than using the Python API. However, as of now,
        ONIX only supports Python API input)

        The _pre_run script also requests OpenMC to calculate the volume of each OpenMC
        Cell with openmc.calculate_volume(). However, the volumes obtained
        via this method can have significant errors, especially for complex or large geometry.
        As of now, it is recommented that the user manually sets the volume of ech BUCell via
        the onix.couple.Couple_openmc.set_volume() command.

        Finally, a very important functionality of import_openmc is to add minute quantities
        of certain nuclides in BUCells' materials so that OpenMC can calculate their cross sections.
        This means that ONIX is going to modify the material.xml and add all nuclides requested by the user
        for cross section calculations. This feature is necessary because OpenMC cannot tally reaction
        rates of nuclides that are not present in the material.

        Parameters
        ----------
        root_cell: openmc.Cell
            Root cell object defined in OpenMC geometry.
        """


        # Printing of ONIX header here since this prerun will start OpenMC simulations
        print (utils.printer.onix_header)

        #Instantiate a system
        system = System(1)
        if self.reac_rank == 'on':
            system.reac_rank_on()
        self.system = system

        # read periodic surfaces  (openmc summary forgets the periodic surfaces coupling)
        # This function stores the periodic coupling of surfaces and it will be used later
        self._periodic_surfaces_dict = read_periodic_surfaces()

        # prerun to access cells and materials objects, to set cell volumes and if chosen
        # add 0 density nuclides
        self._pre_run(root_cell)

        bucell_dict = self._get_bucell_from_cell()
        system.bucell_dict = bucell_dict
        system.bounding_box = self.bounding_box

        # reads, modifies and set settings
        self._read_user_settings()


        # Move input files to input file folder
        self._gen_user_input_folder()
        self._copy_user_input()

        # MOVED TO OPENMC RUN
        # # New material xml file needs to be written with zero dens nuclides
        # self.export_material_to_xml()
        # # New geometry xml file needs to be written with new materials id
        # self.export_geometry_to_xml()
        # # Tallies xml files needs to be writen
        # self.export_tallies_to_xml()
        # # New settings xml files needs to be written
        # self.export_settings_to_xml()

    # Proably for when onix pass dens to OpenMC

    @property
    def bounding_box(self):
        """Returns the bounding box of the geometry.
        """
        return self._bounding_box
     
    def set_bounding_box(self, ll, ur):
        """Sets the bounding box of the geometry.

        Parameters
        ----------
        ll: [x,y,z]
            Extremity of the box located in the negative part of a 3-D plan
        ur: [x,y,z]
            Extremity of the box located in the positive part of a 3-D plan

        """
        self._bounding_box = [ll, ur]

    @property
    def tolerance(self):
        return self._tolerance

    def set_tolerance(self, tolerance):
        """For treating intermediate temperatures,'tolerance' indicates
        a range of temperature within which cross sections may be used.
        """
        self._tolerance = tolerance

    # def _set_material_nuclides(self, cell):

    #   cell_id = cell.id
    #   passlist = cell.passlists
    #   total_dens = cell.total_dens

    #   material = openmc.Material(cell_id)

    #   material.set_density('atom/b-cm', total_dens)

    #   for nuc in passlist:

    #       nuc_name = nuc.name.replace('-', '')
    #       nuc_ao = nuc.dens/total_dens
    #       material.add_nuclide(nuc_name,  nuc_ao)

    # Osbolete
    # def set_settings(self, settings, init_dist):
    # # OpenMC simulation parameters
    #   batches = settings['batches']
    #   inactive = settings['inactive'] 
    #   particles = settings['particles']

    #   # Instantiate a Settings object
    #   settings_file = openmc.Settings()
    #   settings_file.batches = batches
    #   settings_file.inactive = inactive
    #   settings_file.particles = particles
    #   settings_file.output = {'tallies': True}

    #   # Create an initial uniform spatial source distribution over fissionable zones
    #   init_dist = setting.init_dist
    #   shape = init_dist['shape']
    #   low_left_bound = init_dist['low_left']
    #   up_right_bound = init_dist['up_right']
    #   if shape == 'Box':
    #       uniform_dist = openmc.stats.Box(low_left_bound, up_right_bound, only_fissionable=True)
    #   settings_file.source = openmc.source.Source(space=uniform_dist)

    #   # Export to "settings.xml"
    #   settings_file.export_to_xml()

    def _gen_user_input_folder(self):

        utils.gen_folder('user_input')

    def _copy_user_input(self):

        MC_input_path = self.MC_input_path
        user_input_folder_path = os.getcwd() + '/user_input'
        shutil.copyfile(MC_input_path + '/geometry.xml', user_input_folder_path + '/geometry.xml')
        shutil.copyfile(MC_input_path + '/materials.xml', user_input_folder_path + '/materials.xml')
        shutil.copyfile(MC_input_path + '/settings.xml', user_input_folder_path + '/settings.xml')

    def _pre_run(self, root_cell):

        print ('=== OpenMC pre-run ===\n\n\n')

        MC_input_path = self.MC_input_path
        pre_run_path = os.getcwd() +'/pre_run'
        try:
            shutil.rmtree(pre_run_path)
        except OSError:
            pass
        os.mkdir(pre_run_path)

        # Prepare the volume calculation
        #bounding_box = root_cell.bounding_box
        ll = self.bounding_box[0]
        ur = self.bounding_box[1]
        cell_dict = root_cell.get_all_cells()
        cell_list = utils.cell_dict_to_cell_list(cell_dict)
        cell_list.append(root_cell) # Add root_cell so that the total volume is calculated
        vol1 = openmc.VolumeCalculation(cell_list, 100000, lower_left = ll, upper_right = ur)

        settings = openmc.Settings()
        settings.volume_calculations = [vol1]
        settings.temperature = {'method':'nearest', 'tolerance': self.tolerance}
        settings.run_mode='volume'
        settings.export_to_xml(path = pre_run_path + '/settings.xml')

        # Copy the geometry and material file to the new dummy dir
        shutil.copyfile(MC_input_path + '/geometry.xml', pre_run_path + '/geometry.xml')
        shutil.copyfile(MC_input_path + '/materials.xml', pre_run_path + '/materials.xml')

        # By default, the openm_exec is set to 'openmc'
        # For some reasons, this does not work on the cluster (della)
        # On della, we need to explicitly define the absolute path to the bin we want to use

        # If the user has not define the bin path, we assume the simulation is not on a cluster
        if self.openmc_bin_path == None:
            openmc.calculate_volumes(cwd = pre_run_path)
        else:
            openmc.calculate_volumes(cwd = pre_run_path, openmc_exec=self.openmc_bin_path)

        # Read and set initial nuclides dict
        self._set_init_nucl_dict(root_cell)

        # # Read each material object and add 1atm nuclides chosen by the user
        # if self.mode == 'no_const_lib':
        #   self.add_zero_dens_nuclides(self.nucl_list_dict)

        self._set_initial_summary(pre_run_path)
        self._set_cross_sections_path(pre_run_path)
        # Read cross sections xml files, create MC_XS_nucl_list
        self._set_MC_XS_nucl_list()
        self._set_root_universe()
        root_cell_name = 'root cell' #  need to be specified by the user at some point
        self._set_root_cell(root_cell_name)

        # Extract cells from summary, add 1 atm nuclides to their material
        self._change_cell_materials()

        # Read and distribute volumes to cells
        self._set_vol_to_cell(vol1, pre_run_path)

        # pdb.set_trace()

        shutil.rmtree(pre_run_path)

    def _set_vol_to_cell(self, vol1, pre_run_path):

        root_cell = self.root_cell
        cell_dict = root_cell.get_all_cells()
        cell_list = utils.cell_dict_to_cell_list(cell_dict)
        vol2 = vol1.from_hdf5(pre_run_path + '/volume_1.h5')
        for cell in cell_list:
            cell.add_volume_information(vol2)

        # Add volume for the root cell too and set it to system
        root_cell.add_volume_information(vol2)

        # This is now done in pass_vol when user does not define vol manually
        # system = self.system
        # system.total_vol = root_cell.volume



    # About summary
    # There are two type of OpenMC summary stored in onix: initial_summary and updated_summary
    # I have noticed that updating the initial summary was causing some problem in onix, i.e., the material xml file
    # would not update the density. The material was updating densities normaly when relying on the initial summary
    # A temporary solution for that is to separate the initial summary (that will be used to set dens to cells) and the
    # updated summary that will be used to extract densities to get 1g xs for bucells

    @property
    def initial_summary(self):
        """Returns the original summary.h5 file created by OpenMC before ONIX applies change to OpenMC's objects.
        """

        return self._initial_summary
    
    def _set_initial_summary(self, path = os.getcwd()):

        initial_summary = openmc.Summary(path + '/summary.h5')

        ######### OpenMC Summary src does not close the hdf5 file it opens
        ######### When onix tries to shutil.rmtree the pre_run folder, it can't because
        ######### a stream to summary.h5 is still open
        ######### We therefore close it here
        ######### !!!! This should be modified in OpenMC at some points ###########
        initial_summary._f.close()
        ######### !!!! This should be modified in OpenMC at some points ###########
        
        self._initial_summary = initial_summary

    @property
    def updated_summary(self):
        """Returns the updated summary.h5 file created by OpenMC before ONIX applies change to OpenMC's objects.
        """

        return self._updated_summary
    
    def _set_updated_summary(self, path = os.getcwd()):

        updated_summary = openmc.Summary(path + '/summary.h5')

        ######### OpenMC Summary src does not close the hdf5 file it opens
        ######### When onix tries to shutil.rmtree the pre_run folder, it can't because
        ######### a stream to summary.h5 is still open
        ######### We therefore close it here
        ######### !!!! This should be modified in OpenMC at some points ###########
        updated_summary._f.close()
        ######### !!!! This should be modified in OpenMC at some points ###########
        
        self._updated_summary = updated_summary

    @property
    def statepoint(self):

        return self._statepoint

    def _set_statepoint(self, path = os.getcwd()):

        file_list = os.listdir()
        for file in file_list:
            if 'statepoint' in file:
                st_name = file
        statepoint = openmc.StatePoint(path + '/{}'.format(st_name))
        self._statepoint = statepoint

    def _set_kinf(self):

        statepoint = self.statepoint
        kinf = statepoint.k_combined

        system = self.system
        sequence = system.sequence
        sequence._set_macrostep_kinf(kinf)
    

    @property
    def root_cell(self):
        """Returns the root cell object created in the OpenMC geometry."""

        return self._root_cell
    
    def _set_root_universe(self):

        summary = self.initial_summary
        geometry = summary.geometry
        root_universe = geometry.root_universe
        self._root_universe = root_universe

    def _set_root_cell(self, root_cell_name):

        summary = self.initial_summary
        # get_cells_by_name returns a list hence the index 0
        self._root_cell = summary.geometry.get_cells_by_name(root_cell_name)[0]

        add_periodic_surfaces(self._root_cell, self._periodic_surfaces_dict)

        region = self._root_cell.region
        #print (region.get_surfaces())
        for surface_id in region.get_surfaces():
            surface = region.get_surfaces()[surface_id]

    # While the most convenient way would be to extract path from summary.materials. It looks like
    # summary does not store the cross sections path in materials
    def _set_cross_sections_path(self, pre_run_path):

        path_to_materials_xml = pre_run_path + '/materials.xml'
        tree = ET.parse(path_to_materials_xml)
        root = tree.getroot()
        for child in root:
            if child.tag == 'cross_sections':
                self._cross_sections_path = child.text.replace('/cross_sections.xml', '')


    @property
    def materials(self):
        """Returns the materials object created in the OpenMC material input."""

        return self._materials
    
    def _change_cell_materials(self):

        summary = self.initial_summary
        materials = summary.materials
        for bucell_name in self.selected_bucells_name_list:
            # If bucell should only tally initial nuclide, no need to add 1 atm nuclides
            # This caused the cell material branching issue discovered on Nov 2 2019
            # This is now commented out
            # if self.selected_bucells_nucl_list_dict != {}:
            #   if bucell_name in self.selected_bucells_nucl_list_dict:
            #       if self.selected_bucells_nucl_list_dict[bucell_name] == 'initial nuclides':
            #           continue
            cell = summary.geometry.get_cells_by_name(bucell_name)[0]
            self._add_zero_dens_nuclides(cell)


    # def _change_cells_materials(self):

    #   root_cell = self.root_cell
    #   print (root_cell)
    #   cell_dict = root_cell.get_all_cells()
    #   # materials with 1atm nuclides
    #   materials = self.materials
    #   for cell_id in cell_dict:
    #       cell = cell_dict[cell_id]
    #       cell.get_all_materials
    #       material_dict = cell.get_all_materials()
    #       material = material_dict[list(material_dict.keys())[0]]
    #       mat_name = material.name
    #       for new_mat in materials:
    #           if new_mat.name == mat_name:
    #               print(new_mat.get_nuclides(), material.get_nuclides())



    # Add zero dens nuclide for each material
    def _add_zero_dens_nuclides(self, cell):

        material_dict = cell.get_all_materials()

        # A deepcopy of the object material is created and this copy will eventually overwrite the
        # original material that was filling the cell
        # If a deepcopy was not used:
        # -  when a material id is modified (branched out) vis-a-vis the hosting cell
        # the changes would also be applied on the original material object
        # However, this original material object might fill other cells which would mean that
        # when the loop goes to another cell with same material, the change are going to
        # cumulate
        # -  in addition, different cells with same initial material would still only have one common material object
        # (i.e. the branching out fails)
        material = deepcopy(material_dict[list(material_dict.keys())[0]])   

        init_nucl = material.get_nuclides()
        cell_name = cell.name
        mat_name = material.name
        # Not sure if this is necessary
        if self.selected_bucells_nucl_list_dict != {}:
            if cell_name in self.selected_bucells_nucl_list_dict:
                nucl_list_input = self.selected_bucells_nucl_list_dict[cell_name]
                if nucl_list_input == 'initial nuclides':
                    nucl_list = init_nucl
                elif nucl_list_input == 'NAX':
                    NAX_nucl_list_name = utils.zamid_list_to_name_list(data.NAX_nucl_list)
                    NAX_nucl_list_name_new_format = utils.bu_namelist_to_mc_namelist(NAX_nucl_list_name)
                    # I add init_nucl because I believe it is important for init nuclides to be tallied
                    # Their density is usually high enough that their change can influence spectrum etc...
                    nucl_list = [x for x in NAX_nucl_list_name_new_format if x in self.MC_XS_nucl_list] + init_nucl
                    # Here we remove the potential duplicates
                    nucl_list = list(dict.fromkeys(nucl_list))
                else:
                     nucl_list = nucl_list_input
            else:
                nucl_list = self.MC_XS_nucl_list
                #nucl_list = utils.bu_namelist_to_mc_namelist(nucl_list)

        else:
            nucl_list = self.MC_XS_nucl_list

        if not utils.is_lista_in_listb(init_nucl, nucl_list):
            raise Initial_nuclides_not_in_nuclide_list('Some initial nuclides in cell {} material {} are not included in nucl_list'.format(cell_name, mat_name))
        for nucl in nucl_list:
            if nucl not in init_nucl:
                material.add_nuclide(nucl, self.zero_dens_1_atm)

        # Material is rename 'cell name' + 'mat'
        # New material id is 'mat id' + 'cell id'
        material.name = '{} mat'.format(cell_name)
        material.id = int('{}{}'.format(material.id, cell.id))

        # Here, we overwrite the original material object with a copy which id and name have been modified
        # to reflect the fact that this material object is specific to this cell
        cell.fill = material    

    @property
    def sequence(self):
        """Returns the sequence assigned to the Couple_openmc() object."""
        return self._sequence
    
    # @sequence.setter
    # def sequence(self, sequence):

    #   self._sequence = sequence

    def set_sequence(self, sequence):
        """Sets a sequence object to the Couple_openmc object.

        Parameters
        ----------
        sequence: onix.Sequence
            ONIX sequence object defined by the user.


        """

        self._sequence = sequence
        system = self.system
        system.set_sequence(sequence, mode = 'couple')

    # Add zero dens nuclide for each cell
    # Sort of complicated
    # Will deal with that later
    # def add_zero_dens_nuclides(self, root_cell, nucl_list_dict):

    #   cell_dict = root_cell.get_all_cells()

    #   for cell_id in cell_dict:
    #       cell = cell_dict[cell_id]
    #       material_dict = cell.get_all_materials()
    #       material = copy.deepcopy(material_dict[material_dict.key()[0]])
    #       init_nucl = material.get_nuclides()
    #       nucl_list = nucl_list_dict[cell_id]
    #       nucl_list = utils.bu_namelist_to_mc_namelist(nucl_list)
    #       if not is_lista_in_listb(init_nucl, nucl_list):
    #           raise Initial_nuclides_not_in_nuclide_list('Some initial nuclides of material {} are not included in nucl_list'.format(mat_id))
    #       for nucl in nucl_list:
    #           if nucl not in init_nucl:
    #               material.add_nuclide(nucl, 1E-22)

    @property
    def init_nucl_dict(self):

        return self._init_nucl_dict

    # Create a dict with initial nuclides of each cell before
    # cell material is added 1 atm nuclides
    def _set_init_nucl_dict(self, root_cell):

        cell_dict = root_cell.get_all_cells()

        init_nucl_dict = {}
        for cell_id in cell_dict:
            cell = cell_dict[cell_id]
            material_dict = cell.get_all_materials()
            material = material_dict[list(material_dict.keys())[0]]
            init_nucl = material.get_nuclides()
            init_nucl_dict[cell_id] = init_nucl

        self._init_nucl_dict = init_nucl_dict

        # # If no nucl list has been defined, nucl_list_dict = init_nucl
        # if self.nucl_list_dict == None:
        #   self.nucl_list_dict = self.init_nucl_dict

    # Use init_nucl_dict and distribute init_nucl to each bucell
    def _set_init_nucl(self, cell_dict, bucell_dict):

        init_nucl_dict = self.init_nucl_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.init_nucl = init_nucl_dict[bucell_id]

    @property
    def MC_XS_nucl_list(self):
        """Returns the list of nuclides for which cross section data exist in the cross section library.
        """

        return self._MC_XS_nucl_list

    @MC_XS_nucl_list.setter
    def MC_XS_nucl_list(self, MC_XS_nucl_list):

        self._MC_XS_nucl_list = MC_XS_nucl_list

    def _set_MC_XS_nucl_list(self):

        #path_to_xs_xml = os.environ['OPENMC_CROSS_SECTIONS']
        path_to_xs_xml = self._cross_sections_path + '/cross_sections.xml'

        self.MC_XS_nucl_list = []

        tree = ET.parse(path_to_xs_xml)
        root = tree.getroot()

        for child in root:
            if child.attrib['type'] == 'neutron':
                self._MC_XS_nucl_list.append(child.attrib['materials'])

        # Remove trouble makers that appears in JEFF32 cross section    
        # For some reason, OpenMC can't find these nuclides in jeff lib at 800K
        # self._MC_XS_nucl_list.remove('Cu63')
        # self._MC_XS_nucl_list.remove('Cu65')
        # self._MC_XS_nucl_list.remove('Mn55')
        # # THose are not handled by onix
        # self._MC_XS_nucl_list.remove('C0')
        # self._MC_XS_nucl_list.remove('V0')
        # self._MC_XS_nucl_list.remove('Zn0')

        # Remove trouble makers that appears in ENDFVIII cross section  
        # For some reason, OpenMC can't find these nuclides in jeff lib at 800K
        # self._MC_XS_nucl_list.remove('Cu63')
        # self._MC_XS_nucl_list.remove('Cu65')
        # self._MC_XS_nucl_list.remove('Mn55')
        # # THose are not handled by onix
        try:
            self._MC_XS_nucl_list.remove('C0')
        except ValueError:
            pass
        try:
            self._MC_XS_nucl_list.remove('V0')
        except ValueError:
            pass
        try:
            self._MC_XS_nucl_list.remove('Zn0')
        except ValueError:
            pass



    # When volume is passed from OpenMC to ONIX
    def _pass_vol(self, cell_dict, bucell_dict):

        # Need to loop over bucell_dict because there might be more cells than bucells
        for i in bucell_dict:

            cell = cell_dict[i]
            cell_volume = cell.volume
            bucell = bucell_dict[i]
            bucell.vol = cell_volume

        # root_cell is not in bucell_dict but it contains the info on the total volume
        # Here, the total volume is set to system
        system = self.system
        system.total_vol = self.root_cell.volume

    # When volume is set directly by user
    # Right now this should bet set after import openmc as it overwrites volume calculated by openmc
    def set_vol(self, vol_dict):
        """Sets the volume of each BUCell in :math:`cm^{3}`.

        Parameters
        ----------
        vol_dict: dict
            Dictionnary where keys are the names of BUCell (as defined by the user) and entries are volume in :math:`cm^{3}`.

        """

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


    def _pass_nuclide_densities(self, cell_dict, bucell_dict):

        for i in bucell_dict:
            bucell = bucell_dict[i]
            cell = cell_dict[i]
            #init_nucl = self.init_nucl_dict[i]
            init_nucl = bucell.init_nucl
            materials = cell.get_all_materials()
            openmc_dens_dict = materials[list(materials.keys())[0]].get_nuclide_atom_densities()
            onix_dens_dict = {}

            for nucl in openmc_dens_dict:
                #OpenMC xs has cross section for element carbon and Vanadinium (C0, V0) only. onix can't handle that
                if nucl == 'C0' or nucl == 'V0':
                    continue
                onix_nucl = utils.openmc_name_to_onix_name(nucl)
                # if nucl is one of the initial non-zero initial nuclide, pass the density
                if nucl in init_nucl:
                    onix_dens_dict[onix_nucl] = openmc_dens_dict[nucl][1]
                # if nucl is not one of the non-zero initial nuclide, set densiy to zero
                else:
                    onix_dens_dict[onix_nucl] = 0.0

            bucell = bucell_dict[i]
            bucell.set_initial_dens(onix_dens_dict)

    def _get_bucell_from_cell(self):

        root_cell = self.root_cell
        bucell_dict = {}
        cell_dict = root_cell.get_all_cells()

        for i in cell_dict:
            cell = cell_dict[i]
            cell_name = cell.name
            if cell_name in self.selected_bucells_name_list:            
                bucell_dict[i] = Cell(i, cell_name)

        self._set_init_nucl(cell_dict, bucell_dict)
        self._pass_vol(cell_dict, bucell_dict)
        self._pass_nuclide_densities(cell_dict, bucell_dict)

        return bucell_dict


    def _get_flux_tally(self, bucell):

        flux = openmc.Tally(name='{} flux'.format(bucell.name))
        flux.filters = [openmc.CellFilter(bucell.id)]
        flux.filters.append(self.energy_bin)
        flux.scores = ['flux']

        return flux

    def _get_flux_spectrum_tally(self, bucell):

        flux_spectrum = openmc.Tally(name='{} flux spectrum'.format(bucell.name))
        flux_spectrum.filters = [openmc.CellFilter(bucell.id)]
        flux_spectrum.filters.append(self.mg_energy_bin)
        flux_spectrum.scores = ['flux']

        return flux_spectrum

    # Every nuclide presents in cell material will have its tally taken
    # (n,t) is tallied for all nuclides for the moment. It is surely more efficient to only tally it for Lithium
    def _get_all_nucl_rxn_tally(self, bucell):

        nucl_list = self._get_nucl_to_be_tallied(bucell)
        print ('bucell name',bucell.name)
        print ('nucl list when set to tally',nucl_list)
        nucl_list = utils.bu_namelist_to_mc_namelist(nucl_list)
        rxn = openmc.Tally(name='{} rxn rate'.format(bucell.name))
        rxn.filters = [openmc.CellFilter(bucell.id)]
        rxn.filters.append(self.energy_bin)
        rxn.scores = ['fission', '(n,gamma)', '(n,2n)', '(n,3n)', '(n,p)', '(n,a)', '(n,t)']
        #rxn.scores = ['fission', '(n,gamma)']
        #rxn.scores = ['fission', '(n,gamma)', '(n,2n)']
        rxn.nuclides = nucl_list
        
        return rxn

    # Tally dedicated to production reaction rate for tritium
    # Only Lithium6 and Berylium10 are non-negligible tritium producer nuclides through (n,t)
    # def get_H3_rxn_tally(self, bucell):

    #   print ('Creating H3 rxn rate tally')
    #   H3_rxn = openmc.Tally(name='{} H3 rxn rate'.format(bucell.name))
    #   H3_rxn.filters = [openmc.CellFilter(bucell.id)]
    #   H3_rxn.filters.append(self.energy_bin)
    #   H3_rxn.scores = ['(n,t)']
    #   #rxn.scores = ['fission', '(n,gamma)']
    #   #rxn.scores = ['fission', '(n,gamma)', '(n,2n)']
    #   H3_rxn.nuclides = ['Li6', 'Be10']
        
    #   return H3_rxn

    def _export_material_to_xml(self):

        # Collect materials from each cells and put them into a materials object
        materials = openmc.Materials()
        root_cell = self.root_cell
        cell_dict = root_cell.get_all_cells()
        # When different cells have the same material, the material should not be
        # counted twice
        id_list = []
        for cell_id in cell_dict:
            cell = cell_dict[cell_id]
            material_dict = cell.get_all_materials()
            material = material_dict[list(material_dict.keys())[0]]
            if material.id not in id_list:
                materials.append(material)
                id_list.append(material.id)

        # If the input materials.xml file is in cwd, remove it
        try:
            os.remove(os.getcwd() + '/materials.xml')
        except OSError:
            pass

        # Set the cross section path again to materials
        materials.cross_sections = self._cross_sections_path +'/cross_sections.xml'

        materials.export_to_xml()

    def _export_geometry_to_xml(self):

        # Collect each cells and put them into a geometry object
        # Need to re-instantiate a universe that will be filled with the modified root cell
        root_universe = self._root_universe
        geometry = openmc.Geometry(root_universe)

        region = self.root_cell.region
        #print (region.get_surfaces())
        for surface_id in region.get_surfaces():
            surface = region.get_surfaces()[surface_id]
        #quit()

        # If the input materials.xml file is in cwd, remove it
        try:
            os.remove(os.getcwd() + '/geometry.xml')
        except OSError:
            pass

        geometry.export_to_xml()


    def _export_tallies_to_xml(self):

        system = self.system

        bucell_dict = system.bucell_dict

        tallies = openmc.Tallies()

        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            flux = self._get_flux_tally(bucell)
            flux_spectrum = self._get_flux_spectrum_tally(bucell)
            rxn = self._get_all_nucl_rxn_tally(bucell)
            # H3_rxn = self.get_H3_rxn_tally(bucell)
            tallies.append(flux)
            tallies.append(flux_spectrum)
            tallies.append(rxn)
            # tallies.append(H3_rxn)

        # If the input tallies.xml file is in cwd, remove it
        try:
            os.remove(os.getcwd() + '/tallies.xml')
        except OSError:
            pass

        tallies.export_to_xml()

    # settings.xml as provided by the user is read
    # It is then modified and set to couple
    # Since it will not change during the simulation, it is not going
    # to be modified again

    @property
    def settings(self):
        """Returns the setting object defined in OpenMC simulation.
        """

        return self._settings

    def _read_user_settings(self):

        system = self.system
        MC_input_path = self.MC_input_path

        file_path = MC_input_path + '/settings.xml'
        tree = ET.parse(file_path)
        root = tree.getroot()
        settings = openmc.Settings()

        for child in root:
            if child.tag == 'particles':
                settings.particles = int(child.text)
                self.partices = int(child.text)
            if child.tag == 'batches':
                settings.batches = int(child.text)
                self.batches = int(child.text)
            if child.tag == 'inactive':
                settings.inactive = int(child.text)
                self.inactive = int(child.text)
            if child.tag == 'source':
                source = child

        settings.output = {'tallies': False}
        settings.temperature = {'method':'nearest', 'tolerance': self.tolerance}

        ll = self.bounding_box[0]
        ur = self.bounding_box[1]
        #uniform_dist = openmc.stats.Box(ll, ur, only_fissionable=True)
        if source is None:
            point = openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
            settings.source = openmc.source.Source(space=point)
        else:
            settings.source = openmc.source.Source.from_xml_element(source)
        # To reduce the size of the statepoint file
        settings.sourcepoint['write'] = False

        self._settings = settings


    # the setting file created by the user should be stored somewhere
    def _export_settings_to_xml(self):

        settings = self.settings

        # If the input settings.xml file is in cwd, remove it
        try:
            os.remove(os.getcwd() + '/settings.xml')
        except OSError:
            pass

        settings.export_to_xml()

    @property
    def particles(self):

        return self._particles

    @particles.setter
    def particles(self, particles):

        self._particles = particles

    @property
    def batches(self):

        return self._batches

    @batches.setter
    def batches(self, batches):

        self._batches = batches

    @property
    def inactive(self):

        return self._inactive

    @inactive.setter
    def inactive(self, inactive):

        self._inactive = inactive

    @property
    def openmc_bin_path(self):
        """Returns the path to OpenMC executable.
        """

        return self._openmc_bin_path

    @openmc_bin_path.setter
    def openmc_bin_path(self, openmc_bin_path):
        """Sets the path to OpenMC executable. It seems that in certain environement (such as HPC clusters), openmc.run() is unable to find the path to OpenMC executable. Defining the path explicitly here avoid this problem.

        Parameters
        ----------
        openmc_bin_path: str
            The path to OpenMC executable
        """

        self._openmc_bin_path = openmc_bin_path

    def _run_openmc(self):

        # New material xml file needs to be written with zero dens nuclides
        self._export_material_to_xml()
        # New geometry xml file needs to be written with new materials id
        self._export_geometry_to_xml()
        # Tallies xml files needs to be writen
        self._export_tallies_to_xml()
        # Settings xml files needs to be written
        self._export_settings_to_xml()

        #openmc_bin_path = '/tigress/jdtdl/openmc/py3-mpi-190324/bin/openmc'
        #openpc_bin_path = '/tigress/mkutt/openmc/py3-mpi/bin/openmc'

        if self.MPI == 'on':
            if self.openmc_bin_path == None:
                openmc.run(mpi_args=[self._exec, '-n', self._tasks])
            else:
                openmc.run(mpi_args=[self._exec, '-n', self._tasks], openmc_exec = self.openmc_bin_path)
        else:
            if self.openmc_bin_path == None:
                openmc.run()
            else:
                openmc.run(openmc_exec = self.openmc_bin_path)

        self._set_statepoint()
        self._set_updated_summary()
        # Append the new kinf to system sequence
        self._set_kinf()

    # This method set decay_lib_set to yes
    # It can be used when the user does not want the system to be set any decay data
    def no_decay(self):
        """Activate a "no-decay" simulation. This means that the system will not model the decay of nuclides
        """
        self._decay_lib_set = 'yes'

    def set_decay_lib(self, decay_lib_path):
        """Sets the decay library to use for ONIX.

        Parameters
        ----------
        decay_lib_path: str
            Path to a decay library. The library needs to be in ONIX format.
        """

        system = self.system
        self._decay_lib_set = 'yes'
        self._decay_lib_path = decay_lib_path
        system.set_decay_for_all(decay_lib_path)

    def set_default_decay_lib(self):
        """Sets the decay library to the default one (ENDF/B-VIII.0).
        """

        system = self.system
        self._decay_lib_set = 'yes'
        #system.set_default_decay_for_all_no_add()
        self._decay_lib_path = 'default'
        system.set_default_decay_for_all()

    def set_decay_from_object(self, bucell, object):
        """Sets the decay library from an ONIX decay object defined by the user for a specific BUCell.

        Parameters
        ----------
        bucell: onix.Cell
            The BUCell for which the decay data are to be used.
        object: onix.utils.decay_lib
            A decay library object constructed by the user with onix.utils.decay_lib.
        """

        system = self.system
        # This should not set yes since it is only for one bucell
        # Need to be fixed later
        self._decay_lib_set = 'yes'
        bucell = system.get_bucell(bucell)
        bucell.set_decay(object)

    def set_xs_lib(self, xs_lib_path):
        """Sets the cross section library to use for ONIX.

        Parameters
        ----------
        xs_lib_path: str
            Path to a cross section library. The library needs to be in ONIX format.
        """

        system = self.system
        self._xs_lib_set = 'yes'
        system.set_xs_for_all(xs_lib_path)

    def set_default_xs_lib(self):
        """Sets the cross section library to the default one. The default cross section library in ONIX has been obtained by computing middle-burnup one-group cross section for a coupled simumation of a LWR fuel pin-cell. Therefore, this library is only adapted for LWR reactors.
        """

        system = self.system
        self._xs_lib_set = 'yes'
        #system.set_default_decay_for_all_no_add()
        system.set_default_xs_for_all()

    def set_user_fy_lib(self, user_fy_lib_path, complete=False):
        """Sets the fission yield library to use for ONIX. Setting the *complete* parameter to True allows ONIX to complement the data of the provided library with additional data found in ENDF/B-VIII.0.

        Parameters
        ----------
        user_fy_lib_path: str
            Path to a fission yield library. The library needs to be in ONIX format.
        complete: bool
            Indicates to ONIX whether or not to complement the provided library with additional fission yields from ENDF:B-VIII.0 library.
        """

        system = self.system
        self._user_fy_lib_set = 'yes'
        if complete == True:
            self._complete_user_fy_lib = 'yes'
        self._user_fy_lib_path = user_fy_lib_path
        system.set_fy_for_all(user_fy_lib_path, complete)

    def set_default_fy_lib(self):
        """Sets the fission yield library to the default one (ENDF/B-VIII.0).
        """
        system = self.system
        self._default_fy_lib_set = 'yes'
        self._default_fy_lib_path = 'default'
        #system.set_default_fy_for_all_no_add()
        system.set_default_fy_for_all()

    def set_fy_from_object(self, bucell, object):
        """Sets the fission yield library to an ONIX fission yield object defined by the user for a specific BUCell.

        Parameters
        ----------
        bucell: onix.Cell
            The BUCell for which the decay data are to be used.
        object: onix.utils.fy_lib
            A fission yield library object constructed by the user with onix.utils.fy_lib.
        """

        system = self.system
        # This should not set yes since it is only for one bucell
        # Need to be fixed later
        self._fy_lib_set = 'yes'
        bucell = system.get_bucell(bucell)
        bucell.set_fy(object)


    # This method reads, samples the isomeric branching data and xs data
    # using the mg energy bin mid points data and fold them together
    # def fold_sampled_isomeric_xs_data(self):

    #   sampled_iso_data = self.get_sampled_isomeric_branching_data()
    #   sampled_xs_data = self.get_sampled_ng_cross_section_data()

    #   iso_xs_data = {}
    #   for nucl in sampled_iso_data:
    #       iso_data = sampled_iso_data[nucl]
    #       xs_data = sampled_xs_data[nucl]
    #       iso_xs_data[nucl] = {}
    #       iso_xs_data[nucl]['0'] = [x*y for x,y in zip(iso_data['0'], xs_data)]
    #       iso_xs_data[nucl]['1'] = [x*y for x,y in zip(iso_data['1'], xs_data)]
        
    #   self._iso_xs_data = iso_xs_data

    def _set_sampled_isomeric_branching_data(self):

        print ('\n\n\n***********Sampling isomeric branching data***********\n\n\n')
        isomeric_branching_data = data.read_isomeric_data()
        sampled_isomeric_branching_data = {}

        for nucl in isomeric_branching_data:
            nucl_data = isomeric_branching_data[nucl]

            # print (nucl_data['0'])
            sampled_isomeric_branching_data[nucl] = {}
            sampled_isomeric_branching_data[nucl]['0'] = nucl_data['0'](self.mg_energy_mid_points)
            sampled_isomeric_branching_data[nucl]['1'] = nucl_data['1'](self.mg_energy_mid_points)

        self._sampled_isomeric_branching_data = sampled_isomeric_branching_data

    # This method reads and samples the point-wise cross section for ng of the nuclides that
    # have ng isomeric branching data only
    def _set_sampled_ng_cross_section_data(self):

        print ('\n\n\n***********Sampling point-wise cross section data***********\n\n\n')
        sampled_isomeric_branching_data = self._sampled_isomeric_branching_data
        sampled_ng_cross_section_data = {}
        cross_section_path = self._cross_sections_path
        cross_sections_files_name = os.listdir(cross_section_path)

        total_count = 1
        file_name_list = []
        for file_name in cross_sections_files_name:
            nucl_name = file_name.replace('.h5', '')
            if nucl_name in sampled_isomeric_branching_data:
                total_count += 1
                file_name_list.append(file_name)

        start = time.time()
        count = 1
        for file_name in file_name_list:
            nucl_name = file_name.replace('.h5', '')
            if self._few_isomeric == 'on':
                if nucl_name in ['Pm147', 'Am241']: # make it a shorter
                    nucl_path = cross_section_path+'/{}'.format(file_name)
                    xs_data = openmc.data.IncidentNeutron.from_hdf5(nucl_path)
                    ng_xs_data = xs_data[102].xs[xs_data.temperatures[0]]
                    print ('--- Sampling {} (n,gamma) point-wise cross section --- [{}/{}]'.format(nucl_name, count, total_count))
                    sampled_ng_xs_data = ng_xs_data(self.mg_energy_mid_points)
                    sampled_ng_cross_section_data[nucl_name] = sampled_ng_xs_data
                    count += 1
            elif self._few_isomeric == 'off':
                nucl_path = cross_section_path+'/{}'.format(file_name)
                xs_data = openmc.data.IncidentNeutron.from_hdf5(nucl_path)
                ng_xs_data = xs_data[102].xs[xs_data.temperatures[0]]
                print ('--- Sampling {} (n,gamma) point-wise cross section --- [{}/{}]'.format(nucl_name, count, total_count))
                sampled_ng_xs_data = ng_xs_data(self.mg_energy_mid_points)
                sampled_ng_cross_section_data[nucl_name] = sampled_ng_xs_data
                count += 1  
        end = time.time()
        print('\n Time to sample cross sections: {}'.format(end - start))

        self._sampled_ng_cross_section_data = sampled_ng_cross_section_data

    def _set_tallies_to_bucells(self, s):

        system = self.system
        bucell_dict = system.bucell_dict
        sp = self.statepoint
        summary = self.updated_summary
        xs_mode  = self.xs_mode
        sampled_isomeric_branching_data = self._sampled_isomeric_branching_data
        sampled_ng_cross_section_data = self._sampled_ng_cross_section_data

        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]

            # densities from OpenMC are extracted
            # This is for nuclides which are zero in onix but set to 1E-24 in OpenMC
            # 1E-24 fraction percent does not translate into 1E-24 atm/cm3 always
            # Therefore, the code needs to divide the reaction rate by the correct densities
            # taken from summary.geometry.cell.material
            cell = summary.geometry.get_cells_by_name(bucell.name)[0]
            material_dict = cell.get_all_materials()
            material = material_dict[list(material_dict.keys())[0]]
            mc_nuclides_densities = material.get_nuclide_atom_densities()

            flux_tally = sp.get_tally(name = '{} flux'.format(bucell.name))
            flux_spectrum_tally = sp.get_tally(name = '{} flux spectrum'.format(bucell.name))
            rxn_rate_tally = sp.get_tally(name = '{} rxn rate'.format(bucell.name))
            # H3_rxn_rate_tally = sp.get_tally(name = '{} H3 rxn rate'.format(bucell.name))
            bucell._set_MC_tallies(mc_nuclides_densities, flux_tally, flux_spectrum_tally, rxn_rate_tally, sampled_isomeric_branching_data, sampled_ng_cross_section_data, xs_mode, s)

        # YOU NEED TO CREATE LIST TO STORE EACH NEW XS

    # Normalize the flux in each cell with the FMF after each openmc calculation
    # Calculate the power of each cell from the normalized flux
    def _step_normalization(self, s):

        system = self.system
        sequence = self.sequence
        master_bucell = sequence.master_bucell
        if master_bucell == None:
            FMF = sequence._get_system_FMF(system, s)
        else:
            FMF = sequence._get_bucell_FMF(master_bucell, s)

        bucell_list = system.get_bucell_list()
        for bucell in bucell_list:
            bucell_sequence = bucell.sequence
            MC_flux = bucell_sequence.current_MC_flux
            # MC_flux is volume integrated (unit cm.sp
            # FMF is in sp.s
            # to have flux in cm.s you need to divide by volule of the cell
            flux = FMF*MC_flux/bucell.vol
            pow_dens = bucell._update_pow_dens(flux)
            bucell_sequence._set_macrostep_flux(flux)
            bucell_sequence._set_macrostep_pow_dens(pow_dens)

    # Normalize the flux in each cell with the FMF after each openmc calculation
    # Calculate the power of each cell from the normalized flux
    # This method seems obsolete
    def initial_couple_step_normalization(self, norma_mode):

        system = self.system
        sequence = self.sequence
        FMF = sequence._get_FMF1(system, 0)

        bucell_list = system.get_bucell_list()
        for bucell in bucell_list:
            bucell_sequence = bucell.sequence
            MC_flux = bucell_sequence.current_MC_flux
            flux = FMF*MC_flux
            pow_dens = bucell._update_pow_dens(flux)
            bucell_sequence._set_initial_flux(flux)
            bucell_sequence._set_initial_pow_dens(pow_dens)

    def _copy_MC_files(self, s, final='no'):

        if final == 'no':
            step_folder = '/step_{}'.format(s)
        elif final =='yes':
            # If it is the final step, sequence has not been requested to create
            # a step_final folder. This function needs to do it now.
            final_step_path = utils.gen_folder('step_final')
            step_folder = '/step_final'

        openmc_file_path = os.getcwd()+step_folder+'/OpenMC'

        os.mkdir(openmc_file_path)

        all_MC_files = glob.glob('*.xml') + glob.glob('tallies.out') + glob.glob('*.h5')
        for file_name in all_MC_files:
            try:
                shutil.copyfile(os.getcwd()+'/{}'.format(file_name), openmc_file_path + '/{}'.format(file_name))
            except IOError:
                continue

        for file_name in all_MC_files:
            os.remove(os.getcwd() + '/{}'.format(file_name))

    def _change_temperature(self, s):

        # It seems that temperature set in cells trumps the temperature set in material
        # Therefore, this method change temperature in cell

        root_cell = self.root_cell
        cell_dict = root_cell.get_all_cells()
        sequence = self.sequence
        temperature_change_dict = sequence.temperature_change_dict
        if temperature_change_dict != None:
            for cell_id in cell_dict:
                cell = cell_dict[cell_id]
                if cell.name in temperature_change_dict:
                    cell_temperature_change = temperature_change_dict[cell.name]
                    if s in cell_temperature_change:
                        new_temperature = cell_temperature_change[s]
                        cell.temperature = new_temperature


    def _set_dens_to_cells(self):

        summary = self.initial_summary
        system = self.system
        bucell_dict = system.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            tally_nucl_list = utils.mc_namelist_to_bu_namelist(self._get_nucl_to_be_tallied(bucell))
            cell = summary.geometry.get_cells_by_name(bucell.name)[0]
            material_dict = cell.get_all_materials()
            material = material_dict[list(material_dict.keys())[0]]
            ### WARNING Openmc add_nuclide in atom percent is not in absolute atom percent
            ### but in relative ao, i.e., setting U235 and U238 to 0.34 and 0.34 will be the same
            ### as setting them to 50% and 50 %
            ### Therefore the ao value that needs to be updated needs to be normalize by the tot dens of
            ### only the nuclides to be tallied and not by the total dens of all nuclide in bucell
            tally_nucl_list_dens = bucell._get_subtotal_dens_counting_zero_dens(tally_nucl_list)
            material.set_density('atom/b-cm', tally_nucl_list_dens)
            for nucl in tally_nucl_list:
                openmc_nucl_name = utils.onix_name_to_openmc_name(nucl)
                material.remove_nuclide(openmc_nucl_name) 
                # No need to calculate ao, just directly give density
                #nucl_subao = bucell.get_nucl_subao(nucl, tally_nucl_list)
                nucl_dens = bucell._get_nucl_dens_for_openmc(nucl)
                material.add_nuclide(openmc_nucl_name, nucl_dens)

    def _set_MC_XS_nuc_list_to_bucells(self):

        system = self.system
        bucell_dict = system.bucell_dict
        MC_XS_nucl_list = self.MC_XS_nucl_list
        MC_XS_nucl_list_obu_name =  utils.mc_namelist_to_bu_namelist(MC_XS_nucl_list)
        MC_XS_nucl_list_zamid = utils.name_list_to_zamid_list(MC_XS_nucl_list_obu_name)
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.MC_XS_nucl_list = MC_XS_nucl_list_zamid


    def burn(self):
        """Launches the coupled simulation.
        """

        start_time = time.time()

        # If no decay libs have been set, set default libs
        if self._decay_lib_set == 'no':
            self.set_default_decay_lib()
            print ('\n\n\n----  Default decay constants library set for system  ----\n---- {} ----'.format(data.default_decay_b_lib_path))
        else:
            print ('\n\n\n----  User defined path for decay library  ----\n\n')
            print ('----  {}  ----\n\n\n'.format(self._decay_lib_path))
        

        # If user has not set a fission yield library, default library (ENDF/B-VII.O) is set
        if self._user_fy_lib_set == 'no':
            self.set_default_fy_lib()
            print ('\n\n\n----  Default fission yields library set for system  ----\n---- {} ----'.format(data.default_fy_lib_path))
        
        # If user has set a fission library, two options:
        elif self._user_fy_lib_set == 'yes':
            print ('\n\n\n----  User defined path for fission yields library ----\n\n')
            print ('----  {}  ----\n\n\n'.format(self._user_fy_lib_path))

            # 1) Default library used to complete user library
            if self._complete_user_fy_lib == 'yes':
                print ('\n\n\n----  User defined fission yields library completed with default fission yields library ----\n\n')
                print ('----  {}  ----\n\n\n'.format(data.default_fy_lib_path))

            # 2) Default library is not used, ONIX only uses user defined fission yields library

        # self.set_default_fy_lib()
        # print ('\n\n\n----  Default fission yields library set for system  ----\n---- {} ----'.format(data.default_fy_lib_path))
        # # If user has provided own fission yield library, data from this library extand and overwrite default library
        # if self._user_fy_lib_set == 'yes':
        #   print ('\n\n\n----  User defined path for fission yields library ----\n\n')
        #   print ('----  {}  ----\n\n\n'.format(self._fy_lib_path))
        
        #print (self.xs_mode, self._xs_lib_set)
        if self.xs_mode == 'constant lib' and self._xs_lib_set == 'no':
            self.set_default_xs_lib()
            print ('\n\n\n----Default cross section library set for system----\n\n\n')
        else:
            # This method simply passes the MC_XS_nucl_list to each cell so that each cell
            # can then build its own lib_nucl_list
            self._set_MC_XS_nuc_list_to_bucells()

            print ('\n\n\n----  Path for cross sections library ----\n\n')
            print ('----  {}  ----\n\n\n'.format(self._cross_sections_path))

        self._set_sampled_isomeric_branching_data()
        self._set_sampled_ng_cross_section_data()



        system = self.system
        sequence = system.sequence
        norma_mode = sequence.norma_unit

        bucell_list = system.get_bucell_list()
        for bucell in bucell_list:
            # Check if different nuclide list (initial list, lib list and nucl set (user defined set of nuclides to be considered))
            # are consistent with each other
            bucell._check_nucl_list_consistency()

            # Create a list of all nuclides that should be produced, i.e., that belong to the network tree
            #bucell._reduce_nucl_set()


        # This should be somewhere else but for now it is done here
        system.zam_order_passlist()

        #steps_number = sequence.steps_number
        steps_number = sequence.macrosteps_number
        # Shift loop from 1 in order to align loop s and step indexes
        for s in range(1, steps_number+1):

            print ('\n\n\n\n====== STEP {}======\n\n\n\n'.format(s))
            sequence._gen_step_folder(s)
            print (('\n\n\n=== OpenMC Transport {}===\n\n\n'.format(s)))
            self._change_temperature(s)
            self._run_openmc()
            self._set_tallies_to_bucells(s)
            self._step_normalization(s)
            self._copy_MC_files(s)
            print (('\n\n\n=== Salameche Burn {} ===\n\n\n'.format(s)))
            print (utils.printer.salameche_header)
            salameche.burn_step(system, s, 'couple')
            self._set_dens_to_cells()
        
        # This last openmc_run is used to compute the last burnup/time point kinf
        print ('\n\n\n=== OpenMC Transport for Final Point ===\n\n\n')
        self._run_openmc()

        system._gen_output_summary_folder()
        if self.reac_rank == 'on':
            system._print_summary_allreacs_rank()
        system._print_summary_subdens()
        system._print_summary_dens()
        system._print_summary_xs()
        system._print_summary_flux_spectrum(self.mg_energy)
        system._print_summary_kinf()
        system._print_summary_param()
        system._print_summary_isomeric_branching_ratio()
        self._copy_MC_files(s, final='yes')

        run_time = time.time() - start_time
        print ('\n\n\n >>>>>> ONIX burn took {} seconds <<<<<<< \n\n\n'.format(run_time))

class Initial_nuclides_not_in_nuclide_list(Exception):
    """Raise when some initial nuclides are not included in nucl_list """
    pass

class STOP(Exception):
    """Just a way to stop the code"""
    pass
