import numpy as np
import re
import time
import openmc
import onix.compute as compute
from onix.cell import Cell
from onix.system import System
from onix import salameche
from .couple_openmc import *
from .openmc_fix import *
from onix.salameche import mat_builder as mb
from onix.salameche import cram

from onix import utils
from onix import data

class Couple_msr(Couple_openmc):
    """This class is used to execute coupled-mode simulations, for circulating fuel systems

    Through this class, the user can:

    - chose to parallelize the Monte Carlo simulations via the "set_MPI" method

    - set the nuclear data libraries to be used in the simulation

    - set the burnup/time sequence for the simulation
    - define the fuel circulation loop by defining the connetions between OpenMC cells.
    - select which of the Cells defined in the OpenMC input should be depleted (these cells will be known as BUCells)
    - select a list of nuclides in each BUCell which cross sections will be tallied. If the user does not specify any such list of nuclides, ONIX will by
    default tally the cross sections of all
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

    _ELEMENT_SYMBOLS = {symbol.lower(): symbol for symbol in data.nuc_name_dic}
    _ONIX_NUCLIDE_RE = re.compile(r'^([A-Za-z]+)-(\d+)(\*|[mMnN]1?)?$')
    _OPENMC_NUCLIDE_RE = re.compile(r'^([A-Za-z]+)(\d+)(?:_?([mMnN]1?))?$')




    def __init__(self, MC_input_path = None, xs_mode = 'no constant lib', MPI = None, vol_flow_rate = 0):

        super().__init__(MC_input_path, xs_mode, MPI)

        # By default flow rate in the loop is set to 0 m3/s
        self._vol_flow_rate = vol_flow_rate

    @property
    def vol_flow_rate(self):
        return self._vol_flow_rate

    @vol_flow_rate.setter
    def vol_flow_rate(self, value):
        self._vol_flow_rate = value

    def select_bucells(self, loop_dict, loop_nucl_list=None):
        """Selects the cells loop from the OpenMC input that should be depleted.

        If the user does not specify for which nuclides cross sections should be updated,
        only the OpenMC cell objects should be entered. ONIX will by default calculate the cross sections
        for all nuclides available in the cross section library.

        Parameters
        ----------
        loop_dict: a dicotonray of OpenMC cell names. The value for each key is a dicrotary with two keys
                   the 'feed' key is a list of tuples where the 1st input is an openMC cell names that feed materials to the key cell
                   and the second input is the fraction of the volumetric flow rate comeing from this cell.
                   the 'sink' key is a list of tuples where the 1st input is an openMC cell names that recives materials to the key cell
                   and the second input is the fraction of the volumetric flow rate going into this cell.

        loop_nucl_list: the list of nuclides to be tracked
        """
        self.loop_dict = loop_dict
        self.loop_nucl_list = loop_nucl_list

        self.selected_bucells_name_list = list(self.loop_dict.keys())
        #self.selected_bucells_nucl_list_dict = {key: self.loop_nucl_list for key in self.loop_dict}
        
        

        if self.loop_nucl_list is not None:
            self.selected_bucells_nucl_list_dict = {key: self.loop_nucl_list for key in self.loop_dict}
        else:
            self.selected_bucells_nucl_list_dict = {}

    @classmethod
    def _normalize_element_symbol(cls, element_symbol):
        symbol = '{}'.format(element_symbol).strip().replace(' ', '')
        if symbol == '':
            raise ValueError('Element name cannot be empty')

        normalized_symbol = cls._ELEMENT_SYMBOLS.get(symbol.lower())
        if normalized_symbol is None:
            raise ValueError('Unknown element name {}'.format(element_symbol))

        return normalized_symbol

    @classmethod
    def _classify_species_name(cls, species_name):
        species = '{}'.format(species_name).strip().replace(' ', '')
        if species == '':
            raise ValueError('Nuclide or element name cannot be empty')

        if utils.is_zamid(species):
            return 'nuclide', utils.zamid_to_name(species)

        match = cls._ONIX_NUCLIDE_RE.fullmatch(species)
        if match is None:
            match = cls._OPENMC_NUCLIDE_RE.fullmatch(species)

        if match is not None:
            element = cls._normalize_element_symbol(match.group(1))
            mass_number = int(match.group(2))
            state = '*' if match.group(3) is not None else ''
            return 'nuclide', '{}-{}{}'.format(element, mass_number, state)

        if species.lower() in cls._ELEMENT_SYMBOLS:
            return 'element', cls._normalize_element_symbol(species)

        raise ValueError('Unknown nuclide or element name {}'.format(species_name))

    @classmethod
    def _normalize_nuclide_name(cls, nuclide_name):
        species_type, normalized_name = cls._classify_species_name(nuclide_name)
        if species_type != 'nuclide':
            raise ValueError('Unknown nuclide name {}'.format(nuclide_name))

        return normalized_name

    @classmethod
    def _normalize_element_name(cls, element_name):
        species_type, normalized_name = cls._classify_species_name(element_name)
        if species_type == 'element':
            return normalized_name

        return normalized_name.split('-')[0]

    def _build_bucell_lookup(self, bucell_list):
        bucell_index_by_name = {}
        bucell_by_name = {}
        for index, bucell in enumerate(bucell_list):
            bucell_index_by_name[bucell.name] = index
            bucell_by_name[bucell.name] = bucell

        loop_name_set = set(self.loop_dict.keys())
        bucell_name_set = set(bucell_index_by_name.keys())
        if loop_name_set != bucell_name_set:
            missing_from_loop = bucell_name_set - loop_name_set
            missing_from_system = loop_name_set - bucell_name_set
            raise ValueError(
                'loop_dict keys must match the selected BUCells. Missing from loop_dict: {}. '
                'Missing from system.get_bucell_list(): {}.'.format(
                    sorted(missing_from_loop),
                    sorted(missing_from_system),
                )
            )

        return bucell_index_by_name, bucell_by_name

    def _build_passport_lookup(self, passport_list):
        nuclide_index_by_name = {}
        element_indices_by_name = {}

        for index, passport in enumerate(passport_list):
            nuclide_name = self._normalize_nuclide_name(passport.name)
            nuclide_index_by_name[nuclide_name] = index
            element_name = nuclide_name.split('-')[0]
            if element_name not in element_indices_by_name:
                element_indices_by_name[element_name] = []
            element_indices_by_name[element_name].append(index)

        return nuclide_index_by_name, element_indices_by_name

    def _build_element_rate_vector(self, m, element_indices_by_name, element_rate_pairs):
        values = np.zeros(m)

        for element_name, rate in element_rate_pairs:
            normalized_element = self._normalize_element_name(element_name)
            for index in element_indices_by_name.get(normalized_element, []):
                values[index] += rate

        return values

    def _build_nuclide_rate_vector(self, m, nuclide_index_by_name, nuclide_rate_pairs, scale=1.0):
        values = np.zeros(m)

        for nuclide_name, rate in nuclide_rate_pairs:
            normalized_nuclide = self._normalize_nuclide_name(nuclide_name)
            index = nuclide_index_by_name.get(normalized_nuclide)
            if index is not None:
                values[index] += scale * rate

        return values

    def _build_species_rate_vector(
        self,
        m,
        nuclide_index_by_name,
        element_indices_by_name,
        species_rate_pairs,
        scale=1.0,
    ):
        values = np.zeros(m)

        for species_name, rate in species_rate_pairs:
            species_type, normalized_name = self._classify_species_name(species_name)
            if species_type == 'nuclide':
                index = nuclide_index_by_name.get(normalized_name)
                if index is not None:
                    values[index] += scale * rate
            else:
                for index in element_indices_by_name.get(normalized_name, []):
                    values[index] += scale * rate

        return values

    def _get_bucell_index(self, cell_name, bucell_index_by_name):
        if cell_name not in bucell_index_by_name:
            raise KeyError(
                "Loop cell '{}' is referenced in the circulation definition but is not part "
                'of system.get_bucell_list()'.format(cell_name)
            )

        return bucell_index_by_name[cell_name]

    def _get_normalized_passport_order(self, passport_list):
        return [self._normalize_nuclide_name(passport.name) for passport in passport_list]

    @staticmethod
    def _get_order_mismatches(source_order, destination_order):
        mismatch_list = []
        for local_index in range(min(len(source_order), len(destination_order))):
            if source_order[local_index] != destination_order[local_index]:
                mismatch_list.append(
                    (local_index, source_order[local_index], destination_order[local_index])
                )

        return mismatch_list

    def _validate_bucell_block_order(self, bucell_list, reference_order, context):
        for bucell in bucell_list:
            local_order = self._get_normalized_passport_order(bucell.passlist.passport_list)
            if local_order == reference_order:
                continue

            mismatch_list = self._get_order_mismatches(reference_order, local_order)
            mismatch_preview = []
            for local_index, reference_name, local_name in mismatch_list[:10]:
                mismatch_preview.append(
                    'local {}: reference {} != cell {}'.format(
                        local_index,
                        reference_name,
                        local_name,
                    )
                )
            if len(reference_order) != len(local_order):
                mismatch_preview.insert(
                    0,
                    'length mismatch: reference {} != cell {}'.format(
                        len(reference_order),
                        len(local_order),
                    )
                )
            if len(mismatch_list) > 10:
                mismatch_preview.append(
                    '{} additional mismatches omitted'.format(len(mismatch_list) - 10)
                )

            raise ValueError(
                'MSR {} requires every BUCell to use the same passlist order. '
                "Cell '{}' does not match the reference cell. {}".format(
                    context,
                    bucell.name,
                    '; '.join(mismatch_preview),
                )
            )

    def _build_connectivity_matrix(self, system):
        bucell_list = system.get_bucell_list()
        dummy_cell = bucell_list[0]
        passlist = dummy_cell.passlist

        m = len(mb.get_initial_vect(passlist))
        n = m * len(bucell_list)

        connectivity_matrix = np.zeros((n, n))
        fixed_source = np.zeros(n)
        eye_m = np.eye(m)

        bucell_index_by_name, bucell_by_name = self._build_bucell_lookup(bucell_list)
        reference_order = self._get_normalized_passport_order(passlist.passport_list)
        self._validate_bucell_block_order(
            bucell_list,
            reference_order,
            'connectivity matrix construction',
        )

        for i, bucell in enumerate(bucell_list):
            cell_name = bucell.name
            cell_vol = bucell.vol
            passport_list = bucell.passlist.passport_list
            row_slice = slice(i*m, i*m + m)
            loop_data = self.loop_dict[cell_name]

            nuclide_index_by_name, element_indices_by_name = self._build_passport_lookup(passport_list)

            if 'upstream_cells' in loop_data:
                feed_cell_list = loop_data['upstream_cells']
                sink_flow_fraction = 0

                for feed_cell_name, flow_fraction in feed_cell_list:
                    sink_flow_fraction += flow_fraction
                    feed_rate = self.vol_flow_rate * flow_fraction / cell_vol
                    j = self._get_bucell_index(feed_cell_name, bucell_index_by_name)
                    col_slice = slice(j*m, j*m + m)
                    connectivity_matrix[row_slice, col_slice] += feed_rate * eye_m

                sink_rate = self.vol_flow_rate * sink_flow_fraction / cell_vol
                connectivity_matrix[row_slice, row_slice] -= sink_rate * eye_m

            if 'elements_removal' in loop_data:
                removal_vector = self._build_element_rate_vector(
                    m,
                    element_indices_by_name,
                    loop_data['elements_removal'],
                )
                connectivity_matrix[row_slice, row_slice] -= np.diag(removal_vector)

            if 'nuclide_removal' in loop_data:
                removal_vector = self._build_nuclide_rate_vector(
                    m,
                    nuclide_index_by_name,
                    loop_data['nuclide_removal'],
                )
                connectivity_matrix[row_slice, row_slice] -= np.diag(removal_vector)

            if 'external_nuclide_addition' in loop_data:
                addition_vector = self._build_nuclide_rate_vector(
                    m,
                    nuclide_index_by_name,
                    loop_data['external_nuclide_addition'],
                    scale=1e-24/cell_vol,
                )
                fixed_source[row_slice] += addition_vector

            if 'holdup' in loop_data:
                for feed_cell_name, removed_elements in loop_data['holdup']:
                    j = self._get_bucell_index(feed_cell_name, bucell_index_by_name)
                    feed_cell_volume = bucell_by_name[feed_cell_name].vol
                    removal_vector = self._build_element_rate_vector(
                        m,
                        element_indices_by_name,
                        removed_elements,
                    )
                    col_slice = slice(j*m, j*m + m)
                    connectivity_matrix[row_slice, col_slice] += (
                        feed_cell_volume / cell_vol
                    ) * np.diag(removal_vector)

            if 'external_salt_stream' in loop_data:
                flow_rate, feed_stream = loop_data['external_salt_stream']
                connectivity_matrix[row_slice, row_slice] -= flow_rate * eye_m / cell_vol
                density_vector = self._build_nuclide_rate_vector(
                    m,
                    nuclide_index_by_name,
                    feed_stream,
                    scale=flow_rate/cell_vol,
                )
                fixed_source[row_slice] += density_vector

            if 'recycling' in loop_data:
                destination_cell, extracted_nuc = loop_data['recycling']
                j = self._get_bucell_index(destination_cell, bucell_index_by_name)
                destination_cell_volume = bucell_by_name[destination_cell].vol
                transfer_vector = self._build_species_rate_vector(
                    m,
                    nuclide_index_by_name,
                    element_indices_by_name,
                    extracted_nuc,
                    scale=cell_vol/destination_cell_volume,
                )
                row_destination = slice(j*m, j*m + m)
                connectivity_matrix[row_destination, row_slice] += np.diag(transfer_vector)

        return connectivity_matrix, fixed_source


    def burn(self):
        """Launches the coupled simulation.
        """

        start_time = time.time()
        print ('\n\n\n----  Compute backend: {}  ----\n'.format(compute.describe_compute_backend()))

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
        system.zam_order_passlist()
        sequence = system.sequence
        norma_mode = sequence.norma_unit

        bucell_list = system.get_bucell_list()
        connectivity_matrix, fixed_source = self._build_connectivity_matrix(system)

        self.connectivity_matrix = connectivity_matrix
        self.fixed_source = fixed_source
        
        for bucell in bucell_list:
            # Check if different nuclide list (initial list, lib list and nucl set (user defined set of nuclides to be considered))
            # are consistent with each other
            bucell._check_nucl_list_consistency()


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
            norma = sequence.norma_unit
            if norma == 'power':
                self._step_normalization(s)
            else:
                self._step_normalization2(s)
            self._copy_MC_files(s)
            print (('\n\n\n=== Salameche Burn {} ===\n\n\n'.format(s)))
            print (utils.printer.salameche_header)
            self.burn_step(system, s)
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


    def burn_step(self, system, s):

        """Depletes the system for macrostep s.

           Parameters
           ----------
           system: onix.System
                  System to be depleted
           s: int
           Macrostep number
        """

        bucell_list = system.get_bucell_list()
        reac_rank = system.reac_rank

        dummy_cell = bucell_list[0]
        passlist = dummy_cell.passlist
        all_nuc = passlist.nucl_list

        m = len(mb.get_initial_vect(passlist))
        n = m * len(bucell_list)
        reference_order = self._get_normalized_passport_order(passlist.passport_list)
        self._validate_bucell_block_order(
            bucell_list,
            reference_order,
            'burn step assembly',
        )

        giant_B = np.zeros((n,n))
        giant_C = np.zeros((n,n))
        giant_N = np.zeros(n)

        connectivity_matrix = self.connectivity_matrix
        fixed_source = self.fixed_source
        if connectivity_matrix.shape != (n, n):
            raise ValueError(
                'MSR connectivity matrix shape {} does not match the current global '
                'state vector shape {}. A BUCell passlist likely changed after the '
                'connectivity matrix was built.'.format(connectivity_matrix.shape, (n, n))
            )
        if fixed_source.shape != (n,):
            raise ValueError(
                'MSR fixed source shape {} does not match the current global state '
                'vector length {}.'.format(fixed_source.shape, n)
            )

        for i, bucell in enumerate(bucell_list):
            bucell._set_folder()
            passlist = bucell.passlist
            sequence = bucell.sequence

            # microsteps_number is of length s-1
            microsteps_number = sequence.microsteps_number(s-1)

            B = mb.get_xs_mat(passlist)
            C = mb.get_decay_mat(passlist)
            N = mb.get_initial_vect(passlist)
            mb._print_all_mat_to_text(B, C, bucell, s)

            giant_B[i*m:i*m+m, i*m:i*m+m]= B
            giant_C[i*m:i*m+m, i*m:i*m+m]= C
            giant_N[i*m:i*m+m]= N


        A = np.zeros((n,n))
        for i in range(microsteps_number):
            for j, bucell in enumerate(bucell_list):
                sequence = bucell.sequence
                bucell_id = bucell.id
                norma = sequence.norma_unit
                sequence._bucell_time_bu_substep_conversion(bucell, s, i)
                
                time_point = sequence.time_subpoint(s, i)
                bucell_bu_point = sequence.bucell_bu_subpoint(s, i)
                time_substep = sequence.get_time_subintvl(s, i)

                pow_dens = sequence.current_pow_dens
                flux = sequence.current_flux

                if i != 0:
                # Check whether actinides are present in the cell
                # If not, then no flux/power update should be done
                    act = bucell.check_act_presence()
                    if act == 'yes':
                # Now that the density of nuclides is updated, calculate the new substep flux or power density
                       if norma == 'power' and flux != 0.0:
                          flux = bucell._update_flux(pow_dens)

                       elif norma == 'flux':
                          pow_dens = bucell._update_pow_dens(flux)

                sequence._set_substep_flux(flux, s, i)
                sequence._set_substep_pow_dens(pow_dens, s, i)


                A[j*m:j*m+m,:] = (giant_B[j*m:j*m+m,:]*1e-24*flux + giant_C[j*m:j*m+m,:] + connectivity_matrix[j*m:j*m+m,:])

            At = A*time_substep
            
            
            print('micro_step:', i)
            giant_N = cram.CRAM16(At, giant_N)
            if np.any(fixed_source != 0):
                print ('Particular solution')
                tt = time.time()
                integral = cram.CRAM16(At, fixed_source) - fixed_source
                Np = compute.solve(A, integral)
                Np[Np < 0] = 0
                print('Particular solution took:{} s'.format(time.time() - tt))
                giant_N += Np
            
            for j, bucell in enumerate(bucell_list):
                

               print(bucell.name)
               cram.CRAM_density_check(bucell, giant_N[j*m:j*m+m])
               bucell._update_dens(giant_N[j*m:j*m+m], i, microsteps_number)

    # Generate the allreacsdic for each nuclide
               if reac_rank == 'on':

                  bucell._set_allreacs_dic(s, i, microsteps_number)

            
        for bucell in bucell_list:
            sequence = bucell.sequence
            bucell._set_step_dens()
            sequence._set_macrostep_bucell_bu()
            bucell._change_isotope_density(s)
            bucell._change_total_density(s)
            bucell._print_substep_dens(s)

        if reac_rank == 'on':
           system._print_current_allreacs_rank()
        system._copy_cell_folders_to_step_folder(s)

    def _step_normalization2(self, s):

        system = self.system
        sequence = self.sequence

        bucell_list = system.get_bucell_list()
        tot_flux = 0
        for bucell in bucell_list:
            bucell_sequence = bucell.sequence
            MC_flux = bucell_sequence.current_MC_flux
            tot_flux  += MC_flux
        for bucell in bucell_list:
            bucell_sequence = bucell.sequence
            MC_flux = bucell_sequence.current_MC_flux
            avg_flux = sequence.norma_vector[s-1]

            flux = (MC_flux/tot_flux)*avg_flux
            pow_dens = bucell._update_pow_dens(flux)
            bucell_sequence._set_macrostep_flux(flux)
            bucell_sequence._set_macrostep_pow_dens(pow_dens)

