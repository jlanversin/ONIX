import onix.utils as utils
from . import salameche
import os

class System(object):
    """The system oject contains and organizes the ensemble of BUCells.

    When running a standalone simulation, the system object is the the object with which the user can set certain attributes and parameters for the simulation (nuclear data libraries, burnup sequence, new BUCells).
    In a coupled simulation, the class onix.couple.Couple_openmc fullfills that role and passes-on directives to the system object.

    Parameters
    ----------
    system_id: int
        Number id for the system object
    """
    def __init__(self, system_id):

        self._id = system_id
        self._bucell_dict = None
        self._bounding_box = None
        self._sequence = None
        self._output_summary_path = None

        self._reac_rank = 'off'

    @property
    def id(self):
        """Returns the number id of the system object.
        """
        return self._id
    
    @property
    def bucell_dict(self):
        """Returns a dictionnary where keys are BUCell number ids and entries are corresponding BUCells objects"""
        return self._bucell_dict

    @bucell_dict.setter
    def bucell_dict(self, bucell_dict):

        self._bucell_dict = bucell_dict

    # Probably obsolete
    def add_bucell_dict(self, new_bucell_dict):

        # You need to make sure that there are no cells with same ids
        if self.bucell_dict == None:
            self.bucell_dict = new_bucell_dict
        else:
            old_bucell_dict = self.bucell_dict.copy()
            updated_bucell_dict = {**old_bucell_dict, **new_bucell_dict}
            self.bucell_dict = updated_bucell_dict

    def add_bucell(self, new_bucell):
        """Adds a new BUCell to System.

        Parameters
        ----------
        new_BUCell: onix.Cell
        """
        if self.bucell_dict == None:
            self.bucell_dict = {new_bucell.id : new_bucell}
        else:
            # You need to make sure this id is not already taken by other cell
            updated_bucell_dict = self.bucell_dict.copy()
            updated_bucell_dict[new_bucell.id] = new_bucell

            self.bucell_dict = updated_bucell_dict

    def get_bucell(self, name):
        """Returns the BUCell object with corresponding name.

        Parameters
        ----------
        name: str
            Name of the BUCell
        """
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
        """Returns a list of the BUCell objects"""

        bucell_dict = self.bucell_dict
        bucell_list = utils.cell_dict_to_cell_list(bucell_dict)

        return bucell_list

    
    @property
    def bounding_box(self):
        """Returns the bounding box defined by the user for the system"""

        return self._bounding_box

    # This function is only used within the code by Coupe_openmc. No need in standalone.
    @bounding_box.setter
    def bounding_box(self, bounding_box):

        self._bounding_box = bounding_box

    @property
    def reac_rank(self):

        return self._reac_rank
    
    def reac_rank_on(self):
        """Calling this method will tell ONIX to produce production and destrubtion reaction rates ranking for each nuclide and print the data for each BUCells.
        By default ONIX does not produce reaction rates ranking as it takes a lot of memory.
        """
        self._reac_rank = 'on'

    @property
    def total_vol(self):
        """Returns the total volume of the system. The total volume of the system is the volume of all BUCells plus the volume of regions that are not depleted (moderator for instance if no BUCell is created for the moderator).

        In standalone mode, the total volume will simply be the addition of all BUCells' volume as all regions should be depleted.
        """

        return self._total_vol

    @total_vol.setter
    def total_vol(self, total_vol):
        """Sets the total volume of the system.

         The total volume of the system is the volume of all BUCells plus the volume of regions that are not depleted (moderator for instance if no BUCell is created for the moderator).

        In standalone mode, the total volume will simply be the addition of all BUCells' volume as all regions should be depleted.

        Parameters
        ----------
        total_vol: float
            Total volume of the system
        """
        self._total_vol = total_vol

    def get_tot_hm(self):
        """Returns the current total heavy metal mass contained in System (g)."""

        tot_hm = 0
        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            tot_hm += bucell.get_hm()

        return tot_hm

    def get_tot_ihm(self):
        """Returns the total initial heavy metal (IHM) mass contained in System (g)."""

        tot_ihm = 0
        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            tot_ihm += bucell.ihm

        return tot_ihm
    
    def set_default_decay_for_all(self):
        """Sets the decay library to the default one (ENDF/B-VIII.0) for all BUCells.

        This function will add all nuclides for which there are decay data in the library but which are not yet in the onix.Passlist object of the BUCells. This means that the depletion system will be significantly enlarged.
        """

        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.set_default_decay_lib()

    def set_default_decay_for_all_no_add(self):
        """Sets the decay library to the default one (ENDF/B-VIII.0) for all BUCells.

        This function will only take nuclear decay data for nuclides that are already in the onix.Passlist object of the BUCells. This enables to run simulation where the depletion system is not enlarged when adding decay data.
        """
        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.set_default_decay_lib_no_add()

    def set_decay_for_all(self, decay_lib_path):
        """Sets the decay library for all BUCells.

        This function will add all nuclides for which there are decay data in the library but which are not yet in the onix.Passlist object of the BUCells. This means that the depletion system will be significantly enlarged.

        Parameters
        ----------
        decay_lib_path: str
            Path to the decay library provided by the user
        """
        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.set_decay_lib(decay_lib_path)

    def set_default_xs_for_all(self):
        """Sets the cross section library to the default one for all BUCells.

        The default cross section library in ONIX has been obtained by computing middle-burnup one-group cross section for a coupled simumation of a LWR fuel pin-cell. Therefore, this library is only adapted for LWR reactors.
        """

        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.set_default_xs_lib()

    def set_xs_for_all(self, xs_lib_path):
        """Sets the cross section library for all BUCells.

        Parameters
        ----------
        xs_lib_path: str
            Path to the cross section library provided by the user
        """

        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.set_xs_lib(xs_lib_path)

    def set_default_fy_for_all(self):
        """Sets the fission yield library to the default one (ENDF/B-VIII.0) for all BUCells.

        This function will add all fission products for which there are fission yield data in the library but which are not yet in the onix.Passlist object of the BUCells. This means that the depletion system will be significantly enlarged.
        """

        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.set_default_fy_lib()

    def set_default_fy_for_all_no_add(self):
        """Sets the fission yield library to the default one (ENDF/B-VIII.0) for all BUCells.

        This function will only take fission yield data for fission products that are already in the onix.Passlist object of the BUCells. This enables to run simulation where the depletion system is not enlarged when adding fission yield data.
        """

        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.set_default_fy_lib_no_add()

    def set_fy_for_all(self, fy_lib_path, complete):
        """Sets the fission yield library for all BUCells. Setting the *complete* parameter to True allows ONIX to complement the data of the provided library with additional data found in ENDF/B-VIII.0.

        This function will add all fission products for which there are fission yield data in the library but which are not yet in the onix.Passlist object of the BUCells. This means that the depletion system will be significantly enlarged.

        Parameters
        ----------
        decay_lib_path: str
            Path to the decay library provided by the user
        complete: bool
            Indicates to ONIX whether or not to complement the provided library with additional fission yields from ENDF:B-VIII.0 library.
        """

        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell.set_fy_lib(fy_lib_path, complete)

    @property
    def sequence(self):
        """Returns the sequence object associated with System"""

        return self._sequence

    @sequence.setter
    def sequence(self, sequence):

        self._sequence = sequence

    # In addition to set sequence for each cells
    # This method also converts the average power_dens sequence to tot_power sequence if constant_power mode
    # It seems the mode parameters is not used. Will have to investigate if its obsolete
    def set_sequence(self, sequence, mode = 'stand alone'):
        """Sets a sequence object to the System and distributes unique copies of the sequence object to each BUCells.

        Parameters
        ----------
        sequence: onix.Sequence
            ONIX sequence object defined by the user


        """
        # Create tot_pow, av_pow_dens sequence
        # set time seq or bu seq depending on input

        sequence._initial_system_conversion(self)

        self.sequence = sequence
        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            if mode == 'stand alone':
                bucell._set_sequence(sequence)
            elif mode == 'couple':
                bucell._set_sequence(sequence, mode = 'couple')

    @property
    def bu_sec_conv_factor(self):
        
        return self._bu_sec_conv_factor
    

    def _set_bu_sec_conv_factor(self):

        total_vol = self.total_vol
        # Since this get_tot_hm is called before burn, hm is ihm
        total_ihm = self.get_tot_hm()

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

    #   sequence = self.sequence
    #   steps_number = sequence.steps_number
    #   for s in range(steps_number):

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

    def _print_summary_activity(self):

        summary_path = self._output_summary_path
        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell._print_summary_activity(summary_path)

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
        txt += '{:<13.5E}'.format(kinf_seq[0].n)
        for s in range(steps_number):
            txt += '{:<13.5E}'.format(kinf_seq[s+1].n)

        txt += '\n'

        txt += '{:<12}'.format('UNCERTAINTY')
        txt += '{:<13.5E}'.format(kinf_seq[0].s)
        for s in range(steps_number):
            txt += '{:<13.5E}'.format(kinf_seq[s+1].s)

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

    def _copy_cell_folders_to_step_folder(self, s):

        bucell_dict = self.bucell_dict
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell._copy_cell_folders_to_step_folder(s)

    def _gen_output_summary_folder(self):

        name = 'output_summary'
        self._output_summary_path = utils.get_folder_path(name)
        utils.gen_folder(name)

    def burn(self):
        """Launches a standalone simulation."""

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
