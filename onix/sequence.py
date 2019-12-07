import onix.utils as utils
from . import data
import numpy
import uncertainties

class Sequence(object):
    """The sequence object contains information about the burnup sequence for the system,
     including the macrostep vector and its associated units, which are set by the user

    The sequence object also contains information on the flux or power normalization
     """

    def __init__(self, id_number):

        self._id = id_number
        self._macrostep_vector = None
        self._macrostep_unit = None
        self._norma_vector = None
        self._norma_unit = None

        self._macrosteps_number = None
        self._microsteps_number = None


        # self._time_seq = [0]
        # self._time_subseq = [[0]]
        # self._bu_seq = [0]
        # self._bu_subseq = [[0]]

        self._current_time = None
        self._current_bu = None
        self._time_seq = None
        self._bu_seq = None
        self._current_time_subseq = None
        self._current_bu_subseq = None
        self._time_subseq_mat = None
        self._bu_subseq_mat = None

        # Probably unecessary
        self._current_time_intvl = None
        self._current_bu_intvl = None
        self.__time_intvl_seq = None
        self.__bu_intvl_seq = None
        self._current_time_intvl_subseq = None
        self._current_bu_intvl_subseq = None
        self._time_intvl_subseq_mat = None
        self._bu_intvl_subseq_mat = None

        # Average power density and Total power of the system (same for all cells)
        self._current_av_pow_dens = None
        self._current_tot_pow = None
        self._av_pow_dens_seq = None
        self._tot_pow_seq = None
        self._current_av_pow_dens_subseq = None
        self._current_tot_pow_subseq = None
        self._av_pow_dens_subseq_mat = None
        self._tot_pow_subseq_mat = None

        # Average flux of the system
        self._av_flux_seq = None

        # flux and power density for each cell
        self._current_MC_flux = None
        self._current_flux_spectrum = None
        self._current_flux = None
        self._current_pow_dens = None
        self._MC_flux_seq = None
        self._flux_spectrum = None
        self._flux_seq = None
        self._pow_dens_seq = None
        self._current_flux_subseq = None
        self._current_pow_dens_subseq = None
        self._MC_flux_subseq_mat = None
        self._flux_subseq_mat = None
        self._pow_dens_subseq_mat = None

        # System kinf

        self._current_kinf = None
        self._kinf_seq = None

        # Cell isotopic branching
        self._current_isomeric_branching_ratio = None
        self._isomeric_branching_ratio_seq = None

        # Isotopic change
        self._isotopic_change_dict = None

        # Density change
        self._density_change_dict = None

        # Temperature change
        self._temperature_change_dict = None

    # def _set_from_input(self, sequence_dict, passlist,  bu_sec_conv_factor):

    #     sequence = sequence_dict
    #     self._bu_sec_conv_factor = bu_sec_conv_factor

    #     self.set_sequence(sequence['unit_vector'], sequence['unit'])

    #     self.set_norma(sequence['norma_vector'], sequence['normalization'] )

    #     self.microstep_vector = sequence['microstep_vector']

    #     self.flux_approximation = sequence['flux_approximation']

    #     #self.set_flux_approximation(sequence['flux_approximation'])

    #     self._cell_conversion(passlist, bu_sec_conv_factor, mode)


    # I AM NOT SO SURE THIS FUNCTION IS STILL NEDDED
    # SYSTEM CONVERSION ALREADY CREATE TOT_POW, AV_POW_DENS AND TIME/BU SEQ
    # THE REMAINING (FLUX, POW_DENS, MC_FLUX, TIME/BU), WILL BE DYNAMICALLY
    # SET DURING THE BURN LOOP
    # THE ONLY THING NEEDED IS TO SET THE BU_SEC_CONV_FACTOR
    # def _cell_conversion(self, passlist, bu_sec_conv_factor, mode):

    #     #This method will convert the user input sequence into time, burnup, flux and power sequence if possible

    #     # First, we set the bu_sec_conv_factor to the sequence
    #     self._set_initial_bucell_bu()

    #     # if self._macrostep_unit in ['s', 'm', 'y', 'd']:

        #     # Create the time sequence
        #     # += because self._time_seq is already initialised
        #     self._time_seq += self._macrostep_vector
        #     for s in range(len(self._time_seq)-1):
        #         substeps = self._microstep_vector[s]
        #         #time_intvl += [self._time_seq[i+1]- self._time_seq[i]]
        #         time_intvl = self.get_time_intvl(s+1)
        #         _time_subintvl_val = time_intvl/substeps
        #         self._time_subseq_mat.append([_time_subintvl_val*(j+1) + self._time_seq[s] for j in range(substeps)])
        #         # self._time_subintvl.append([_time_subintvl_val for j in range(substeps)])

        #     if mode == 'stand alone':
        #     # If mode is stand alone and power is the normalization, you can directly convert time sequence to bu sequence
        #         if self._norma_unit == 'power':

        #             # Create the bu sequence
        #             self._bu_seq += [a*(b*bu_sec_conv_factor) for a, b in zip(self._time_seq[1:], self._norma_vector)]
        #             for i in range(len(self._bu_seq)-1):
        #                 substeps = self._microstep_vector[i]
        #                 self._bu_intvl += [self._bu_seq[i+1]- self._bu_seq[i]]
        #                 bu_substep_val = self._bu_intvl[i]/substeps
        #                 self._bu_subseq.append([bu_substep_val*(j+1) + self._bu_seq[i] for j in range(substeps)])
        #                 self._bu_subintvl.append([bu_substep_val for j in range(substeps)])

        # elif self._macrostep_unit == 'MWd/kg':

        #     self._bu_seq += self._macrostep_vector
        #     for i in range(len(self._bu_seq)-1):
        #         substeps = self._substeps[i]
        #         self._bu_intvl += [self._bu_seq[i+1]- self._bu_seq[i]]
        #         bu_substep_val = self._bu_intvl[i]/substeps
        #         self._bu_subseq.append([bu_substep_val*(j+1) + self._bu_seq[i] for j in range(substeps)])
        #         #self._bu_subintvl.append([bu_substep_val for j in range(substeps)])
        
        #     if mode == 'stand alone':
        #     # If mode is stand alone and power is the normalization, you can directly convert  bu sequence to time sequence
        #         if self._norma_unit == 'power':

        #             self._time_seq += [a/(b*bu_sec_conv_factor) for a, b in zip(self._bu_seq[1:], self._norma_vector)]
        #             for i in range(len(self._time_seq)-1):
        #                 substeps = self._substeps[i]
        #                 self._time_intvl += [self._time_seq[i+1]- self._time_seq[i]]
        #                 _time_subintvl_val = self._time_intvl[i]/substeps
        #                 self._time_subseq.append([_time_subintvl_val*(j+1) + self._time_seq[i] for j in range(substeps)])
        #                 self._time_subintvl.append([_time_subintvl_val for j in range(substeps)])

        # # If norma_value is flux then flux_seq is norma_vector
        # # I HAVE NO IDEA WHAT THIS MEANS IN COUPLE MODE SO FAR
        # if self._norma_unit == 'flux':
        #     # Flux
        #     self._flux_seq = self._norma_vector
        #     for s in range(self._steps_number):
        #         flux_value = self._flux_seq[s]
        #         flux_list = [flux_value]*self._microstep_vector[s]
        #         self._flux_subseq.append(flux_list)

        #     # Initial Power
        #     initial_flux = self._flux_seq[0]
        #     initial_power = self._update_pow_dens(passlist, initial_flux)
        #     self._power_seq.append(initial_power)


        # # WARNING THIS IS PROBANLY OBSOLETE
        # # If mode is stand alone and norma_value is power then power_seq is norma_vector
        # if mode == 'stand alone':
        #     if self._norma_unit == 'power':
        #         # Power
        #         self._power_seq = self._norma_vector
        #         for step in range(self._steps_number):
        #             power_value = self._power_seq[step]
        #             power_list = [power_value]*self._microstep_vector[step]
        #             self._power_subseq.append(power_list)

        #         #Initial flux
        #         initial_power = self._power_seq[0]
        #         initial_flux = self.update_flux(passlist, initial_power)
        #         self.flux = initial_flux
        #         self._flux_seq.append(initial_flux)



    # THIS SHOULD ONLY BE IMPLEMENTED IN STAND ALONE
    # COUPLED WILL USE FLUX FROM OPENMC
    # def _set_initial_flux(self, passlist):

    #     initial_power = self._power_seq[0]
    #     initial_flux = self.update_flux(passlist, initial_power)
    #     self.flux = initial_flux
    #     self._flux_seq.append(initial_flux)


    def _initial_system_conversion(self, system):

        total_vol = system.total_vol
        norma_vector = self.norma_vector
        norma_mode = self.norma_unit

        system._set_bu_sec_conv_factor()
        self._set_initial_time()
        self._set_initial_system_bu()

        # Directly create the whole time seq and subseq from macrostep_vector and substeps
        if self._macrostep_unit in ['s', 'm', 'y', 'd']:
            # Create the time sequence
            # += because self._time_seq is already initialised
            conv_to_sec = data.time_dic[self._macrostep_unit]
            self._time_seq += [x*conv_to_sec for x in self._macrostep_vector]
            for s in range(len(self._time_seq)-1):
                microsteps_number = self._microstep_vector[s]
                #time_intvl += [self._time_seq[i+1]- self._time_seq[i]]
                time_intvl = self.get_time_intvl(s+1)
                time_subintvl_val = time_intvl/microsteps_number
                self._time_subseq_mat.append([time_subintvl_val*(j+1) + self._time_seq[s] for j in range(microsteps_number)])

            if norma_mode == 'power':
                # Only when norma_unit is power
                # The initial average power density is used both for s=0 and s=1
                self._av_pow_dens_seq = [norma_vector[0]]+norma_vector
                cm3_to_liter = 1E-3 # av_pow is in kW/l
                self._tot_pow_seq = [x*total_vol*cm3_to_liter for x in self._av_pow_dens_seq]
                self._initial_system_time_bu_conversion(system)

            if norma_mode == 'flux':

                self._av_flux_seq = [norma_vector[0]]+norma_vector


                # Right now, I set burnup seq and sub_seq to zeros
                # so that the print function does not crash
                # Later improvements will need to develop module to dynmically update
                # bu. Look at "dynamic_system_time_bu_conversion"

                self._system_bu_seq = [0]*len(self._time_seq)
                for s in range(len(self._system_bu_seq)-1):
                    microsteps_number = self._microstep_vector[s]
                    self._system_bu_subseq_mat.append([0 for j in range(microsteps_number)])

        elif self._macrostep_unit == 'MWd/kg':
        # Directly create the whole system_bu seq and subseq from macrostep_vector and microsteps_number
            self._system_bu_seq += self._macrostep_vector
            for s in range(len(self._system_bu_seq)-1):
                microsteps_number = self._microstep_vector[s]
                system_bu_intvl = self.get_system_bu_intvl(s+1)
                system_bu_subintvl_val = system_bu_intvl/microsteps_number
                #self._system_bu_subseq.append([system_bu_substep_val*(j+1) + self._system_bu_seq[i] for j in range(microsteps_number)])
                self._system_bu_subseq_mat.append([system_bu_subintvl_val*(j+1) + self._system_bu_seq[s] for j in range(microsteps_number)])
                #self._system_bu_subintvl.append([system_bu_substep_val for j in range(microsteps_number)])

            if norma_mode == 'power':
                # Only when norma_unit is power
                # The initial average power density is used both for s=0 and s=1
                self._av_pow_dens_seq = [norma_vector[0]]+norma_vector
                cm3_to_liter = 1E-3 # av_pow is in kW/l
                self._tot_pow_seq = [x*total_vol*cm3_to_liter for x in self._av_pow_dens_seq]

                self._initial_system_bu_time_conversion(system)

            if norma_mode == 'flux':

                self._av_flux_seq = [norma_vector[0]]+norma_vector

                # Right now, I set burnup seq and sub_seq to zeros
                # so that the print function does not crash
                # Later improvements will need to develop module to dynmically update
                # bu. Look at "dynamic_system_time_bu_conversion"

                self._system_bu_seq = [0]*len(self._time_seq)
                for s in range(len(self._system_bu_seq)-1):
                    microsteps_number = self._microstep_vector[s]
                    self._system_bu_subseq_mat.append([0 for j in range(microsteps_number)])


        

    # Old version
    # def _time_bu_substep_conversion(self, powe_dens, s, ss):

    #     if self._macrostep_unit in ['s', 'm', 'y', 'd']:

    #         microsteps_number = self.microsteps_number(s)
    #         _time_subintvl = self._time_subintvl[s][0] # Since d is the user default unit, each time substep has the same length
    #         bu_substep_val = _time_subintvl*self._bu_sec_conv_factor*power
    #         if ss == 0:
    #             self._bu_subintvl.append([])
    #             self._bu_subseq.append([])
    #             self._bu_subseq[s+1].append(bu_substep_val + self._bu_seq[-1]) #If this is a new step, the last value is the last step value
    #         else:
    #             self._bu_subseq[s+1].append(bu_substep_val + self._bu_subseq[s+1][-1]) #If this the same  step, the last value is the last substep value

    #         self._bu_subintvl[s].append(bu_substep_val)

    #         if ss == microsteps_number - 1: # Time to  update bu_step and bu_seq
    #             bu_step = sum(self._bu_subintvl[s])
    #             self._bu_intvl += [bu_step]
    #             self._bu_seq += [self._bu_seq[-1] + bu_step]


    #     elif self._macrostep_unit == 'MWd/kg':

    #         microsteps_number = self._microsteps_number[s]
    #         bu_substep = self._bu_subintvl[s][0] # Since MWd/kg is the user default unit, each bu substep has the same length
    #         _time_subintvl_val = bu_substep/(self._bu_sec_conv_factor*power)
    #         if i == 0:
    #             self._time_subintvl.append([])
    #             self._time_subseq.append([])
    #             self._time_subseq[s+1].append(_time_subintvl_val + self._time_seq[-1]) #If this is a new step, the last value is th
    #         else:
    #             self._time_subseq[s+1].append(_time_subintvl_val + self._time_subseq[s][-1]) #If this the same  step, the last value is the last substep value

    #         self._time_subintvl[s].append(_time_subintvl_val)

    #         if i == self._microsteps_number - 1: # Time to  update time_step and time_seq
    #             time_step = sum(self._time_subintvl[s])
    #             self._time_intvl += [time_step]
    #             self._time_seq += [self._time_seq[-1] + time_step]

    # This convert time_seq and subseq to system bu seq and subseq using av_pow_dens
    def _initial_system_time_bu_conversion(self, system):

        av_pow_dens_seq = self._av_pow_dens_seq
        bu_sec_conv_factor = system.bu_sec_conv_factor
        macrosteps_number = self.macrosteps_number
        time_seq = self.time_seq
        time_subseq_mat = self.time_subseq_mat
        microstep_vector = self._microstep_vector

        system_bu = 0
        system_bu_seq = []
        for s in range(1, macrosteps_number+1):
            time_intvl = self.get_time_intvl(s)
            system_bu_intvl = time_intvl*av_pow_dens_seq[s]*bu_sec_conv_factor # the initial power density is is set to both av_pow_dens_seq[0] and av_pow_dens_seq[1]
            system_bu += system_bu_intvl
            system_bu_seq.append(system_bu)
        self._system_bu_seq += system_bu_seq

        # Incorrect way to calculate system bu from time
        #self._system_bu_seq += [x*y*bu_sec_conv_factor for x, y in zip(time_seq[1:],av_pow_dens_seq)]

        for s in range(len(self._system_bu_seq)-1):
            microsteps_number = self._microstep_vector[s]
            system_bu_intvl = self.get_system_bu_intvl(s+1)
            system_bu_subintvl_val = system_bu_intvl/microsteps_number
            self._system_bu_subseq_mat.append([system_bu_subintvl_val*(j+1) + self._system_bu_seq[s] for j in range(microsteps_number)])

        # Bucell can't convert bu to time. Time is converted from BU by system
        # elif self._macrostep_unit == 'MWd/kg':
        #     microsteps_number = self.microsteps_number(s-1)
        #     bu_subintvl = self.get_bu_subintvl(s,0) # Since MWd/kg is the user default unit, each bu substep has the same length
        #     time_subintvl_val = bu_subintvl/(self._bu_sec_conv_factor*pow_dens)
        #     time_substep_val = time_subintvl_val + self.current_time
        #     self._set_substep_time(time_substep_val, ss)

    # This convert bu_seq and subseq to time seq and subseq using av_pow_dens
    def _initial_system_bu_time_conversion(self, system):

        av_pow_dens_seq = self._av_pow_dens_seq
        bu_sec_conv_factor = system.bu_sec_conv_factor
        macrosteps_number = self.macrosteps_number
        system_bu_seq = self.system_bu_seq
        system_bu_subseq_mat = self.system_bu_subseq_mat
        microstep_vector = self._microstep_vector

        time = 0
        time_seq = []
        for s in range(1, macrosteps_number+1):
            system_bu_intvl = self.get_system_bu_intvl(s)
            time_intvl = system_bu_intvl/(av_pow_dens_seq[s]*bu_sec_conv_factor) # the initial power density is is set to both av_pow_dens_seq[0] and av_pow_dens_seq[1]
            time += time_intvl
            time_seq.append(time)
        self._time_seq += time_seq

        # Incorrect way to calculate time from system bu
        #self._time_seq += [x/(y*bu_sec_conv_factor) for x, y in zip(system_bu_seq[1:],av_pow_dens_seq)]
       
        for s in range(len(self._time_seq)-1):
            microsteps_number = self._microstep_vector[s]
            #time_intvl += [self._time_seq[i+1]- self._time_seq[i]]
            time_intvl = self.get_time_intvl(s+1)
            time_subintvl_val = time_intvl/microsteps_number
            self._time_subseq_mat.append([time_subintvl_val*(j+1) + self._time_seq[s] for j in range(microsteps_number)])

    # Looks obsolete
    # def dynamic_system_time_bu_conversion(self, system, s):

    #     # calculate total power
    #     tot_pow = 0
    #     for bucell_name in system.bucell_dict:
    #         bucell = bucell_dict[bucell_name]
    #         bucell_vol = bucell.vol
    #         bucell_pow_dens = bucell.current_pow_dens
    #         bucell_pow = bucell_vol*bucell_pow_dens
    #         tot_pow += bucell_pow

    #     system._append_tot_pow_seq(tot_pow)

    #     # total ihm
    #     tot_ihm = system.get_tot_ihm()

    #     #current time
    #     current_time = self.time









    # This convert time to bucell bu for each substep
    def _bucell_time_bu_substep_conversion(self, bucell, s, ss):

        pow_dens = self.current_pow_dens
        bu_sec_conv_factor = bucell.bu_sec_conv_factor
        # microsteps_number length is s-1
        microsteps_number = self.microsteps_number(s-1)
        time_subintvl = self.get_time_subintvl(s, 0) # Since d is the user default unit, each time substep has the same length
        print (time_subintvl)
        print (bu_sec_conv_factor)
        print (pow_dens)
        bucell_bu_subintvl_val = time_subintvl*bu_sec_conv_factor*pow_dens
        bucell_bu_substep_val = bucell_bu_subintvl_val + self.current_bucell_bu
        self._set_substep_bucell_bu(bucell_bu_substep_val, ss)

        # Bucell can't convert bu to time. Time is converted from BU by system
        # elif self._macrostep_unit == 'MWd/kg':
        #     microsteps_number = self.microsteps_number(s-1)
        #     bu_subintvl = self.get_bu_subintvl(s,0) # Since MWd/kg is the user default unit, each bu substep has the same length
        #     time_subintvl_val = bu_subintvl/(self._bu_sec_conv_factor*pow_dens)
        #     time_substep_val = time_subintvl_val + self.current_time
        #     self._set_substep_time(time_substep_val, ss)

    # Update the flux sequence
    def _set_flux(self, new_flux, s, ss):

        # Update current flux
        self.flux = new_flux



        # If i = 0 (i.e. if we begin a new step), we need to create a new sublist in flux_subseq
        if ss == 0:
            self._flux_subseq.append([])
            if s > 0: # the first value is already initiated by conversion of openmc
                self._flux_seq.append(new_flux)

        self._flux_subseq[-1].append(new_flux)

        # # In the case of power norma, the flux seq must be updated too
        # if self._norma_unit == 'power':
        #     microsteps_number = self.microsteps_number(s)
        #     if i == microsteps_number - 1: # Time to  update flux_seq
        #         self._flux_seq.append(flux)

    # Update the flux subsequence (done by )
    def _set_subflux(self, new_subflux, ss):

        # Update current flux
        self.flux = new_flux

        # Update current flux subseq
        self._append_current_flux_subseq(self, new_flux, ss)

        # # If i = 0 (i.e. if we begin a new step), we need to create a new sublist in flux_subseq
        # if ss == 0:
        #     self._flux_subseq.append([])
        #     if s > 0: # the first value is already initiated by conversion of openmc
        #         self._flux_seq.append(new_flux)

        # self._flux_subseq[-1].append(new_flux)

        # # In the case of power norma, the flux seq must be updated too
        # if self._norma_unit == 'power':
        #     microsteps_number = self.microsteps_number(s)
        #     if i == microsteps_number - 1: # Time to  update flux_seq
        #         self._flux_seq.append(flux)

    def _set_subpower(self, new_subpower, ss):


        # Update current power density
        self.pow_dens = new_pow_dens

        # Update current power density subsequence
        self._append_current_pow_dens_subseq(self, new_subpower, ss)
        # # If i = 0 (i.e. if we begin a new step), we need to create a new sublist in power_subseq
        # if ss == 0:
        #     self._power_subseq.append([])
        #     if s > 0: # the first value is already initiated in conversion
        #         self._power_seq.append(power)

        # self._power_subseq[s].append(power)

        # # In the case of flux norma, the power seq must be updated too
        # if self._norma_unit == 'flux':
        #     microsteps_number = self.microsteps_number(s)
        #     if i == microsteps_number - 1: # Time to  update power_seq
        #         self._power_seq.append(power)

    def _set_power(self, power, s, ss):

        # If i = 0 (i.e. if we begin a new step), we need to create a new sublist in power_subseq
        if ss == 0:
            self._power_subseq.append([])
            if s > 0: # the first value is already initiated in conversion
                self._power_seq.append(power)

        self._power_subseq[s].append(power)

        # # In the case of flux norma, the power seq must be updated too
        # if self._norma_unit == 'flux':
        #     microsteps_number = self.microsteps_number(s)
        #     if i == microsteps_number - 1: # Time to  update power_seq
        #         self._power_seq.append(power)


    # Method provided by MCODE
    def _get_FMF1(self, system, s):

        tot_pow = self.tot_pow_seq[s] # in kW
        bucell_dict = system.bucell_dict
        deno = 0
        fission_rate_list = []
        for bucell_id in bucell_dict:
            bucell = bucell_dict[bucell_id]
            bucell_sequence = bucell.sequence
            MC_flux = bucell_sequence.current_MC_flux # MC_flux unit of neutron.cm per source particle
            vol = bucell.vol # probably useless, MC_flux is already volume integrated
            passlist = bucell.passlist
            passport_list = passlist.passport_list
            conv_Mev_J = 1.60218e-13
            conv_J_kJ = 1E-3
            fission_energy = 0
            for nucl in passport_list:
                if nucl.get_FAM() == 'ACT' and nucl.current_xs != None:
                    fission_E = nucl.fission_E
                    # xs is obtained by dividing reaction rates (rate per source neutron) by flux (neutron.cm per source neutron)
                    # and density (cm-3), therefore, xs are just in cm2
                    fission_xs = nucl.current_xs['fission'][0]
                    dens = nucl.current_dens
                    fission_energy += dens*fission_xs*MC_flux*fission_E*conv_Mev_J*conv_J_kJ #kW


            deno += fission_energy

        return tot_pow/deno


##### steps and norma info ######


    @property
    def macrostep_vector(self):
        """Returns the macrostep vector"""
        return self._macrostep_vector

    @macrostep_vector.setter
    def macrostep_vector(self, macrostep_vector):
        """Sets the macrostep vector"""
        self._macrostep_vector = macrostep_vector

    @property
    def macrostep_unit(self):
        """Returns the units of the macrosteps

        Units can be time in second or burnup in MWd/kg
        """
        return self._macrostep_unit

    @macrostep_unit.setter
    def macrostep_unit(self, macrostep_unit):
        """Sets the units of the macrosteps

        Units can be time in second or burnup in MWd/kg"""
        self._macrostep_unit = macrostep_unit

    @property
    def macrosteps_number(self):
        """Returns the number of macrosteps"""
        return self._macrosteps_number

    @macrosteps_number.setter
    def macrosteps_number(self, macrosteps_number):
        """Sets the number of macrosteps"""
        self._macrosteps_number = macrosteps_number


    def set_norma(self, norma_vector, norma_unit):
        """Sets information on normalization

        In version 0.10, the couple mode can only accept normalizations of the flux
        where the total power of the system is held constant
        The standalone mode can only accept normalization of the total power
        where the flux of the system if held constant

        In version 0.10, in standalone mode, the flux set by the user is going to be
        used for all BUCells of the system. In other words, all BUCells will have the
        same neutron flux"""
        self._norma_vector = norma_vector
        self._norma_unit = norma_unit

    @property
    def norma_vector(self):
        """Returns the normalization vector"""
        return self._norma_vector

    @norma_vector.setter
    def norma_vector(self, norma_vector):
        """Sets the normalization vector

        One value for normalization should be set per macrostep
        Hence, the size of the vector norma_vector should be the same
        as the size of macrostep_vector
        """
        self._norma_vector = norma_vector

    @property
    def norma_unit(self):
        """Returns the normalization unit"""
        return self._norma_unit

    @norma_unit.setter
    def norma_unit(self, norma_unit):
        """Sets the normalization unit

        Power density: kW/l
        In version 0.10, the user should specify the power density as the total power
        divided by the total volume of the system. The total volume of the system should includes
        all regions of the system, even regions that are not BUCells (for example, it should include
        the water region around the fuel BUCell)

        Neutron flux: cm^-2 s^-1
        """

        self._norma_unit = norma_unit


    @property
    def isotopic_change_dict(self):
        return self._isotopic_change_dict

    def set_isotopic_change(self, cell, cell_isotopic_change, unit='number density'):
        """Manually changes the isotopic densities of user-specified nuclides in the BUCell for user-specified
        macrosteps

        Parameter:

        unit: specifies the unit of the new isotopic density
        'number density' (default) for density in atm per cm3
        'atom fraction' for density as atomic fraction
        
        IMPORTANT: If the method "set_density_change" is also used for the same BUCell, the isotopic change set by the user
        must be in atom fraction"""
        # If this is the first cell which is set isotopic change, create the dict
        if self.isotopic_change_dict == None:
            self._isotopic_change_dict = {}

        self._isotopic_change_dict[cell.name] = cell_isotopic_change
        self._isotopic_change_dict[cell.name]['unit'] = unit




    @property
    def density_change_dict(self):
        return self._density_change_dict

    def set_density_change(self, cell, cell_density_change):
        """Manually changes the total density (in atm per cm3) of the material of the BUCell for user-specified
        macrosteps
        
        IMPORTANT: This method trumps the other method "set_isotopic_change", i.e., the new isotopic densities set in
        set_isotopic_change will be renormalized by the new value from set_density_change"""

        # If this is the first cell which is set density change, create the dict
        if self.density_change_dict == None:
            self._density_change_dict = {}

        self._density_change_dict[cell.name] = cell_density_change

    @property
    def temperature_change_dict(self):
        return self._temperature_change_dict

    def set_temperature_change(self, cell, cell_temperature_change):
        """Manually changes the temperature (Kelvin) of the material of the BUCell for user-specified
        macrosteps"""
        # If this is the first cell which is set density change, create the dict
        if self.temperature_change_dict == None:
            self._temperature_change_dict = {}

        self._temperature_change_dict[cell.name] = cell_temperature_change    

    @property
    def flux_approximation(self):
        """Returns the method used for approximating the flux between two microsteps

        Warning: This method should not be used by the user in version 0.10"""
        return self._flux_approximation

    @flux_approximation.setter
    def flux_approximation(self, flux_approximation):
        """Sets the method used for approximating the flux between two microsteps

        Warning: This method should not be used by the user in version 0.10"""
        self._flux_approximation = flux_approximation

    @property
    def microstep_vector(self):
        """Returns the microstep vector
        
        Time: seconds

        Burnup: MWd/kg
        """
        return self._microstep_vector

    @microstep_vector.setter
    def microstep_vector(self, microstep_vector):
        """Sets the microstep vector
        
        Time: seconds

        Burnup: MWd/kg
        """
        self._microstep_vector = microstep_vector

    # bu_sec_conv_factor should be set to each cell or system, not sequence
    # def set_bu_sec_conv_factor(self, bu_sec_conv_factor):

    #     self._bu_sec_conv_factor = bu_sec_conv_factor


    # Looks obsolete
    # def gen_initial_step_folder(self):

    #     name = 'step_0'
    #     utils.gen_folder(name)

    def _gen_step_folder(self, s):

        # Folder numbering starts with 1 not 0
        name = 'step_{}'.format(s)
        utils.gen_folder(name)

    def microsteps_number(self, s):
        """Returns the number of microsteps for macrostep s"""
        return self._microstep_vector[s]

    def set_macrostep(self, macrostep_vector, macrostep_unit):
        """Sets the macrostep vector and the macrostep units"""
        self.macrostep_vector = macrostep_vector
        self.macrostep_unit = macrostep_unit
        self.macrosteps_number = len(macrostep_vector)


######## Set initial values


    def _set_initial_time(self):

        self.current_time = 0.0
        self.time_seq = [0.0]
        self._time_subseq_mat = [[0.0]]

    def _set_initial_bucell_bu(self):

        self.current_bucell_bu = 0.0
        self.bucell_bu_seq = [0.0]
        self._bucell_bu_subseq_mat = [[0.0]]

    def _set_initial_system_bu(self):

        self.current_system_bu = 0.0
        self.system_bu_seq = [0.0]
        self._system_bu_subseq_mat = [[0.0]]

    def _set_initial_MC_flux(self, new_MC_flux):

        self.current_MC_flux = new_MC_flux
        self.MC_flux_seq = [new_MC_flux]
      #  self._MC_flux_subseq_mat = [[new_MC_flux]]

    def _set_initial_kinf(self, new_kinf):

        self.current_kinf = new_kinf
        self.kinf_seq = [new_kinf]

    def _set_initial_flux(self, new_flux):

        self.current_flux = new_flux
        self.flux_seq = [new_flux]
        self._flux_subseq_mat = [[new_flux]]
   
    def _set_initial_pow_dens(self, new_pow_dens):

        self.current_pow_dens = new_pow_dens
        self.pow_dens_seq = [new_pow_dens]
        self._pow_dens_subseq_mat = [[new_pow_dens]]



######## Set Step values


    # After the completion of a burn step
    # Add time to time_seq and add the current time_subseq to subseq_mat
    def _set_macrostep_time(self):

        time = self.current_time
        self._append_time_seq(time)
       # self._append_time_subseq_mat()

    # After the completion of a burn step
    # Add bu to bu_seq and add the current bu_subseq to subseq_mat
    def _set_macrostep_bucell_bu(self):

        bucell_bu = self.current_bucell_bu
        # If this is the fist step
        if self.bucell_bu_seq == None:
            self.bucell_bu_seq = [bucell_bu]
        else:
            self._append_bucell_bu_seq(bucell_bu)   
       # self._append_flux_subseq_mat()

    # After an MC simulation
    # There are no substep MC_flux, therefore, set_macrostep_MC_flux should
    # be the one to update the current MC_flux
    def _set_macrostep_MC_flux(self, MC_flux):

        # If this is the fist step
        if self.current_MC_flux == None:
            self.MC_flux_seq = [MC_flux]
        else:
            self._append_MC_flux_seq(MC_flux)

        self.current_MC_flux = MC_flux
       # self._append_time_subseq_mat()

    # After an MC simulation
    # There are no substep MC_flux, therefore, set_macrostep_MC_flux should
    # be the one to update the current MC_flux
    def _set_macrostep_flux_spectrum(self, flux_spectrum):

        # If this is the fist step
        if self.current_flux_spectrum == None:
            self.flux_spectrum_seq = [flux_spectrum]
        else:
            self._append_flux_spectrum_seq(flux_spectrum)

        self.current_flux_spectrum = flux_spectrum
       # self._append_time_subseq_mat()

    # After an MC simulation
    # There are no substep kinf (as least calculated by MC), therefore, set_macrostep_kinf should
    # be the one to update the current kinf
    def _set_macrostep_kinf(self, kinf):

        # If kinf from previous cycles have already been stored
        if isinstance(self.current_kinf, uncertainties.UFloat):
            self._append_kinf_seq(kinf)
        # If this is the fist step
        elif self.current_kinf == None:
            self.kinf_seq = [kinf]

        self.current_kinf = kinf

    # After an MC simulation
    def _set_macrostep_flux(self, flux):

        # If this is the fist step
        if self.current_flux == None:
            self.flux_seq = [flux]
        else:
            self._append_flux_seq(flux)

        self.current_flux = flux
       # self._append_flux_subseq_mat()

    # After an MC simulation
    def _set_macrostep_pow_dens(self, pow_dens):

        # If this is the fist step
        if self.current_pow_dens == None:
            self.pow_dens_seq = [pow_dens]
        else:
            self._append_pow_dens_seq(pow_dens)
            
        self.current_pow_dens = pow_dens
       # self._append_pow_dens_subseq_mat()

    def _set_macrostep_isomeric_branching_ratio(self, isomeric_branching_ratio):

        # If this is the fist step
        if self.current_isomeric_branching_ratio == None:
            self.isomeric_branching_ratio_seq = [isomeric_branching_ratio]
        else:
            self._append_isomeric_branching_ratio_seq(isomeric_branching_ratio)

        self.current_isomeric_branching_ratio = isomeric_branching_ratio


######## Set Substep values


    # Add time to time_seq and add the current time_subseq to subseq_mat
    def _set_substep_time(self, time, ss):

        self.current_time = time
        self._append_time_subseq_mat(time, ss)
      #  self._append_current_time_subseq(time, ss)

    # Add bu to bu_subseq and add the current bu_subseq to subseq_mat
    def _set_substep_system_bu(self, bu, ss):

        self.current_system_bu = bu
        self._append_system_bu_subseq_mat(bu, ss)

    # Add bu to bu_seq and add the current bu_subseq to subseq_mat
    def _set_substep_bucell_bu(self, bu, ss):

        self.current_bucell_bu = bu
        self._append_bucell_bu_subseq_mat(bu, ss)
       # self._append_current_bu_subseq(bu, ss)

    # # After the completion of a burn step
    # # Add time to time_seq and add the current time_subseq to subseq_mat
    # def _set_substep_MC_flux(self, MC_flux, ss):

    #     self.current_MC_flux = MC_flux
    #     self._append_MC_flux_subseq_mat(MC_flux, ss)
    #   #  self._append_current_MC_flux_subseq(MC_flux, ss)

    # Add time to time_seq and add the current time_subseq to subseq_mat
    def _set_substep_flux(self, flux, s, ss):

        # If this is the first step, first substep
        if s == 1 and ss == 0:
            # First bracket is for flux0,0
            # append_flux_subseq_mat will put flux 1,0 (but they are the same)
            self._flux_subseq_mat = [[flux]]

        self._append_flux_subseq_mat(flux, ss)
        self.current_flux = flux
       # self._append_current_flux_subseq(flux, ss)

    # After the completion of a burn step
    # Add time to time_seq and add the current time_subseq to subseq_mat
    def _set_substep_pow_dens(self, pow_dens, s, ss):

        # If this is the first step
        if s == 1 and ss == 0:
            # First bracket is for flux0,0
            # append_flux_subseq_mat will put flux 1,0 (but they are the same)
            self._pow_dens_subseq_mat = [[pow_dens]]

        self._append_pow_dens_subseq_mat(pow_dens, ss)
        self.current_pow_dens = pow_dens
       # self._append_current_pow_dens_subseq(pow_dens, ss)


##### average power density info ######


    # Looks obsolete
    # @property
    # def av_pow_dens_seq(self):

    #     return self._av_pow_dens_seq


##### total power info ######

    # So fat I have decided that total power will only get updated every step

    # Power density that is currently set
    @property
    def current_tot_pow(self):
        """Returns the current total power of the system
        """
        if self._current_tot_pow is None:
            pass  # define exception for undefined variable
        return self._current_tot_pow

    @current_tot_pow.setter
    def current_tot_pow(self, new_tot_pow):
        """Sets the current total power of the system
        """
        self._current_tot_pow = new_tot_pow

    @property
    def tot_pow_seq(self):
        """Returns the sequence of total power of the system
        """
        return self._tot_pow_seq

    @tot_pow_seq.setter
    def tot_pow_seq(self, new_tot_pow_seq):
        """Sets the sequence of total power of the system
        """
        self._tot_pow_seq = new_tot_pow_seq

    def _append_tot_pow_seq(self, new_tot_pow):

        self.current_tot_pow = new_tot_pow
        self._tot_pow_seq.append(new_tot_pow)
    
    def tot_pow_point(self, s):
        """Returns the total power value for macrostep s
        """
        return self._tot_pow_seq[s]



##### time info ######

    @property
    def current_time(self):
        """Returns the time in second corresponding to the current macrostep or microstep
        """    
        if self._current_time is None:
            pass  # define exception for undefined variable
        return self._current_time

    @current_time.setter
    def current_time(self, new_time):
        """Sets the time in second corresponding to the current microstep (or macrostep)
        """ 
        self._current_time = new_time

    @property
    def time_seq(self):
        """Returns the time macrosequence
        """       
        return self._time_seq

    @time_seq.setter
    def time_seq(self, new_time_seq):
        """Sets the time macrosequence
        """ 
        self._time_seq = new_time_seq

    def _append_time_seq(self, new_time):

        self._time_seq.append(new_time)
    
    @property
    def current_time_subseq(self):
        """Returns the time microsequence in second corresponding to the current macrostep
        """ 
        return self._current_time_subseq

    def _append_current_time_subseq(self, new_time, ss):

        # Start of a new step
        if ss == 0:
            self._current_time_subseq = []
        self.current_time_subseq.append(new_time)

    @property
    def time_subseq_mat(self):
        """Returns a list of time microsequences in second where each microsequence correspond 
        to one macrostep
        """ 
        return self._time_subseq_mat

    @time_subseq_mat.setter
    def time_subseq_mat(self, time_subseq_mat):
        """Sets a list of time microsequences in second where each microsequence correspond 
        to one macrostep
        """
        self._time_subseq_mat = time_subseq_mat

    def _append_time_subseq_mat(self, time, ss):

        #self._time_subseq_mat.append(self.current_time_subseq)
        if ss == 0:
            self._time_subseq_mat.append([])
        self._time_subseq_mat[-1].append(time)

    def time_point(self, s):
        """Returns the time in second for macrostep s
        """
        return self._time_seq[s]

    def time_subpoint(self, s, ss):
        """Returns the time in second for macrostep s and microstep ss
        """
        return self._time_subseq_mat[s][ss]

    def get_time_intvl(self, s):
        """Gets the time interval in second between macrostep s nd macrostep s-1
        """
        if s == 0:
            raise Step_0 ('Step 0 has no interval')

        time_intvl = self.time_point(s) - self.time_point(s-1)

        return time_intvl

    def get_time_subintvl(self, s, ss):
        """Gets the time interval in second between macrostep s and macrostep s-1
        """
        if s == 0:
            raise Step_0 ('Step 0 has no subinterval')

        if ss == 0:
            time_subintvl = self.time_subpoint(s, ss) - self.time_point(s-1)

        else:
            time_subintvl = self.time_subpoint(s, ss) - self.time_subpoint(s, ss-1)

        return time_subintvl



##### system bu info ######

    @property
    def current_system_bu(self):
        """Returns burnup level of current macro or microstep
        """
        if self._current_system_bu is None:
            pass  # define exception for undefined variable
        return self._current_system_bu

    @current_system_bu.setter
    def current_system_bu(self, new_system_bu):
        """Sets burnup level of current macro or microstep
        """
        self._current_system_bu = new_system_bu

    @property
    def system_bu_seq(self):
        """Returns the burnup macrosequence
        """  
        return self._system_bu_seq

    @system_bu_seq.setter
    def system_bu_seq(self, new_system_bu_seq):
        """Sets the burnup macrosequence
        """ 
        self._system_bu_seq = new_system_bu_seq

    def _append_system_bu_seq(self, new_system_bu):

        self._system_bu_seq.append(new_system_bu)
    
    @property
    def current_system_bu_subseq(self):
        """Returns the burnup microsequence of current macrostep
        """ 
        return self._current_system_bu_subseq

    def _append_current_system_bu_subseq(self, new_system_bu, ss):

        # Start of a new step
        if ss == 0:
            self._current_system_bu_subseq = []
        self.current_system_bu_subseq.append(new_system_bu)

    @property
    def system_bu_subseq_mat(self):
        """Returns the list of burnup microsequences with one microsequences per macrostep
        """ 
        return self._system_bu_subseq_mat

    @system_bu_subseq_mat.setter
    def system_bu_subseq_mat(self, system_bu_subseq_mat):
        """Sets the list of burnup microsequences with one microsequences per macrostep
        """ 
        self._bu_subseq_mat = bu_subseq_mat 

    def _append_system_bu_subseq_mat(self, new_system_bu, ss):

        if ss == 0:
            self._system_bu_subseq_mat.append([])
        self._system_bu_subseq_mat[-1].append(new_system_bu)

    def system_bu_point(self, s):
        """Returns the burnup level at macrostep s"""
        return self._system_bu_seq[s]

    def system_bu_subpoint(self, s, ss):
        """Returns the burnup level at macrostep s and microstep ss"""
        return self._system_bu_subseq_mat[s][ss]

    def get_system_bu_intvl(self, s):
        """Gets the burnup interval between macrostep s-1 and macrostep s"""
        if s == 0:
            raise Step_0 ('Step 0 has no interval')

        system_bu_intvl = self.system_bu_point(s) - self.system_bu_point(s-1)

        return system_bu_intvl

    def get_system_bu_subintvl(self, s, ss):
        """Gets the burnup interval between microstep ss-1 and microstep ss for macrostep s"""
        if s == 0:
            raise Step_0 ('Step 0 has no subinterval')

        if ss == 0:
            system_bu_subintvl = self.system_bu_subpoint(s, ss) - self.system_bu_point(s-1)

        else:
            system_bu_subintvl = self.system_bu_subpoint(s, ss) - self.system_bu_subpoint(s, ss-1)

        return system_bu_subintvl



##### bucell bu info ######

    # Power density that is currently set
    @property
    def current_bucell_bu(self):
        """Returns burnup level of current macro or microstep
        for BUCell"""
        if self._current_bucell_bu is None:
            pass  # define exception for undefined variable
        return self._current_bucell_bu

    @current_bucell_bu.setter
    def current_bucell_bu(self, new_bucell_bu):
        """Returns burnup level of current macro or microstep
        for BUCell"""
        self._current_bucell_bu = new_bucell_bu

    @property
    def bucell_bu_seq(self):
        """Returns the burnup macrosequence for BUCell
        """ 
        return self._bucell_bu_seq

    @bucell_bu_seq.setter
    def bucell_bu_seq(self, new_bucell_bu_seq):
        """Sets the burnup macrosequence for BUCell
        """ 
        self._bucell_bu_seq = new_bucell_bu_seq

    def _append_bucell_bu_seq(self, new_bucell_bu):

        self._bucell_bu_seq.append(new_bucell_bu)
    
    @property
    def current_bucell_bu_subseq(self):
        """Returns the burnup microsequence of current macrostep
        for BUCell""" 
        return self._current_bucell_bu_subseq

    def _append_current_bucell_bu_subseq(self, new_bucell_bu, ss):

        # Start of a new step
        if ss == 0:
            self._current_bucell_bu_subseq = []
        self.current_bucell_bu_subseq.append(new_bucell_bu)

    @property
    def bucell_bu_subseq_mat(self):
        """Returns the list of burnup microsequences with one microsequences per macrostep
        for BUCell"""
        return self._bucell_bu_subseq_mat

    @bucell_bu_subseq_mat.setter
    def bucell_bu_subseq_mat(self, bucell_bu_subseq_mat):
        """Sets the list of burnup microsequences with one microsequences per macrostep
        for BUCell""" 
        self._bu_subseq_mat = bu_subseq_mat 

    def _append_bucell_bu_subseq_mat(self, new_bucell_bu, ss):

        if ss == 0:
            self._bucell_bu_subseq_mat.append([])
        self._bucell_bu_subseq_mat[-1].append(new_bucell_bu)

    def bucell_bu_point(self, s):
        """Returns the burnup level at macrostep s for BUCell"""
        return self._bucell_bu_seq[s]

    def bucell_bu_subpoint(self, s, ss):
        """Returns the burnup level at macrostep s and microstep ss for BUCell"""
        return self._bucell_bu_subseq_mat[s][ss]

    def get_bucell_bu_intvl(self, s):
        """Gets the burnup interval between macrostep s-1 and macrostep s for BUCell"""
        if s == 0:
            raise Step_0 ('Step 0 has no interval')

        bucell_bu_intvl = self.bucell_bu_point(s) - self.bucell_bu_point(s-1)

        return bucell_bu_intvl

    def get_bucell_bu_subintvl(self, s, ss):
        """Gets the burnup interval between microstep ss-1 and microstep ss for macrostep s for BUCell"""
        if s == 0:
            raise Step_0 ('Step 0 has no subinterval')

        if ss == 0:
            bucell_bu_subintvl = self.bucell_bu_subpoint(s, ss) - self.bucell_bu_point(s-1)

        else:
            bucell_bu_subintvl = self.bucell_bu_subpoint(s, ss) - self.bucell_bu_subpoint(s, ss-1)

        return bucell_bu_subintvl


##### flux info ######

    # Power density that is currently set
    @property
    def current_flux(self):
        """Returns neutron flux of current macro or microstep
        for BUCell"""
        if self._current_flux is None:
            pass  # define exception for undefined variable
        return self._current_flux

    @current_flux.setter
    def current_flux(self, new_flux):
        """Sets neutron flux of current macro or microstep
        for BUCell"""
        self._current_flux = new_flux

    @property
    def flux_seq(self):
        """Returns the neutron flux macrosequence for BUCell
        """ 
        return self._flux_seq

    @flux_seq.setter
    def flux_seq(self, new_flux_seq):
        """Sets the neutron flux macrosequence for BUCell
        """ 
        self._flux_seq = new_flux_seq

    def _append_flux_seq(self, new_flux):

        self._flux_seq.append(new_flux)
    
    @property
    def current_flux_subseq(self):
        """Returns the neutron flux microsequence of current macrostep
        for BUCell""" 
        return self._current_flux_subseq

    def _append_current_flux_subseq(self, new_flux, ss):

        # Start of a new step
        if ss == 0:
            self._current_flux_subseq = []
        self.current_flux_subseq.append(new_flux)

    @property
    def flux_subseq_mat(self):
        """Returns the list of neuron flux microsequences with one microsequences per macrostep
        for BUCell"""
        return self._flux_subseq_mat  

    def _append_flux_subseq_mat(self, flux, ss):

      #  self.flux_subseq_mat.append(self.current_flux_subseq)
        if ss == 0:
            self._flux_subseq_mat.append([])
        self._flux_subseq_mat[-1].append(flux)

    def flux_point(self, s):
        """Returns the neutron flux at macrostep s for BUCell"""
        return self._flux_seq[s]

    def flux_subpoint(self, s, i):
        """Returns the neutron flux at macrostep s and microstep ss for BUCell"""
        return self._flux_subseq_mat[s][i]



##### pow_dens info ######

    # Power density that is currently set
    @property
    def current_pow_dens(self):
        """Returns power density of current macro or microstep
        for BUCell"""
        if self._current_pow_dens is None:
            pass  # define exception for undefined variable
        return self._current_pow_dens

    @current_pow_dens.setter
    def current_pow_dens(self, new_pow_dens):
        """Sets power density of current macro or microstep
        for BUCell"""
        self._current_pow_dens = new_pow_dens

    @property
    def pow_dens_seq(self):
        """Returns the power density macrosequence for BUCell
        """ 
        return self._pow_dens_seq

    @pow_dens_seq.setter
    def pow_dens_seq(self, new_pow_dens_seq):
        """Sets the power density macrosequence for BUCell
        """
        self._pow_dens_seq = new_pow_dens_seq

    def _append_pow_dens_seq(self, new_pow_dens):

        self._pow_dens_seq.append(new_pow_dens)
    
    @property
    def current_pow_dens_subseq(self):
        """Returns the power density microsequence of current macrostep
        for BUCell"""
        return self._current_pow_dens_subseq

    def _append_current_pow_dens_subseq(self, new_pow_dens, ss):

        # Start of a new step
        if ss == 0:
            self._current_pow_dens_subseq = []
        self.current_pow_dens_subseq.append(new_pow_dens)

    @property
    def pow_dens_subseq_mat(self):
        """Returns the list of neuron flux microsequences with one microsequences per macrostep
        for BUCell"""
        return self._pow_dens_subseq_mat  

    def _append_pow_dens_subseq_mat(self, pow_dens, ss):

      #  self.pow_dens_subseq_mat.append(self.current_pow_dens_subseq)
        if ss == 0:
            self._pow_dens_subseq_mat.append([])
        self._pow_dens_subseq_mat[-1].append(pow_dens)

    def pow_dens_point(self, s):
        """Returns the power density at macrostep s for BUCell"""
        return self._pow_dens_seq[s]

    def pow_dens_subpoint(self, s, i):
        """Returns the power density at macrostep s and microstep ss for BUCell"""
        return self._pow_dens_subseq_mat[s][i]



##### MC_flux info ######

    # Power density that is currently set
    @property
    def current_MC_flux(self):
        """Returns normalized neutron flux of current macro or microstep
        for BUCell"""
        if self._current_MC_flux is None:
            pass  # define exception for undefined variable
        return self._current_MC_flux

    @current_MC_flux.setter
    def current_MC_flux(self, new_MC_flux):
        """Sets normalized neutron flux of current macro or microstep
        for BUCell"""
        self._current_MC_flux = new_MC_flux

    @property
    def MC_flux_seq(self):
        """Returns the normalized neutron flux macrosequence for BUCell
        """ 
        return self._MC_flux_seq

    @MC_flux_seq.setter
    def MC_flux_seq(self, new_MC_flux_seq):
        """Sets the normalized neutron flux macrosequence for BUCell
        """ 
        self._MC_flux_seq = new_MC_flux_seq

    def _append_MC_flux_seq(self, new_MC_flux):

        self._MC_flux_seq.append(new_MC_flux)

    def MC_flux_point(self, s):
        """Returns the normalized neutron flux at macrostep s for BUCell"""
        return self._MC_flux_seq[s]


##### flux_spectrum info ######

    @property
    def current_flux_spectrum(self):
        """Returns flux spectrum of current macrostep
        for BUCell"""
        if self._current_flux_spectrum is None:
            pass  # define exception for undefined variable
        return self._current_flux_spectrum

    @current_flux_spectrum.setter
    def current_flux_spectrum(self, new_flux_spectrum):
        """Sets flux spectrum of current macro or microstep
        for BUCell"""
        self._current_flux_spectrum = new_flux_spectrum

    @property
    def flux_spectrum_seq(self):
        """Returns the flux spectrum macrosequence for BUCell
        """ 
        return self._flux_spectrum_seq

    @flux_spectrum_seq.setter
    def flux_spectrum_seq(self, new_flux_spectrum_seq):
        """Sets the flux spectrum macrosequence for BUCell
        """
        self._flux_spectrum_seq = new_flux_spectrum_seq

    def _append_flux_spectrum_seq(self, new_flux_spectrum):

        self._flux_spectrum_seq.append(new_flux_spectrum)


##### kinf info ######

    # Power density that is currently set
    @property
    def current_kinf(self):
        """Returns kinf of current macrostep of system"""
        if self._current_kinf is None:
            pass  # define exception for undefined variable
        return self._current_kinf

    @current_kinf.setter
    def current_kinf(self, new_kinf):
        """Sets kinf of current macrostep of system"""
        self._current_kinf = new_kinf

    @property
    def kinf_seq(self):
        """Returns the kinf macrosequence
        """ 
        return self._kinf_seq

    @kinf_seq.setter
    def kinf_seq(self, new_kinf_seq):
        """Sets the kinf macrosequence
        """ 
        self._kinf_seq = new_kinf_seq

    def _append_kinf_seq(self, new_kinf):

        self._kinf_seq.append(new_kinf)
    
    def kinf_point(self, s):
        """Returns kinf at macrostep s"""
        return self._kinf_seq[s]

##### branching ratio info ######

    @property
    def current_isomeric_branching_ratio(self):
        """Returns isomeric branching ratios of current macro or microstep
        for BUCell"""
        if self._current_isomeric_branching_ratio is None:
            pass  # define exception for undefined variable
        return self._current_isomeric_branching_ratio

    @current_isomeric_branching_ratio.setter
    def current_isomeric_branching_ratio(self, new_isomeric_branching_ratio):
        """Sets isomeric branching ratios of current macro or microstep
        for BUCell"""
        self._current_isomeric_branching_ratio = new_isomeric_branching_ratio

    @property
    def isomeric_branching_ratio_seq(self):
        """Returns the isomeric branching ratios macrosequence for BUCell
        """ 
        return self._isomeric_branching_ratio_seq

    @isomeric_branching_ratio_seq.setter
    def isomeric_branching_ratio_seq(self, new_isomeric_branching_ratio_seq):
        """Sets the isomeric branching ratios macrosequence for BUCell
        """ 
        self._isomeric_branching_ratio_seq = new_isomeric_branching_ratio_seq

    def _append_isomeric_branching_ratio_seq(self, new_isomeric_branching_ratio_seq):

        self._isomeric_branching_ratio_seq.append(new_isomeric_branching_ratio_seq)

class Step_0(Exception):
    """Raise when the user try to access subinterval for the first step"""
    pass









