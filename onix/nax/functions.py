import onix.utils as utils
import onix.data as d
import math as m
import numpy as np
import time
import matplotlib.pyplot as plt
import openmc
import itertools


 # This class stores the path of the simulation of the batch
 # and also extract corresponding cross sections
class Batch(object):

    """Batch objects store neutronics data from a coupled simulation modelling the irradiation of a fuel batch.

    Parameters
    ----------
    path: str
        Path to a simulation's directory 
    """

    def __init__(self, path):

        self._path = path
        self._path_output = path +'/output_summary'

    @property
    def path(self):
        """Returns the path of the simulation's directory.
        """
        return self._path

    @property
    def path_output(self):
        """Returns the path of the simulation's output summary directory.
        """
        return self._path_output

    @property
    def ng_xs_seq_dict(self):
        """Returns a dictionnary where keys are nuclides names and entries are one-group (n,gamma) cross sections evolution during the fuel batch exposition.
        """
        return self._ng_xs_seq_dict
        
    @property
    def tot_xs_seq_dict(self):
        """Returns a dictionnary where keys are nuclides names and entries are one-group total absorption cross sections evolution during the fuel batch exposition.
        """
        return self._tot_xs_seq_dict

    # Reads the n,gamma and total abs xs for nuclide list in cell
    def read_nuclide_list_xs(self, nuclide_list, cell):
        """Reads (n,gamma) and total absorption one-group cross sections evolution from coupled simulation's output for a given list of nuclides and for a given BUCell. This data is stored in cross section dictionnaries.

        Parameters
        ----------
        nuclide_list

        """
        ng_xs_seq_dict = {}
        tot_xs_seq_dict= {}
        for nuclide in nuclide_list:
            #print (self.path)
            ng_xs_seq_dict[nuclide] = utils.read_xs_seq(nuclide, '(n,gamma)', self.path, cell)
            tot_xs_seq_dict[nuclide] = utils.get_tot_xs(nuclide, self.path, cell)

        self._ng_xs_seq_dict =  ng_xs_seq_dict
        self._tot_xs_seq_dict =  tot_xs_seq_dict


# This function go over NAX nucl list and list the chain of isotopes connected through (n,gamma) reactions
# It builds a dict when each key is the z number and entry is list of tuble with first element nuclide name, second element natural abundance
# This function does not take into account if ngamma cross sections are non-zero (this might reduce chain)
# In addition, this function will remove single nuclides and nuclide at start of chains and which have zero abundance

# This function is probably obsolete. Might consider removing it.
def list_NAX_ng_chain():

    NAX_nucl_list = d.NAX_nucl_list

    z_chain_dict = {}

    for zamid in NAX_nucl_list:
        z = utils.get_zamid_z(zamid)
        name_old_format = utils.zamid_to_name(zamid)
        name_new_format = utils.onix_name_to_openmc_name(name_old_format)
        nat_abun = utils.get_zamid_natural_abundance(zamid)

        if z not in z_chain_dict:
            # If nuclide is first AND has zero abun, the code does not add it to the dict
            if nat_abun == 0.0:
                continue
            else:
                z_chain_dict[z] = [(name_new_format,nat_abun)]
        else:
            z_chain_dict[z].append((name_new_format,nat_abun))

    # Remove single nuclide
    for z in list(z_chain_dict):
        if len(z_chain_dict[z]) == 1:
            del z_chain_dict[z]

    return z_chain_dict





def review_all_ratio_candidates(NAX_cell, operation_history, path, ratio_uncertainty):

    """This function plots a series of graphs containing isotopic ratios evolution and their associated relative errors on fluence for a given operation history. It can be used to identify potential good fluence indicators.

    **Detail description**: the function builds (n,gamma) chains from all isotopes which are available in the one-group cross section library with onix.nax.list_NAX_ng_chain_from_output. Then, it depletes each chain with a simple and fast analytical integrator (based on Bateman solution) according to the operation history provided. The obtained isotopic evolutions are used to construct ratios evolution as well as relative error on fluence associated with these ratios. Both ratios evolution and relative errors on fluence are displayed on graphs for the user to assess good fluence indicators.

    **Note**: Since the relative error on fluence depends on the relative error on the measured ratios, the user must provide an assumed relative error that affect ratios measurements. These uncertainties affecting measurements depend on the isotopes densities, the ratio values and the technology used to measure the ratios. ONIX makes the approximation that the uncertainties affecting all ratios are the same.

    Parameters
    ----------
    NAX_cell: str
        Name of the BUCell for which fluence indicators are to be identified
    operation_history: List of Tuples
        Irradiation history of the reactor. This history can include successive fuel batches of different types. Each tuples represent the operation of a batch. The first element of the Tuple is the Batch object (onix.nax.Batch) corresponding to the batch type and the second element is the number of days for which this fuel batch has been run.
    path: str
        Path to the simulation's directory
    ratio_uncertainty: float
        Assumed uncertainty on the measured ratios provided by the user
    """
    chain_dict = list_NAX_ng_chain_from_output(path ,NAX_cell, 2)

    for key in chain_dict:
        # if key != 68:
        #   continue
        print ('\n\n\n chain dict\n\n\n', chain_dict, '\n\n\n')

        # It seems that list_NAX_ng_chain_from_output finds that there are only one chain per element
        # Hence taking only the first element of the list of chains
        # This is not satisfactory but will do for the moment
        chain = chain_dict[key][0]
        init_abun = get_nat_abun_list_from_chain(chain)
        nuclide_list = get_nuclide_list_from_chain(chain)

        history_matrix_list = get_history_matrix_list(operation_history, chain, NAX_cell)
        step_break_indexes = locate_step_break(history_matrix_list)
        batch_break_indexes = locate_batch_break_end_start(history_matrix_list)
        combine_indexes = get_combine_indexes(step_break_indexes, batch_break_indexes)

        history_matrix = concatenate_history_matrix(history_matrix_list)
        history_fluence = concatenate_history_fluence(history_matrix_list)
        ratio_evolution = get_ratio_evolution(chain, history_matrix)

        ratio_derivative_dict = get_ratio_derivative_dict(ratio_evolution, history_fluence, batch_break_indexes, step_break_indexes)
        fluence_derivative_dict = get_fluence_derivative_dict(ratio_derivative_dict)
        history_mid_fluence = get_history_mid_fluence(history_fluence, step_break_indexes)

        plot_fluence_relative_error_with_ratio_history(ratio_evolution, fluence_derivative_dict, history_mid_fluence, combine_indexes, chain, history_matrix, history_fluence, ratio_uncertainty)
        # plot_ng_chain_densities_history(chain, history_matrix, history_fluence, 'no')
        # plot_ng_chain_ratio_history(chain, ratio_evolution, history_fluence)
        # plot_chain_fluence_relative_error_history(ratio_evolution, fluence_derivative_dict, history_mid_fluence, combine_indexes)
        # plot_ng_chain_ratio_derivative_history(ratio_evolution, ratio_derivative_dict, history_mid_fluence, combine_indexes)


def review_selected_ratio_candidates(NAX_cell, operation_history, path, selected_list, ratio_uncertainty, cut_off=None):

    """This function plots a series of graphs containing isotopic ratios evolution and their associated relative errors on fluence for a selected list of provided ratios  and for a given operation history. It can be used to identify potential good fluence indicators.

    **Detail description**: the function builds (n,gamma) chains from all isotopes which are available in the one-group cross section library with onix.nax.list_NAX_ng_chain_from_output. Then, it depletes chains which contain the selected ratios provided by the user with a simple and fast analytical integrator (based on Bateman solution) according to the operation history provided. The obtained isotopic evolutions are used to construct ratios evolution as well as relative error on fluence associated with these ratios. Both ratios evolution and relative errors on fluence are displayed on graphs for the user to assess good fluence indicators.

    **Note**: Since the relative error on fluence depends on the relative error on the measured ratios, the user must provide an assumed relative error that affect ratios measurements. These uncertainties affecting measurements depend on the isotopes densities, the ratio values and the technology used to measure the ratios. ONIX makes the approximation that the uncertainties affecting all ratios are the same.

    Parameters
    ----------
    NAX_cell: str
        Name of the BUCell for which fluence indicators are to be identified
    operation_history: List of Tuples
        Irradiation history of the reactor. This history can include successive fuel batches of different types. Each tuples represent the operation of a batch. The first element of the Tuple is the Batch object (onix.nax.Batch) corresponding to the batch type and the second element is the number of days for which this fuel batch has been run.
    path: str
        Path to the simulation's directory
    selected_list: List of str
        List of nuclide's names selected by the user
    ratio_uncertainty: float
        Assumed uncertainty on the measured ratios provided by the user
    """
    chain_dict = list_NAX_ng_chain_from_output(path ,NAX_cell, 2)
    selected_fluence_derivative_dict = {}
    selected_ratio_dict = {}
    z_to_name_dict = d.nuc_zz_dic

    for z in chain_dict:
        name = z_to_name_dict[z]

        skip = 'yes'
        for ratio in selected_list:
            if name in ratio:
                skip = 'no'

        if skip == 'yes':
            continue


        chain = chain_dict[z][0]
        init_abun = get_nat_abun_list_from_chain(chain)
        nuclide_list = get_nuclide_list_from_chain(chain)

        history_matrix_list = get_history_matrix_list(operation_history, chain, NAX_cell)
        step_break_indexes = locate_step_break(history_matrix_list)
        batch_break_indexes = locate_batch_break_end_start(history_matrix_list)
        batch_break_indexes2 = locate_batch_break(history_matrix_list)
        combine_indexes = get_combine_indexes(step_break_indexes, batch_break_indexes)

        history_matrix = concatenate_history_matrix(history_matrix_list)
        history_fluence = concatenate_history_fluence(history_matrix_list)
        ratio_evolution = get_ratio_evolution(chain, history_matrix)

        ratio_name_list = ratio_evolution[1]
        print (ratio_name_list)
        ratios = []
        for ratio in selected_list:
            if ratio in ratio_name_list:
                ratios.append(ratio)

        if ratios == []:
            continue


        for ratio in ratios:
            print (ratio)
            if ratio == 'Dy-163/Dy-161':
                print (ratio)
                print (ratio_evolution[0][ratio])
                print (history_fluence)
            if ratio == 'Fe-57/Fe-56':
                print (ratio)
                print (ratio_evolution[0][ratio])
                #print (history_fluence)
            # if ratio == 'Fe-57/e-56':
            #   print (ratio)
            #   print (ratio_evolution[0][ratio])
            #   print (history_fluence)
            if ratio == 'Dy-164/Dy-163':
                print (ratio)
                print (ratio_evolution[0][ratio])
                #print (history_fluence)


        ratio_derivative_dict = get_ratio_derivative_dict(ratio_evolution, history_fluence, batch_break_indexes, step_break_indexes)
        fluence_derivative_dict = get_fluence_derivative_dict(ratio_derivative_dict)
        history_mid_fluence = get_history_mid_fluence(history_fluence, step_break_indexes)

        for ratio in ratios:
            selected_fluence_derivative_dict[ratio] = fluence_derivative_dict[ratio]
            selected_ratio_dict[ratio] = ratio_evolution[0][ratio]

    #quit()

    plot_selected_ratio_history(selected_list, selected_ratio_dict, history_fluence, batch_break_indexes = batch_break_indexes2)
    plot_selected_ratio_history_sup1(selected_list, selected_ratio_dict, history_fluence, batch_break_indexes = batch_break_indexes2)

    plot_selected_fluence_relative_error_history(selected_list, selected_ratio_dict, selected_fluence_derivative_dict, history_mid_fluence, combine_indexes, ratio_uncertainty, cut_off = cut_off, batch_break_indexes = batch_break_indexes2)
        #plot_fluence_relative_error_with_ratio_history(ratio_evolution, fluence_derivative_dict, history_mid_fluence, combine_indexes, chain, history_matrix, history_fluence)
    #plot_ng_chain_densities_history(chain, history_matrix, history_fluence, 'no')
    #plot_selected_ratio_history(selected_list, selected_ratio_dict, history_fluence)
        # plot_selected_fluence_relative_error_history(ratio_evolution, fluence_derivative_dict, history_mid_fluence, combine_indexes)
        # plot_ng_chain_ratio_derivative_history(ratio_evolution, ratio_derivative_dict, history_mid_fluence, combine_indexes)


def plot_pu_prod(fuel_cell, NAX_cell, operation_history, path, scale_up_factor = None):

    """This function plots the plutonium production in a reactor against fluence measured in a specified NAX BUCell (typically, a region where fluence indicators are measured) and according to a provided operation history.

    **Detail description**: The function goes over the density output of all the batches' simulation that constitute the operation history and extracts the plutonium evolution. Likewise, it extracts the fluence evolution for each batches in the NAX BUCell (the BUCell where fluence indicators are measured). It then constructs a plutonium evolution history against fluence according to the operation history provided. Finally, the function scales up the plutonium production to the entire reactor using the provided scale_up_factor. Graphs of plutonium evolution are displayed where each batch refueling is indicated.

    Parameters
    ----------
    fuel_cell: str
        Name of the BUCell where plutonium is produced
    NAX_cell: str
        Name of the BUCell where fluence indicators are measured
    operation_history: List of Tuples
        Irradiation history of the reactor. This history can include successive fuel batches of different types. Each tuples represent the operation of a batch. The first element of the Tuple is the Batch object (onix.nax.Batch) corresponding to the batch type and the second element is the number of days for which this fuel batch has been run.
    path: str
        Path to the simulation's directory
    scale_up_factor: float
        Factor to be used to scale up the plutonium production from the simulation's system (a pin-cell or an assembly for instance) to the whole reactor.

    """
    pu_prod_history_matrix_list = get_pu_prod_history_matrix_list(operation_history, fuel_cell, NAX_cell)
    concatenated_pu_prod_history_matrix = concatenate_pu_prod_history_matrix(pu_prod_history_matrix_list)
    concatenated_pu_cum_prod_history_matrix = concatenate_pu_cum_prod_history_matrix(pu_prod_history_matrix_list)

    concatenated_history_fluence = concatenate_history_fluence_from_pu_prod_matrix_list(pu_prod_history_matrix_list)
    batch_break_indexes = locate_batch_break_from_pu_prod_matrix_list(pu_prod_history_matrix_list)

    mass_pu_prod_history_matrix = convert_density_to_mass(concatenated_pu_prod_history_matrix, path, fuel_cell, scale_up_factor = scale_up_factor)
    mass_pu_cum_prod_history_matrix = convert_density_to_mass(concatenated_pu_cum_prod_history_matrix, path, fuel_cell, scale_up_factor = scale_up_factor)

    plot_mass_pu_prod_against_fluence(mass_pu_prod_history_matrix, concatenated_history_fluence, batch_break_indexes = batch_break_indexes)
    plot_mass_pu_cum_prod_against_fluence(mass_pu_cum_prod_history_matrix, concatenated_history_fluence, batch_break_indexes = batch_break_indexes)


# The fact that it depends on step is troulbesome
# If the xs of a nuclide is non-zero initially, I dont think it is possible that it just goes to zero during the operation
def list_NAX_ng_chain_from_output(path, cell, step):

    """This function builds chains of nuclides connected via (n,gamma) reactions. These chains will be used to find the best isotopic ratios for fluence estimations.

    **Detail description**: The function builds dictionnaries where keys are atomic number (z) and entries are list of (n,gamma) chains (there can be multiple chains per z number). The elements in the list are tuples where the first element is the name of the nuclide and the second element is its natural abundance. The way chains are constructed is as follows:
        - The algorithm gathers all isotopes of same z which have data in the one-group cross section library of the simulation (nuclides with no cross section data are not considered)
        - Then it removes all isotopes which comes before the first naturally occuring isotopes (in the (n,gamma) chain order)
        - Then it removes all isotopes which are not produced by their (n,gamma) precursor in the (n,gamma) chain
        - Then it removes all isotopes which have half life smaller than 10,000 years (changes due to decay should be negligible compared to changes due to (n,gamma) reactions)

    **Note 1**: Isotopes from a (n,gamma) chain which are produced by other neutron-induced reactions at a rate similar or higher than the (n,gamma) reaction should also be removed. However, since the production rate depends not only on the one-group cross section but also on the density of the precursor, it requires knowing the actual density of the elements present in the material studied. The objective of the NAX module is to identify potential good fluence indicators among all possible isotopes without knowing the actual isotopics of the material studied. Therefore, this criteria is not implemented in the NAX module.

    **Note 2**: All isotopes tested are assumed to be present in the material according to their natural abundance. Absolute density of each element (families of isotopes) does not play a role here as ratios are taken between isotopes of the same elements.

    **Note 3**: In the current version, the user needs to specify a macrostep which will be used by the code to verify the values of the one-group (n,gamma) cross sections. It is reasonably assumed that if a nuclide's one-group cross section is non-zero at a specific macrostep, it is non-zero all along the simulation.

    Parameters
    ----------
    path: str
        Path to a simulation's directory
    cell: str
        Name of the BUCell to be studied
    step: int
        Macrostep at which tests on one-group (n,gamma) cross section are to be done

    """

    #This function will build dict (entries are z) of list (one list per chain, there can be multiple chain per z if ng xs is zero/decay too big/produced by other element somewhere and breaks chain)
# of tuples (first element is name and second elt is nat abun)
# It is build from passlist (which will have cross sections) and the NAX_z_list

# First the algorithm gather all isotopes of z existing from output xs lib (nuclides with no xs are not considered)
# Then it eliminates all those which are before the first naturally occuring isotope
# Then it eliminates those which have decay greater than threshold
# Then eliminates those which are not produced by ng of precursor
# Then eliminates those for which another reaction than ng from another element produced them at a rate that is bigger than 1E-3 of the ng production rate
# Verification that depends on cross section are repeated at each new step
# All nuclides used in ratio must either be directly linked through other valid isotopes or linked through a chain of valid isotopes
# If one isotopes is linked through other via an invalid isotope, then it should be removed

# Remark: both the criteria that nuclide needs to be produced by non-zero ng by precursor and criteria that
# ng rate must be 1000 bigger than next production rate can be done by looking reac rate.


    NAX_z_list = d.NAX_z_list
    z_chain_dict = {}

    # This method should detect which decay lib was used in the run
    # This is not impplemented yet and here ENDFVIII decay is set
    #decay_path ='/home/julien/Open-Burnup.dev/onix/data/other_libs/ENDFVIII/decay_lib'
    #decay_dict = d.read_decay_lib(decay_path)
    decay_dict = d.default_decay_lib_b

    xs_lib_nucl_list = utils.read_xs_nucl(path, cell)

    # This part list the nuclides with z that exist in the xs lib and create lists of chain starting with the first non-zero isotopes
    # It also remove nuclide with decay constant that are too big

    z_list = []
    for z in NAX_z_list:

        new_chain = 'yes'
        z_chains = []

        #Only takes nuclide which have cross section
        for name in xs_lib_nucl_list:
            zamid = utils.name_to_zamid(name)
            nucl_z = utils.get_name_z(name)

            if nucl_z == z:

                # Eliminate nuclide with half-life smaller than 10,000 years
                # This create a hole and thus a new chain
                if zamid in decay_dict:
                    if 'total decay' in decay_dict[zamid]:
                        if decay_dict[zamid]['total decay']>2.2e-12:
                            new_chain = 'yes'
                            continue

                # Eliminate isomeric state, it creates non linear chains
                if utils.get_zamid_s(zamid) == 1:
                    continue

                nucl_nat_abun = utils.get_name_natural_abundance(name)
                #if z not in z_chain_dict:
                #if z not in z_list:
                if new_chain == 'yes':

                    # This guy is candidate to be first nuclide in chain

                    # Here we test if its ng cross section is non-zero for step
                    #path_to_xs = path_to_output+'/{}_xs_lib'.format(cell)
                    nucl_ng_seq = utils.read_xs_seq(name,'(n,gamma)',path, cell)
                    step_ng = nucl_ng_seq[step]

                    # If this isotope has a zero ng xs, it can't be first nuclide in chain
                    if step_ng == 0.0:
                        continue

                    # If this isotope has zero natural abundance, it can't be the first nuclide in chain
                    elif nucl_nat_abun == 0.0:
                        continue

                    else:
                        #z_chain = [(name, nucl_nat_abun)]
                        z_chains.append([(name, nucl_nat_abun)])
                        z_list.append(z)
                        new_chain = 'no'
                        #z_chain_dict[z] = [(name, nucl_nat_abun)]
                elif new_chain == 'no':

                    # This guy is candidate to be in the chain somewhere (but not as first)

                    # There is an important subtlety here.
                    # You want to remove a nuclide from the chain if it is produced significantly by another process than ng
                    # This depends obviously on all the precursors density in the system (if the wrong precursors is way larger than the ng precusros for example)
                    # However, this would only be doable when you know the initial density of your material or when you do an actual NAX work on a material you know the impurities in
                    # WWhat we would like to do is comparing multiple ratio usefulness without any information on the density situation of the material (when you want to 
                    # compare all nuclides, what density should you put in the material? The default 1E-22? Their natural abundance? A real case application will never 
                    # be like that)
                    # I could compare cross section instead of reaction rate. But when undertaking theoritical evaluation of ratio usefulness, I dont think looking at xs relative magnitude is relevant
                    # Beause in the practical case, the magnitude of the rate might not be the same as xs
                    # I thus propose not to discard nuclide on the basis on the relative cross section of their precursor

                    # It is at least necessary to verify that ng precursor has non-zero cross section for ng
            
                    # prod_rank = utils.read_nuclide_reac_rank(name, step+1, path_to_output+'/cell_2_reacs_rank')[1]
                    ng_precursor_zamid = utils.find_zamid_precursor(zamid, '(n,gamma)')
                    ng_precursor_name = utils.zamid_to_name(ng_precursor_zamid)

                    # If ng precursor is not in xs lib, remove nuclide
                    # This create a hole and thus a new chain
                    if ng_precursor_name not in xs_lib_nucl_list:
                        new_chain = 'yes'
                        continue


                    precursor_xs_seq = utils.read_xs_seq(ng_precursor_name,'(n,gamma)',path,cell)
                    precursor_step_xs = precursor_xs_seq[step]
                    # If ng precursor has 0 ng xs, remove nuclide
                    # This create a hole and thus a new chain
                    if precursor_step_xs == 0.0:
                        new_chain = 'yes'
                        continue

                    else:
                        z_chains[-1].append((name, nucl_nat_abun))

        # This part remove isolated nuclide
        z_chains_no_single = [chain for chain in z_chains if len(chain)>1]
        # for chain in z_chains:
        #   if len(chain) == 1:
        #       print (chain)
        #       z_chains.remove(chain)

        # Certain z chains will be just empty ([]) because all chains were single chains
        if len(z_chains_no_single) != 0: 
            z_chain_dict[z] = z_chains_no_single

    # The chain dict obtain will contain holes
    # If a nuclide is isolated at the beginning or at the end of the chain, it can be removed
    # However, when removing a nuclide that is isolated at the beginning, the code must be careful that the
    # chain lead has to have abundance non zero
    #Otherwise, if the chain has hole in its middle, the code has to divide chain in severail subchains, again, making sure that
    # each chain lead has non-zero abundance

    return z_chain_dict

def get_chain_nuclide_index(chain, nuclide):

    index = 0
    for data in chain:
        if data[0] == nuclide:
            break
        index += 1

    return index

def get_chain_nuclide_nat_abun(chain, nuclide):

    for data in chain:
        if data[0] == nuclide:
            nat_abun = data[1]

    return nat_abun

def get_chain_nuclide_name(chain, index):

    data = chain[index]
    name = data[0]

    return name

def get_nat_abun_list_from_chain(chain):

    nat_abun_list = []
    for data in chain:
        nat_abun_list.append(data[1])

    return nat_abun_list

def get_nuclide_list_from_chain(chain):

    nuclide_list = []
    for data in chain:
        nuclide_list.append(data[0])

    return nuclide_list


def bateman_term(chain, init_abun, nuclide, step, path, cell, fluence, ng_xs_seq_dict, tot_xs_seq_dict):

    i = get_chain_nuclide_index(chain, nuclide)
    #nat_abun = get_chain_nuclide_nat_abun(chain, nuclide)

    res = 0
    for j in range(i+1):

        f = 1
        j_name = get_chain_nuclide_name(chain, j)

        j_tot_xs = tot_xs_seq_dict[j_name][step]

        for k in range(i+1):
            k_name = get_chain_nuclide_name(chain, k)
            if k < i:
                k_ng_xs = ng_xs_seq_dict[k_name][step]
                f *= k_ng_xs
            if k != j:
                k_tot_xs = tot_xs_seq_dict[k_name][step]
                f /= (k_tot_xs - j_tot_xs)
        # print(f, f*N0)
        # While in the coefficient, we can keep xs  in barn unit (coefficients are unit-less), we need to convert xs to cm unit in the exp
        res += f * init_abun * np.exp(-j_tot_xs*1E-24*fluence)
    return res


# This build the Bateman solution for a specific step (with specific xs values) for a specific nuclide
def bateman_step_solution(chain, init_abun_list, nuclide, step, path, cell, fluence, ng_xs_seq_dict, tot_xs_seq_dict):
    i = get_chain_nuclide_index(chain, nuclide)
    res = 0

    for j in range(i+1):

        res += bateman_term(chain[j:], init_abun_list[j], nuclide, step, path, cell, fluence, ng_xs_seq_dict, tot_xs_seq_dict)
    return res

# This funtion builds a matrix of bateman step solution where:
# Each line corresponds to one isotope ratio, each element of the line is a list of data generated by the Bateman step solutions for step s with 
# associated cross sections at different fluence point (separated by dphi)
def bateman_step_solution_matrix(chain, abun, batch, cell, dphi, EFPD):

    """This function calculates the density evolution of a (n,gamma) chain over the length of one fuel batch irradiation. The Bateman solution uses the one-group cross sections obtained from the neutronics simulation of the batch.

    Parameters
    ----------
    chain: dict
         Dictionnary where keys are z number and entries are (n,gamma) chains
    abun: List
         List of natural abundance of each isotope of the (n,gamma) chain
    batch: onix.nax.Batch
        Batch object from which one-group cross section evolution are extracted to construct the Bateman solution for the density evlution of the (n,gamma) chain
    cell: str
        Name of the BUCell where fluence indicators are to be measured
    dphi: float
        The fluence interval that separates each fluence point where Bateman solution is constructed

    """


    path = batch.path
    ng_xs_seq_dict = batch.ng_xs_seq_dict
    tot_xs_seq_dict = batch.tot_xs_seq_dict

    fluence_seq = utils.get_fluence_seq_until_time(path, cell, EFPD)


    #abun = get_nat_abun_list_from_chain(chain)
    nuclide_list = get_nuclide_list_from_chain(chain)

    mat = [[None for x in range(len(fluence_seq)-1)] for y in range(len(chain))]
    fluence_points_seq = []

    for step in range(len(fluence_seq)-1):

        bos_fluence = fluence_seq[step]
        eos_fluence = fluence_seq[step+1]
        #step_fluence_length = step_fluence_length_seq[step]
        step_fluence_length = eos_fluence - bos_fluence

        # List of cumulated fluence interval
        fluence_int = [dphi*(i+1) for i in range(int(step_fluence_length/dphi))] + [step_fluence_length]
        # List of fluence points for each substep
        #fluence_points = [dphi*(i+1)+bos_fluence for i in range(int(step_fluence_length/dphi))] + [eos_fluence]
        fluence_points = [dphi*(i+1)+bos_fluence for i in range(int(step_fluence_length/dphi))] + [eos_fluence]

        fluence_points_seq.append(fluence_points)

        for index in range(len(chain)):
            nuclide = get_chain_nuclide_name(chain, index)
            step_point_wise_density = []
            for fluence in fluence_int:
                step_point_wise_density.append(bateman_step_solution(chain, abun, nuclide, step, path, cell, fluence, ng_xs_seq_dict, tot_xs_seq_dict))
            #step_point_wise_density = [bateman_step_solution(chain, abun, nuclide, step, path, cell, fluence) for fluence in fluence_points]
            mat[index][step] = step_point_wise_density

        # Update abun (abun list is replaced by last time point abun)
        for index in range(len(chain)):
            abun[index] = mat[index][step][-1]

        #quit()

    return mat, fluence_points_seq

# The fluence in pu_prod_matrix should be the fluence of the NAX material, not of the fissile material
def pu_prod_matrix(batch, fuel_cell, NAX_cell, EFPD):

    path = batch.path_output

    fluence_subseq = utils.get_fluence_subseq_until_time(path, NAX_cell, EFPD)

    pu_subseq_mat = utils.get_pu_subseq_mat(path, fuel_cell, EFPD)

    return pu_subseq_mat, fluence_subseq


# This function builds a list of bateman_step_solution_matrix according to the operation_history set by the user
def get_history_matrix_list(operation_history, chain, cell):

    nuclide_list = get_nuclide_list_from_chain(chain)

    # Look for all different batches in operation history
    # and make them read cross section of nuclide list in cell
    batch_list = []
    for data in operation_history:
        batch = data[0]
        if batch not in batch_list:
            batch_list.append(batch)
            batch.read_nuclide_list_xs(nuclide_list, cell)

    abun = get_nat_abun_list_from_chain(chain)

    dphi = 1E19

    history_matrix_list = []
    for data in operation_history:
        batch = data[0]
        EFPD = data[1]
        step_solution_matrix = bateman_step_solution_matrix(chain, abun, batch, cell, dphi, EFPD)
        history_matrix_list.append(step_solution_matrix)

    return history_matrix_list

def get_pu_prod_history_matrix_list(operation_history, fuel_cell, NAX_cell):

    pu_prod_history_matrix_list = []

    for data in operation_history:
        batch = data[0]
        EFPD = data[1]
        pu_prod_history = pu_prod_matrix(batch, fuel_cell, NAX_cell, EFPD)
        pu_prod_history_matrix_list.append(pu_prod_history)

    return pu_prod_history_matrix_list

# This concatenation will NOT add previous batch density to current batch density
def concatenate_pu_prod_history_matrix(pu_prod_history_matrix_list):

    pu_prod_history_matrix = []

    matrix_size = len(pu_prod_history_matrix_list[0][0])

    for i in range(matrix_size):
        concatenated_dens = []
        for matrix in pu_prod_history_matrix_list:
            pu_prod_matrix = matrix[0]
            for dens in pu_prod_matrix[i]:
                concatenated_dens.append(dens)

        pu_prod_history_matrix.append(concatenated_dens)

    return pu_prod_history_matrix


# This concatenation will add previous batch density to current batch density
def concatenate_pu_cum_prod_history_matrix(pu_prod_history_matrix_list):

    pu_cum_prod_history_matrix = []

    matrix_size = len(pu_prod_history_matrix_list[0][0])

    for i in range(matrix_size):
        concatenated_dens = []
        eos_dens = 0
        for matrix in pu_prod_history_matrix_list:
            pu_prod_matrix = matrix[0]
            for dens in pu_prod_matrix[i]:
                concatenated_dens.append(dens + eos_dens)
            eos_dens = concatenated_dens[-1]

        pu_cum_prod_history_matrix.append(concatenated_dens)

    return pu_cum_prod_history_matrix

# Convert number density to kg. Can be scaled up to the whole system
def convert_density_to_mass(concatenated_pu_prod_history_matrix, path, fuel_cell, scale_up_factor = None):

    NA = d.NA
    path = path + '/output_summary'
    vol = utils.read_BUCell_vol(path, fuel_cell)

    mass_pu_prod_history_matrix = []

    Pu_isotopes_name = d.Pu_isotopes_name
    Pu_isotopes_zamid = d.Pu_isotopes_zamid
    default_atm_mass_lib = d.default_atm_mass_lib
    Pu_isotopes_mass = []
    for zamid in Pu_isotopes_zamid:
        Pu_isotopes_mass.append(default_atm_mass_lib[zamid])

    for i in range(len(Pu_isotopes_name)):

        atm_mass = Pu_isotopes_mass[i]
        dens_seq = concatenated_pu_prod_history_matrix[i]
        if scale_up_factor != None:
            mass_seq = [x*1E24*vol*atm_mass*scale_up_factor*1E-3/NA for x in dens_seq]
        else:
            mass_seq = [x*1E24*vol*atm_mass/NA for x in dens_seq]
        mass_pu_prod_history_matrix.append(mass_seq)


    return mass_pu_prod_history_matrix


def cumulate_pu_prod_history_matrix(concatenate_pu_prod_history_matrix):

    cumulate_pu_prod_history_matrix = []

    for dens_seq in concatenate_pu_prod_history_matrix:
        cum_dens_seq = []
        cum_dens_seq.append(dens_seq[0])
        for dens in dens_seq[1:]:
            cum_dens = dens + cum_dens_seq[-1]
            cum_dens_seq.append(cum_dens)

        cumulate_pu_prod_history_matrix.append(cum_dens_seq)

    return cumulate_pu_prod_history_matrix

# This method concatenate each step_solution_matrix in history_matrix_list together so that each nuclide evolution is represented by one array
def concatenate_history_matrix(history_matrix_list):

    matrix_size = len(history_matrix_list[0][0])

    history_matrix = []

    for nuclide_index in range(matrix_size):
        nuclide_concatenated_evolution = []
        for step_solution_matrix in history_matrix_list:
            for nuclide_solution in step_solution_matrix[0][nuclide_index]:
                for density in nuclide_solution:
                    nuclide_concatenated_evolution.append(density)

        history_matrix.append(nuclide_concatenated_evolution)

    return history_matrix


# This method concatenate each flux_points sequence in history_matrix_list together so that fluence evolution is rerepsented by one array
def concatenate_history_fluence(history_matrix_list):

    # The eos_fluence at the end of each step_solution is used to cumulate the fluence
    eos_fluence = 0
    fluence_concatenated_evolution = []
    for matrix in history_matrix_list:
        fluence_seq = matrix[1]
        for fluence_points in fluence_seq :
            for fluence in fluence_points:
                fluence_concatenated_evolution.append(fluence+eos_fluence)
        eos_fluence = fluence_concatenated_evolution[-1]


    return fluence_concatenated_evolution

def fraction_derivative(concatenate_history_matrix, history_fluence, chain):

    fraction_derivative_dict = {}
    for j in range(len(chain)):
        data = chain[j]
        name = data[0]
        fraction_derivative = []
        fraction_evolution = concatenate_history_matrix[j]
        index = 0
        for i in range(len(fraction_evolution)-1):
            # if index in step_break_indexes:
            #   index += 1
            #   continue
            #derivative = abs((fraction_evolution[i+1]-fraction_evolution[i])/(history_fluence[i+1]-history_fluence[i]))
            derivative = (fraction_evolution[i+1]-fraction_evolution[i])/(history_fluence[i+1]-history_fluence[i])
            # print ('ratio diff',ratio_evolution[i+1]-ratio_evolution[i])
            # print ('fluence diff',(history_fluence[i+1]-history_fluence[i]))
            fraction_derivative.append(derivative)
            index += 1

        fraction_derivative_dict[name] = fraction_derivative

    return fraction_derivative_dict


    # This method concatenate each flux_points sequence in history_matrix_list together so that fluence evolution is rerepsented by one array
def concatenate_history_fluence_from_pu_prod_matrix_list(history_matrix_list):

    # The eos_fluence at the end of each step_solution is used to cumulate the fluence
    eos_fluence = 0
    fluence_concatenated_evolution = []
    for matrix in history_matrix_list:
        fluence_seq = matrix[1]
        for fluence in fluence_seq :
            fluence_concatenated_evolution.append(fluence+eos_fluence)
        eos_fluence = fluence_concatenated_evolution[-1]


    return fluence_concatenated_evolution

# This method stores the indexes of each new steps
# It is used in derivative to locate the index where discontinuity appears and bypass them 
def locate_step_break(history_matrix_list):

    index_list = []
    index = -1
    for step_solution_matrix in history_matrix_list:
        for fluence_points in step_solution_matrix[1]:
            if index != -1:
                index_list.append(index)
            for fluence in fluence_points:
                index  += 1

    return index_list

def locate_batch_break(history_matrix_list):

    index_list = []
    index = -1
    for batch_solution_matrix in history_matrix_list:
        if index != -1:
            #index_list.append(index-1)
            index_list.append(index)
        for fluence_points in batch_solution_matrix[1]:
            for fluence in fluence_points:
                index  += 1

    # add the before last index 
    index_list = index_list + [index-1]

    return index_list

def locate_batch_break_end_start(history_matrix_list):

    index_list = []
    index = -1
    for batch_solution_matrix in history_matrix_list:
        if index != -1:
            index_list.append(index-1)
            index_list.append(index)
        for fluence_points in batch_solution_matrix[1]:
            for fluence in fluence_points:
                index  += 1

    # add the before last index 
    index_list = index_list + [index-1]

    return index_list

def locate_batch_break_from_pu_prod_matrix_list(history_matrix_list):

    index_list = []
    index = -1
    for batch_solution_matrix in history_matrix_list:
        if index != -1:
            #index_list.append(index-1)
            index_list.append(index)
        for fluence in batch_solution_matrix[1]:
            index  += 1

    # add the before last index 
    index_list = index_list + [index]

    return index_list

def get_ratio_evolution(chain, history_matrix):

    ratio_evolution_dict = {}
    ratio_name_list = []

    # Lower A as denominator
    for i in range(len(chain)-1):
        data = chain[i]
        nucl_name = data[0]
        nucl_nat_abun = data[1]
        if nucl_nat_abun == 0:
            # We turn the ratio upside down
            num_name = data[0]
            num_nat_abun = data[1]
            num_evolution = history_matrix[i]
            for j in range(i+1, len(chain)):
                data = chain[j]
                deno_name = data[0]
                deno_evolution = history_matrix[j]
                ratio_evolution = [x/y for x,y in zip(num_evolution, deno_evolution)]
                ratio_name = '{}/{}'.format(num_name, deno_name)
                ratio_evolution_dict[ratio_name] = ratio_evolution
                ratio_name_list.append(ratio_name)

        elif nucl_nat_abun != 0:
            deno_name = nucl_name
            deno_nat_abun = nucl_nat_abun
            deno_evolution = history_matrix[i]
            for j in range(i+1, len(chain)):
                data = chain[j]
                num_name = data[0]
                num_evolution = history_matrix[j]
                ratio_evolution = [x/y for x,y in zip(num_evolution, deno_evolution)]
                ratio_name = '{}/{}'.format(num_name, deno_name)
                ratio_evolution_dict[ratio_name] = ratio_evolution
                ratio_name_list.append(ratio_name)

    # #Higher A as denominator
    # for i in range(len(chain)-1):
    #   data = chain[i-1]
    #   nucl_name = data[0]
    #   nucl_nat_abun = data[1]
    #   if nucl_nat_abun == 0:
    #       # We turn the ratio upside down
    #       num_name = data[0]
    #       num_nat_abun = data[1]
    #       num_evolution = history_matrix[i-1]
    #       for j in range(i+1, len(chain)):
    #           data = chain[j-1]
    #           deno_name = data[0]
    #           deno_evolution = history_matrix[j-1]
    #           ratio_evolution = [x/y for x,y in zip(num_evolution, deno_evolution)]
    #           ratio_name = '{}/{}'.format(num_name, deno_name)
    #           if ratio_name in ratio_evolution_dict:
    #               continue
    #           ratio_evolution_dict[ratio_name] = ratio_evolution
    #           ratio_name_list.append(ratio_name)

    #   elif nucl_nat_abun != 0:
    #       deno_name = nucl_name
    #       deno_nat_abun = nucl_nat_abun
    #       deno_evolution = history_matrix[i-1]
    #       for j in range(i+1, len(chain)):
    #           data = chain[j-1]
    #           num_name = data[0]
    #           num_evolution = history_matrix[j-1]
    #           ratio_evolution = [x/y for x,y in zip(num_evolution, deno_evolution)]
    #           ratio_name = '{}/{}'.format(num_name, deno_name)
    #           if ratio_name in ratio_evolution_dict:
    #               continue
    #           ratio_evolution_dict[ratio_name] = ratio_evolution
    #           ratio_name_list.append(ratio_name)

    # onix always take the ratio of the higher A over lower A. However, this configuration fails sometime because
    # the lower A isotope drops to zero and the ratio presents a singulatiry
    # In this situation, onix invert the ratio
    #invert_ratio(ratio_evolution_dict, ratio_name_list)

    return ratio_evolution_dict, ratio_name_list

def invert_ratio(ratio_evolution_dict, ratio_name_list):

    for i in range(len(ratio_name_list)):
        ratio_name = ratio_name_list[i]
        ratio_evolution = ratio_evolution_dict[ratio_name]
        invert = 'no'
        for ratio in ratio_evolution:
            if ratio > 1000:
                invert = 'yes'
                break

        if invert == 'yes':
            inverted_ratio_evolution = [1/x for x in ratio_evolution]
            inverted_ratio_name = invert_ratio_name(ratio_name)
            ratio_name_list[i] = inverted_ratio_name
            ratio_evolution_dict[inverted_ratio_name] = inverted_ratio_evolution
            del ratio_evolution_dict[ratio_name]

def invert_ratio_name(ratio_name):

    separate_name = ratio_name.split('/')
    inverted_ratio_name = '{}/{}'.format(separate_name[1], separate_name[0])

    return inverted_ratio_name


# Derivative of ratio over fluence
def get_ratio_derivative_dict(ratio_evolution, history_fluence, batch_break_indexes, step_break_indexes):

    ratio_derivative_dict = {}
    ratio_evolution_dict = ratio_evolution[0]

    for ratio_name in ratio_evolution_dict:
        ratio_derivative = []
        smoother_batch_derivative_list = []
        ratio_evolution = ratio_evolution_dict[ratio_name]
        index = 0
        for i in range(len(history_fluence)-1):
            # if index in step_break_indexes:
            #   index += 1
            #   continue
            #derivative = abs((ratio_evolution[i+1]-ratio_evolution[i])/(history_fluence[i+1]-history_fluence[i]))
            mid_ratio = (ratio_evolution[i+1]+ratio_evolution[i])/2
            derivative = abs((ratio_evolution[i+1]-ratio_evolution[i])/(history_fluence[i+1]-history_fluence[i]))
            # print ('ratio diff',ratio_evolution[i+1]-ratio_evolution[i])
            # print ('fluence diff',(history_fluence[i+1]-history_fluence[i]))
            ratio_derivative.append(derivative)
            index += 1

        ratio_derivative_dict[ratio_name] = ratio_derivative
        #ratio_derivative_dict[ratio_name] = smoother_ratio_derivative
    return ratio_derivative_dict

# Derivative of fluence over ratio
def get_fluence_derivative_dict(ratio_derivative_dict):

    fluence_derivative_dict = {}
    for ratio_name in ratio_derivative_dict:
        ratio_derivative = ratio_derivative_dict[ratio_name]
        fluence_derivative = [1/x for x in ratio_derivative]
        fluence_derivative_dict[ratio_name] = fluence_derivative

    return fluence_derivative_dict



# since the derivative of the ratio is found at midpoints between two ratio points,
# we need to compute the equivalient mid point fluence
# mid fluence also bypass step break points
def get_history_mid_fluence(history_fluence, step_break_indexes):

    history_mid_fluence = []
    for i in range(len(history_fluence)-1):
        # if i in step_break_indexes:
        #   continue
        history_mid_fluence.append((history_fluence[i+1]+history_fluence[i])/2)

    return history_mid_fluence

# Probably useless
# Only ratio where the denon is initially non-zero is usable
# def get_usable_ratio(chain):

#   ratio_name_list = []
#   for i in range(len(chain)-1):
#       data = chain[i]
#       deno_name = data[0]
#       deno_nat_abun = data[1]
#       if deno_nat_abun == 0:
#           continue
#       else:
#           for j in range(i+1, len(chain)):
#               data = chain[j]
#               num_name = data[0]
#               ratio_name_list.append('{}/{}'.format(num_name, deno_name))

#   return ratio_name_list

def sample_ratio_evolution_on_fluence_grid(ratio_evolution, old_fluence_grid, new_fluence_grid):

    tabulated_ratio = openmc.data.Tabulated1D(ratio_evolution, old_fluence_grid)

    sampled_ratio = tabulated_ratio(new_fluence_grid)

    return sampled_ratio

def get_combine_indexes(index_list1, index_list2):

    combined_indexes = list(set(index_list1+index_list2))
    combined_indexes.sort()

    return combined_indexes
    
def sampled_index(batch_break_indexes):

    point_per_batch = 5

    previous_break_index = 0
    sampled_index_list = []
    for i in range(len(batch_break_indexes)):
        break_index = batch_break_indexes[i]
        int_length = break_index - previous_break_index
        step = m.floor(int_length/point_per_batch)
        batch_sampled_index = [j*step + previous_break_index for j in range(point_per_batch)] +[break_index-1]
        previous_break_index = break_index
        sampled_index_list.append(batch_sampled_index)

    sampled_index = list(itertools.chain.from_iterable(sampled_index_list))


    return sampled_index

def sample_data_with_sample_indexes(sampled_indexes, data):


    sampled_data = []
    for i in sampled_indexes:
        sampled_data.append(data[i])

    return sampled_data
    

def plot_ng_chain_ratio_derivative_history(ratio_evolution, ratio_derivative_dict, history_mid_fluence, sampled_index):

    ratio_name_list = ratio_evolution[1]

    sampled_history_fluence = sample_data_with_sample_indexes(sampled_index, history_mid_fluence)

    fig, ax = plt.subplots()
    count = 0
    linestyle = '-'
    for i in range(len(ratio_name_list)):
        if count > 8:
            linestyle = ':'
        ratio_name = ratio_name_list[i]
        ratio_derivative = ratio_derivative_dict[ratio_name]
        sampled_ratio_derivative = sample_data_with_sample_indexes(sampled_index, ratio_derivative)
        plt.plot(sampled_history_fluence, sampled_ratio_derivative, linestyle = linestyle, label = ratio_name)
        count += 1
    plt.ylabel('Ratio derivative', fontsize=16)
    plt.grid()
    ax.yaxis.get_offset_text().set_fontsize(16)
    plt.legend(prop={'size': 15})
    plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    plt.tick_params(labelsize=15)
    plt.xlabel('Fluence', fontsize=16)
    plt.show()

def plot_chain_fluence_relative_error_history(ratio_evolution, fluence_derivative_dict, history_mid_fluence, sampled_index):

    ratio_name_list = ratio_evolution[1]

    sampled_history_fluence = sample_data_with_sample_indexes(sampled_index, history_mid_fluence)

    fig, ax = plt.subplots()
    count = 0
    linestyle = '-'
    for i in range(len(ratio_name_list)):
        if count > 8:
            linestyle = ':'
        ratio_name = ratio_name_list[i]
        fluence_derivative = fluence_derivative_dict[ratio_name]
        sampled_fluence_derivative = sample_data_with_sample_indexes(sampled_index, fluence_derivative)
        relative_error_seq = [x/y for x, y in zip(sampled_fluence_derivative, sampled_history_fluence)]
        plt.plot(sampled_history_fluence, relative_error_seq, linestyle = linestyle, label = ratio_name)
        count +=1
    plt.ylabel('Fluence relative error', fontsize=16)
    plt.grid()
    ax.yaxis.get_offset_text().set_fontsize(16)
    plt.legend(prop={'size': 15})
    plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    plt.tick_params(labelsize=15)
    plt.xlabel('Fluence', fontsize=16)
    plt.yscale('log')
    plt.show()


def plot_selected_fluence_relative_error_history(selected_list, selected_ratio_dict, selected_fluence_derivative_dict, history_mid_fluence, sampled_index, ratio_uncertainty, cut_off=None, batch_break_indexes= None):

    plt.style.use('dark_background')

    ratio_name_list = selected_list
    sampled_history_fluence = sample_data_with_sample_indexes(sampled_index, history_mid_fluence)

    marker_list = ['.', 'x']
    linestyle_list = [':', '--']

    fig, ax = plt.subplots()
    count = 0
    linestyle = '-'
    for i in range(len(ratio_name_list)):
        if count > 8:
            linestyle = ':'
        if count > 16:
            linestyle = '--'
        ratio_name = ratio_name_list[i]
        fluence_derivative = selected_fluence_derivative_dict[ratio_name]

        ratio_evolution = selected_ratio_dict[ratio_name]
        mid_ratio_evolution = [(x+y)/2 for x,y in zip(ratio_evolution[:-1], ratio_evolution[1:])]

        cut = None
        if cut_off != None:
            if ratio_name in cut_off:
                cut = cut_off[ratio_name]

        sample = 'no'

        if sample == 'no':
        ########### NON SAMPLED

            relative_error_seq = [x*ratio_uncertainty*k/y for x,k, y in zip(fluence_derivative, mid_ratio_evolution, history_mid_fluence)]
            #relative_error_seq = [x*ratio_uncertainty/y for x, y in zip(fluence_derivative, history_mid_fluence)]
            smooth_relative_error_seq = utils.moving_average(relative_error_seq, 100)
            x_seq = [i for i in range(len(history_mid_fluence))]
            # Non sampled
            if cut != None:
                plt.plot(history_mid_fluence[:cut], relative_error_seq[:cut], linestyle = linestyle, label = ratio_name)
            else:
                plt.plot(history_mid_fluence, relative_error_seq, linestyle = linestyle, label = ratio_name)
                #plt.plot(x_seq, relative_error_seq, linestyle = linestyle, label = ratio_name)

            # Threeshold line
            line = [1E-2 for x in range(len(history_mid_fluence))]
            plt.plot(history_mid_fluence, line, 'r', linestyle = '--')


        if sample == 'yes':
        ########### NON SAMPLED

            ########## SAMPLED

            sampled_mid_ratio_evolution = sample_data_with_sample_indexes(sampled_index, mid_ratio_evolution)
            sampled_fluence_derivative = sample_data_with_sample_indexes(sampled_index, fluence_derivative)
            relative_error_seq = [x*ratio_uncertainty*k/y for x,k, y in zip(sampled_fluence_derivative,sampled_mid_ratio_evolution, sampled_history_fluence)]
            smooth_relative_error_seq = utils.moving_average(relative_error_seq, 100)
            #ax1.plot(sampled_history_fluence, relative_error_seq, linestyle = linestyle, label = ratio_name)   

            x_seq = [i for i in range(len(sampled_history_fluence))]
            if cut != None:
                plt.plot(sampled_history_fluence[:cut], relative_error_seq[:cut], linestyle = linestyle, label = ratio_name)
            else:
                #plt.plot(x_seq, relative_error_seq, linestyle = linestyle, label = ratio_name)
                plt.plot(sampled_history_fluence, relative_error_seq, linestyle = linestyle, label = ratio_name)

            # Threeshold line
            line = [1E-2 for x in range(len(sampled_history_fluence))]
            plt.plot(sampled_history_fluence, line, 'r', linestyle = '--')




        ### Old cut off
        # if cut_off != None:
        #   if ratio_name in cut_off:
        #       cut = cut_off[ratio_name]
        #       plt.plot(history_mid_fluence[:cut], relative_error_seq[:cut], linestyle = linestyle, label = ratio_name)
        #       #plt.plot(sampled_history_fluence[:cut], relative_error_seq[:cut], linestyle = linestyle, label = ratio_name)
        #   else:
        #       plt.plot(history_mid_fluence, relative_error_seq, linestyle = linestyle, label = ratio_name)
        # else:
        #   x_seq = [i for i in range(len(history_mid_fluence))]
        #   #plt.plot(history_mid_fluence, relative_error_seq, linestyle = linestyle, label = ratio_name)
        #   plt.plot(x_seq, relative_error_seq, linestyle = linestyle, label = ratio_name)
        #   #plt.plot(sampled_history_mid_fluence[:cut], relative_error_seq[:cut], linestyle = linestyle, label = ratio_name)
        
        count +=1

    #Fluence at each batch break point
    batch_break_fluence = [history_mid_fluence[i] for i in batch_break_indexes]
    #plt.axvline(x=0, linestyle = '--', color ='grey')
    for fluence in batch_break_fluence:
        plt.axvline(x=fluence, linestyle = '--', color ='grey')



    plt.ylabel('Fluence relative error', fontsize=16)
    ax.grid(color = 'dimgray')
    ax.yaxis.get_offset_text().set_fontsize(16)
    ax.xaxis.get_offset_text().set_fontsize(16)
    plt.legend(prop={'size': 12})
    plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    plt.tick_params(labelsize=15)
    plt.xlabel('Fluence [cm/cm$^{3}$]', fontsize=16)
    plt.yscale('log')


    # Put these freaking Shutdown dates on the top
    ax3 = ax.twiny()
    ax3.set_xlim(ax.get_xlim())
    ax3.set_xticks(batch_break_fluence)
    ax3.tick_params(labelsize=11)
    ax3.set_xticklabels(['Shutdown\nApril 1994', 'Shutdown\nApril 2005', 'Shutdown\nJuly 2007', 'Shutdown\nOctober 2015', 'Shutdown\nMarch 2018'])

    #plt.xscale('log')
    plt.show()

def plot_fluence_relative_error_with_ratio_history(ratio_evolution, fluence_derivative_dict, history_mid_fluence, sampled_index, chain, history_matrix, history_fluence, ratio_uncertainty):

    #ratio_evoluton_dict = ratio_evolution[0]
    ratio_name_list = ratio_evolution[1]
    ratio_evolution_dict = ratio_evolution[0]
    nuclide_list = get_nuclide_list_from_chain(chain)
    sampled_history_fluence = sample_data_with_sample_indexes(sampled_index, history_mid_fluence)

    fig, (ax1, ax2) = plt.subplots(2, sharex = True)
    count = 0
    linestyle = '-'
    for i in range(len(ratio_name_list)):
        if count > 8:
            linestyle = ':'
        ratio_name = ratio_name_list[i]
        ratio_evolution = ratio_evolution_dict[ratio_name]
        mid_ratio_evolution = [(x+y)/2 for x,y in zip(ratio_evolution[:-1], ratio_evolution[1:])]
        fluence_derivative = fluence_derivative_dict[ratio_name]

        sampled_mid_ratio_evolution = sample_data_with_sample_indexes(sampled_index, mid_ratio_evolution)
        sampled_fluence_derivative = sample_data_with_sample_indexes(sampled_index, fluence_derivative)
        relative_error_seq = [x*ratio_uncertainty*k/y for x,k, y in zip(sampled_fluence_derivative,sampled_mid_ratio_evolution, sampled_history_fluence)]
        ax1.plot(sampled_history_fluence, relative_error_seq, linestyle = linestyle, label = ratio_name)    

        #relative_error_seq = [x*ratio_uncertainty*k/y for x,k, y in zip(fluence_derivative, mid_ratio_evolution, history_fluence)]
        #ax1.plot(history_fluence[:-1], relative_error_seq, linestyle = linestyle, label = ratio_name)
        
        count +=1
    # Threeshold line
    line = [1E-2 for x in range(len(sampled_history_fluence))]
    ax1.plot(sampled_history_fluence, line, 'r', linestyle = '--')

    #line = [1E-2 for x in range(len(history_fluence))]
    #ax1.plot(sampled_history_fluence, line, 'r', linestyle = '--')
    #ax1.plot(history_fluence, line, 'r', linestyle = '--')

    ax1.set_ylabel('Fluence relative error', fontsize=16)
    ax1.grid()
    ax1.yaxis.get_offset_text().set_fontsize(16)
    ax1.legend(prop={'size': 6})
    ax1.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    ax1.tick_params(labelsize=15)
    ax1.set_yscale('log')
    #ax1.set_xscale('log')

    for i in range(len(nuclide_list)):
        ax2.plot(history_fluence, history_matrix[i], label = nuclide_list[i])
    ax2.set_xlabel('Fluence', fontsize=16)
    ax2.set_ylabel('Atomic fraction', fontsize=16)
    ax2.grid()
    ax2.yaxis.get_offset_text().set_fontsize(16)
    ax2.xaxis.get_offset_text().set_fontsize(16)
    ax2.legend(prop={'size': 15})
    ax2.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    ax2.tick_params(labelsize=15)

    plt.show()






def plot_chain_fraction_derivative_history(fraction_derivative_dict, history_mid_fluence):

    fig, ax = plt.subplots()
    for name in fraction_derivative_dict:
        fraction_derivative = fraction_derivative_dict[name]
        plt.plot(history_mid_fluence, fraction_derivative, label = name)
    plt.ylabel('Fraction derivative', fontsize=16)
    plt.grid()
    ax.yaxis.get_offset_text().set_fontsize(16)
    plt.legend(prop={'size': 15})
    plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    plt.tick_params(labelsize=15)
    plt.xlabel('Fluence', fontsize=16)
    plt.show()



def plot_ng_chain_ratio_history(chain, ratio_evolution, history_fluence):

    ratio_evoluton_dict = ratio_evolution[0]
    ratio_name_list = ratio_evolution[1]

    if len(ratio_name_list) == 1:
        plt.figure(1)
        ratio_name = ratio_name_list[0]
        ratio = ratio_evoluton_dict[ratio_name]
        plt.plot(history_fluence, ratio, label = ratio_name)
        plt.ylabel('Ratio')
        plt.grid()
        plt.legend()
        plt.xlabel('Fluence')
        plt.show()

    elif len(ratio_name_list) > 1:
        f, ax_list = plt.subplots(len(ratio_name_list), sharex=True)
        for i in range(len(ax_list)):
            ax = ax_list[i]
            ratio_name = ratio_name_list[i]
            ratio = ratio_evoluton_dict[ratio_name]
            ax.plot(history_fluence, ratio, label = ratio_name)
            ax.set_ylabel('Ratio')
            ax.grid()
            ax.legend()

        ax.set_xlabel('Fluence')
        plt.show()

def plot_selected_ratio_history(selected_list, selected_ratio_dict, history_fluence, batch_break_indexes=None):

    ratio_name_list = selected_list
    plt.style.use('dark_background')

    f, ax = plt.subplots()
    for i in range(len(ratio_name_list)):
        ratio_name = ratio_name_list[i]
        ratio = selected_ratio_dict[ratio_name]
        ax.plot(history_fluence, ratio, label = ratio_name)

    #Fluence at each batch break point
    batch_break_fluence = [history_fluence[i] for i in batch_break_indexes]
    #plt.axvline(x=0, linestyle = '--', color ='grey')
    for fluence in batch_break_fluence:
        print (fluence)
        plt.axvline(x=fluence, linestyle = '--', color ='grey')

    #quit()



    # plt.axvspan(batch_break_fluence[0], batch_break_fluence[2], alpha=0.2, color='grey')
    # plt.axvspan(batch_break_fluence[9], batch_break_fluence[12], alpha=0.2, color='grey')


    ax.set_ylabel('Ratio', fontsize=16)
    ax.set_xlabel('Fluence [cm/cm$^{3}$]', fontsize=16)
    ax.yaxis.get_offset_text().set_fontsize(16)
    ax.xaxis.get_offset_text().set_fontsize(16)
    plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    ax.grid(color = 'dimgray')
    plt.tick_params(labelsize=15)
    ax.set_ylim(bottom = 0, top=3.1)
    plt.legend(prop={'size': 12})

    # Put these freaking Shutdown dates on the top
    ax3 = ax.twiny()
    ax3.set_xlim(ax.get_xlim())
    ax3.set_xticks(batch_break_fluence)
    ax3.tick_params(labelsize=11)
    ax3.set_xticklabels(['Shutdown\nApril 1994', 'Shutdown\nApril 2005', 'Shutdown\nJuly 2007', 'Shutdown\nOctober 2015', 'Shutdown\nMarch 2018'])

    #ax3.set_ylim(bottom = 0, top=3.1)

    plt.show()


def plot_selected_ratio_history_sup1(selected_list, selected_ratio_dict, history_fluence, batch_break_indexes=None):

    ratio_name_list = selected_list
    #plt.style.use('dark_background')

    f, ax = plt.subplots()
    for i in range(len(ratio_name_list)):
        ratio_name = ratio_name_list[i]
        ratio = selected_ratio_dict[ratio_name]
        ratio_sup1 = convert_ratio_to_sup1(ratio)
        ax.plot(history_fluence, ratio_sup1, label = ratio_name)

    #Fluence at each batch break point
    batch_break_fluence = [history_fluence[i] for i in batch_break_indexes]
    #plt.axvline(x=0, linestyle = '--', color ='grey')
    for fluence in batch_break_fluence:
        print (fluence)
        plt.axvline(x=fluence, linestyle = '--', color ='grey')

    #quit()

    # plt.axvspan(batch_break_fluence[0], batch_break_fluence[2], alpha=0.2, color='grey')
    # plt.axvspan(batch_break_fluence[9], batch_break_fluence[12], alpha=0.2, color='grey')


    ax.set_ylabel('Ratio', fontsize=16)
    ax.set_xlabel('Fluence [cm/cm$^{3}$]', fontsize=16)
    ax.yaxis.get_offset_text().set_fontsize(16)
    ax.xaxis.get_offset_text().set_fontsize(16)
    plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    ax.grid(color = 'dimgray')
    plt.tick_params(labelsize=15)
    ax.set_ylim(bottom = 0, top=100)
    plt.legend(prop={'size': 12})

    # Put these freaking Shutdown dates on the top
    ax3 = ax.twiny()
    ax3.set_xlim(ax.get_xlim())
    ax3.set_xticks(batch_break_fluence)
    ax3.tick_params(labelsize=11)
    ax3.set_xticklabels(['Shutdown\nApril 1994', 'Shutdown\nApril 2005', 'Shutdown\nJuly 2007', 'Shutdown\nOctober 2015', 'Shutdown\nMarch 2018'])

    #ax3.set_ylim(bottom = 0, top=3.1)

    plt.show()

def convert_ratio_to_sup1(ratio):

    for i in range(len(ratio)):
        ratio_value = ratio[i]
        if ratio_value < 1:
            ratio[i] = 1/ratio_value

    return ratio

def plot_ng_chain_densities_history(chain, history_matrix, history_fluence, different_axes):

    nuclide_list = get_nuclide_list_from_chain(chain)

    #### This is the data from onix ####
    #### Just to compare evolution of Bateman for whole history to Salameche density evolution
    #### for one simulation

    cell = 'Clad'
    path = '/home/julien/Open-Burnup.dev/test_cluster/della/IVB/case5/cycle1/activation/test12_Ti-B-La_solution3/output_summary'

    fluence_seq = utils.get_fluence_seq(path, cell)
    salameche_densities = []
    path = path + '/{}_dens'.format(cell)
    for nuclide in nuclide_list:
        nuclide_density = utils.read_dens(nuclide, path)
        normalized_density = [x*1E8 for x in nuclide_density]
        salameche_densities.append(normalized_density)

    #########

    if different_axes == 'yes':

        f, ax_list = plt.subplots(len(nuclide_list), sharex=True)
        for i in range(len(ax_list)):
            ax = ax_list[i]
            ax.plot(history_fluence, history_matrix[i], label = nuclide_list[i])
            ax.plot(fluence_seq, salameche_densities[i], marker = 'x')
            ax.set_ylabel('Density [atm/cm3]')
            ax.grid()
            ax.legend()

        ax.set_xlabel('Fluence')

    if different_axes == 'no':

        fig, ax = plt.subplots()
        count = 0
        linestyle = '-'
        for i in range(len(nuclide_list)):
            if count > 8:
                linestyle = ':'
            plt.plot(history_fluence, history_matrix[i], linestyle = linestyle, label = nuclide_list[i])
        ax.yaxis.get_offset_text().set_fontsize(16)
        plt.ylabel('Atomic fraction', fontsize = 16)
        plt.xlabel('Fluence', fontsize = 16)
        plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
        plt.tick_params(labelsize=15)
        plt.legend(prop={'size': 15})
        plt.grid()
        count += 1

    plt.show()


def get_eos_abun_from_matrix(step_solution_matrix):

    mat = step_solution_matrix[0]
    for index in range(len(mat)):
        abun[index] = mat[index][-1][-1]

    return abun

def plot_ng_chain_densities(chain, fluence_points_seq, mat):

    nuclide_list = get_nuclide_list_from_chain(chain)

    chain_densities = []

    # concatenate densities for each nuclide
    for line in mat:
        joint_densities = []
        for densities_array in line:
            for density in densities_array:
                joint_densities.append(density)
        chain_densities.append(joint_densities)

    # concatenate fluence points
    joint_fluence = []
    for fluence_points in fluence_points_seq:
        for fluence in fluence_points:
            joint_fluence.append(fluence)

    f, ax_list = plt.subplots(len(nuclide_list), sharex=True)
    for i in range(len(ax_list)):
        ax = ax_list[i]
        ax.plot(joint_fluence, chain_densities[i], marker = '.', label = nuclide_list[i])
        ax.set_ylabel('Density [atm/cm3]')
        ax.grid()
        ax.legend()

    ax.set_xlabel('Fluence')
    plt.show()

def plot_compare_ng_chain_densities_with_salameche(chain, fluence_points_seq, mat, path, cell):

    nuclide_list = get_nuclide_list_from_chain(chain)

    chain_densities = []

    fluence_seq = utils.get_fluence_seq(path, cell)
    salameche_densities = []
    path = path + '/{}_dens'.format(cell)
    for nuclide in nuclide_list:
        nuclide_density = utils.read_dens(nuclide, path)
        normalized_density = [x*1E8 for x in nuclide_density]
        salameche_densities.append(normalized_density)

    # concatenate densities for each nuclide
    for line in mat:
        joint_densities = []
        for densities_array in line:
            for density in densities_array:
                joint_densities.append(density)
        chain_densities.append(joint_densities)

    # concatenate fluence points
    joint_fluence = []
    for fluence_points in fluence_points_seq:
        for fluence in fluence_points:
            joint_fluence.append(fluence)

    f, ax_list = plt.subplots(len(nuclide_list), sharex=True)
    for i in range(len(ax_list)):
        ax = ax_list[i]
        ax.plot(joint_fluence, chain_densities[i], marker = '.', label = nuclide_list[i])
        ax.plot(fluence_seq, salameche_densities[i], marker = 'x')
        ax.set_ylabel('Density [atm/cm3]')
        ax.grid()
        ax.legend()

    ax.set_xlabel('Fluence')
    plt.show()

def plot_cum_pu_prod_against_fluence(cumulate_pu_prod_history_matrix, concatenate_history_fluence):

    Pu_isotopes = d.Pu_isotopes

    plt.figure(1)
    for i in range(len(Pu_isotopes)):
        plt.plot(concatenate_history_fluence, cumulate_pu_prod_history_matrix[i], label = Pu_isotopes[i])
    plt.ylabel('Cumulative Density [atm/cm3]')
    plt.xlabel('Fluence')
    plt.grid()
    plt.legend()

    plt.show()

def plot_pu_prod_against_fluence(concatenate_pu_prod_history_matrix, concatenate_history_fluence):

    Pu_isotopes = d.Pu_isotopes_name

    plt.figure(1)
    for i in range(len(Pu_isotopes)):
        plt.plot(concatenate_history_fluence, concatenate_pu_prod_history_matrix[i], label = Pu_isotopes[i])
    plt.ylabel('Cumulative Density [atm/cm3]')
    plt.xlabel('Fluence')
    plt.grid()
    plt.legend()

    plt.show()


def plot_mass_pu_prod_against_fluence(mass_pu_prod_history_matrix, concatenate_history_fluence, batch_break_indexes = None):

    Pu_isotopes = d.Pu_isotopes_name

    # Tot mass of Pu
    tot_pu_seq = []

    for i in range(len(mass_pu_prod_history_matrix[0])):
        tot_pu = 0
        for j in range(len(mass_pu_prod_history_matrix)):
            tot_pu += mass_pu_prod_history_matrix[j][i]
        tot_pu_seq.append(tot_pu)

    f, ax1 = plt.subplots()

    # Fluence at each batch break point
    batch_break_fluence = [concatenate_history_fluence[i] for i in batch_break_indexes]
    for fluence in batch_break_fluence:
        plt.axvline(x=fluence, linestyle = '--', color ='k')

    # tot Pu mass at each batch break point
    batch_break_pu = [round(tot_pu_seq[i], 1) for i in batch_break_indexes]
    print (batch_break_pu)

    #ax2 = ax1.twinx()

    # for i in range(len(batch_break_indexes)):
    #   index = batch_break_indexes[i]
    #   # if i < 2:
    #   #   batch_fluence_seq = concatenate_history_fluence[:index +1]
    #   #   batch_tot_pu_line = [round(tot_pu_seq[index ], 1) for x in range(len(batch_fluence_seq))]
    #   #   ax1.plot(batch_fluence_seq,batch_tot_pu_line, linestyle = '--', color ='darkgrey')
    #   # else:
    #   #   batch_fluence_seq = concatenate_history_fluence[index:] + [2.37E22]
    #   #   batch_tot_pu_line = [round(tot_pu_seq[index ], 1) for x in range(len(batch_fluence_seq))]
    #   #   ax2.plot(batch_fluence_seq,batch_tot_pu_line, linestyle = '--', color ='darkgrey')
    #   batch_fluence_seq = concatenate_history_fluence[:index +1]
    #   batch_tot_pu_line = [round(tot_pu_seq[index ], 1) for x in range(len(batch_fluence_seq))]
    #   ax1.plot(batch_fluence_seq,batch_tot_pu_line, linestyle = '--', color ='darkgrey')

    previous_index = 0
    for i in range(len(batch_break_indexes)):
        index = batch_break_indexes[i]
        batch_fluence_seq = concatenate_history_fluence[previous_index:index+1]
        batch_tot_pu = tot_pu_seq[previous_index:index+1]
        # if i<2:
        #   ax1.plot(batch_fluence_seq, batch_tot_pu, linestyle = '-', color ='orange')
        # else:
        #   ax2.plot(batch_fluence_seq, batch_tot_pu, linestyle = '-', color ='orange')
        ax1.plot(batch_fluence_seq, batch_tot_pu, linestyle = '-', color ='orange')
        previous_index = index+1

    max_fluence = max(concatenate_history_fluence)

    ax1.set_ylabel('Mass [kg]', fontsize=16)
    #ax2.set_ylabel('Mass [kg]', fontsize=16)
    #ax1.set_yticks([0, batch_break_pu[0], batch_break_pu[1]])
    ax1.set_xlabel('Fluence [neutrons cm$^{-2}$]', fontsize=16)
    ax1.set_xlim(left = 0, right = max_fluence)
    # ax1.set_ylim(bottom = 0,top = 31)
    # ax2.set_ylim(bottom = 0, top = 31)
    #ax2.yaxis.tick_right()
    #ax2.set_yticks([0, batch_break_pu[2], batch_break_pu[3]])
    #ax.yaxis.get_offset_text().set_fontsize(15)
    ax1.xaxis.get_offset_text().set_fontsize(15)
    #plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    # ax1.grid()
    # ax2.grid()
    ax1.tick_params(labelsize=15)
    #ax2.tick_params(labelsize=15)
    #plt.legend(prop={'size': 12})

    # Put these freaking Shutdown dates on the top
    ax3 = ax1.twiny()
    ax3.set_xlim(left = 0, right= max_fluence)
    ax3.set_xticks(batch_break_fluence)
    ax3.tick_params(labelsize=11)
    #ax3.set_xticklabels(['Shutdown\nApril 1994', 'Shutdown\nApril 2005', 'Shutdown\nJuly 2007', 'Shutdown\nOctober 2015'])
    #ax3.set_xticklabels(['Shutdown\nApril 1994', 'Shutdown\nApril 2005', 'Shutdown\nJuly 2007', 'Shutdown\nOctober 2015'])

    plt.show()

def plot_mass_pu_cum_prod_against_fluence(mass_pu_prod_history_matrix, concatenate_history_fluence, batch_break_indexes = None):

    Pu_isotopes = d.Pu_isotopes_name

    # Tot mass of Pu
    tot_pu_seq = []

    for i in range(len(mass_pu_prod_history_matrix[0])):
        tot_pu = 0
        for j in range(len(mass_pu_prod_history_matrix)):
            tot_pu += mass_pu_prod_history_matrix[j][i]
        tot_pu_seq.append(tot_pu)

    f, ax1 = plt.subplots()
    # for i in range(len(Pu_isotopes)):
    #   ax.plot(concatenate_history_fluence, mass_pu_prod_history_matrix[i], label = Pu_isotopes[i])
    
    ax1.plot(concatenate_history_fluence, tot_pu_seq, 'orange')

    # Fluence at each batch break point
    batch_break_fluence = [concatenate_history_fluence[i] for i in batch_break_indexes]
    for fluence in batch_break_fluence:
        plt.axvline(x=fluence, linestyle = '--', color ='k')

    # tot Pu mass at each batch break point
    batch_break_cum_pu = [round(tot_pu_seq[i],1) for i in batch_break_indexes]


    # for i in batch_break_indexes:
    #   batch_fluence_seq = concatenate_history_fluence[:i]
    #   batch_tot_pu_line = [round(tot_pu_seq[i],1) for x in range(len(batch_fluence_seq))]
    #   ax1.plot(batch_fluence_seq,batch_tot_pu_line, linestyle = '--', color ='darkgrey')

    max_fluence = max(concatenate_history_fluence)

    ax1.set_ylabel('Mass [kg]', fontsize=16)
    ax1.set_xlabel('Fluence [neutrons cm$^{-2}$]', fontsize=16)
    ax1.set_xlim(left = 0, right= max_fluence)
    #ax1.set_ylim(bottom = 0,top = 31)
    #ax.yaxis.get_offset_text().set_fontsize(15)
    ax1.set_yticks([0] + batch_break_cum_pu)
    ax1.xaxis.get_offset_text().set_fontsize(15)
    #plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
    ax1.grid()
    plt.tick_params(labelsize=15)

    # Put these freaking Shutdown dates on the top
    ax3 = ax1.twiny()
    ax3.set_xlim(left = 0, right= max_fluence)
    #ax3.set_xlim(ax1.get_xlim())
    ax3.set_xticks(batch_break_fluence)
    ax3.tick_params(labelsize=11)
    #ax3.set_xticklabels(['Shutdown\nApril 1994', 'Shutdown\nApril 2005', 'Shutdown\nJuly 2007', 'Shutdown\nOctober 2015'])

    #plt.legend(prop={'size': 12})
    plt.show()


def plot_mass_pu_prod_against_ratio(ratio_name, sampled_ratio, mass_pu_prod):

    Pu_isotopes = d.Pu_isotopes_name

    plt.figure(1)
    for i in range(len(Pu_isotopes)):
        plt.plot(sampled_ratio, mass_pu_prod[i], label = Pu_isotopes[i])
    plt.ylabel('Mass [g]')
    plt.xlabel('{} evolution'.format(ratio_name))
    plt.grid()
    plt.legend()

    plt.show()

#N1972367430

            









        










