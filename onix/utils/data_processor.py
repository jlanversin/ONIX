from .functions import *
import ast
import re
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
import onix.data as d
#import pylab

def read_nuclide_reac_rank(nuclide, step, path):

    file = open(path, 'r')
    lines = file.readlines()

    zamid = name_to_zamid(nuclide)

    search = 'nuclide'
    for i in range(0,len(lines)):
        l = lines[i]
        if l == ' ==={}({}) ===\n'.format(nuclide, zamid):
            search = 'step'
        elif l == 'STEP {}\n'.format(step) and search == 'step':
            data = lines[i+2]
            break

    data = ast.literal_eval(data)


    destruction = {}

    dest_total = 0
    for tuples in data:
        if '-' not in tuples[0]:
            if tuples[1] == 0.0:
                continue
            dest_total += tuples[1]


    for tuples in data:
        if '-' not in tuples[0]:
            if tuples[1] == 0.0:
                continue
            destruction[tuples[0]] = [tuples[1], tuples[1]/dest_total] # val and percent


    production = {}

    prod_total = 0
    for tuples in data:
        if '-' in tuples[0]:
            reaction_val = tuples[1]
            if reaction_val == 0.0:
                continue    
            prod_total += reaction_val

    for tuples in data:
        if '-' in tuples[0]:
            parent = tuples[0].split()[0]
            reaction = tuples[0].split()[1]
            reaction_val = tuples[1]
            if reaction_val == 0.0:
                continue    
            production[parent] = [reaction, reaction_val, reaction_val/prod_total]

    return [destruction, production]

def plot_bucell_nuclide_network(nuclide, step, path, cell, threshold):
    
    """Plots a network diagram of the destruction and production reaction rates of a specified nuclide at a given macrostep for a given BUCell. Reaction rates are in :math:`barn^{-2}cm^{-1}s^{-1}`. Production channels are indicated with the name of the parent nuclide, destruction channels are indicated with the name of the reaction.

    Parameters
    ----------
    nuclide: str
        Name of the nuclide
    step: int
        Macrostep for which the diagram should be plotted
    path: str
        Path of the simulation directory
    cell: str
        Name of the BUCell for which the diagram should be plotted
    threshold: float
        Value under which reaction rates are not shown on diagram

    """
    path_to_rank = path +'/output_summary/cell_{}_reacs_rank'.format(cell)

    file = open(path_to_rank, 'r')
    lines = file.readlines()

    zamid = name_to_zamid(nuclide)

    print (nuclide, zamid)
    search = 'nuclide'
    for i in range(0,len(lines)):
        l = lines[i]
        if l == ' ==={}({}) ===\n'.format(nuclide, zamid):
            search = 'step'
        elif l == 'STEP {}\n'.format(step) and search == 'step':
            data = lines[i+2]
            break

    data = ast.literal_eval(data)


    destruction = {}

    dest_total = 0
    for tuples in data:
        if '-' not in tuples[0]:
            # if tuples[1] == 0.0:
            #   continue
            if tuples[1] < threshold:
                continue
            dest_total += tuples[1]


    for tuples in data:
        if '-' not in tuples[0]:
            # if tuples[1] == 0.0:
            #   continue
            if tuples[1] < threshold:
                continue
            destruction[tuples[0]] = [tuples[1], tuples[1]/dest_total] # val and percent


    production = {}

    prod_total = 0
    for tuples in data:
        if '-' in tuples[0]:
            reaction_val = tuples[1]
            # if reaction_val == 0.0:
            #   continue    
            if reaction_val < threshold:
                continue    
            prod_total += reaction_val

    for tuples in data:
        if '-' in tuples[0]:
            parent = tuples[0].split()[0]
            reaction = tuples[0].split()[1]
            reaction_val = tuples[1]
            # if reaction_val == 0.0:
            #   continue    
            if reaction_val < threshold:
                continue
            production[parent] = [reaction, reaction_val, reaction_val/prod_total]


    G = nx.MultiDiGraph()

    for parent in production:
        label = '{}\n{:.2E}[{:.2%}]'.format(production[parent][0], production[parent][1], production[parent][2])
        G.add_edge(parent, nuclide, label = label, length = 10)

    for edge in destruction:
        G.add_edge(nuclide, edge, label = '{:.2E}[{:.2%}]'.format(destruction[edge][0],destruction[edge][1] ), length = 10)

    # Get target nuclide index
    index = 0
    for node in G.nodes():
        if node == nuclide:
            break
        index += 1


    edges = G.edges()
    edge_labels = []

    edge_labels=dict([((u,v,),d['label'])
                     for u,v,d in G.edges(data=True)])
    node_color = []
    for node in G.nodes():
        if node == nuclide:
            node_color.append('mediumaquamarine')
        if node in production:
            node_color.append('lightskyblue')
        if node in destruction:
            node_color.append('darksalmon')

    edges = G.edges()
    # edge_weights = [G[u][v]['weight'] for u,v in edges]
    # red_edges = [('C','D'),('D','A')]
    # edge_colors = ['black' if not edge in red_edges else 'red' for edge in G.edges()]
    pos=nx.circular_layout(G, scale = 2)
    pos[nuclide] = np.array([0, 0])
    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_labels)
    nx.draw_networkx_labels(G,pos)
    #nx.draw_networkx_edges(G,pos, edges = edges)
    nx.draw_networkx_edges(G,pos)
    nx.draw(G,pos, node_size=4000, node_color = node_color, font_size = 6)

    plt.show()

def plot_nuclide_dens_from_passport(bucell, nuclide):
    """Plots the density evolution of a given nuclide (in :math:`atm barn^{-2}cm{-1}`) in a given BUCell against time (in days). This method requires that the corresonding Passport and BUCell objects are already defined in the Python environment.

    Parameters
    ----------
    bucell: onix.Cell
    nuclide: onix.Passport
    """
    sequence = bucell.sequence
    time_seq = sequence.time_seq.copy()
    dens_seq = nuclide.dens_seq 

    time_seq = [i/(24*3600) for i in time_seq]# time_seq is in seconds

    plt.figure(1)
    plt.plot(time_seq, dens_seq, color = 'orange', marker = 'o')
    plt.xlabel('Time [day]')
    plt.ylabel('Density [atm/cm3]')
    plt.grid()
    plt.title('{} {} density evolution'.format(bucell.name, nuclide.name))

    plt.show()

def plot_nuclide_dens(bucell, nuclide):

    path = os.getcwd() +'/{}_dens'.format(bucell)

    time_seq = read_time_seq(path)
    dens_seq = read_dens(nuclide, path)

    plt.figure(1)
    plt.plot(time_seq, dens_seq, color = 'orange', marker = 'o')
    plt.xlabel('Time [day]')
    plt.ylabel('Density [atm/cm3]')
    plt.grid()
    plt.title('{} {} density evolution'.format(bucell, nuclide))

    plt.show()

def plot_nuclide_dens_from_path(bucell, nuclide, path_to_simulation):

    """Plots the density evolution of a given nuclide (in :math:`atm barn^{-2}cm{-1}`) in a given BUCell against time (in days).

    Parameters
    ----------
    bucell: str
        Name of the BUCell
    nuclide: str
        Name of the nuclide
    path_to_simulation: str
        Path to simulation directory
    """

    path = path_to_simulation + '/output_summary/{}_dens'.format(bucell)

    time_seq = read_time_seq(path)
    dens_seq = read_dens(nuclide, path)

    plt.figure(1)
    plt.plot(time_seq, dens_seq, color = 'orange', marker = 'o')
    plt.xlabel('Time [day]')
    plt.ylabel('Density [atm/cm3]')
    plt.grid()
    plt.title('{} {} density evolution'.format(bucell, nuclide))

    plt.show()

def plot_nuclide_group_dens_from_path(bucell, nuclide_list, path_to_simulation):

    """Plots the total density evolution of a group of nuclides (in :math:`atm barn^{-2}cm{-1}`) in a given BUCell against time (in days).

    Parameters
    ----------
    bucell: str
        Name of the BUCell
    nuclide_list: str
        List of the names of the nuclides
    path_to_simulation: str
        Path to simulation directory
    """
    path = path_to_simulation + '/output_summary/{}_dens'.format(bucell)

    time_seq = read_time_seq(path)

    group_dens_seq = [0.0]*len(time_seq)
    for nuclide in nuclide_list:
        dens_seq = read_dens(nuclide, path)
        for i in range(len(group_dens_seq)):
            group_dens_seq[i] = group_dens_seq[i] + dens_seq[i]


    plt.figure(1)
    plt.plot(time_seq, group_dens_seq, color = 'orange', marker = 'o')
    plt.xlabel('Time [day]')
    plt.ylabel('Density [atm/cm3]')
    plt.grid()
    plt.title('{} {} density evolution'.format(bucell, nuclide_list))

    plt.show()


def plot_xs_time_evolution(bucell, nuclide, xs_name):

    path = os.getcwd() +'/{}_xs_lib'.format(bucell)

    time_seq = read_time_seq(path)
    xs_seq = read_xs(nuclide, xs_name, path)

    marker_list = ['x', '+', 'o', '*', '^', 's']
    color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

    plt.figure(1)
    plt.plot(time_seq, xs_seq, color = 'teal', marker = 'o')

    plt.xlabel('Time in Days')
    plt.ylabel('Eff. XS in barn')
    plt.legend()
    plt.grid()
    plt.title('{} {} effective cross sections evolution'.format(bucell, nuclide))

    plt.show()

def plot_xs_bu_evolution(bucell_list, nuclide, xs_name):

    index = 0
    for bucell in bucell_list:
        path = os.getcwd() +'/{}_xs_lib'.format(bucell)

        xs_seq = read_xs_seq(nuclide, xs_name, path)

        marker_list = ['x', '+', 'o', '*', '^', 's']
        color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

        plt.figure(1)
        plt.plot(bu_seq, xs_seq, color = color_list[index], marker = marker_list[index], label = bucell)

        index += 1

    bu_seq = read_bu_seq(path)

    plt.xlabel('BU in MWd/kg')
    plt.ylabel('Eff. XS in barn')
    plt.legend()
    plt.grid()
    plt.title('{} {} effective cross sections evolution'.format(nuclide, xs_name))

    plt.show()


def plot_xs_time_evolution_from_path(bucell, nuclide, xs_name, path):

    """Plots the one-group cross section evolution (in barn) of a given reaction for a given nuclide in a given BUCell against time (in days).

    Parameters
    ----------
    bucell: str
        Name of the BUCell
    nuclide_list: str
        Name of the nuclide
    xs_name: str
        Name of the reaction
    path_to_simulation: str
        Path to simulation directory
    """

    time_seq = read_time_seq()
    xs_seq = read_xs_seq(nuclide, xs_name, path)

    marker_list = ['x', '+', 'o', '*', '^', 's']
    color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

    plt.figure(1)
    plt.plot(time_seq, xs_seq, color = 'teal', marker = 'o')

    plt.xlabel('Time in days')
    plt.ylabel('Eff. XS in barn')
    plt.legend()
    plt.grid()
    plt.title('{} {} effective cross sections evolution'.format(bucell, nuclide))

    plt.show()

def plot_xs_bu_evolution_from_path(bucell_list, nuclide, xs_name, path):

    """Plots the one-group cross section evolution (in barn) of a given reaction for a given nuclide in a given BUCell against burnup (MWd/kg).

    Parameters
    ----------
    bucell: str
        Name of the BUCell
    nuclide_list: str
        Name of the nuclide
    xs_name: str
        Name of the reaction
    path_to_simulation: str
        Path to simulation directory
    """

    index = 0
    for bucell in bucell_list:
        path_xs = path +'/output_summary/{}_xs_lib'.format(bucell)

        xs_seq = read_xs_seq(nuclide, xs_name, path, bucell)

        marker_list = ['x', '+', 'o', '*', '^', 's']
        color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

        plt.figure(1)
        bu_seq = read_bu_seq(path_xs)
        plt.plot(bu_seq, xs_seq, color = color_list[index], marker = marker_list[index], label = bucell)

        index += 1


    plt.xlabel('BU in MWd/kg')
    plt.ylabel('Eff. XS in barn')
    plt.legend()
    plt.grid()
    plt.title('{} {} effective cross sections evolution'.format(nuclide, xs_name))

    plt.show()

# Compare xs evolution for the same nuclide, for various xs for different runs
def compare_xs_bu_evolution_from_path(bucell, nuclide, xs_name_list, path_list, name_list):

    """Plots multiple cross sections from multiple simulation directories in a series of subplots with the same burnup y-axis. Evolution of the same cross section for different simulations are plotted in the same subplots.

    Parameters
    ----------
    bucell: str
        Name of the BUCell
    nuclide: str
        Name of the nuclide
    xs_name_list: List
        List of the names of the cross sections
    path_list: List
        List of the path to the different simulations' directories
    name_list:
        List of the names of the different simulations
    """

    bu_seq_list = []
    xs_seq_list = []
    for path in path_list:

        bu_seq_list.append(read_bu_seq(path))
        xs_seq_list.append([])

        for xs in xs_name_list:
            xs_seq_list[-1].append(read_xs_seq(nuclide, xs, path))


    marker_list = ['x', '+', 'o', '*', '^', 's']
    color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']
    
    # plt.figure(1)

    # for i in range(len(path_list)):

    #   for j in range(len(xs_name_list)):

    #       plt.plot(bu_seq_list[i], xs_seq_list[i][j], label = '{} {}'.format(name_list[i], xs_name_list[j]))

    # Three subplots sharing both x/y axes
    f, ax_tuple = plt.subplots(len(xs_name_list), sharex=True)
    ax_tuple[0].set_title('{} {} {} effective cross sections evolution'.format(bucell, nuclide, xs_name_list))
    for i in range(len(xs_name_list)):
        for j in range(len(path_list)):
            ax_tuple[i].plot(bu_seq_list[j], xs_seq_list[j][i], label = name_list[j])
        ax_tuple[i].set_ylabel('{} Eff. XS [barn]'.format(xs_name_list[i]))
        ax_tuple[i].grid()
        ax_tuple[i].set_xlabel('BU [MWd/kg]')
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.legend()
    plt.show()

def plot_kinf_from_path(path_to_simulation):

    """Plots the multiplication factor against time (in days).

    Parameters
    ----------
    path_to_simulation: str
        Path to the simulation's directory
    """

    path = path_to_simulation + '/output_summary/kinf'

    time_seq = read_time_seq(path)
    kinf_seq = read_kinf_seq(path)

    plt.figure(1)
    plt.plot(time_seq, kinf_seq)
    plt.xlabel('Time [day]')
    plt.ylabel('kinf')
    plt.grid()
    plt.title('kinf evolution')

    plt.show()  

def plot_flux(bucell):

    path = os.getcwd() +'/{}_xs_lib'.format(bucell)

    time_seq = read_time_seq(path)
    flux_seq = read_flux(path)

    plt.figure(1)
    plt.step(time_seq, flux_seq, where = 'pre', color = 'blue')
    plt.xlabel('Time [day]')
    plt.ylabel('Flux [neutron/cm3.s-1]')
    plt.grid()
    plt.title('{} neutron flux evolution'.format(bucell))

    plt.show()

def plot_flux_from_path(bucell, path_to_simulation):

    """Plots the neutron flux against time (in days).

    Parameters
    ----------
    bucell: str
        Name of the BUCell
    path_to_simulation: str
        Path to the simulation's directory
    """

    path = path_to_simulation + '/output_summary/{}_dens'.format(bucell)

    time_seq = read_time_seq(path)
    flux_seq = read_flux(path)

    plt.figure(1)
    plt.step(time_seq, flux_seq, where = 'pre', color = 'blue')
    plt.xlabel('Time [day]')
    plt.ylabel('Flux [neutron/cm3.s-1]')
    plt.grid()
    plt.title('{} neutron flux evolution'.format(bucell))

    plt.show()

def plot_flux_spectrum_bu_evolution_from_path(bucell_list, steps_list, path):

    """Plots the neutron spectrum for different macrosteps for multiple BUCells. Each plot represents the spectrum evolution in one BUCell.

    Parameters
    ----------
    bucell_list: List
        List of the names of the BUCells
    steps_list: sList
        List of the macrosteps for which neutron spectrum should be plotted
    path: str
        Path to the simulation's directory
    """
    bucell_index = 0
    for bucell in bucell_list:
        path_flux_spectrum = path +'/{}_flux_spectrum'.format(bucell)

        flux_spectrum_list = read_flux_spectrum(path_flux_spectrum, steps_list)

        energy_mid_points = read_energy_mid_points(path_flux_spectrum)

        marker_list = ['x', '+', 'o', '*', '^', 's']
        color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

        plt.figure(bucell_index)
        index = 0
        for flux_spectrum in flux_spectrum_list:
            plt.plot(energy_mid_points, flux_spectrum, color = color_list[index], marker = marker_list[index], label=steps_list[index])
            index += 1


        plt.xlabel('eV')
        plt.ylabel('Neutron spectrum')
        plt.legend()
        plt.grid()
        plt.title('neutron spectrum in cell {} for steps {} '.format(bucell, steps_list))
        
    plt.show()

def plot_lethargy_spectrum_bu_evolution_from_path(bucell_list, steps_list, path):

    """Plots the neutron spectrum in lethargy units for different macrosteps for multiple BUCells. Each plot represents the spectrum evolution in one BUCell.

    Parameters
    ----------
    bucell_list: List
        List of the names of the BUCells
    steps_list: List
        List of the macrosteps for which neutron spectrum should be plotted
    path: str
        Path to the simulation's directory
    """
    bucell_index = 0
    for bucell in bucell_list:
        path_flux_spectrum = path +'/output_summary/{}_flux_spectrum'.format(bucell)

        flux_spectrum_list = read_flux_spectrum(path_flux_spectrum, steps_list)

        energy_mid_points = read_energy_mid_points(path_flux_spectrum)
        energy_bin_length = read_energy_bin_length(path_flux_spectrum)

        marker_list = ['x', '+', 'o', '*', '^', 's']
        color_list = ['r', 'b', 'g', 'k', 'brown', 'orange']

        plt.figure(bucell_index)
        index = 0
        for flux_spectrum in flux_spectrum_list:
            lethargy_spectrum = [x*y/z for x,y,z in zip(flux_spectrum, energy_mid_points, energy_bin_length)]
            plt.plot(energy_mid_points, lethargy_spectrum, color = color_list[index], marker = marker_list[index], label=steps_list[index])
            index += 1


        plt.xlabel('eV')
        plt.ylabel('Neutron spectrum [lethargy unit]')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.grid()
        plt.title('neutron spectrum in cell {} for steps {} '.format(bucell, steps_list))
        
        bucell_index += 1

    plt.show()

def plot_xs_dens_flux(bucell, xs_nuclide, xs_name, dens_nuclide, xs_path, dens_path):

    """Plots the evolution of a cross section, a nuclide density and the neutron flux in a specified BUCell in three subplots witht the same y-axis in days. This plot can be useful to understand the dynamic between the neutron flux and the absorption rates of certain nuclides.

    Parameters
    ----------
    bucell: str
        Name of the BUCell
    xs_nuclide: str
        Name of the nuclide for which the cross section should be plotted
    xs_name: str
        Name of the cross section to plot
    dens_nuclide: str
        Name of the nuclide for which the density should be plotted
    xs_path: str
        Path to the simulation's directory from which the cross section data should be extracted
    dens_path: str
        Path to the simulation's directory from which the nuclides density data should be extracted. The neutron flux will be taken from the same simulation.
    """

    time_subseq = read_time_seq(dens_path)
    time_seq = read_time_seq(xs_path)
    xs_seq = read_xs_seq(xs_nuclide, xs_name, xs_path)
    dens_subseq = read_dens(dens_nuclide, dens_path)
    flux_subseq = read_flux(dens_path)

    # Three subplots sharing both x/y axes
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    ax1.plot(time_seq, xs_seq, color = 'teal', marker = 'o')
    ax1.set_ylabel('{}\n{} Eff. XS [barn]'.format(xs_nuclide, xs_name))
    ax1.set_title('{} cell'.format(bucell))
    ax1.grid()
    ax2.plot(time_subseq, dens_subseq, color = 'orange', marker = 'o')
    ax2.set_ylabel('{}\nDensity [atm/barn-cm]'.format(dens_nuclide))
    ax2.grid()
    ax3.step(time_subseq, flux_subseq, color = 'blue')
    ax3.set_ylabel('Neutron flux [neutron/cm2-s]')
    ax3.set_xlabel('Time [day]')
    ax3.grid()
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.show()


# This function reads the old format for densities
# Should be removed
def read_time_seq_old_version(path):

    time_file = open(path, 'r')

    lines = time_file.readlines()

    # Find and store time
    for line in lines:
        if line != '\n':
            if line.split()[0] == 'TIME':
                time_seq = [float(x) for x in line.split()[1:]]
                break

    return time_seq

def read_time_seq(path):

    """Reads the time sequence from a simulation's directory and returns a list of the time sequence.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    """
    time_file = open(path, 'r')

    lines = time_file.readlines()

    # Find and store time
    for line in lines:
        print (line.split())
        if line != '\n':
            if line.split()[2] == 'TIME':
                time_seq = [float(x) for x in line.split()[4:]]
                break

    return time_seq

def get_step_time_length_seq(path):

    """Creates a list with time length between each macrosteps.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    """

    time_seq = read_time_seq(path)

    if time_seq[0] == 0:
        step_time_length = [x-y for x,y in zip(time_seq[1:], time_seq[:-1])]
    else:
        step_time_length = [time_seq[0]] +  [x-y for x,y in zip(time_seq[1:], time_seq[:-1])]
    
    return step_time_length


def read_bu_seq(path):

    """Reads the burnup sequence from a simulation's directory and returns a list of the burnup sequence.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    """
    bu_file = open(path, 'r')

    lines = bu_file.readlines()

    # Find and store time
    for line in lines:
        if line != '\n':
            if line.split()[0] == 'SYSTEM-BU':
                bu_seq = [float(x) for x in line.split()[1:]]
                break

    return bu_seq

def read_kinf_seq(path):

    """Reads the multiplication factor evolution from a simulation's directory and returns a list of the multiplication factor evolution.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    """
    kinf_file = open(path, 'r')

    lines = kinf_file.readlines()

    for line in lines:
        if line != '\n':
            if line.split()[0] == 'K-INF':
                kinf_seq = [float(x) for x in line.split()[1:]]
                break

    return kinf_seq


def read_flux(path):

    """Reads the macrostep-wise neutron flux evolution from a simulation's directory and returns a list with the neutron flux evolution.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    """
    flux_file = open(path, 'r')

    lines = flux_file.readlines()

    for line in lines:
        if line != '\n':
            if line.split()[0] == 'FLUX':
                flux_seq = [float(x) for x in line.split()[1:]]
                break
    # Add an initial 0 for the flux to have same len for flux_seq and time_seq
    flux_seq = [0] + flux_seq

    return flux_seq

def read_flux_subseq(path):

    """Reads the microstep-wise neutron flux evolution from a simulation's directory and returns a list with the neutron flux evolution.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    """

    flux_file = open(path, 'r')

    lines = flux_file.readlines()

    for line in lines:
        if line != '\n':
            if line.split()[0] == 'FLUX':
                flux_seq = [float(x) for x in line.split()[1:]]
                break

    return flux_seq

def get_fluence_seq(path, cell):

    """Returns a list with the macrostep-wise fluence evolution in a specificed BUCell

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    """

    dens_file = path +'/output_summary/{}_dens'.format(cell)

    flux_seq = read_flux(dens_file)
    time_seq = read_time_seq(dens_file)

    #Reminder: flux starts with 0

    fluence_seq = [0]
    pre_time = 0
    for i in range(1, len(time_seq)):
        time = time_seq[i]
        flux = flux_seq[i]
        time_int = (time-pre_time)*24*3600
        fluence_seq.append(time_int*flux + fluence_seq[i-1])
        pre_time = time

    return fluence_seq

def get_fluence_subseq(path, cell):

    """Returns a list with the microstep-wise fluence evolution in a specificed BUCell

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    """

    subdens_file = path +'/{}_subdens'.format(cell)

    flux_subseq = read_flux(subdens_file)
    # Subdens output not yet under new format
    time_subseq = read_time_seq_old_version(subdens_file)

    #Reminder: flux starts with 0

    fluence_subseq = [0]
    pre_time = 0
    for i in range(1, len(time_subseq)):
        time = time_subseq[i]
        flux = flux_subseq[i]
        time_int = (time-pre_time)*24*3600
        fluence_subseq.append(time_int*flux + fluence_subseq[i-1])
        pre_time = time

    return fluence_subseq





def read_flux_spectrum(path, steps_list):
    """Reads the neutron flux spectra for multiple macrosteps from a flux spectrum output file and returns it in a list.

    Parameters
    ----------
    path: str
        Path to the neutron flux spectrum output file
    step_list: List
        List of macrostep for which spectra should be read
    """
    flux_spectrum_file = open(path, 'r')

    lines = flux_spectrum_file.readlines()

    flux_spectrum_list = []

    step_count = 0
    for line in lines[6:]: # The flux spectrum data always start at the 6th line
        if step_count in steps_list:
            flux_spectrum_list.append([float(x) for x in line.split()[3:]])
        step_count += 1

    return flux_spectrum_list

def read_energy_mid_points(path):

    """Reads the mid-points of the energy groups against which the multi-group neutron flux is stored and returns it in a list.

    Parameters
    ----------
    path: str
        Path to the neutron flux spectrum output file
    """

    flux_spectrum_file = open(path, 'r')

    lines = flux_spectrum_file.readlines()

    energy_mid_points = [float(x) for x in lines[2].split()[1:]]

    return energy_mid_points

def read_energy_bin_length(path):
    
    """Reads the energy interval of each energy bin group used to store the multi-group neutron flux and returns it in a list.

    Parameters
    ----------
    path: str
        Path to the neutron flux spectrum output file
    """
    flux_spectrum_file = open(path, 'r')

    lines = flux_spectrum_file.readlines()

    energy_bin_length = [float(x) for x in lines[1].split()[1:]]

    return energy_bin_length

def read_dens(nuclide, path):
    
    """Reads the density (in :math:`atm barn^{-2}cm^{-1}`) of a specified nuclide from the provided density file.

    Parameters
    ----------
    nuclide: str
        Name of the nuclide
    path: str
        Path to the density file
    """
    zamid = name_to_zamid(nuclide)

    dens_seq = []

    dens_file = open(path, 'r')

    lines = dens_file.readlines()

    for line in lines:
        if line != '\n':
            if line.split()[1] == zamid:
                dens_seq = [float(x) for x in line.split()[2:]]
                break

    return dens_seq

# For simulations using old output format
def read_dens_old_version(nuclide, path):

    zamid = name_to_zamid(nuclide)

    dens_seq = []

    dens_file = open(path, 'r')

    lines = dens_file.readlines()

    for line in lines:
        if line != '\n':
            if line.split()[0] == zamid:
                dens_seq = [float(x) for x in line.split()[1:]]
                break

    return dens_seq

# cumulative dens
# This is not correct and should be deleted
def get_cum_dens(nuclide, path):

    dens_seq = read_dens(nuclide, path)

    cum_dens_seq = []
    cum_dens_seq.append(dens_seq[0])
    for i in range(1, len(dens_seq)):
        cum_dens_seq.append(dens_seq[i] + cum_dens_seq[i-1])

    return cum_dens_seq

# This is not correct and should be deleted
def convert_dens_seq_to_cum_dens_seq(dens_seq):

    cum_dens_seq = []
    cum_dens_seq.append(dens_seq[0])
    for i in range(1, len(dens_seq)):
        cum_dens_seq.append(dens_seq[i] + cum_dens_seq[i-1])

    return cum_dens_seq

# This should be in utils.functions
# It has been moved to functions
# def get_nucl_atomic_mass(nucl):

#   zamid = name_to_zamid(nucl)
#   zaid = zamid[:-1]
#   if zaid in d.default_atm_mass_lib:
#       M = d.default_atm_mass_lib[zaid]
#   else:
#       M = int(get_zamid_a(zamid))

#   return M

# calculate total density at certain step
def get_total_density(path, cell, step):

    """Returns the total density of the material (in :math:`atm barn^{-2}cm^{-1}`) in a BUCell at a given macrostep.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    step: int
        Macrostep number

    """
    nucl_name_list = read_dens_nucl(path, cell)
    dens_path = path+'/{}_dens'.format(cell)
    dens_file = open(dens_path )

    total_density = 0

    for nucl in nucl_name_list:
        dens = read_dens(nucl, dens_path)[step]
        total_density += dens

    return total_density    

# calculate total mass density at certain step
def get_total_mass_density(path, cell, step):

    """Returns the total mass density of the material (in :math:`g cm^{-3}`) in a BUCell at a given macrostep.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    step: int
        Macrostep number
    """

    nucl_name_list = read_dens_nucl(path, cell)
    dens_path = path+'/{}_dens'.format(cell)
    dens_file = open(dens_path )

    NA = d.NA

    total_mass_density = 0

    for nucl in nucl_name_list:
        dens = read_dens(nucl, dens_path)[step]
        M = get_nucl_atomic_mass(nucl)
        mass_density = dens*(M/NA)
        total_mass_density += mass_density

    return total_mass_density


# cumulative plutonium production
# This is not correct and should be deleted
def get_cum_pu_subseq_mat(path, cell, EFPD):

    path = path +'/{}_subdens'.format(cell)

    name_list = ['Pu-238', 'Pu-239', 'Pu-240', 'Pu-241', 'Pu-242', 'Pu-243']

    final_substep =find_substep_from_time(path, cell, EFPD)
    
    time_subseq = read_time_seq(path)

    t_before = time_subseq[final_substep]
    t_after = time_subseq[final_substep+1]

    cum_pu_subseq_mat = []
    for name in name_list:
        dens_subseq = read_dens(name, path)
        dens_subseq_until_time = dens_subseq[:final_substep+1]
        dens_before = dens_subseq[final_substep]
        dens_after = dens_subseq[final_substep+1]
        pair1 = [t_before, dens_before]
        pair2 = [t_after, dens_after]
        interpolated_dens = interpolation_between_two_points(pair1, pair2, EFPD)
        dens_subseq = dens_subseq_until_time + [interpolated_dens]
        cum_dens_subseq = convert_dens_seq_to_cum_dens_seq(dens_sbuseq)

        cum_pu_subseq_mat.append(cum_dens_subseq)

    return cum_pu_subseq_mat


# linear interpolation between two points
# This should be in utils.function
# # It has been moved to utils.function
# def interpolation_between_two_points(pair1, pair2, x):

#   a = (pair2[1] - pair1[1])/(pair2[0] - pair1[0])
#   b = (pair1[1]*pair2[0] - pair1[0]*pair2[1])/(pair2[0] - pair1[0])

#   y = a*x+b

#   return y


def read_xs_seq(nuclide, xs_name, path, cell):

    """Reads the cross section evolution (in barn) for a specified reactions for a given nuclide in a given BUCell and returns it in a list.

    Parameters
    ----------
    nuclide: str
        Name of the nuclide
    xs_name: str
        Name of the cross section
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    """
    path = path + '/output_summary/{}_xs_lib'.format(cell)

    zamid = name_to_zamid(nuclide)
    xs_name_found = 'no'

    xs_file = open(path, 'r')

    lines = xs_file.readlines()

    # Search for the line
    line_index = 0
    for line in lines:
        if line != '\n':
            if line.split()[0] == nuclide:
                break
        line_index += 1

    # First line needs to be treated differently
    if lines[line_index].split()[2] == xs_name:
        xs_name_found = 'yes'
        xs_seq= [float(x) for x in lines[line_index].split()[3:]]

    xs_loop = line_index+1

    while lines[xs_loop].split()[0] == zamid:
        if lines[xs_loop].split()[1] == xs_name:
            xs_name_found = 'yes'
            xs_seq = [float(x) for x in lines[xs_loop].split()[2:]]
        xs_loop += 1

    if xs_name_found == 'no':
        raise xs_name_not_found("nuclide {} has no data for cross section {}".format(nuclide, xs_name))
    else:
        return xs_seq

def get_time_averaged_xs(nuclide, xs_name, path, cell):
    
    """Computes time-averaged cross section (in barn) for a specified reactions for a given nuclide in a given BUCell.

    Parameters
    ----------
    nuclide: str
        Name of the nuclide
    xs_name: str
        Name of the cross section
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    """

    xs_lib_path = path + '/{}_dens'.format(cell)
    xs_seq = read_xs_seq(nuclide, xs_name, xs_lib_path, cell)
    dens_path = path + '/{}_dens'.format(cell)
    time_seq = read_time_seq(dens_path)
    tot_time = time_seq[-1]

    av_xs = 0

    for i in range(len(xs_seq)):

        xs = xs_seq[i]
        time_bos = time_seq[i]
        time_eos = time_seq[i+1]
        time_coeff = (time_eos - time_bos)/tot_time
        av_xs += xs*time_coeff

    return av_xs

def get_time_averaged_flux(path, cell):
    
    """Computes the time-averaged neutron flux (in :math:`cm^{-2}s^{-1}`) in a given BUCell.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    """

    xs_lib_path = path + '/{}_dens'.format(cell)
    flux_seq = read_flux(xs_lib_path)
    time_seq = read_time_seq(xs_lib_path)
    tot_time = time_seq[-1]

    av_flux = 0

    for i in range(len(flux_seq)-1):

        flux = flux_seq[i+1]
        time_bos = time_seq[i]
        time_eos = time_seq[i+1]
        time_coeff = (time_eos - time_bos)/tot_time
        av_flux += flux*time_coeff

    return av_flux


def get_tot_xs(nuclide, path, cell):

    """Computes the total one-group cross section evolution (in barn) for a specified nuclide in a given BUCell.

    Parameters
    ----------
    nuclide: str
        Name of the nuclide
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    """
    xs_name_list = ['fission','(n,gamma)','(n,2n)','(n,3n)','(n,p)','(n,a)','(n,gamma)X']

    tot_xs_seq = []
    i = 0
    for xs_name in xs_name_list:
        try:
            xs_seq = read_xs_seq(nuclide, xs_name, path, cell)
        except xs_name_not_found:
            continue
        if i == 0:
            tot_xs_seq = xs_seq
        else:
            tot_xs_seq = [x+y for x,y in zip(tot_xs_seq,xs_seq)]

        i += 1

    return tot_xs_seq




# This method list all nuclides that are present in xs lib
def read_xs_nucl(path, bucell):

    """Returns the list of nuclides for which one-group cross sections have been calculated during a simulation in a given BUCell.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    """

    path = path +'/output_summary/{}_xs_lib'.format(bucell)

    xs_lib_file = open(path)

    lines = xs_lib_file.readlines()

    nucl_name_list = []
    for line in lines:
        if line == '\n':
            continue
        line = line.split()
        if line[0].split('-')[0] in d.nuc_name_dic:
            nucl_name_list.append(line[0])

    xs_lib_file.close()


    return nucl_name_list

# make a list of all nuclide present in the dens_file
def read_dens_nucl(path, cell):

    """Returns the list of nuclides for which densities have been calculated during a simulation in a given BUCell.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    cell: str
        Name of the BUCell
    """

    dens_path = path + '/{}_dens'.format(cell)
    dens_file = open(dens_path, 'r')

    lines = dens_file.readlines()

    nucl_zamid_list = []
    for line in lines[7:]:
        line = line.split()
        nucl_zamid_list.append(line[1])

    dens_file.close()

    nucl_name_list = zamid_list_to_name_list(nucl_zamid_list)

    return nucl_name_list


def rank_nuclide_per_dens(bucell, step_list, path):
    """Print a ranking of nuclides in a BUCell according to their densities for multiple macrosteps in a text file.

    Parameters
    ----------
    bucell: str
       Name of the BUCell
    step_list: List
       List of macrostep numbers
    path: str
       Path to the simulation's directory
    """
    dens_path = path +'/output_summary' + '/{}_dens'.format(bucell)
    dens_file = open(dens_path, 'r')

    lines = dens_file.readlines()

    dens_list_per_step = []

    for step in step_list:

        dens_dict = {}
        # Data starts at 8th line
        for line in lines[7:]:
            line = line.split()
            dens_dict[line[0]] = float(line[step+1])

        # for key, value in sorted(dens_dict.iteritems(), key=lambda (k,v): (v,k)):
        #     print ("%s: %s" % (key, value))

        sorted_dens_dict = sorted(dens_dict.items(), key=lambda kv: kv[1], reverse=True)
        dens_list_per_step.append(sorted_dens_dict)

    dens_file.close()

    cwd = os.getcwd()
    sorted_dens_file = open('ranked dens', 'w')

    txt = ''
    for step in step_list:
        txt += '{:<20}'.format(step)
    txt += '\n\n'

    for i in range(len(dens_list_per_step[0])):
        for step in step_list:
            dens_list = dens_list_per_step[step]
            txt += '{:<8}{:<12.2E}'.format(dens_list[i][0], dens_list[i][1])
        txt += '\n' 

    sorted_dens_file.write(txt)
    sorted_dens_file.close()

# This function is broken
def rank_nuclide_per_reac_rate(bucell, step_list, path, file_name):

    dens_path = path +'/output_summary' + '/{}_dens'.format(bucell)
    dens_file = open(dens_path, 'r')
    xs_path = path +'/output_summary' + '/{}_xs_lib'.format(bucell)
    xs_file = open(xs_path, 'r')

    # Read densities
    lines = dens_file.readlines()
    dens_dict_per_step = []
    for step in step_list:
        dens_dict = {}
        # Data starts at 8th line
        for line in lines[7:]:
            line = line.split()
            name = zamid_to_name(line[0])
            dens_dict[name] = float(line[step+1])

        dens_dict_per_step.append(dens_dict)

    # Read xs
    lines = xs_file.readlines()
    xs_dict_per_step = []
    for step in step_list:
        xs_dict = {}
        count = 0
        # Data starts at 8th line
        for line in lines[7:]:
            if line == '\n':
                continue
            line = line.split()

            if line[0].split('-')[0] in d.nuc_name_dic:
                # For first data, you need to read data first
                if count == 0:
                    nucl_name = line[0]
                    abs_xs = float(line[step+3])

                # Reached new nuclide, need to store data
                else:
                    xs_dict[nucl_name] = abs_xs
                    nucl_name = line[0]
                    abs_xs = float(line[step+3])

            else:

                abs_xs += float(line[step+2])
                count += 1

        xs_dict_per_step.append(xs_dict)

    flux_per_step = []
    flux_seq = read_flux(dens_path)
    for step in step_list:
        flux_per_step.append(flux_seq[step+1]) # flux starts with 0


    # Now we will go over each nuclide in the dict per step and multiply xs with density
    #print (xs_dict_per_step)

    sorted_reac_dict_per_step = []
    total_abs_per_step = []
    for step in range(len(step_list)):
        dens_dict = dens_dict_per_step[step]
        xs_dict = xs_dict_per_step[step]
        reac_dict = {}
        for nucl in xs_dict:
            if nucl in dens_dict:
                nucl_dens = dens_dict[nucl]
                nucl_xs = xs_dict[nucl]
                reac_rate = nucl_dens*nucl_xs
                reac_dict[nucl] = reac_rate

        # Create a sorted list of tuples
        sorted_reac_tuple = sorted(reac_dict.items(), key=lambda kv: kv[1], reverse=True)
        sorted_reac_dict_per_step.append(sorted_reac_tuple)
        total_abs = 0
        for i in sorted_reac_tuple:
            total_abs += i[1]
        total_abs_per_step.append(total_abs)

    dens_file.close()
    xs_file.close()


    cwd = os.getcwd()
    sorted_reac_file = open('{} ranked react'.format(file_name), 'w')

    txt = ''
    for step in step_list:
        txt += '{:<20}'.format(step)
    txt += '\n\n'

    for step in range(len(step_list)):
        txt += 'flux={:<10.5E}'.format(flux_per_step[step])

    txt += '\n'

    for step in range(len(step_list)):
        txt += 'tot-abs={:<10.5E}'.format(total_abs_per_step[step])

    txt += '\n\n'

    for i in range(len(sorted_reac_dict_per_step[0])):
        for step in range(len(step_list)):
            reac_tuple_list = sorted_reac_dict_per_step[step]
            txt += '{:<8}{:<12.2E}'.format(reac_tuple_list[i][0], reac_tuple_list[i][1])
        txt += '\n' 

    sorted_reac_file.write(txt)
    sorted_reac_file.close()

def plot_matrix_from_compressed_matrix(path, step, cell):

    path_to_xs = path +'/step_{}'.format(step) +'/{}_cell'.format(cell) +'/matrix/xs_mat'
    path_to_decay = path +'/step_{}'.format(step) +'/{}_cell'.format(cell) +'/matrix/decay_mat'
    file_xs = open(path_to_xs, 'r')
    file_decay = open(path_to_decay, 'r')
    lines_xs = file_xs.readlines()
    lines_decay = file_decay.readlines()

    plt.figure(1)

    count = 0
    # x_vect = [i for i in range(len(lines))]
    for i in range(len(lines_xs)):
        line_xs = lines_xs[i]
        line_decay = lines_decay[i]
        zamid = line_xs.split('|')[0]
        line_elt_xs = line_xs.split(':')[1]
        elts_xs = line_elt_xs.split(',')[:-1] # Last element empty because of last coma in each line
        current_line_xs = []
        current_x_vect_xs = []
        elt_val_xs = len(lines_xs) - count  # The value assigned is the index of the nuclide starting from the last
        for elt_xs in elts_xs:
            elt_index = int(elt_xs.split()[0])
            current_x_vect_xs .append(elt_index)
            current_line_xs.append(elt_val_xs)
        if i == len(lines_xs)-1:
            plt.scatter(current_x_vect_xs , current_line_xs , marker='s', color = 'k', s = 4, label = 'cross section')
        else:
            plt.scatter(current_x_vect_xs , current_line_xs , marker='s', color = 'k', s = 4)
        line_elt_decay = line_decay.split(':')[1]
        elts_decay = line_elt_decay.split(',')[:-1] # Last element empty because of last coma in each line
        current_line_decay = []
        current_x_vect_decay = []
        elt_val_decay = len(lines_decay) - count  # The value assigned is the index of the nuclide starting from the last
        for elt_decay in elts_decay:
            elt_index = int(elt_decay.split()[0])
            current_x_vect_decay.append(elt_index)
            current_line_decay.append(elt_val_decay)
        if i == len(lines_xs)-1:    
            plt.scatter(current_x_vect_decay , current_line_decay , marker='+', color = 'r', s = 4, label = 'decay')
        else:
            plt.scatter(current_x_vect_decay , current_line_decay , marker='+', color = 'r', s = 4)
        #mat.append(current_line)
        count += 1

    file_xs.close()
    file_decay.close()
    plt.legend()
    plt.show()

def plot_matrix_bysign_from_compressed_matrix(path, step, cell):
    """Plots the transmuation matrix for a given BUCell at a given macrostep. Negative elements are in black while positive elements are in red.

    Parameters
    ----------
    path: str
        Path to the simulation's directory
    step: int
        Macrostep number
    cell: str
        Name of the BUCell
    """   
    #plt.style.use('dark_background')
    path_to_xs = path +'/step_{}'.format(step) +'/{}_cell'.format(cell) +'/matrix/xs_mat'
    path_to_decay = path +'/step_{}'.format(step) +'/{}_cell'.format(cell) +'/matrix/decay_mat'
    file_xs = open(path_to_xs, 'r')
    file_decay = open(path_to_decay, 'r')
    lines_xs = file_xs.readlines()
    lines_decay = file_decay.readlines()

    size = 0.8

    plt.figure(1, figsize=(7.3,5.5))

    count = 0
    # x_vect = [i for i in range(len(lines))]
    for i in range(len(lines_xs)):
        line_xs = lines_xs[i]
        line_decay = lines_decay[i]
        zamid = line_xs.split('|')[0]
        line_elt_xs = line_xs.split(':')[1]
        elts_xs = line_elt_xs.split(',')[:-1] # Last element empty because of last coma in each line
        current_line_xs = []
        current_x_vect_xs = []
        #elt_val_xs = len(lines_xs) - count  # The value assigned is the index of the nuclide starting from the last
        elt_val_xs = count
        for elt_xs in elts_xs[1:]:
            elt_index = int(elt_xs.split()[0])
            current_x_vect_xs .append(elt_index)
            current_line_xs.append(elt_val_xs)
        if i == len(lines_xs)-1:
            plt.scatter(current_x_vect_xs , current_line_xs , marker='s', color = 'red', s = size, label = 'Production terms')
        else:
            plt.scatter(current_x_vect_xs , current_line_xs , marker='s', color = 'red', s = size)
        line_elt_decay = line_decay.split(':')[1]
        elts_decay = line_elt_decay.split(',')[:-1] # Last element empty because of last coma in each line
        current_line_decay = []
        current_x_vect_decay = []
        #elt_val_decay = len(lines_decay) - count  # The value assigned is the index of the nuclide starting from the last
        elt_val_decay = count
        for elt_decay in elts_decay:
            elt_index = int(elt_decay.split()[0])
            current_x_vect_decay.append(elt_index)
            current_line_decay.append(elt_val_decay)
        if i == len(lines_xs)-1:    
            plt.scatter(current_x_vect_decay , current_line_decay , marker='s', color = 'red', s = size)
        else:
            plt.scatter(current_x_vect_decay , current_line_decay , marker='s', color = 'red', s = size)
        #mat.append(current_line)
        count += 1

    # cover diag elt
    for i in range(len(lines_xs)):
        if i == 0:
            plt.scatter([i], [i],  marker='s', color = 'black', s = size, label = 'Destruction terms')
        else:
            plt.scatter([i], [i],  marker='s', color = 'black', s = size)
    file_xs.close()
    file_decay.close()
    plt.tick_params(labelsize=14)
    plt.gca().invert_yaxis()
    plt.ylabel('Row index', fontsize=16)
    plt.xlabel('Column index', fontsize=16)
    plt.legend(prop={'size': 15}, markerscale = 6)
    plt.grid(color = 'darkgray')
    plt.show()






def plot_nuclide_chart_color_per_nuclear_data(decay_path, fy_path, path_to_xs_xml):
    """Plots a chart of the nuclide where nuclides are colored only if present in the decay, cross section or fission yield library. Different colors are used to indicate that a nuclide has data in the decay, cross section or fission yield library. The 

    Parameters
    ----------
    decay_path: str
        Path to a decay library
    fy_path: str
        Path to a fission yield library
    """
    fy_dict = d.read_fy_lib(fy_path)
    fy_unordered_keys = get_keylist_from_dict(fy_dict)
    fy_nucl_list = order_nuclide_per_z(fy_unordered_keys)

    decay_dict = d.read_decay_lib(decay_path)
    decay_unordered_keys = get_keylist_from_dict(decay_dict)
    decay_nucl_list = order_nuclide_per_z(decay_unordered_keys)

    MC_xs_nucl_name_list = get_openmc_xs_nucl_list(path_to_xs_xml)
    xs_nucl_list_name_list = mc_namelist_to_bu_namelist(MC_xs_nucl_name_list)
    xs_nucl_list = name_list_to_zamid_list(xs_nucl_list_name_list)
    xs_nucl_list = order_nuclide_per_z(xs_nucl_list)

    NAX_nucl_list = d.NAX_nucl_list

    orphan_NAX_z = []
    orphan_NAX_n = []

    only_decay_z = []
    only_decay_n = []
    only_decay_NAX_z = []
    only_decay_NAX_n = []   

    only_fy_z = []
    only_fy_n = []
    only_fy_NAX_z = []
    only_fy_NAX_n = []

    only_xs_z = []
    only_xs_n = []
    only_xs_NAX_z = []
    only_xs_NAX_n = []
    only_xs_zamid = []

    decay_xs_z = []
    decay_xs_n = []
    decay_xs_NAX_z = []
    decay_xs_NAX_n = []

    decay_fy_z = []
    decay_fy_n = []
    decay_fy_NAX_z = []
    decay_fy_NAX_n = []

    xs_fy_z = []
    xs_fy_n = []
    xs_fy_NAX_z = []
    xs_fy_NAX_n = []

    decay_xs_fy_z = []
    decay_xs_fy_n = []
    decay_xs_fy_NAX_z = []
    decay_xs_fy_NAX_n = []

    for decay_zamid in decay_nucl_list:

        if get_zamid_s(decay_zamid) == 1:
            continue

        in_fy = 'no'
        for fy_zamid in fy_nucl_list:
            if fy_zamid == decay_zamid:
                in_fy = 'yes'
                break

        in_xs = 'no'
        for xs_zamid in xs_nucl_list:
            if xs_zamid == decay_zamid:
                in_xs = 'yes'
                break

        in_NAX = 'no'
        for NAX_zamid in NAX_nucl_list:
            if NAX_zamid == decay_zamid:
                in_NAX = 'yes'
                break

        if in_fy =='no' and in_xs =='no' and in_NAX =='no':
            only_decay_z.append(get_zamid_z(decay_zamid))
            only_decay_n.append(get_zamid_n(decay_zamid))

        if in_fy =='no' and in_xs =='no' and in_NAX =='yes':
            only_decay_NAX_z.append(get_zamid_z(decay_zamid))
            only_decay_NAX_n.append(get_zamid_n(decay_zamid))

        if in_fy =='yes' and in_xs =='no' and in_NAX =='no':
            decay_fy_z.append(get_zamid_z(decay_zamid))
            decay_fy_n.append(get_zamid_n(decay_zamid))

        if in_fy =='yes' and in_xs =='no' and in_NAX =='yes':
            decay_fy_NAX_z.append(get_zamid_z(decay_zamid))
            decay_fy_NAX_n.append(get_zamid_n(decay_zamid))

        if in_fy =='no' and in_xs =='yes' and in_NAX =='no':
            decay_xs_z.append(get_zamid_z(decay_zamid))
            decay_xs_n.append(get_zamid_n(decay_zamid))

        if in_fy =='no' and in_xs =='yes' and in_NAX =='yes':
            decay_xs_NAX_z.append(get_zamid_z(decay_zamid))
            decay_xs_NAX_n.append(get_zamid_n(decay_zamid))

        if in_fy =='yes' and in_xs =='yes' and in_NAX =='no':
            decay_xs_fy_z.append(get_zamid_z(decay_zamid))
            decay_xs_fy_n.append(get_zamid_n(decay_zamid))

        if in_fy =='yes' and in_xs =='yes' and in_NAX =='yes':
            decay_xs_fy_NAX_z.append(get_zamid_z(decay_zamid))
            decay_xs_fy_NAX_n.append(get_zamid_n(decay_zamid))

    for fy_zamid in fy_nucl_list:

        if get_zamid_s(fy_zamid) == '1':
            continue

        in_xs = 'no'
        for xs_zamid in xs_nucl_list:
            if fy_zamid == xs_zamid:
                in_xs = 'yes'
                break

        in_decay = 'no'
        for decay_zamid in decay_nucl_list:
            if fy_zamid == decay_zamid:
                in_decay = 'yes'
                break

        in_NAX = 'no'
        for NAX_zamid in NAX_nucl_list:
            if fy_zamid == NAX_zamid:
                in_NAX = 'yes'
                break

        if in_xs =='no' and in_decay =='no' and in_NAX =='no':
            only_fy_z.append(get_zamid_z(fy_zamid))
            only_fy_n.append(get_zamid_n(fy_zamid))

        if in_xs =='no' and in_decay =='no' and in_NAX =='yes':
            only_fy_NAX_z.append(get_zamid_z(fy_zamid))
            only_fy_NAX_n.append(get_zamid_n(fy_zamid))

        if in_xs =='yes' and in_decay =='no' and in_NAX =='no':
            xs_fy_z.append(get_zamid_z(fy_zamid))
            xs_fy_n.append(get_zamid_n(fy_zamid))

        if in_xs =='yes' and in_decay =='no' and in_NAX =='yes':
            xs_fy_NAX_z.append(get_zamid_z(fy_zamid))
            xs_fy_NAX_n.append(get_zamid_n(fy_zamid))

    for xs_zamid in xs_nucl_list:

        if get_zamid_s(xs_zamid) == '1':
            continue

        in_decay = 'no'
        for decay_zamid in decay_nucl_list:
            if xs_zamid == decay_zamid:
                in_decay = 'yes'
                break

        in_fy = 'no'
        for fy_zamid in fy_nucl_list:
            if xs_zamid == fy_zamid:
                in_fy = 'yes'
                break

        in_NAX = 'no'
        for NAX_zamid in NAX_nucl_list:
            if xs_zamid == NAX_zamid:
                in_NAX = 'yes'
                break

        if in_fy =='no' and in_decay =='no' and in_NAX =='no':
            only_xs_z.append(get_zamid_z(xs_zamid))
            only_xs_n.append(get_zamid_n(xs_zamid))

        if in_fy =='no' and in_decay =='no' and in_NAX =='yes':
            only_xs_NAX_z.append(get_zamid_z(xs_zamid))
            only_xs_NAX_n.append(get_zamid_n(xs_zamid))

    for NAX_zamid in NAX_nucl_list:

        if get_zamid_s(NAX_zamid) == '1':
            continue

        in_decay = 'no'
        for decay_zamid in decay_nucl_list:
            if NAX_zamid == decay_zamid:
                in_decay = 'yes'
                break

        in_fy = 'no'
        for fy_zamid in fy_nucl_list:
            if NAX_zamid == fy_zamid:
                in_fy = 'yes'
                break

        in_xs = 'no'
        for xs_zamid in xs_nucl_list:
            if NAX_zamid == xs_zamid:
                in_xs = 'yes'
                break

        if in_fy =='no' and in_decay =='no' and in_xs =='no':
            orphan_NAX_z.append(get_zamid_z(NAX_zamid))
            orphan_NAX_n.append(get_zamid_n(NAX_zamid))




    # for i in range(len(decay_fy_n)):
    #   if decay_fy_n[i] == 52 and decay_fy_z[i] ==41:
    #       print (i)
    #       print ('410920 in decay_fy')

    # for i in range(len(decay_xs_fy_NAX_n)):
    #   if decay_xs_fy_NAX_n[i] == 52 and decay_xs_fy_NAX_z[i] ==41:
    #       print (i)
    #       print ('410920 in decay_xs_fy_NAX')


    size = 13

    # setup the plot
    plt.figure(1)

    plt.scatter(only_decay_n, only_decay_z, marker = 's', color='k', label = 'decay only', s = size)
    plt.scatter(only_decay_NAX_n, only_decay_NAX_z, marker = '.', color='k', label = 'decay only NAX', s = size)
    
    plt.scatter(only_xs_n, only_xs_z, marker = 's', color='gold', label = 'xs only', s = size)
    plt.scatter(only_xs_NAX_n, only_xs_NAX_z, marker = '.', color='gold', label = 'xs only NAX', s = size)
    
    plt.scatter(only_fy_n, only_fy_z, marker = 's', color='pink', label = 'fy only', s = size)
    plt.scatter(only_fy_NAX_n, only_fy_NAX_z, marker = '.', color='pink', label = 'fy only NAX', s = size)


    plt.scatter(decay_fy_n, decay_fy_z, marker = 's', color='b', label = 'decay & fy', s = size)
    plt.scatter(decay_fy_NAX_n, decay_fy_NAX_z, marker = '.', color='b', label = 'decay & fy NAX', s = size)


    plt.scatter(decay_xs_n, decay_xs_z, marker = 's', color='lime', label = 'decay & xs', s = size)
    plt.scatter(decay_xs_NAX_n, decay_xs_NAX_z, marker = '.', color='lime', label = 'decay & xs NAX', s = size)


    plt.scatter(xs_fy_n, xs_fy_z, marker = 's', color='fuchsia', label = 'xs & fy', s = size)
    plt.scatter(xs_fy_NAX_n, xs_fy_NAX_z, marker = '.', color='fuchsia', label = 'xs & fy NAX', s = size)


    plt.scatter(decay_xs_fy_n, decay_xs_fy_z, marker = 's', color='r', label = 'decay xs fy', s = size)
    plt.scatter(decay_xs_fy_NAX_n, decay_xs_fy_NAX_z, marker = '.', color='r', label = 'decay xs fy NAX', s = size)

    plt.scatter(orphan_NAX_n, orphan_NAX_z, marker = '.', color='grey', label = 'orphan NAX', s = size)


    plt.title('Nuclide colored per data availability\nNAX nuclides included')
    plt.xlabel('N')
    plt.ylabel('Z')
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2, markevery=50)
    plt.minorticks_on()
    # plt.grid('on')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    #plt.set_minor_frequency(2)

    # for i, txt in enumerate(label_list):
    #   plt.annotate(txt, (a_list1[i], z_list1[i]))

    plt.show()

def plot_compare_two_nuclear_data_on_nuclide_chart(decay_path1, fy_path1,decay_path2, fy_path2, path_to_xs_xml):

    """Plots a chart of the nuclide where nuclides are colored only if present in the decay, cross section or fission yield library. Black color indicates that the nuclide has data in decay_path1 or fy_path1. Red color indicates that nuclide has data in decay_path2 or fy_path2 (or both in libraries 1 and 2).

    Parameters
    ----------
    decay_path1: str
        Path to decay library 1
    fy_path1: str
        Path to fission yield library 1
    decay_path2: str
        Path to decay library 2
    fy_path2: str
        Path to fission yield library 2
    path_to_xs_xml: str
        Path to the cross_sections.xml file in a cross section library
    """

    fy_dict1 = d.read_fy_lib(fy_path1)
    fy_unordered_keys1 = get_keylist_from_dict(fy_dict1)
    fy_nucl_list1 = order_nuclide_per_z(fy_unordered_keys1)

    decay_dict1 = d.read_decay_lib(decay_path1)
    decay_unordered_keys1 = get_keylist_from_dict(decay_dict1)
    decay_nucl_list1 = order_nuclide_per_z(decay_unordered_keys1)

    fy_dict2 = d.read_fy_lib(fy_path2)
    fy_unordered_keys2 = get_keylist_from_dict(fy_dict2)
    fy_nucl_list2 = order_nuclide_per_z(fy_unordered_keys2)

    decay_dict2 = d.read_decay_lib(decay_path2)
    decay_unordered_keys2 = get_keylist_from_dict(decay_dict2)
    decay_nucl_list2 = order_nuclide_per_z(decay_unordered_keys2)

    MC_xs_nucl_name_list = get_openmc_xs_nucl_list(path_to_xs_xml)
    xs_nucl_list_name_list = mc_namelist_to_bu_namelist(MC_xs_nucl_name_list)
    xs_nucl_list = name_list_to_zamid_list(xs_nucl_list_name_list)
    xs_nucl_list = order_nuclide_per_z(xs_nucl_list)

    NAX_nucl_list = d.NAX_nucl_list

    decay_z1 = []
    decay_n1 = []
    fy_z1 = []
    fy_n1 = []
    decay_z2 = []
    decay_n2 = []
    fy_z2 = []
    fy_n2 = []
    xs_z2 = []
    xs_n2 = []
    NAX_z2 = []
    NAX_n2 = []

    for decay_zamid1 in decay_nucl_list1:

        if get_zamid_s(decay_zamid1) == 1:
            continue

        decay_z1.append(get_zamid_z(decay_zamid1))
        decay_n1.append(get_zamid_n(decay_zamid1))


    for fy_zamid1 in fy_nucl_list1:

        if get_zamid_s(fy_zamid1) == '1':
            continue

        fy_z1.append(get_zamid_z(fy_zamid1))
        fy_n1.append(get_zamid_n(fy_zamid1))

    for decay_zamid2 in decay_nucl_list2:

        if get_zamid_s(decay_zamid2) == 1:
            continue

        decay_z2.append(get_zamid_z(decay_zamid2))
        decay_n2.append(get_zamid_n(decay_zamid2))


    for fy_zamid2 in fy_nucl_list2:

        if get_zamid_s(fy_zamid2) == '1':
            continue

        fy_z2.append(get_zamid_z(fy_zamid2))
        fy_n2.append(get_zamid_n(fy_zamid2))

    for xs_zamid2 in xs_nucl_list:

        if get_zamid_s(xs_zamid2) == '1':
            continue

        xs_z2.append(get_zamid_z(xs_zamid2))
        xs_n2.append(get_zamid_n(xs_zamid2))

    for NAX_zamid2 in NAX_nucl_list:

        if get_zamid_s(NAX_zamid2) == '1':
            continue

        NAX_z2.append(get_zamid_z(NAX_zamid2))
        NAX_n2.append(get_zamid_n(NAX_zamid2))



    size = 13

    # setup the plot
    plt.figure(1,figsize=(7.3,5.5))

    # Plot nuclear data 1
    plt.scatter(decay_n1, decay_z1, marker = 's', color='k', label = 'Full', s = size)
    plt.scatter(fy_n1, fy_z1, marker = 's', color='k', s = size)
    plt.scatter(xs_n2, xs_z2, marker = 's', color='k',s = size)
    #plt.scatter(NAX_n2, NAX_z2, marker = 's', color='k', s = size)

    # Plot nuclear data 2   
    plt.scatter(decay_n2, decay_z2, marker = 's', color='r', label = 'Reduced', s = size)
    plt.scatter(fy_n2, fy_z2, marker = 's', color='r', s = size)
    plt.scatter(xs_n2, xs_z2, marker = 's', color='r',s = size)
    #plt.scatter(NAX_n2, NAX_z2, marker = 's', color='b', s = size, label = 'Archaeology')


    plt.xlabel('Number of neutrons', fontsize=16)
    plt.ylabel('Number of protons', fontsize=16)
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='r', linestyle='-', alpha=0.2, markevery=50)
    plt.minorticks_on()
    plt.tick_params(labelsize=14)
    # plt.grid('on')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend(prop={'size': 15})
    #plt.set_minor_frequency(2)
    plt.subplots_adjust(left=0.21, bottom=0.110, right=0.98, top=0.98)

    # for i, txt in enumerate(label_list):
    #   plt.annotate(txt, (a_list1[i], z_list1[i]))

    plt.show()


def plot_nuclide_chart_compare_fy(lib1_path, lib2_path, fissile_parent):

    fy_dict1 = d.read_fy_lib(lib1_path)
    fy_dict2 = d.read_fy_lib(lib2_path)

    unordered_keys1 = get_keylist_from_dict(fy_dict1)
    unordered_keys2 = get_keylist_from_dict(fy_dict2)

    ordered_keys1 = order_nuclide_per_z(unordered_keys1)
    ordered_keys2 = order_nuclide_per_z(unordered_keys2)

    rel_diff = []
    label_list = []
    nucl_list = []
    fy_list1 = []
    fy_list2 = []

    # sum over fy in list 1
    for nuclide in ordered_keys1:
        fy_1 = fy_dict1[nuclide][fissile_parent][0]
        fy_list1.append(fy_1)

    # sum over fy in list 2
    for nuclide in ordered_keys2:
        fy_2 = fy_dict2[nuclide][fissile_parent][0]
        fy_list2.append(fy_2)


    for nuclide in ordered_keys1:
        if nuclide in ordered_keys2:
            fy_1 = fy_dict1[nuclide][fissile_parent][0]
            fy_2 = fy_dict2[nuclide][fissile_parent][0]

            if fy_1 == 0 and fy_2 == 0:
                label = 'both = 0'
                rel_diff_val = 0
            elif fy_1 == 0 and fy_2 != 0: # if fy1 = 0 but fy2 != 1, set label to 100%. The inverse situation yield - 100%
                label = 'ori = 0'
                rel_diff_val = 0
            elif fy_1 != 0 and fy_2 == 0: # if fy1 = 0 but fy2 != 1, set label to 100%. The inverse situation yield - 100%
                label = 'jef33 = 0'
                rel_diff_val = 0
            else:
                rel_diff_val = (fy_2 - fy_1)*100/fy_1

                if abs(rel_diff_val) > 500:
                    rel_diff_val = 500*abs(rel_diff_val)/rel_diff_val # set value to 500 or -500

                label = '{:4.2E}'.format((fy_2 - fy_1)*100/fy_1)

            # if rel_diff_val == 0:
            #   rel_diff_log.append(abs(rel_diff_val))
            #   rel_diff.append(rel_diff_val)
            # else:
            rel_diff.append(rel_diff_val)
            label_list.append(label)

            nucl_list.append(nuclide)



    not_in_jeff33_zamid = []
    for zamid in ordered_keys1:
        if zamid not in ordered_keys2:
            not_in_jeff33_zamid.append(zamid)

    not_in_endf4_zamid = []
    for zamid in ordered_keys2:
        if zamid not in ordered_keys1:
            not_in_endf4_zamid.append(zamid)

    # test if not_in list overlap
    for zamid in nuclide:
        if zamid in not_in_endf4_zamid:
            print (zamid)
        if zamid in not_in_jeff33_zamid:
            print (zamid)


    # loop over the z
    # z list for lib 1
    z_list1 = []
    for zamid in ordered_keys1:
        if zamid in ordered_keys2:
            z = get_zamid_z(zamid)
            z_list1.append(z)


    # loop over the z for the nuclide in endf4 but not jeff33
    not_in_jeff33_z_list = []
    for zamid in not_in_jeff33_zamid:
        z = get_zamid_z(zamid)
        not_in_jeff33_z_list.append(z)

    # loop over the z for the nuclide in jeff33 but not endf4
    not_in_endf4_z_list = []
    for zamid in not_in_endf4_zamid:
        z = get_zamid_z(zamid)
        not_in_endf4_z_list.append(z)

    # # z list for lib 2
    # z_list2 = []
    # for zamid in ordered_keys2:
    #   z = get_zamid_z(zamid)
    #   z_list2.append(z)

    # loop over the a
    # a list for lib 1
    a_list1 = []
    for zamid in ordered_keys1:
        if zamid in ordered_keys2:
            a = get_zamid_a(zamid)
            a_list1.append(a)

    # loop over the a for the nuclide in endf4 but not jeff33
    not_in_jeff33_a_list = []
    for zamid in not_in_jeff33_zamid:
        a = get_zamid_a(zamid)
        not_in_jeff33_a_list.append(a)

    # loop over the z for the nuclide in jeff33 but not endf4
    not_in_endf4_a_list = []
    for zamid in not_in_endf4_zamid:
        a = get_zamid_a(zamid)
        not_in_endf4_a_list.append(a)

    # loop over the n
    # a list for lib 1
    n_list1 = []
    for zamid in ordered_keys1:
        if zamid in ordered_keys2:
            n = get_zamid_n(zamid)
            n_list1.append(n)

    # loop over the n for the nuclide in endf4 but not jeff33
    not_in_jeff33_n_list = []
    for zamid in not_in_jeff33_zamid:
        n = get_zamid_n(zamid)
        not_in_jeff33_n_list.append(n)

    # loop over the n for the nuclide in jeff33 but not endf4
    not_in_endf4_n_list = []
    for zamid in not_in_endf4_zamid:
        n = get_zamid_n(zamid)
        not_in_endf4_n_list.append(n)


    # check if list overlap

    # # z list for lib 2
    # a_list2 = []
    # for zamid in ordered_keys2:
    #   a = get_zamid_a(zamid)
    #   a_list2.append(a)


    N = len(rel_diff) # Number of labels

    # setup the plot
    plt.figure(1)

    tag = rel_diff# Tag each point with a corresponding label    

    # define the colormap
    cmap = plt.cm.Spectral
    # extract all colors from the .jet map
    # cmaplist = [cmap(i) for i in range(cmap.N)]
    # # create the new map
    # cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)     

    # define the bins and normalize
    bounds = np.linspace(0,N,N+1)
    #norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # make the scatter
    scat_a_z = plt.scatter(n_list1,z_list1,c=tag,cmap=cmap, marker='s', s = [20]*len(z_list1), norm=MidpointNormalize(midpoint=0,vmin=-100, vmax=max(rel_diff)))
    #scat_a_z = plt.scatter(a_list1,z_list1,c=tag,cmap=cmap, marker='s', s = [20]*len(z_list1), norm=MidpointNormalize(midpoint=0,vmin=-100, vmax=max(rel_diff)))
    
    scat2 = plt.scatter(not_in_jeff33_a_list,not_in_jeff33_z_list, marker = 'o', color='k', s = [20]*len(not_in_jeff33_z_list), label='Not in JEF33')
    scat3 = plt.scatter(not_in_endf4_a_list,not_in_endf4_z_list, marker = '*', color='k', s = [20]*len(not_in_endf4_z_list), label='Not in ENDF4')
    
    #cb = plt.colorbar(scat, spacing='proportional',ticks=bounds)
    #plt.imshow(ras, cmap=cmap, clim=(elev_min, elev_max), norm=MidpointNormalize(midpoint=mid_val,vmin=elev_min, vmax=elev_max))

    cb = plt.colorbar(scat_a_z, spacing='proportional')
    cb.set_label('JEF33 < ENDF4              JEF33 > ENDF4')
    # create the colorbar
    plt.title('Relative Difference in Percent between ENDF4 and JEF33 Fission Yields')
    plt.xlabel('A')
    plt.ylabel('Z')
    plt.grid('on')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()

    # for i, txt in enumerate(label_list):
    #   plt.annotate(txt, (a_list1[i], z_list1[i]))

    plt.show()

def plot_compare_libs(lib1_path, lib2_path, fissile_parent):

    fy_dict1 = d.read_fy_lib(lib1_path)
    fy_dict2 = d.read_fy_lib(lib2_path)

    unordered_keys1 = get_keylist_from_dict(fy_dict1)
    unordered_keys2 = get_keylist_from_dict(fy_dict2)

    ordered_keys1 = order_nuclide_per_z(unordered_keys1)
    ordered_keys2 = order_nuclide_per_z(unordered_keys2)

    rel_diff = []
    abs_diff = []
    nucl_list = []
    fy_list1 = []
    fy_list2 = []

    # sum over fy in list 1
    for nuclide in ordered_keys1:
        fy_1 = fy_dict1[nuclide][fissile_parent][0]
        fy_list1.append(fy_1)

    # sum over fy in list 2
    for nuclide in ordered_keys2:
        fy_2 = fy_dict2[nuclide][fissile_parent][0]
        fy_list2.append(fy_2)


    for nuclide in ordered_keys1:
        if nuclide in ordered_keys2:
            fy_1 = fy_dict1[nuclide][fissile_parent][0]
            fy_2 = fy_dict2[nuclide][fissile_parent][0]

            if fy_1 == 0 and fy_2 == 0:
                rel_diff_val = 0
            elif fy_1 == 0 and fy_2 != 0: # if fy1 = 0 but fy2 != 1, set label to 100%. The inverse situation yield - 100%
                rel_diff_val = 100
            else:
                rel_diff_val = (fy_2 - fy_1)*100/fy_1

            abs_diff_val = fy_2 -fy_1

            rel_diff.append(rel_diff_val)
            abs_diff.append(abs_diff_val)


            nucl_list.append(nuclide)

    x_vect = [i for i in range(len(rel_diff))]

    z_vect = []
    x_z_vect = []
    count = 0
    for i in range(len(nucl_list)):
        if count == 0:
            z_vect.append(nucl_list[i][:-4])
            x_z_vect.append(i)
        count += 1

        if count == 49:
            count = 0

    plt.figure(1)
    plt.title('Relative percent difference between fission yield data')
    plt.xticks(x_z_vect, z_vect)
    plt.plot(x_vect, rel_diff, 'ro')
    plt.plot(x_vect, [100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')
    plt.ylabel('Relative percent difference')
    plt.xlabel('Atomic number')
    #plt.yscale('log')

    plt.figure(2)
    plt.xticks(x_vect, nucl_list)
    plt.plot(x_vect, rel_diff, 'ro')
    plt.plot(x_vect, [100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')
    #plt.yscale('log')

    plt.figure(3)
    plt.xticks(x_vect, nucl_list)
    plt.plot(x_vect, abs_diff, 'ro')
    # plt.plot(x_vect, [100]*len(rel_diff), 'k')
    # plt.plot(x_vect, [-100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')



    plt.show()


def plot_compare_libs_sum_over_parents(lib1_path, lib2_path, parent_list):

    fy_dict1 = d.read_fy_lib(lib1_path)
    fy_dict2 = d.read_fy_lib(lib2_path)

    unordered_keys1 = get_keylist_from_dict(fy_dict1)
    unordered_keys2 = get_keylist_from_dict(fy_dict2)

    ordered_keys1 = order_nuclide_per_z(unordered_keys1)
    ordered_keys2 = order_nuclide_per_z(unordered_keys2)

    rel_diff = []
    abs_diff = []
    nucl_list = []
    # fy_list1 = []
    # fy_list2 = []

    # # sum over fy in list 1
    # for nuclide in ordered_keys1:
    #   sum_fy = 0
    #   for parent in parent_list:
    #       sum_fy += fy_dict1[nuclide][parent][0]
    #   fy_list1.append(sum_fy)

    # # sum over fy in list 2
    # for nuclide in ordered_keys2:
    #   sum_fy = 0
    #   for parent in parent_list:
    #       sum_fy += fy_dict2[nuclide][parent][0]
    #   fy_list2.append(sum_fy)


    for nuclide in ordered_keys1:
        if nuclide in ordered_keys2:
            sum_fy1 = 0
            for parent in parent_list:
                sum_fy1 += fy_dict1[nuclide][parent][0]
            sum_fy2 = 0
            for parent in parent_list:
                sum_fy2 += fy_dict2[nuclide][parent][0]


            if sum_fy1 == 0 and sum_fy2 == 0:
                rel_diff_val = 0
            elif sum_fy1 == 0 and sum_fy2!= 0: # if fy1 = 0 but fy2 != 1, set label to 100%. The inverse situation yield - 100%
                rel_diff_val = 100
            else:
                rel_diff_val = (sum_fy2 - sum_fy1)*100/sum_fy1

            abs_diff_val = sum_fy2 -sum_fy1

            rel_diff.append(rel_diff_val)
            abs_diff.append(abs_diff_val)


            nucl_list.append(nuclide)

    x_vect = [i for i in range(len(rel_diff))]

    z_vect = []
    x_z_vect = []
    count = 0
    for i in range(len(nucl_list)):
        if count == 0:
            z_vect.append(nucl_list[i][:-4])
            x_z_vect.append(i)
        count += 1

        if count == 49:
            count = 0

    plt.figure(1)
    plt.title('Relative percent difference between fission yield data')
    plt.xticks(x_z_vect, z_vect)
    plt.plot(x_vect, rel_diff, 'ro')
    plt.plot(x_vect, [100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')
    plt.ylabel('Relative percent difference')
    plt.xlabel('Atomic number')
    #plt.yscale('log')

    plt.figure(2)
    plt.xticks(x_vect, nucl_list)
    plt.plot(x_vect, rel_diff, 'ro')
    plt.plot(x_vect, [100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')
    #plt.yscale('log')

    plt.figure(3)
    plt.xticks(x_vect, nucl_list)
    plt.plot(x_vect, abs_diff, 'ro')
    # plt.plot(x_vect, [100]*len(rel_diff), 'k')
    # plt.plot(x_vect, [-100]*len(rel_diff), 'k')
    plt.plot(x_vect, [-0.0]*len(rel_diff), 'k--')

    plt.show()
















#### The following functions are mainly used by the NAX module


# This method build the fluence sequence until a final time
# Warning this function is only valid for non-fissile material as the flux
# is the same within each step
def get_fluence_seq_until_time(path, cell, final_time):

    final_step = find_step_from_time(path, cell, final_time)
    extra_fluence = get_extra_fluence_from_time(path, cell, final_time)
    fluence_seq = get_fluence_seq(path, cell)

    fluence_seq_until_time = fluence_seq[:final_step+1] + [fluence_seq[final_step] + extra_fluence]

    return fluence_seq_until_time

def get_fluence_subseq_until_time(path, cell, final_time):

    final_substep =find_substep_from_time(path, cell, final_time)
    extra_fluence = get_extra_subfluence_from_time(path, cell, final_time)
    fluence_subseq = get_fluence_subseq(path, cell)

    fluence_subseq_until_time = fluence_subseq[:final_substep+1] + [fluence_subseq[final_substep] + extra_fluence]

    return fluence_subseq_until_time

# This function calculate the additional fluence from the previous time point
# to where the final time is set
def get_extra_fluence_from_time(path, cell, time):

    step = find_step_from_time(path, cell, time)
    dens_file = path +'/output_summary/{}_dens'.format(cell)
    # flux seq has an added zero at the beginning of the arry
    flux = read_flux(dens_file)[step+1]
    previous_time_point = read_time_seq(dens_file)[step]

    time_int = time-previous_time_point

    return flux*time_int*24*3600

# This function calculate the additional fluence from the previous time point
# to where the final time is set
# From subdens file
def get_extra_subfluence_from_time(path, cell, time):

    substep = find_substep_from_time(path, cell, time)
    subdens_file = path +'/{}_subdens'.format(cell)
    # flux seq has an added zero at the beginning of the arry

    flux_subseq = read_flux_subseq(subdens_file)
    flux= flux_subseq[substep]
    # Subdens output not yet under new format
    previous_time_point = read_time_seq_old_version(subdens_file)[substep]

    time_int = time-previous_time_point

    return flux*time_int*24*3600


def find_step_from_time(path, cell, time):

    dens_file = path +'/output_summary/{}_dens'.format(cell)
    time_seq = read_time_seq(dens_file)

    step = 0
    for t in time_seq[1:]:
        if time <= t:
            break
        step += 1

    return step

def find_substep_from_time(path, cell, time):

    subdens_file = path +'/{}_subdens'.format(cell)
    # Subdens output not yet under new format
    time_seq = read_time_seq_old_version(subdens_file)

    substep = 0
    for t in time_seq[1:]:
        if time <= t:
            break
        substep += 1

    return substep


def get_step_fluence_length(path, cell):

    fluence_seq = get_fluence_seq(path, cell)

    step_fluence_length = [x-y for x,y in zip(fluence_seq[1:], fluence_seq[:-1])]
    
    return step_fluence_length


def get_pu_subseq_mat(path, cell, EFPD):

    final_substep =find_substep_from_time(path, cell, EFPD)

    path = path +'/{}_subdens'.format(cell)

    name_list = d.Pu_isotopes_name
    
    # Subdens output not yet under new format
    time_subseq = read_time_seq_old_version(path)

    t_before = time_subseq[final_substep]
    t_after = time_subseq[final_substep+1]

    pu_subseq_mat = []
    for name in name_list:
        dens_subseq = read_dens_old_version(name, path)
        dens_subseq_until_time = dens_subseq[:final_substep+1]
        dens_before = dens_subseq[final_substep]
        dens_after = dens_subseq[final_substep+1]
        pair1 = [t_before, dens_before]
        pair2 = [t_after, dens_after]
        interpolated_dens = interpolation_between_two_points(pair1, pair2, EFPD)
        dens_subseq = dens_subseq_until_time + [interpolated_dens]

        pu_subseq_mat.append(dens_subseq)

    return pu_subseq_mat



class xs_name_not_found(Exception):
    """Raise when the user tries to access fission XS for a nuclide which fission XS have not been set yet """
    pass





# def bucells_average_density(bucell_list, step_list):

#   bucell_vol_list = []

#   system_parameters_file = open(os.getcwd() + '/output_summary/system_parameters')
    







 
