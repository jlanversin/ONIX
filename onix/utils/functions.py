# Small tricks and calculations to make the life of a nuclear engineer easier
from math import log
import os
import shutil
from onix.data import time_dic
import onix.data as d
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.colors as colors
import xml.etree.ElementTree as ET

NA = 6.02214086e+23

def decay_to_halflife(decay_constant, unit):
    """Converts a decay constant into a half life in specified units.

    Parameters
    ----------
    decay_constant: float
        Decay constant of the nuclide in :math:`s^{-1}`
    unit: str
        Units in which half life is returned
        Possible unit entries:
            - 's' for seconds
            - 'm' for minutes
            - 'h' for hours
            - 'd' for days
            - 'y' for years
            - '1e3y' for 1000 years
            - '1e6y'
            - '1e9y'
    """

    half_life_s = log(2)/decay_constant
    half_life = half_life_s/time_dic[unit]

    return half_life


def halflife_to_decay(half_life, unit):

    """Converts a half life into a decay constant. User can choose half life input units with the unit parameter.

    Parameters
    ----------
    half life: float
        Half life of the nuclide
    unit: str
        Units in which half life is entered
        Possible unit entries:
            - 's' for seconds
            - 'm' for minutes
            - 'h' for hours
            - 'd' for days
            - 'y' for years
            - '1e3y' for 1000 years
            - '1e6y'
            - '1e9y'
    """
    half_life_s = half_life*time_dic[unit]
    decay_constant = log(2)/half_life_s

    return decay_constant

def halflife_to_second(half_life, unit):

    """Converts a half life from specified units into seconds.

    Parameters
    ----------
    half life: float
        Half life of the nuclide
    unit: str
        Units in which half life is returned
        Possible unit entries:
            - 's' for seconds
            - 'm' for minutes
            - 'h' for hours
            - 'd' for days
            - 'y' for years
            - '1e3y' for 1000 years
            - '1e6y'
            - '1e9y'
    """

    half_life_s = half_life*time_dic[unit]

    return half_life_s


def is_int(s):

    """Check whether an object is an integer.
    
    Parameters
    ----------
    s: NA
        Object to be checked
    """
    try:
        int(s)
        return True
    except ValueError:
        return False

def is_zamid(string):

    """Check whether a string is a nuclide's z-a-m id.
    
    Parameters
    ----------
    string: str
    """
    statement = True

    for number in string:
        if not is_number(number):
            statement = False

    return statement

def is_name(string):

    """Check whether a string is a nuclide's name.
    
    Parameters
    ----------
    string: str
    """

    if '-' in string:
        return True
    else:
        return False

def is_list_redundant(l):

    """Check whether a list is redundant.
    
    Parameters
    ----------
    l: List
        List to check
    """
    result = False
    l_set = set(l)
    if len(l) > len(l_set):
        result = True

def get_list_redundant_elt(l):

    """Returns the redundant elements in a list.
    
    Parameters
    ----------
    l: List
        List to check
    """

    count_elt = []
    for elt in l:
        if elt in count_elt:
            redundant_elt.append(elt)
        else:
            count_elt.append(elt)

    return redundant_elt

def zamid_list_to_name_list(zamid_list):

    """Converts a list of nuclides' z-a-m ids into a list of names.
    
    Parameters
    ----------
    zamid_list: List of str
        List of z-a-m ids
    """

    name_list = []
    for zamid in zamid_list:
        name = zamid_to_name(zamid)
        name_list.append(name)

    return name_list

def name_list_to_zamid_list(name_list):

    """Converts a list of nuclides' names into a list of z-a-m ids.
    
    Parameters
    ----------
    name_list: List of str
        List of names
    """

    zamid_list = []
    for name in name_list:
        zamid = name_to_zamid(name)
        zamid_list.append(zamid)

    return zamid_list

def get_zamid_z(zamid):

    """Gets the atomic number (z) from a nuclide's z-a-m id.
    
    Parameters
    ----------
    zamid: str
        z-a-m id of a nuclide
    """

    z = int(zamid[:-4])

    return z

def get_name_z(name):

    """Gets the atomic number (z) from a nuclide's name.
    
    Parameters
    ----------
    name: str
        Name of a nuclide
    """
    zamid = name_to_zamid(name)
    z = int(zamid[:-4])

    return z

def get_zamid_a(zamid):

    """Gets the mass number (a) from a nuclide's z-a-m id.
    
    Parameters
    ----------
    zamid: str
        z-a-m id of a nuclide
    """
    a = int(zamid[-4:-1])

    return a

def get_zamid_n(zamid):

    """Gets the number of neutrons (n) from a nuclide's z-a-m id.
    
    Parameters
    ----------
    zamid: str
        z-a-m id of a nuclide
    """
    a = int(zamid[-4:-1])
    z = int(zamid[:-4])

    return a - z

def get_zamid_s(zamid):

    """Gets the state of a nuclide from a its z-a-m id. 1 = first excited state, 0 = ground state.
    
    Parameters
    ----------
    zamid: str
        z-a-m id of a nuclide
    """
    s = int(zamid[-1])

    return s

def zamid_to_name(zamid):

    """Converts a nuclide's z-a-m id into the nuclide's name.
    
    Parameters
    ----------
    zamid: str
        z-a-m id of a nuclide
    """
    dic = d.nuc_name_dic

    if len(zamid) == 5:
        nz = int(zamid[0:1])
        na = int(zamid[1:4])
        state = int(zamid[4])
    if len(zamid) == 6:
        nz = int(zamid[0:2])
        na = int(zamid[2:5])
        state = int(zamid[5])
    if len(zamid) == 7:
        nz = int(zamid[0:3])
        na = int(zamid[3:6])
        state = int(zamid[6])

    if state == 0:
        nuc_name = '{}-{}'.format(d.nuc_zz_dic[nz], na) 
    else:
        nuc_name = '{}-{}*'.format(d.nuc_zz_dic[nz], na)

    return nuc_name

def name_to_zamid(name):

    """Converts a nuclide's name into the nuclide's z-a-m id.
    
    Parameters
    ----------
    name: str
        Name of a nuclide
    """
    dic = d.nuc_name_dic

    elt_name = name.split('-')[0]
    na = int(name.split('-')[1].replace('*',''))
    if '*' in name:
        state = 1
    else:
        state = 0
    zzaaam = 10000*d.nuc_name_dic[elt_name] + na*10 + state
    zamid = str(zzaaam)

    return zamid


def get_hm(passlist, hm_vol):

    """Gets the mass of heavy metal in a Passlist object.
    
    Parameters
    ----------
    passlist: onix.Passlist
    hm_vol: float
        Volume of the region containing heavy metal
    """
    hmmd = 0 # Heavy Metal Mass Density
    for i in passlist.passport_list:
        if int(i.zamid[:-4]) > 90:
            hmmd += i.current_dens*1E+24*i.mass/NA

    hm = hmmd*hm_vol

    return hm

def get_nucl_atomic_mass(nucl):

    """Gets the atomic mass of a nuclide (in grams).
    
    Parameters
    ----------
    nucl: str
        Name of a nuclide
    """
    zamid = name_to_zamid(nucl)
    zaid = zamid[:-1]
    if zaid in d.default_atm_mass_lib:
        M = d.default_atm_mass_lib[zaid]
    else:
        M = int(get_zamid_a(zamid))

    return M

def convert_mass_to_atom(mass, nuclide):

    """Converts the mass quantity (in grams) of a given nuclide type into number of atoms.
    
    Parameters
    ----------
    mass: float
        Mass quantity of the nuclide species
    nuclide: str
        Name of a nuclide
    """
    if is_name(nuclide):
        zamid = name_to_zamid(nuclide)
        zaid =zamid[:-1]
    else:
        zamid = nuclide
        zaid =zamid[:-1]

    molar_mass = d.default_atm_mass_lib[zaid]

    atom = mass*NA/molar_mass

    return atom

def convert_atom_to_mass(atom, nuclide):

    """Converts the quantity of a given nuclide from number of atoms into mass (in grams).
    
    Parameters
    ----------
    atom: float
        Number of atoms of the nuclide
    nuclide: str
        Name of a nuclide
    """
    if is_name(nuclide):
        zamid = name_to_zamid(nuclide)
        zaid =zamid[:-1]
    else:
        zamid = nuclide
        zaid =zamid[:-1]

    molar_mass = d.default_atm_mass_lib[zaid]

    mass = atom*molar_mass/NA

    return mass

def get_bu_sec_conv_factor(vol, ihm):
    """Computes the factor that converts seconds into burnup units (MWd/kg) from a given volume and a Initial Heavy Metal mass (IHM).

    Parameters
    ----------
    vol: float
        Volume of the region
    ihm: float
        Initial Heavy Metal mass of the region
    """
    bu_sec_conv_factor = vol*1e-3/(ihm*24*3600) # Unit in L/g

    return bu_sec_conv_factor

def get_keylist_from_dict(dict):

    keylist = list(dict.keys())

    return keylist

def get_decay_nucl(decay_a_lib):
    """Gets the list of nuclides from a decay dictionnary.

    Parameters
    ----------
    decay_a_lib: dict
        A decay dictionnary
    """
    decay_nucl = []
    for zamid in decay_a_lib:
        decay_nucl.append(zamid)

    return decay_nucl

def get_xs_nucl(xs_lib):
    """Gets the list of nuclides from a cross section dictionnary.

    Parameters
    ----------
    xs_lib: dict
        A cross section dictionnary
    """
    xs_nucl = []
    for zamid in xs_lib:
        xs_nucl.append(zamid)

    return xs_nucl

def get_fy_nucl(fy_lib):
    """Gets the list of fission products from a fission yield dictionnary.

    Parameters
    ----------
    fy_lib: dict
        A fission yield dictionnary
    """
    fy_nucl = []
    for zamid in fy_lib:
        fy_nucl.append(zamid)

    return fy_nucl

def get_all_nucl(list_of_dict):

    unfiltered_list = []
    for dictionary in list_of_dict:
        unfiltered_list += get_keylist_from_dict(dictionary)

    all_nucl = list(set(unfiltered_list))

    return all_nucl

def is_lista_in_listb(lista, listb):

    """Check whether elements from a list (lista) are all contained in another list (listb).

    Parameters
    ----------
    lista: List
    listb: List
    """
    result = all(elem in listb for elem in lista)

    return result

def get_fy_parent_nucl(fy_lib):

    """Gets the list of fission parents from a fission yield dictionnary.

    Parameters
    ----------
    fy_lib: dict
        A fission yield dictionnary
    """
    fy_nucl = get_fy_nucl(fy_lib)
    fy_parent = []
    sample_zamid = fy_nucl[0]
    sample = fy_lib[sample_zamid]

    for fission_parent in sample:
        fy_parent.append(fission_parent)

    return fy_parent


def get_cell_folder_path(file_name, *dir_path):

    folder_name = '{}_cell'.format(file_name)

    if not dir_path:
        dir_path = os.getcwd()

    folder_path = dir_path + '/' + folder_name

    return folder_path

def gen_cell_folder(name, *dir_path):

    if dir_path:
        folder_path = get_cell_folder_path(name, dir_path)
    elif not dir_path:
        folder_path = get_cell_folder_path(name)

    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)

    os.makedirs(folder_path)

def get_folder_path(folder_name, *dir_path):

    if not dir_path:
        dir_path = os.getcwd()

    folder_path = dir_path + '/' + folder_name

    return folder_path

def gen_folder(folder_name, *dir_path):

    if dir_path:
        folder_path = get_folder_path(folder_name, dir_path)
    elif not dir_path:
        folder_path = get_folder_path(folder_name)

    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)

    os.makedirs(folder_path)

# Convert a cell dictionary into a cell list.
# Dict keys must be cell IDs
# List is ordered according to cells ID
def cell_dict_to_cell_list(cell_dict):

    cell_list = []
    for ID in cell_dict:
        cell_list.append(cell_dict[ID])

    changed = True
    while changed:
        changed = False
        for i in range(len(cell_list) - 1):
            cell_id0 = cell_list[i].id
            cell_id1 = cell_list[i+1].id
            if cell_id0 > cell_id1:
                cell_list[i], cell_list[i+1] = cell_list[i+1], cell_list[i]
                changed = True


    return cell_list

def is_number(s):
    """Check whether s is a float.

    Parameters
    ----------
    s: NA
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def openmc_name_to_onix_name(name):

    """Converts a nuclide's name written with OpenMC format ('U235_m1') into ONIX format ('U-235*').

    **Note**: ONIX will adopt OpenMC name format in the next release.

    Parameters
    ----------
    name: str
        Name of a nuclide in OpenMC format ('U235_m1')
    """
    i = 0
    while not is_number(name[i]):
        i += 1

    am = name[i:]
    am = am.replace('_m1', '*')
    am = am.replace('m', '*') # used in jeff3.3
    am = am.replace('n', '*') # used in jeff3.3
    onix_name = name[:i] + '-' + am



    return onix_name

def onix_name_to_openmc_name(name):

    """Converts a nuclide's name written with ONIX format ('U-235*') into OpenMC format ('U235_m1').

    **Note**: ONIX will adopt OpenMC name format in the next release.

    Parameters
    ----------
    name: str
        Name of a nuclide in ONIX format ('U-235*')
    """
    openmc_name = name.replace('-', '')
    openmc_name = openmc_name.replace('*','_m1')

    return openmc_name

def bu_namelist_to_mc_namelist(name_list):

    """Converts a list of nuclides' names written with ONIX format ('U-235*') into OpenMC format ('U235_m1').

    **Note**: ONIX will adopt OpenMC name format in the next release.

    Parameters
    ----------
    name_list: List of str
        List of nuclides' names in ONIX format ('U-235*')
    """
    new_name_list = []
    for name in name_list:
        new_name = onix_name_to_openmc_name(name)
        new_name_list.append(new_name)

    return new_name_list

def mc_namelist_to_bu_namelist(name_list):

    """Converts a list of nuclides' names written with OpenMC format ('U235_m1') into ONIX format ('U-235*').

    **Note**: ONIX will adopt OpenMC name format in the next release.

    Parameters
    ----------
    name: List of str
        List of nuclides' names in OpenMC format ('U235_m1')
    """
    new_name_list = []
    for name in name_list:
        new_name = openmc_name_to_onix_name(name)
        new_name_list.append(new_name)

    return new_name_list

def order_nuclide_per_z(nucl_list):

    """Orders a list of nuclides' z-a-m ids according to their atomic number (z).

    Parameters
    ----------
    nucl_list: List of str
        List of nuclides' z-a-m ids
    """
    changed = True
    while changed:
        changed = False
        for i in range(len(nucl_list) - 1):
            zamid0 = int(nucl_list[i])
            zamid1 = int(nucl_list[i+1])

            if zamid0 > zamid1:
                nucl_list[i], nucl_list[i+1] = nucl_list[i+1], nucl_list[i]
                changed = True

    return nucl_list

def order_nuclide_name_per_z(nucl_name_list):

    """Orders a list of nuclides' names (in OpenMC format) according to their atomic number (z).

    Parameters
    ----------
    nucl_list: List of str
        List of nuclides' names
    """
    nucl_name_list_old_format = mc_namelist_to_bu_namelist(nucl_name_list)
    zamid_list = name_list_to_zamid_list(nucl_name_list_old_format)
    ordered_zamid_list = order_nuclide_per_z(zamid_list)
    ordered_nucl_name_list_old_format = zamid_list_to_name_list(ordered_zamid_list)
    ordered_nucl_name_list = bu_namelist_to_mc_namelist(ordered_nucl_name_list_old_format)

    return ordered_nucl_name_list

def order_nuclide_per_a(nucl_list):

    """Orders a list of nuclides' z-a-m ids according to their mass number (a).

    Parameters
    ----------
    nucl_list: List of str
        List of nuclides' z-a-m ids
    """
    changed = True
    while changed:
        changed = False
        for i in range(len(nucl_list) - 1):
            zamid0 = nucl_list[i]
            zamid1 = nucl_list[i+1]


            if zamid0 > zamid1:
                nucl_list[i], nucl_list[i+1] = nucl_list[i+1], nucl_list[i]
                changed = True

    return nucl_list


# set the colormap and centre the colorbar
class MidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))



def get_openmc_xs_nucl_list(path_to_xs_xml):
    """Returns the list of nuclides for which there are cross section data in a specified HDF5 cross section library directory (produced with OpenMC).

    Parameters
    ----------
    path_to_xs_xml: str
        Path to the cross_sections.xml file in a cross section library
    """
    #path_to_xs_xml = os.environ['OPENMC_CROSS_SECTIONS']
    #path_to_xs_xml = '/home/julien/Open-Burnup.dev/ENDFVIII_cross_sections/cross_sections.xml'
    MC_XS_nucl_list = []

    tree = ET.parse(path_to_xs_xml)
    root = tree.getroot()

    for child in root:
        if child.attrib['type'] == 'neutron':
            MC_XS_nucl_list.append(child.attrib['materials'])

    # Remove trouble makers     
    # For some reason, OpenMC can't find these nuclides in jeff lib at 800K
    MC_XS_nucl_list.remove('Cu63')
    MC_XS_nucl_list.remove('Cu65')
    MC_XS_nucl_list.remove('Mn55')

    try:
        MC_XS_nucl_list.remove('C0')
    except ValueError:
        pass
    try:
        MC_XS_nucl_list.remove('V0')
    except ValueError:
        pass
    try:
        MC_XS_nucl_list.remove('Zn0')
    except ValueError:
        pass

    return MC_XS_nucl_list


def convert_spectrum_to_janis_weighting_format(path_to_simulation, bucell, BU):
    """Converts a neutron spectrum output file computed in an ONIX simulations to a format suitable to be used to fold cross sections in JANIS4.1.

    Parameters
    ----------
    path_to_simulation: str
        Path to an ONIX simulation's directory with a neutron spectrum file
    bucell: str
        Name of a BUCell
    BU: float
        Burnup level at which the spectrum is to be taken
    """
    path = path_to_simulation +'/output_summary/{}_flux_spectrum'.format(bucell)
    spectrum_file = open(path)

    lines = spectrum_file.readlines()

    #Energy bins
    energy_bins = lines[0].split()[1:]

    #Energy mid points
    energy_mid_points = lines[2].split()[1:]

    # The data starts at the 7th line
    for line in lines[6:]:
        line = line.split()
        if float(line[1]) == BU:
            spectrum_lethargy = line[3:]

    spectrum = [float(x)/float(y) for x,y in zip(spectrum_lethargy, energy_mid_points)]

    energy_bin_file = open('energy_bin.gst', 'w')

    txt_bin = 'neutron group structure......anl 299 group\n'
    for i in range(len(energy_bins)-1):
        Emin = energy_bins[i]
        Emax = energy_bins[i+1]
        txt_bin += '{} {} {}\n'.format(i+1, Emin, Emax)

    energy_bin_file.write(txt_bin)
    energy_bin_file.close()

    spectrum_file = open('spectrum.txt', 'w')
    txt_spectrum = ''
    for i in range(len(spectrum)):
        E = energy_mid_points[i]
        flux = spectrum[i]
        txt_spectrum += '{} {}\n'.format(E, flux)

    spectrum_file.write(txt_spectrum)
    spectrum_file.close()


def get_zamid_natural_abundance(zamid):

    """Gets the natural abundance of a nuclide (values from 0 to 1).

    Parameters
    ----------
    zamid: str
        z-a-m id of a nuclide
    """
    name_old_format = zamid_to_name(zamid)
    name_new_format = onix_name_to_openmc_name(name_old_format)
    nat_abun_dict = d.NATURAL_ABUNDANCE
    if name_new_format in nat_abun_dict:
        nat_abun = d.NATURAL_ABUNDANCE[name_new_format]
    else:
        nat_abun = 0.0

    return nat_abun

def get_name_natural_abundance(name):

    """Gets the natural abundance of a nuclide (values from 0 to 1).

    Parameters
    ----------
    name: str
        Name of a nuclide
    """
    name_new_format = onix_name_to_openmc_name(name)
    nat_abun_dict = d.NATURAL_ABUNDANCE
    if name_new_format in nat_abun_dict:
        nat_abun = d.NATURAL_ABUNDANCE[name_new_format]
    else:
        nat_abun = 0.0

    return nat_abun

def find_zamid_precursor(zamid, reaction):

    """Finds the precursor of a nuclide via a specified reaction (except fission reactions).
    
    Parameters
    ----------
    zamid: str
        z-a-m id of a nuclide
    reaction: str
        Name of the reaction
        Possible names:
            - (n,gamma)
            - (n,2n)
            - (n,3n)
            - (n,p)
            - (n,a)
            - (n,t)

    """
    xs_prod_fromS_toS = d.xs_prod_fromS_toS
    zamid_shift = xs_prod_fromS_toS[reaction]
    precursor_zamid = int(zamid) - 10000*zamid_shift[0] - 10*zamid_shift[1] - zamid_shift[2]

    return str(precursor_zamid)


def smooth_triangle(data, degree, dropVals=False):
    triangle=np.array(list(range(degree)) + [degree] + list(range(degree)[::-1])) + 1
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(sum(point)/sum(triangle))
    if dropVals:
        return smoothed
    #smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    # while len(smoothed) < len(data):
    #     smoothed.append(smoothed[-1])

        
    return smoothed

def moving_average(data, window):

    weights = np.repeat(1.0, window)/window
    smas = np.convolve(data, weights, 'valid')

    print (smas)
    print (data[:int(window/2)])

    return list(data[:int(window/2)]) + list(smas) + list(data[-int(window/2):-1])

def read_BUCell_vol(path, cell):

    path_to_parameters = path +'/system_parameters'

    parameter_file = open(path_to_parameters)

    lines = parameter_file.readlines()

    search = 'BuCell'

    for line in lines:
        print (line)
        if line == 'BuCell {}\n'.format(cell):
            search = 'volume'

        if search == 'volume':
            if line.split()[0] == 'Volume':
                vol = line.split()[3]
                break

    return float(vol)

# linear interpolation between two points
def interpolation_between_two_points(pair1, pair2, x):
    """Linearly interpolates between two points to find the ordinate to a given abscissa value (x).

    Parameters
    ----------
    pair1: List of two floats
    pair2: List of two floats
    x: float
        Abscissa of the point for which the ordinate is calculated
    """
    a = (pair2[1] - pair1[1])/(pair2[0] - pair1[0])
    b = (pair1[1]*pair2[0] - pair1[0]*pair2[1])/(pair2[0] - pair1[0])

    y = a*x+b

    return y


class Empty_argument(Exception):
    """Raise when the user calls decay_halflife_conv without entering any argument """
    pass