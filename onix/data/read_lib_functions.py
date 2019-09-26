import os
import math as m
import numpy as np
from .list_and_dict import *
import openmc

def read_mass_lib(mass_lib_path):

    mass_list = {}
    with open(mass_lib_path, 'r') as atm_mass_file:
        line = atm_mass_file.readlines()
        for i in range(40, 3352):
            if line[i][0] == ' ':
                nz, na = int(line[i].split()[2]), int(line[i].split()[3])
            elif line[i][0] == '0':
                nz, na = int(line[i].split()[3]), int(line[i].split()[4])

            zaid = str(1000*nz + na)
            _mass = line[i].split()[-3:-1]
            # This part is to convert the weird number format of mass.mass12 into a normal format
            mass_deci = _mass[1].replace('.','').replace('#','')
            #  mass_deci = '{}{}'.format(_mass[1].split('.')[0], _mass[1].split('.')[1])
            mass = float('{}.{}'.format(_mass[0], mass_deci))
            mass_list[zaid] = mass

    return mass_list

def read_decay_lib(decay_lib_path):

    with open(decay_lib_path, 'r') as decay_file:
        lines = decay_file.readlines()
        dict_list = {}

        previous_zamid = lines[3].split()[1]
        new_nuclide = 'yes'
        # Data starts at line 4
        count = 0
        for i in range(len(lines)-3):
        #for line in lines[3:]:
            line = lines[i+3]

            line = line.split()

            # This part check if the new line is about a new nuclide
            if count != 0:
                col_1_val = line[0]
                if col_1_val != previous_zamid:
                    new_nuclide = 'yes'


            if new_nuclide == 'yes':
                if count != 0:
                    dict_list[previous_zamid] = data_dict
                data_dict = {}
                new_nuclide = 'no'
                zamid = line[1]
                previous_zamid = zamid
                if line[2] == 'stable':
                    data_dict = 'stable'
                else:
                    data_dict['half-life'] = float(line[3])
                    total_decay = m.log(2)/float(line[3])
                    data_dict['total decay'] = total_decay

            else:
                data_dict[line[1]] = float(line[2])

            # If this is the last line, we need to add the last dictionrary
            if i == len(lines)-4:
                dict_list[previous_zamid] = data_dict


            count += 1      

    return dict_list


def conv_decay_b_a(decay_b):

    decay_a = {}
    for i in decay_b:
        if decay_b[i] == 'stable':
            decay_a[i] = 'stable'
        else:  
            decay_a[i] = {}
            decay_a[i]['half-life'] = decay_b[i]['half-life']
            decay_a[i]['total decay'] = decay_b[i]['total decay']
            for reac in decay_b[i]:
                if reac not in ['half-life', 'total decay']:
                    decay_a[i][reac] = decay_b[i][reac]*decay_b[i]['total decay']
        
    return decay_a

# Outdated version for decay. Adapted for ORIGEN lib format

# def read_decay_lib(decay_lib_path):

#     with open(decay_lib_path, 'r') as decay_file:
#         line = decay_file.readlines()
#         dic_list = {}
#         r = 0
#         for l in line:
#             if r == 0:
#                 if len(l.split()) > 1 and l.split()[0].split('-')[0] in nuc_name_dic:
#                     zamid = l.split()[1]
#                     unit = l.split()[2]
#                     val = l.split()[3]
#                     dic_list[zamid] = {'unit':val}
#                     r = r + 1
#             elif r > 0 and r < 6:
#                 decay = l.split()[1]
#                 val = l.split()[2]
#                 dic_list[zamid][decay] = val
#                 r = r + 1
#             elif r == 6 :
#                 decay = l.split()[1]
#                 val = l.split()[2]
#                 dic_list[zamid][decay] = val
#                 r = 0
        
#     decay_b_dic = {}
#     for i in dic_list:
#         dic = dic_list[i]
#         decay_b_dic[i] = {}
#         if dic['unit'] == 'stable':
#             for j in range(9):
#                 decay_b_dic[i][decay_key_b[j]] = stable_dic_b[j]
                
#         else:        
#             hl_s = float(time_dic[dic['unit']]*float(dic['half-life']))
#             total_decay = m.log(2)/hl_s

#             # WARNING, while the X in origen are fraction of their respective decay (betanegX is shown as fraction of betaneg), here it is shown as fraction of total decay

#             decay_b_dic[i]['unit'] = dic['unit']
#             decay_b_dic[i]['half-life'] = float(dic['half-life'])
#             decay_b_dic[i]['total'] = total_decay
#             rest = (1 - float(dic['betapos']) -float(dic['alpha']) - float(dic['gamma'])) # The rest, i.e., what should be betaneg and betanegX
#             if rest < 0: # Because ORIGEN is messed up and sometimes, the addition of the other reactions is bigger than one
#                 rest = 0
#             decay_b_dic[i]['betaneg'] = rest*(1 - float(dic['betanegX']))
#             decay_b_dic[i]['betanegX'] = rest*float(dic['betanegX'])
#             decay_b_dic[i]['betapos'] = float(dic['betapos'])
#             decay_b_dic[i]['betaposX'] = float(dic['betapos'])*float(dic['betaposX'])
#             decay_b_dic[i]['alpha'] = float(dic['alpha'])
#             decay_b_dic[i]['gamma'] = float(dic['gamma'])
        

#     return decay_b_dic

# def conv_decay_b_a(decay_b):

#     decay_a = {}
#     for i in decay_b:
#         decay_a[i] = {}
        
#         if decay_b[i]['half-life'] == 'stable':
#             for j in range(8):
#                 decay_a[i][decay_key_a[j]] = stable_dic_a[j]

#         else:
#             # decay_a[i]['half-life'] = float(decay_b[i]['half-life']*time_dic[decay_b[i]['unit']])
#             # decay_a[i]['total'] = decay_b[i]['total']
#             # decay_a[i]['betaneg'] = (1 - decay_b[i]['betapos'] -decay_b[i]['alpha'] - decay_b[i]['gamma'])*decay_b[i]['total']*(1 - decay_b[i]['betanegX'])
#             # decay_a[i]['betanegX'] = (1 - decay_b[i]['betapos'] - decay_b[i]['alpha'] - decay_b[i]['gamma'])*decay_b[i]['total']*decay_b[i]['betanegX']
#             # decay_a[i]['betapos'] = decay_b[i]['betapos']*decay_b[i]['total']
#             # decay_a[i]['betaposX'] = decay_b[i]['betaposX']*decay_a[i]['betapos']
#             # decay_a[i]['alpha'] = decay_b[i]['alpha']*decay_b[i]['total']
#             # decay_a[i]['gamma'] = decay_b[i]['gamma']*decay_b[i]['total']

#             decay_a[i]['half-life'] = float(decay_b[i]['half-life']*time_dic[decay_b[i]['unit']])
#             decay_a[i]['total'] = decay_b[i]['total']
#             decay_a[i]['betaneg'] = decay_b[i]['betaneg']*decay_b[i]['total']
#             decay_a[i]['betanegX'] = decay_b[i]['betanegX']*decay_b[i]['total']
#             decay_a[i]['betapos'] = decay_b[i]['betapos']*decay_b[i]['total']
#             decay_a[i]['betaposX'] = decay_b[i]['betaposX']*decay_b[i]['total']
#             decay_a[i]['alpha'] = decay_b[i]['alpha']*decay_b[i]['total']
#             decay_a[i]['gamma'] = decay_b[i]['gamma']*decay_b[i]['total']

#             # total = decay_a[i]['total']
#             # expo = np.floor(m.log10(total))
#             # if i == '621460':
#             #     print(expo)
#             #     print((decay_a[i]))
        
#     return decay_a


def read_xs_lib(xs_lib_path):

    with open(xs_lib_path, 'r') as xs_file:
        line = xs_file.readlines()
        r = 0
        xs_dic = {}
        for l in line:
            # if r == 0:
            #     if len(l.split()) > 1 and l.split()[0].split('-')[0] in nuc_name_dic:
            #         zamid = l.split()[1]
            #         xs = l.split()[2]
            #         val = [float(k) for k in l.split()[3:5]]
            #         removal_val += removal_val + val[0]
            #         xs_dic[zamid] = {xs:val}
            #         summ[0] = summ[0] + val[0]
            #         r = 1
            # elif r > 0:
            #     if l.split()[1] == 'removal':
            #         r = 0
            #     xs = l.split()[1]
            #     val = [float(k) for k in l.split()[2:4]]
            #     xs_dic[zamid][xs] = val

            if len(l.split()) > 1 and l.split()[0].split('-')[0] in nuc_name_dic:
                zamid = l.split()[1]
                xs = l.split()[2]
                val = [float(k) for k in l.split()[3:5]]
                xs_dic[zamid] = {xs:val}
                r = 1
            elif len(l.split()) > 1 and r == 1:
                xs = l.split()[1]
                val = [float(k) for k in l.split()[2:4]]
                xs_dic[zamid][xs] = val
            else:
                r=0

    # Add the removal entry and values to the dic
    for zamid in xs_dic:
        removal_val = 0
        for xs in xs_dic[zamid]:
            removal_val += xs_dic[zamid][xs][0]
        xs_dic[zamid]['removal'] = [0, 0]
        xs_dic[zamid]['removal'][0] = removal_val

    # Add the removal uncertainty values to the dic
    for zamid in xs_dic:
        removal_unc = 0
        if xs_dic[zamid]['removal'][0] != 0: # if removal = 0 then no need to go through the following
            for xs in xs_dic[zamid]:
                val_weight = xs_dic[zamid][xs][0]/xs_dic[zamid]['removal'][0]
                removal_unc += xs_dic[zamid][xs][1]*val_weight
            xs_dic[zamid]['removal'][1] = removal_unc

    return xs_dic


def read_fy_lib(fy_lib_path):

    with open(fy_lib_path, 'r') as fy_file:
        line = fy_file.readlines()
        r = 0
        read =  0
        fy_dic = {}
        for i in range(0,len(line)):
            l = line[i]
            if l == '--- Fission Products Yields ---\n':
                read = 1
            if read == 1:
                # This is the first line for the nuclide
                if r == 0:
                    if len(l.split()) > 1 and l.split()[0].split('-')[0] in nuc_name_dic:
                        zamid = l.split()[1]
                        father_nuc = l.split()[2]
                        fy = [float(k) for k in l.split()[3:5]]
                        fy_dic[zamid] = {father_nuc:fy}
                        r = 1
                # This is not the first line for the nuclide
                elif r > 0:
                    father_nuc = l.split()[1]
                    fy = [float(k) for k in l.split()[2:4]]
                    fy_dic[zamid][father_nuc] = fy
                    # If this the last line of the file, exit of the loop
                    if i == len(line) - 1:
                        break
                    # If next line is a new nuclide and is not empty, reset r to 0
                    elif len(l.split()) > 1 and line[i+1].split()[0].split('-')[0] in nuc_name_dic:
                        r = 0

    return fy_dic


def default_xs_mat_from_Btxt():

    Btxt_rel_path = '/default_libs/Btxt'
    dir_path = os.path.dirname(__file__)
    Btxt_path = dir_path + Btxt_rel_path
    Btxt = open(Btxt_path, 'r')
    B_line = Btxt.readlines()

    N = len(B_line)
    default_xs_mat = np.zeros((N,N))
    for line in B_line:
        data = line.split('|')[1]
        row = int(data.split(':')[0])
        col_data = data.split(':')[1].split(',')
        col_data.remove('\n') #The format of Ctxt and Btxt contains a last coma that needs to be removed
        for data in col_data:
            col = int(data.split()[0])
            val = float(data.split()[1])
            default_xs_mat[row][col] =  val

    return default_xs_mat

def default_decay_mat_from_Ctxt():

    Ctxt_rel_path = '/default_libs/Ctxt'
    dir_path = os.path.dirname(__file__)
    Ctxt_path = dir_path + Ctxt_rel_path
    Ctxt = open(Ctxt_path, 'r')
    C_line = Ctxt.readlines()

    N = len(C_line)
    default_decay_mat = np.zeros((N,N))
    for line in C_line:
        data = line.split('|')[1]
        row = int(data.split(':')[0])
        col_data = data.split(':')[1].split(',')
        col_data.remove('\n') #The format of Ctxt and Btxt contains a last coma that needs to be removed
        for data in col_data:
            col = int(data.split()[0])
            val = float(data.split()[1])
            default_decay_mat[row][col] =  val

    return default_decay_mat

def default_nucl_list_from_txt():

    mattxt_rel_path = '/default_libs/Ctxt'
    dir_path = os.path.dirname(__file__)
    mattxt_path = dir_path + mattxt_rel_path
    mattxt = open(mattxt_path, 'r')
    mat_line = mattxt.readlines()

    N = len(mat_line)
    default_nucl_list = []
    for line in mat_line:
        zamid = line.split('|')[0]
        default_nucl_list.append(zamid)

    return default_nucl_list

# This might not be appropriate to have this function here as it needs OpenMC
def read_isomeric_data():

    # This command will find the absolute path of read_lib_fuctions.py
    # Since read_lib_fuctions.py is located in data, the default isomeric data is just in __file__path + '/isomeric_data/eaf-2010-multiplicities'
    __file__path = os.path.abspath(os.path.dirname(__file__))

    isomeric_data_path = __file__path+ '/isomeric_data/eaf-2010-multiplicities'

    # onix does not consider branching that goes to 2nd excited state
    file_name_list = [x for x in os.listdir(isomeric_data_path) if '_2.' not in x]

    isomeric_branching_dict = {}
    for file_name in file_name_list:
        name = file_name.replace('.csv', '')
        target_name = name.split('_')[0]
        target_state = name.split('_')[1]
        daughter_state = name.split('_')[3]
        if target_state == '1':
            nucl_name = target_name + '_m1'
        else:
            nucl_name = target_name
        if nucl_name not in isomeric_branching_dict:
            isomeric_branching_dict[nucl_name] = {}
        file = open(isomeric_data_path + '/' + file_name, 'r')
        lines = file.readlines()
        energy_grid = []
        data = []
        for line in lines[1:]:
            line = line.split(',')
            energy_grid.append(float(line[1]))
            data.append(float(line[2]))
        tabulated_data = openmc.data.Tabulated1D(energy_grid, data)
        isomeric_branching_dict[nucl_name][daughter_state] = tabulated_data
    
    # for some reason some parasitycal element are added to the dict (ex:.~lock.I129)
    # These lines remove them
    wrong_entries = []
    for key in isomeric_branching_dict:
        if 'lock' in key:
            wrong_entries.append(key)
    for key in wrong_entries:
        del isomeric_branching_dict[key]

    return isomeric_branching_dict