
from onix.passport import Passport
import onix.utils as utils
import onix.data as data

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def find_E_index(s):
    length = len(s)
    j = 0
    index = []
    for i in range(length):
        if s[i] == 'E':
            index.append(i)
            j = j + 1
    return index

time_dic = {}

#fathers_zamid = ['902320','922330','922350','922380', '942390','942410','962450','982490']
#fathers_name = ['Th232', 'U233', 'U234', 'U235', 'U236', 'U238', 'Np237', 'Np238', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'Am241', 'Am243', 'Cm243', 'Cm244', 'Cm245']
fathers_name = ['Th227',
'Th229',
'Th232',
'Pa231',
'U232',
'U233',
'U234',
'U235',
'U236',
'U237',
'U238',
'Np237',
'Np238',
'Pu238',
'Pu239',
'Pu240',
'Pu241',
'Pu242',
'Am241',
'Am242_m1',
'Am243',
'Cm242',
'Cm243',
'Cm244',
'Cm245',
'Cm246',
'Cm248',
'Cf249',
'Cf251',
'Es254',
'Fm255']

father_name_fast = ['Th232','Pa231','U233','U233','U234','U235','U236','U237','U238','Np237','Np238','Pu238','Pu239','Pu240','Pu241','Pu242','Am241','Am243','Cm242','Cm243','Cm244','Cm246','Cm248'
]

fathers_zamid = [utils.name_to_zamid(utils.openmc_name_to_onix_name(x)) for x in father_name_fast]
reduced_nucl_set = data.reduced_nucl_set





# obu_name_list = utils.mc_namelist_to_bu_namelist(fathers_name)
# zamid_list = utils.name_list_to_zamid_list(obu_name_list)

#path_to_jeff33_lib = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/ENDFVIII_ind_fy_noheader.csv'
path_to_ENDF8_fast = '/home/julien/ONIX/ONIX/onix/data/other_libs/ENDFVIII/original/ENDFVIII_ind_fy_fast.csv'

# jeff33_lib = open(path_to_jeff33_lib)
# jeff33_line = jeff33_lib.readlines()

ENDF8_fast_lib = open(path_to_ENDF8_fast)
ENDF8_fast_line = ENDF8_fast_lib.readlines()

#obu_fy_lib = open('/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/obu_endfVIII_fy_lib', 'w')
onix_fy_fast_lib = open('/home/julien/ONIX/ONIX/onix/data/other_libs/ENDFVIII/original/onix_endfVIII_fy_fast_lib', 'w')


#fathers_zamid = ['902320','922330','922350','922380', '942390','942410','962450','982490']
#fathers_name = ['Th232', 'U233', 'U234', 'U235', 'U236', 'U238', 'Np237', 'Np238', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'Am241', 'Am243', 'Cm243', 'Cm244', 'Cm245']



fp_data = []


fp_index = 0
for i in range(len(ENDF8_fast_line)):
    line = ENDF8_fast_line[i]
    #fp_data.append([])
    line = line.split(';')
    new_line = ['0.0' if x in ['', '\n'] else x for x in line]
    fp_data.append(new_line)

    # zamid = line.split()[1]
    # fp_data[fp_index][0] = zamid
    # line2 = ori_line[i+1]
    # for j in range(1,9):
    #     fp_data[fp_index][j] = float(line2.split()[j])
    # fp_index += 1


fp_txt = '\n\n--- Fission Products Yields ---\n\n'
fp_txt += '===================================================================================\n'
fp_txt += '||  Nuclide Name  ||    ID    || Father ||     Value     ||  Uncertainty  ||\n'
fp_txt += '===================================================================================\n'
for i in range(len(fp_data)):
    print (i)
    name = fp_data[i][0]
    obu_name = utils.openmc_name_to_onix_name(name)
    print (fp_data[i])
    zamid = utils.name_to_zamid(obu_name)
    # if zamid not in reduced_nucl_set:
    #     continue
    count = 0
    for j in range(int((len(fp_data[i])-1)/2)):
        val = fp_data[i][2*j+1]
        uncertainty = fp_data[i][2*j+2].rstrip()

        if count == 0:
            fp_txt += '{:^19}'.format(obu_name)
            count = 1
        elif count == 1:
            fp_txt += '{:^19}'.format('')

        fp_txt += '{:^12}'.format(zamid)
        fp_txt += '{:^10}'.format(fathers_zamid[j])
        fp_txt += '{:^17.10E}'.format(float(val)*100)
        fp_txt += '{:^17.10E}'.format(float(uncertainty))
        fp_txt += '\n'


onix_fy_fast_lib.write(fp_txt)
