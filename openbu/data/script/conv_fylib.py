# Convert the Origen xs file into an OpenBU fy file
from openbu.passport import passport

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



path_to_orilib = '/home/julien/origen22/libs/pwrue.lib'

ori_lib = open(path_to_orilib)
ori_line = ori_lib.readlines()

fy_lib = open('/home/julien/Open-Burnup.dev/openbu/data/fy_lib', 'w')

fathers = ['902320','922330','922350','922380', '942390','942410','962450','982490']

fp_data = []


fp_index = 0
for i in range(len(ori_line)):
    line = ori_line[i]
    if float(line.split()[0]) == 606:
        if is_int(line.split()[1]) == True and line.split()[8] == '1.0':
            fp_data.append([0.0]*9)
            zamid = line.split()[1]
            fp_data[fp_index][0] = zamid
            line2 = ori_line[i+1]
            for j in range(1,9):
                fp_data[fp_index][j] = float(line2.split()[j])
            fp_index += 1



fp_txt = '\n\n--- Fission Products Yields ---\n\n'
fp_txt += '===================================================================================\n'
fp_txt += '||  Nuclide Name  ||    ID    || Father ||     Value     ||  Uncertainty  ||\n'
fp_txt += '===================================================================================\n'
for i in range(len(fp_data)):
    zamid = fp_data[i][0]
    passpt = passport(zamid)
    name = passpt.nuc_name
    for j in range(8):
        if j == 0:
            fp_txt +='        {}         {}    {}           {}            0\n'.format(name, zamid, fathers[j], fp_data[i][1])
        else:
            fp_txt +='                       {}    {}           {}            0\n'.format(zamid, fathers[j], fp_data[i][j+1])

txt = fp_txt

fy_lib.write(txt)
