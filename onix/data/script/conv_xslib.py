# Convert the Origen xs file into an OpenBU decay file
from openbu.passport import Passport

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



#path_to_orilib = '/home/julien/origen22/libs/pwrue.lib'
path_to_orilib = '/home/julien/origen/libs/pwrpupu.lib'
#path_to_orilib = '/home/julien/princeton/origen/libs/pwrue.lib'

ori_lib = open(path_to_orilib)
ori_line = ori_lib.readlines()

#xs_lib = open('/home/julien/Open-Burnup.dev/openbu/data/xs_lib4', 'w')
xs_lib = open('/home/julien/Open-Burnup.dev/openbu/data/xs_lib_pupu', 'w')



#   2  922350  1  (time unit)   2.221E 16 (half-time) 0.0 (FBX)      0.0  (FPEC)      0.0   (FPECX)    1.000E 00 (FA) 0.0   (FIT)    
#   2                2.600E-09 0.0       4.418E 00 7.200E-01 4.000E-12 3.000E-05 

# 3 -> unit time
# 4 -> ht
# 5 -> FBX
# 6 -> FPEC
# 7 -> FPECX
# 8 -> FA

# alpha -> FA
# beta+ -> FPEC (FPECX % leads to excited state)
# beta- -> 1 - FPEC - FA (FBX % leads to excited state)



# Actinides

# act_data = {}
# fp_data = {}

# act_xs_key = ['(n,gamma)','(n,2n)','(n,3n)','fission','(n,gamma)X','(n,2n)X', 'removal']
# fp_xs_key = ['(n,gamma)','(n,2n)','(n,a)','(n,p)','(n,gamma)X','(n,2n)X', 'removal']

# for i in ori_line:
#     # Actinides
#     if float(i.split()[0]) == 605:
#         if is_int(i.split()[1]) == True:
#             zamid = i.split()[1]
#             act_data[zamid] = [0.0]*7
#             for j in range(6):
#                 act_data[zamid][j] = float(i.split()[j+2])
#             act_data[zamid][6] = sum(act_data[zamid][0:6])
#     # FPs
#     elif float(i.split()[0]) == 606:
#         if is_int(i.split()[1]) == True:
#             zamid = i.split()[1]
#             fp_data[zamid] = [0.0]*7
#             for j in range(6)
#:                 fp_data[zamid][j] = float(i.split()[j+2])
#             fp_data[zamid][6] = sum(fp_data[zamid][0:6])

avt_data = []
act_data = []
fp_data = []

avt_xs_key = ['(n,gamma)','(n,2n)','(n,a)','(n,p)','(n,gamma)X','(n,2n)X', 'removal']
act_xs_key = ['(n,gamma)','(n,2n)','(n,3n)','fission','(n,gamma)X','(n,2n)X', 'removal']
fp_xs_key = ['(n,gamma)','(n,2n)','(n,a)','(n,p)','(n,gamma)X','(n,2n)X', 'removal']

avt_index = 0
act_index = 0
fp_index = 0
for i in ori_line:
    # Actinides
    if float(i.split()[0]) == 211:
        if is_int(i.split()[1]) == True:
            act_data.append([0.0]*8)
            zamid = i.split()[1]
            act_data[act_index][0] = zamid
            for j in range(1,7):
                act_data[act_index][j] = float(i.split()[j+1])
            act_data[act_index][7] = sum(act_data[act_index][1:7])
            act_index += 1
    # FPs
    elif float(i.split()[0]) == 212:
        if is_int(i.split()[1]) == True:
            fp_data.append([0.0]*8)
            zamid = i.split()[1]
            fp_data[fp_index][0] = zamid
            for j in range(1,7):
                fp_data[fp_index][j] = float(i.split()[j+1])
            fp_data[fp_index][7] = sum(fp_data[fp_index][1:7])
            fp_index += 1

# A list that stores all the nuclides zamid of the FPs
fp_nucl = [None]*len(fp_data)
for i in range(len(fp_data)):
    fp_nucl[i] = fp_data[i][0]

print(fp_nucl)

for i in ori_line:
    # Activation products
    if float(i.split()[0]) == 210:
        if is_int(i.split()[1]) == True:
            if i.split()[1] in fp_nucl:
                continue
            avt_data.append([0.0]*8)
            zamid = i.split()[1]
            avt_data[avt_index][0] = zamid
            for j in range(1,7):
                avt_data[avt_index][j] = float(i.split()[j+1])
            avt_data[avt_index][7] = sum(avt_data[avt_index][1:7])
            avt_index += 1


# Writing actinide part

act_txt = '\n\n--- Actinides ---\n\n'
act_txt += '===================================================================================\n'
act_txt += '||  Nuclide Name  ||    ID    || Reaction Type ||     Value     ||  Uncertainty  ||\n'
act_txt += '===================================================================================\n'
for i in range(len(act_data)):
    zamid = act_data[i][0]
    passpt = Passport(zamid)
    name = passpt.name
    act_txt +='        {}         {}    (n,gamma)           {}            0\n'.format(name, zamid, act_data[i][1])
    act_txt +='                      {}     (n,2n)                  {}            0\n'.format(zamid,act_data[i][2])
    act_txt +='                      {}     (n,3n)                  {}            0\n'.format(zamid,act_data[i][3])
    act_txt +='                      {}     fission              {}            0\n'.format(zamid,act_data[i][4])
    act_txt +='                      {}     (n,gamma)X          {}            0\n'.format(zamid,act_data[i][5])
    act_txt +='                      {}     (n,2n)X                 {}            0\n'.format(zamid,act_data[i][6])
 #   act_txt +='                      {}     removal              {}            0\n'.format(zamid,act_data[i][7])

# Writing FP part

fp_txt = '\n\n--- Fission Products ---\n\n'
fp_txt += '===================================================================================\n'
fp_txt += '||  Nuclide Name  ||    ID    || Reaction Type ||     Value     ||  Uncertainty  ||\n'
fp_txt += '===================================================================================\n'
for i in range(len(fp_data)):
    zamid = fp_data[i][0]
    passpt = Passport(zamid)
    name = passpt.name
    fp_txt +='        {}         {}    (n,gamma)           {}            0\n'.format(name, zamid, fp_data[i][1])
    fp_txt +='                      {}     (n,2n)                  {}            0\n'.format(zamid,fp_data[i][2])
    fp_txt +='                      {}     (n,a)               {}            0\n'.format(zamid,fp_data[i][3])
    fp_txt +='                      {}     (n,p)                   {}            0\n'.format(zamid,fp_data[i][4])
    fp_txt +='                      {}     (n,gamma)X          {}            0\n'.format(zamid,fp_data[i][5])
    fp_txt +='                      {}     (n,2n)X                 {}            0\n'.format(zamid,fp_data[i][6])
 #   fp_txt +='                      {}     removal              {}            0\n'.format(zamid,fp_data[i][7])



avt_txt = '\n\n--- Activation Products ---\n\n'
avt_txt += '===================================================================================\n'
avt_txt += '||  Nuclide Name  ||    ID    || Reaction Type ||     Value     ||  Uncertainty  ||\n'
avt_txt += '===================================================================================\n'
for i in range(len(avt_data)):
    zamid = avt_data[i][0]
    passpt = Passport(zamid)
    name = passpt.name
    avt_txt +='        {}         {}    (n,gamma)           {}            0\n'.format(name, zamid, avt_data[i][1])
    avt_txt +='                      {}     (n,2n)                  {}            0\n'.format(zamid,avt_data[i][2])
    avt_txt +='                      {}     (n,a)               {}            0\n'.format(zamid,avt_data[i][3])
    avt_txt +='                      {}     (n,p)                   {}            0\n'.format(zamid,avt_data[i][4])
    avt_txt +='                      {}     (n,gamma)X          {}            0\n'.format(zamid,avt_data[i][5])
    avt_txt +='                      {}     (n,2n)X                 {}            0\n'.format(zamid,avt_data[i][6])
 #   avt_txt +='                      {}     removal              {}            0\n'.format(zamid,avt_data[i][7])


txt = act_txt + fp_txt + avt_txt

xs_lib.write(txt)

