# Convert the Origen decay file into an OpenBU decay file
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
    for i in range(length-1):
        if s[i] == 'E' and s[i+1] == ' ':
            index.append(i)
            j = j + 1
    return index

time_dic = {1.0:'s', 2.0:'m', 3.0:'h', 4.0:'d', 5.0:'y', 6.0:'stable', 7.0:'1e3y', 8.0:'1e6y', 9.0:'1e9y'}

# In Princeton
#path_to_orilib = '/home/julien/origen22/libs/decay.lib'
# While in CHina
path_to_orilib = '/home/julien/princeton/origen/libs/decay.lib'

ori_lib = open(path_to_orilib)
ori_line = ori_lib.readlines()
decay_lib = open('/home/julien/Open-Burnup.dev/openbu/data/decay_lib', 'w')


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

avt_data = []
act_data = []
fp_data = []

avt_index = 0
act_index = 0
fp_index = 0

for i in ori_line:
    # Actinides
    if float(i.split()[0]) == 2:
        if is_int(i.split()[1]) == True:
            index = find_E_index(i)
            h = 0
            for k in range(len(index)):
                i = i[:index[k]+1-h] + i[index[k]+1+1-h:]
                h = h + 1

            act_data.append([0.0]*8)
            zamid = i.split()[1]
            act_data[act_index][0] = zamid
            for j in range(1,8):
                act_data[act_index][j] = float(i.split()[j+1])
            act_index += 1
    # FPs
    elif float(i.split()[0]) == 3:
        if is_int(i.split()[1]) == True:
            index = find_E_index(i)
            h = 0
            for k in range(len(index)):
                i = i[:index[k]+1-h] + i[index[k]+1+1-h:]
                h = h + 1

            fp_data.append([0.0]*8)
            zamid = i.split()[1]
            fp_data[fp_index][0] = zamid
            for j in range(1,8):
                fp_data[fp_index][j] = float(i.split()[j+1])
            fp_index += 1

# A list that stores all the nuclides zamid of the FPs and Act to check redonduncy in avt
fp_nucl = [None]*len(fp_data)
act_nucl = [None]*len(act_data)
for i in range(len(fp_data)):
    fp_nucl[i] = fp_data[i][0]
for i in range(len(act_data)):
    act_nucl[i] = act_data[i][0]

redundant_count = 0
for i in ori_line:
    # Activation products
    if float(i.split()[0]) == 1:
        if is_int(i.split()[1]) == True:
            if i.split()[1] in fp_nucl+act_nucl:
                redundant_count += 1
                continue
            index = find_E_index(i)
            h = 0
            for k in range(len(index)):
                i = i[:index[k]+1-h] + i[index[k]+1+1-h:]
                h = h + 1

            avt_data.append([0.0]*8)
            zamid = i.split()[1]
            avt_data[avt_index][0] = zamid
            for j in range(1,8):
                avt_data[avt_index][j] = float(i.split()[j+1])
            avt_index += 1

print(redundant_count)
print(avt_index)


act_txt = '\n\n--- Actinides ---\n\n'
act_txt += '===================================================================================\n'
act_txt += '||  Nuclide Name  ||    ID    || Reaction Type ||     Value     ||  Uncertainty  ||\n'
act_txt += '===================================================================================\n'
for i in range(len(act_data)):
    zamid = act_data[i][0]
    passpt = passport(zamid)
    name = passpt.name
    act_txt +='        {}         {}    unit           {}            0\n'.format(name, zamid, time_dic[act_data[i][1]])
    act_txt +='                      {}     half-life                  {}            0\n'.format(zamid,act_data[i][2])
    act_txt +='                      {}     betanegX                  {}            0\n'.format(zamid,act_data[i][3])
    act_txt +='                      {}     betapos              {}            0\n'.format(zamid,act_data[i][4])
    act_txt +='                      {}     betaposX          {}            0\n'.format(zamid,act_data[i][5])
    act_txt +='                      {}     alpha                 {}            0\n'.format(zamid,act_data[i][6])
    act_txt +='                      {}     gamma              {}            0\n'.format(zamid,act_data[i][7])

# Writing FP part

fp_txt = '\n\n--- Fission Products ---\n\n'
fp_txt += '===================================================================================\n'
fp_txt += '||  Nuclide Name  ||    ID    || Reaction Type ||     Value     ||  Uncertainty  ||\n'
fp_txt += '===================================================================================\n'
for i in range(len(fp_data)):
    zamid = fp_data[i][0]
    passpt = passport(zamid)
    name = passpt.name
    fp_txt +='        {}         {}    unit           {}            0\n'.format(name, zamid, time_dic[fp_data[i][1]])
    fp_txt +='                      {}     half-life                  {}            0\n'.format(zamid,fp_data[i][2])
    fp_txt +='                      {}     betanegX               {}            0\n'.format(zamid,fp_data[i][3])
    fp_txt +='                      {}     betapos                   {}            0\n'.format(zamid,fp_data[i][4])
    fp_txt +='                      {}     betaposX          {}            0\n'.format(zamid,fp_data[i][5])
    fp_txt +='                      {}     alpha                 {}            0\n'.format(zamid,fp_data[i][6])
    fp_txt +='                      {}     gamma              {}            0\n'.format(zamid,fp_data[i][7])



    avt_txt = '\n\n--- Activation Products ---\n\n'
avt_txt += '===================================================================================\n'
avt_txt += '||  Nuclide Name  ||    ID    || Reaction Type ||     Value     ||  Uncertainty  ||\n'
avt_txt += '===================================================================================\n'
for i in range(len(avt_data)):
    zamid = avt_data[i][0]
    passpt = passport(zamid)
    name = passpt.name
    avt_txt +='        {}         {}    unit           {}            0\n'.format(name, zamid, time_dic[avt_data[i][1]])
    avt_txt +='                      {}     half-life                  {}            0\n'.format(zamid,avt_data[i][2])
    avt_txt +='                      {}     betanegX               {}            0\n'.format(zamid,avt_data[i][3])
    avt_txt +='                      {}     betapos                   {}            0\n'.format(zamid,avt_data[i][4])
    avt_txt +='                      {}     betaposX          {}            0\n'.format(zamid,avt_data[i][5])
    avt_txt +='                      {}     alpha                 {}            0\n'.format(zamid,avt_data[i][6])
    avt_txt +='                      {}     gamma              {}            0\n'.format(zamid,avt_data[i][7])


txt = act_txt + fp_txt + avt_txt

decay_lib.write(txt)

