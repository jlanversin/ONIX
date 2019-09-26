from openbu.passport import Passport
import openbu.utils as utils
import openbu.data as d
import math as m



# Read old format origen decay lib

#ori_lib = '/home/julien/Open-Burnup.dev/openbu/data/other_libs/origen2.2/decay_lib'
ori_lib ='/home/julien/Open-Burnup.dev/test/argonne/decay_lib_argonne'

with open(ori_lib, 'r') as decay_file:
    line = decay_file.readlines()
    dic_list = {}
    r = 0
    for l in line:
        if r == 0:
            if len(l.split()) > 1 and l.split()[0].split('-')[0] in d.nuc_name_dic:
                zamid = l.split()[1]
                unit = l.split()[2]
                val = l.split()[3]
                dic_list[zamid] = {'unit':val}
                r = r + 1
        elif r > 0 and r < 6:
            decay = l.split()[1]
            val = l.split()[2]
            dic_list[zamid][decay] = val
            r = r + 1
        elif r == 6 :
            decay = l.split()[1]
            val = l.split()[2]
            dic_list[zamid][decay] = val
            r = 0
    
decay_b_dic = {}
for i in dic_list:
    dic = dic_list[i]
    decay_b_dic[i] = {}
    if dic['unit'] == 'stable':
        for j in range(9):
            decay_b_dic[i][d.decay_key_b[j]] = d.stable_dic_b[j]
            
    else:        
        hl_s = float(d.time_dic[dic['unit']]*float(dic['half-life']))
        total_decay = m.log(2)/hl_s

        # WARNING, while the X in origen are fraction of their respective decay (betanegX is shown as fraction of betaneg), here it is shown as fraction of total decay

        decay_b_dic[i]['unit'] = dic['unit']
        decay_b_dic[i]['half-life'] = float(dic['half-life'])
        decay_b_dic[i]['total'] = total_decay
        rest = (1 - float(dic['betapos']) -float(dic['alpha']) - float(dic['gamma'])) # The rest, i.e., what should be betaneg and betanegX
        if rest < 0: # Because ORIGEN is messed up and sometimes, the addition of the other reactions is bigger than one
            rest = 0
        decay_b_dic[i]['betaneg'] = rest*(1 - float(dic['betanegX']))
        decay_b_dic[i]['betanegX'] = rest*float(dic['betanegX'])
        decay_b_dic[i]['betapos'] = float(dic['betapos'])
        decay_b_dic[i]['betaposX'] = float(dic['betapos'])*float(dic['betaposX'])
        decay_b_dic[i]['alpha'] = float(dic['alpha'])
        decay_b_dic[i]['gamma'] = float(dic['gamma'])




nucl_list = []
for nucl in decay_b_dic:
	print (nucl)
	nucl_list.append(nucl)

order_nucl_list = utils.order_nuclide_per_z(nucl_list)

new_format = open('/home/julien/Open-Burnup.dev/openbu/data/other_libs/argonne/decay_lib_new_format', 'w')








#time_dic = {'s': 1, 'm': 60, 'h':3600, 'd': 24*3600, 'y': 24*3600*365.25, '1e3y':1e3*24*3600*365.25, '1e6y':1e6*24*3600*365.25, '1e9y':1e9*24*3600*365.25}


txt = '===================================================================================\n'
txt += '||  Nuclide Name  ||    ID    || Reaction Type ||   Branching   ||  Uncertainty  ||\n'
txt += '===================================================================================\n'
for i in range(len(order_nucl_list)):
	zamid = order_nucl_list[i]
	passpt = Passport(zamid)
	name = passpt.name
	decay_dict = decay_b_dic[zamid]
	#unit = decay_dict[zamid][1]['unit']
	unit = decay_dict['unit']
	#half_life = decay_dict[zamid][1]['half_life']
	if unit != 'n/a':
		half_life = decay_dict['half-life']*d.time_dic[unit]
	#half_life_error = decay_dict[zamid][1]['half_life_error']

	if unit == 'n/a':
		txt += '{:^19}{:^12}   {:<17}\n'.format(name, zamid, 'stable')
	# else:
	# 	#txt += '{:^19}{:^12}   {:<17}{:<17}\n'.format(name, zamid, 'unit',unit)
	# 	txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format(name, zamid, 'half-life', float(half_life), float(half_life_error))
	# 	for reac in reac_dict:
	# 		txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format('', zamid, reac, float(reac_dict[reac]['br']), float(reac_dict[reac]['br_error']))

	# txt += {}
	# txt += ''
	# txt +='        {}         {}    unit           {}            0\n'.format(name, zamid, time_dic[act_data[i][1]])
	else:
		txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format(name, zamid, 'half-life', float(half_life), 0.0)
		if float(decay_dict['betaneg']) != 0.0:
			txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format('', zamid, 'betaneg', float(decay_dict['betaneg']), 0.0)
		if float(decay_dict['betanegX']) != 0.0:
			txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format('', zamid, 'betanegX', float(decay_dict['betanegX']), 0.0)
		if float(decay_dict['betapos']) != 0.0:
			txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format('', zamid, 'betapos', float(decay_dict['betapos']), 0.0)
		if float(decay_dict['betaposX']) != 0.0:
			txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format('', zamid, 'betaposX', float(decay_dict['betaposX']), 0.0)
		if float(decay_dict['alpha']) != 0.0:
			txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format('', zamid, 'alpha', float(decay_dict['alpha']), 0.0)
		if float(decay_dict['gamma']) != 0.0:
			txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format('', zamid, 'gamma', float(decay_dict['gamma']), 0.0)
	# txt +='                      {}     betanegX                  {}            0\n'.format(zamid,act_data[i][3])
	# txt +='                      {}     betapos              {}            0\n'.format(zamid,act_data[i][4])
	# txt +='                      {}     betaposX          {}            0\n'.format(zamid,act_data[i][5])
	# txt +='                      {}     alpha                 {}            0\n'.format(zamid,act_data[i][6])
	# txt +='                      {}     gamma              {}            0\n'.format(zamid,act_data[i][7])



new_format.write(txt)