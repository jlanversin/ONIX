from openbu.passport import Passport
import openbu.data as data


column_title = ['line', 'zaid', 'br', 'br error', 'daughter s', 'energy', 'energy error', 'half life', 'half life error', 's', 'name', 'type']
reac_name_dict = {'IT': 'gamma', 'beta-': 'betaneg', 'ec/beta+': 'betapos', 'alpha': 'alpha', 'neutron': 'neutron', 'proton': 'proton',
'n':'neutron', 'p':'proton'}
state_dict = {'0.0':'', '0':'', '1.0':'X', '1':'X'}

reduced_nucl_set= data.reduced_nucl_set

# jeff33_file = open('/home/julien/Open-Burnup.dev/openbu/data/other_libs/jeff33/jeff3.3-decay.csv')
# OpenBU_format = open('/home/julien/Open-Burnup.dev/openbu/data/other_libs/jeff33/obu_jeff33_decay', 'w')

ENDF8_file = open('/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/endf-b-viii.0-decay-openmc.csv', 'r')
OpenBU_format = open('/home/julien/Open-Burnup.dev/openbu/data/other_libs/ENDFVIII/obu_ENDFVIII_decay_reduced', 'w')

lines = ENDF8_file.readlines()


decay_dict = {}
nucl_index = []
nucl_name = []
from_excited_zamid = []
to_excited_zamid = []

i = 0
for line in lines:
	data_dict = {}
	line = line.split(',')

	s = line[9]
	reac_type = line[11].replace('\n', '')
	daughter_s = line[4]

	# istopes in their second excited state are not considered
	if s == '2' or s == '3':
		i += 1
		continue

	# reaction to isotopes to lvl 2 and 3 are ignored
	if daughter_s in ['2', '2.0', '3', '3.0']:
		i += 1
		continue	

	# double reactions are not consideed
	if '(' in reac_type:
		i += 1
		continue

		# in ENDFVIII file format, double reactions are separated by a space
	if ' ' in reac_type:
		i += 1
		continue

	# Spontaneous fission is not considered
	if reac_type ==  'sf':
		i += 1
		continue

	#endf = 'no'
	endf = 'yes'
	# The format for ENDFVIII that Moritz gave me is weird, zamid is written as z0amid...
	if endf == 'yes':
		zamid = line[1][:-4] + line[1][-3:] + '{}'.format(s)
	else:
		zamid = line[1] + '{}'.format(s)

	name = line[10]



	# Need to initialize nucl_index and reac_dict 
	if i == 0:
		nucl_index.append(zamid)
		nucl_name.append(name)
		reac_dict = {}


	# If next line is about a new nuclide
	# Add the former reac_dic created to decay_dict
	# Create a new one
	if i != 0:
		if zamid != nucl_index[-1]:
			decay_dict[nucl_index[-1]] = [reac_dict, hl_dict]
			nucl_index.append(zamid)
			reac_dict = {}
			hl_dict = {}

	br = line[2]
	br_error = line[3]
	half_life = line[7]
	half_life_error = line[8]

	# ENDFVIII has nuclide that are not labelled as stable but that have decay 0 seconds
	if reac_type != 'stable' and half_life != '0':
		reac_type = reac_name_dict[reac_type] + state_dict[daughter_s]
		unit = 's'
	if reac_type == 'stable' or half_life == '0':
		half_life = 'stable'
		half_life_error = 0
		unit = 'stable'

	reac_dict[reac_type] = {'br':br, 'br_error':br_error, 'daughter_s':daughter_s}
	hl_dict = {'half_life': half_life, 'half_life_error': half_life_error, 'unit':unit}

	# if this is the last line
	print (i)
	print (len(lines)-1)
	if i == len(lines)-1:
		print (zamid)
		decay_dict[zamid] = [reac_dict, hl_dict]

	i += 1

# print (decay_dict['420900'])

# quit()



txt = '===================================================================================\n'
txt += '||  Nuclide Name  ||    ID    || Reaction Type ||   Branching   ||  Uncertainty  ||\n'
txt += '===================================================================================\n'
for i in range(len(nucl_index)):
	zamid = nucl_index[i]
	if zamid not in reduced_nucl_set:
		continue
	passpt = Passport(zamid)
	name = passpt.name
	unit = decay_dict[zamid][1]['unit']
	half_life = decay_dict[zamid][1]['half_life']
	half_life_error = decay_dict[zamid][1]['half_life_error']
	reac_dict = decay_dict[zamid][0]

	if half_life == 'stable':
		txt += '{:^19}{:^12}   {:<17}\n'.format(name, zamid, 'stable')
	else:
		#txt += '{:^19}{:^12}   {:<17}{:<17}\n'.format(name, zamid, 'unit',unit)
		txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format(name, zamid, 'half-life', float(half_life), float(half_life_error))
		for reac in reac_dict:
			txt += '{:^19}{:^12}   {:<17}{:<17.5E}{:<17.5E}\n'.format('', zamid, reac, float(reac_dict[reac]['br']), float(reac_dict[reac]['br_error']))

	# txt += {}
	# txt += ''
	# txt +='        {}         {}    unit           {}            0\n'.format(name, zamid, time_dic[act_data[i][1]])
	# txt +='                      {}     half-life                  {}            0\n'.format(zamid,act_data[i][2])
	# txt +='                      {}     betanegX                  {}            0\n'.format(zamid,act_data[i][3])
	# txt +='                      {}     betapos              {}            0\n'.format(zamid,act_data[i][4])
	# txt +='                      {}     betaposX          {}            0\n'.format(zamid,act_data[i][5])
	# txt +='                      {}     alpha                 {}            0\n'.format(zamid,act_data[i][6])
	# txt +='                      {}     gamma              {}            0\n'.format(zamid,act_data[i][7])



OpenBU_format.write(txt)

