#This file convert the libraries into a text matrix. The text matrix is not supposed to be read by
# human eyes but by the code to directly build the matrix in the python environment
import os
from openbu.passport import passport
from openbu import data as d

def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def azm_sort(seq):
    """Sort the nuclides by mass number"""
    changed = True
    while changed:
        changed = False
        for i in range(len(seq) - 1):
            a0 = int(seq[i][-4:-1])
            a1 = int(seq[i+1][-4:-1])
            z0 = int(seq[i][:-4])
            z1 = int(seq[i+1][:-4])
            state0 = int(seq[i][-1])
            state1 = int(seq[i+1][-1])
            azm0 = 1000*a0 + 10*z0 + state0
            azm1 = 1000*a1 + 10*z1 + state1
            if azm0 > azm1:
                seq[i], seq[i+1] = seq[i+1], seq[i]
                changed = True
    return None

def pointer_dic(s):
	dic = {}
	for i in range(len(s)):
		dic[s[i]] = i

	return dic

path_to_decay = os.path.join(os.path.dirname(__file__), '../default_libs/decay_lib')
path_to_xs = os.path.join(os.path.dirname(__file__), '../default_libs/xs_lib')
path_to_fy = os.path.join(os.path.dirname(__file__), '../default_libs/fy_lib')

decay_file = open(path_to_decay, 'r')
xs_file = open(path_to_xs, 'r')
fy_file = open(path_to_xs, 'r')

decay_line = decay_file.readlines()
xs_line = xs_file.readlines()
fy_line = fy_file.readlines()

decay_nucl = []

for line in decay_line:
	if len(line.split()) > 1:
		if is_int(line.split()[1]):
			if line.split()[1] in decay_nucl:
 				continue
			decay_nucl.append(line.split()[1])


xs_nucl = []

for line in xs_line:
	if len(line.split()) > 1:
		if is_int(line.split()[1]):
			if line.split()[1] in xs_nucl:
				continue
			xs_nucl.append(line.split()[1])

nucl_list = list(set(decay_nucl + xs_nucl))

print((len(nucl_list)))

azm_sort(nucl_list)

pointer_dic = pointer_dic(nucl_list)

passlist = []
for i in nucl_list:
	passlist.append(passport(i))

# Read the xs data and pass them to the list
for i in range(len(passlist)):

    nuc_pass = passlist[i]
    zamid = nuc_pass.nuc_zzaaam
    if zamid in xs_nucl:

        xs = d.default_xs_lib[zamid]
        nuc_pass.set_xs(xs)

 # Read the decay data and pass them to the list
for i in range(len(passlist)):

    nuc_pass = passlist[i]
    zamid = nuc_pass.nuc_zzaaam
    
    if zamid in decay_nucl:
        decay_a = d.default_decay_lib_a[zamid]
        decay_b = d.default_decay_lib_b[zamid]
    
        nuc_pass.set_decay(decay_a, decay_b)

# Read the fission yields data and pass them to the list
for i in range(len(passlist)):

    nuc_pass = passlist[i]
    zamid = nuc_pass.nuc_zzaaam
    FAM = nuc_pass.check_FAM()
    if FAM == 'FP' and zamid in xs_nucl:
    	# I consder FP as z<89 but there are non FP for z<89. Since I don't want to follow ORIGEN classification of
    	# FP, I use a try statement
    	try:
            fy = d.default_fy_lib[zamid]
            nuc_pass.set_fy(fy)
        except KeyError:
        	pass
#pass_decay(passlist)

# U5_pointer = pointer_dic['922350']
# U5_pass = passlist[U5_pointer]
# print U5_pass.xs
# print U5_pass.decay_a


# X5_pointer = pointer_dic['541350']
# X5_pass = passlist[X5_pointer]
# print X5_pass.fy


# Writting the cross section matrix file
Btxt = ''
i = 0
for child in passlist:
	Btxt += '{}|'.format(child.nuc_zzaaam)
	Btxt += '{}:'.format(i)
	if child.nuc_zzaaam in xs_nucl:
		child_removal = -child.xs['removal'][0]
		Btxt += '{} {},'.format(i, child_removal)

	# Cross sections (non fission)
	xs_parent = child.xs_parent
	for j in xs_parent:

		if xs_parent[j] not in xs_nucl:
			continue
		index = pointer_dic[xs_parent[j]]
		parent_pass = passlist[index]
		if j not in parent_pass.xs:
			continue
		xs_val = parent_pass.xs[j][0]
		if xs_val == 0.0:
			continue
		Btxt += ' {} {},'.format(index, xs_val)

	# Fission yields
	if child.fy != None:
		fy = child.fy 
		for j in fy:
			if fy[j] == 0.0:
				continue
			index = pointer_dic[j]
			parent_pass = passlist[index]
			fission_val = parent_pass.xs['fission'][0]
			val = fission_val*fy[j][0]*1e-2
			Btxt += ' {} {},'.format(index, val)

	Btxt += '\n'
	i += 1

B = open('Btxt', 'w')
B.write(Btxt)

Ctxt = ''
i = 0
for child in passlist:
	Ctxt += '{}|'.format(child.nuc_zzaaam)
	Ctxt += '{}:'.format(i)
	if child.nuc_zzaaam in decay_nucl:
		child_removal = -child.decay_a['total']
		Ctxt += '{} {},'.format(i, child_removal)

	d_parent = child.d_parent
	for j in d_parent:

		if d_parent[j] not in decay_nucl:
			continue
		index = pointer_dic[d_parent[j]]
		parent_pass = passlist[index]
		d_val = parent_pass.decay_a[j]
		if d_val == 0.0:
			continue
		Ctxt += ' {} {},'.format(index, d_val)
	Ctxt += '\n'
	i += 1

C = open('Ctxt', 'w')
C.write(Ctxt)



