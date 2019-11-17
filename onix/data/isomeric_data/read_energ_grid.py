from os import listdir

#file_name = [x.split('. |_') for x in listdir('eaf-2010-multiplicities')]
file_name_list = [x for x in listdir('eaf-2010-multiplicities') if '_2.' not in x]


energy_grid_mat = []
for file_name in file_name_list:
	file = open('eaf-2010-multiplicities/' + file_name, 'r')
	lines = file.readlines()
	energy_grid = []
	for line in lines[1:]:
		line = line.split(',')
		energy_grid.append(line[1])

	energy_grid_mat.append(energy_grid)

#print (energy_grid_mat)


intersection = set(energy_grid_mat[0]).intersection(*energy_grid_mat)

union = list(set().union(*energy_grid_mat)) 

# print (intersection)
# print (len(union))


