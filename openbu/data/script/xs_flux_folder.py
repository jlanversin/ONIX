# This first part of the script should find the index of the energy array for the value that are closest
# openbu flux spectrum energy bin

xs_energy_bin = [i for i in range(100)]

spectrum_energy_bin = [1.2, 5.3,  20, 53.5, 89.01, 100]

xs_energy_index = []
for spect_point in spectrum_energy_bin:

	for j in range(len(xs_energy_bin)):
		xs_point = xs_energy_bin[j]
		if j == 0:
			diff = abs(spect_point-xs_point)

		else:
			if diff < abs(spect_point-xs_point):

				break

	xs_energy_index.append(j)

print (xs_energy_bin)
print (xs_energy_index)