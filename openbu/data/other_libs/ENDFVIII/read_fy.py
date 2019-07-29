import os
import pandas as pd
import openmc.data
from openmc.data import ATOMIC_SYMBOL

directory = 'ENDF-B-VIII.0_nfy'
files = [file for file in os.listdir(directory) if ".endf" in file]
files.sort()

f = files[3] # just pick one for random

for file in files:

	file_name = file.split('_')
	nucl_name = file_name[1]+file_name[2]
	nucl_name = nucl_name.replace('.endf', '')
	print (nucl_name)

#f = files[3] # just pick one for random

# if f[-5:] == ".endf":
#     fpdata = openmc.data.FissionProductYields(os.path.join(directory, f))
#     print(fpdata.nuclide)
#     #print(fpdata.cumulative)
#     print(fpdata.independent)