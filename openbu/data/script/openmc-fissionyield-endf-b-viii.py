import os
import pandas as pd
import openmc.data
from openmc.data import ATOMIC_SYMBOL

directory = 'endf-b-viii.0-nfylib-single'
files = os.listdir(directory)
files.sort()

f = files[3] # just pick one for random

for f in files:
    if f[-5:] == ".endf":
        fpdata = openmc.data.FissionProductYields(os.path.join(directory, f))
        print(fpdata.nuclide, fpdata.energies)
