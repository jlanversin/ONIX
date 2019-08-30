import os
import pandas as pd
import openmc.data
from openmc.data import ATOMIC_SYMBOL

directory = 'endf-b-viii.0-decaylib-single'
files = os.listdir(directory)
files.sort()

data = {'name': [], 'ZA': [], 'isomeric_state': [],
        'half_life': [], 'half_life_error': [],
        'type': [],
        'daughter_isomeric_state': [],
        'energy': [], 'energy_error': [],
        'br': [], 'br_error': []}
stableskip = ['half_life', 'half_life_error', 'daughter_isomeric_state', 'energy', 'energy_error', 'br', 'br_error']

nolist = []
# not checked...

for f in files:
    print(f)
    decaydata = openmc.data.Decay(os.path.join(directory, f))
    jname = "{}{}".format(ATOMIC_SYMBOL[decaydata.nuclide['atomic_number']], decaydata.nuclide['mass_number'])
    ZA = decaydata.nuclide['atomic_number'] * 1000 + decaydata.nuclide['mass_number']
    if decaydata.nuclide['isomeric_state'] != 0:
        jname += "*"
    if decaydata.nuclide['stable']:
        data['name'].append(jname)
        data['ZA'].append(ZA)
        data['isomeric_state'].append(decaydata.nuclide['isomeric_state'])
        data['type'].append('stable')
        for i in stableskip:
            data[i].append('NA')
    else:    
        for m in decaydata.modes:
            data['name'].append(jname)
            data['ZA'].append(ZA)
            data['isomeric_state'].append(decaydata.nuclide['isomeric_state'])
            hl = decaydata.half_life
            data['half_life'].append(hl.n)
            data['half_life_error'].append(hl.s)
            mode = ' '.join(m.modes)
            data['type'].append(mode)
            data['daughter_isomeric_state'].append(m._daughter_state)
            e = m.energy
            data['energy'].append(e.n)
            data['energy_error'].append(e.s)
            br = m.branching_ratio
            data['br'].append(br.n)
            data['br_error'].append(br.s)
    
decaydf = pd.DataFrame(data)

decaydf.to_csv("endf-b-viii.0-decay-openmc.csv")

##%%
# files that could not be read do to error
for f in nolist:
    a = Evaluation(os.path.join(directory, f), verbose = False)
    print(a.target['zsymam'])

