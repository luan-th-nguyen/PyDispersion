#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from src import readers, writers

reader = getattr(readers, 'su')
writer = getattr(writers, 'su') 

def read(path, filename, endian):
    traces = reader(path, filename, endian)
    return traces



#****************************************************************************
# Data description
#****************************************************************************
NR = 50 
tsteps = 10000 
dt = 1.0e-4 #ms
dx = 50*(0.001e3) #mm
time = np.arange(0,tsteps)*dt


#****************************************************************************
# input
#****************************************************************************
path_in = '../../data/viscoelastic_waves/'
path_out = './'
if not os.path.exists(path_out):
    os.makedirs(path_out)

comp = 'Uz'
fname = comp + '_file_single.su'
traces = read(path_in, fname, 'big')

fig, ax = plt.subplots(figsize=(6.0,5))

ampl_max = np.max(traces[0].data)
for ri in range(NR):
    plt.plot(time, ri*dx + (100.0/ampl_max)*traces[ri].data, 'b-', linewidth=2.0) # obs. radial component
    print 'max of data plot: %e' % np.max(traces[ri].data)

    print 'rec.num. ' + str(ri)

ax.set_xlim([time[0], time[-1]])
ax.set_ylim([0-dx, 0+NR*dx])
plt.gca().invert_yaxis()
ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
ax.tick_params(axis = 'both', which = 'minor', labelsize = 14)
plt.ylabel('Offset [mm]', fontsize = 14)
plt.xlabel('Time [ms]', fontsize = 14)
#plt.show()

plt.savefig(path_out + 'traces_dispersion_analysis'  +  '.png',dpi=300)
