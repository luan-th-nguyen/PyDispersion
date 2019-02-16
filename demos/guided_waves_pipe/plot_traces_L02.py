#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, axes, sci
from matplotlib import cm, colors
from matplotlib.font_manager import FontProperties
from numpy import amin, amax, ravel
from numpy.random import rand
import scipy.io as sio


path_in = '../../data/guided_waves_pipe/data_L02/'
path_out = './'
#if not os.path.exists(path_out):
#    os.makedirs(path_out)


#****************************************************************************
# Data description
#****************************************************************************
nrecs_x = 24 #48 #24 #240 #360
delta_deg = 0.0 #7.5 #15 #1.5 # 1
nrecs_y = 1
nrecs = nrecs_x * nrecs_y

rinc_x = 1
rinc_y = 1

tsteps = 6000
dt = 0.000035


#****************************************************************************
# input
#****************************************************************************
data_plot = np.zeros ((nrecs_x, tsteps))
dat_type = 'semd'
#dat_type = 'adj'
comp1 = 'FXZ'
comp2 = 'FXZ'
sn = 1
prec = 'single'

rec_y = 1
rec_sta = (rec_y-1)*nrecs_x
rec_fin = rec_y*nrecs_x
rec_buf = rec_fin - rec_sta


fig, ax = plt.subplots(figsize=(6.0,5))
time = np.arange(0,tsteps)*dt
#offset = np.arange(-delta_deg*nrecs_x//2,delta_deg*nrecs_x//2,delta_deg)
offset = 40.0 + np.arange(0,nrecs_x*10,10)
index = 0
for rec_x in range(rec_sta+1,rec_fin+1):

    file_name_in1 = path_in + 'P.R' + str(int(rec_x)) + '.' + comp1 + '.' + dat_type
    file_name_in2 = path_in + 'P.R' + str(int(rec_x)) + '.' + comp2 + '.' + dat_type
    xz1 = np.genfromtxt(file_name_in1)
    xz2 = np.genfromtxt(file_name_in2)


    # projection
    alpha = np.arctan(xz2[:tsteps,1]/(1.0e-40 + xz1[:tsteps,1])) # angle of projection
    #alpha = np.angle(1.0j*xz2[:,1] + xz1[:,1]) - deg*np.pi/180.0 # angle of projection
    # negative or positive of the projected vector on [cos(phi) sin(phi)] 
    deg = 0.0
    direction = np.sign(np.cos(deg*np.pi/180.0)*xz1[:tsteps,1]*np.cos(alpha) + np.sin(deg*np.pi/180.0)*xz2[:tsteps,1]*np.cos(alpha))    
    data_plot[index,:] = direction*np.sqrt(xz1[:tsteps,1]**2 + xz2[:tsteps,1]**2)*np.cos(alpha) # scalar radial component
    #data_plot[index,:] = xz1[:,1]

    # plot
    #plt.plot(deg + 5.0e4*xz1[:,1],time, 'b-') # x component
    #print np.max(xz1[:,1])
    #plt.plot(deg + 5.0e4*xz2[:,1],time, 'b-') # y component
    #print np.max(xz2[:,1])
    plt.plot(offset[index] + 0.9e3*data_plot[index,:], time, 'b-', linewidth=2.0) # obs. radial component
    #plt.plot(offset[index] + 0.4e4*np.flip(data_plot[index,:],0), time, 'b-', linewidth=2.0) # diff. radial component
    print 'max of data plot: %e' % np.max(data_plot[index,:])

    # next trace
    index += 1
    print 'rec.num. ' + str(rec_x) + '/\tdeg. ' + str(deg)

ax.set_ylim([time[-1], time[0]])
#ax.set_xlim([-10, 240])
ax.set_xlim([30, 280])
ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
ax.tick_params(axis = 'both', which = 'minor', labelsize = 14)
plt.xlabel('Offset [mm]', fontsize = 14)
plt.ylabel('Time [ms]', fontsize = 14)
#plt.show()

plt.savefig(path_out + 'traces_supershot_dispersion_analysis_L02'  +  '.png',dpi=300)
