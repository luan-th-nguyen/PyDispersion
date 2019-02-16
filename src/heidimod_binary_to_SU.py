#!/usr/bin/env python
import os
import shutil
import sys
import numpy as np
from obspy import read
from subprocess import call
from obspy.core import Stream, Stats, Trace
#from obspy import read
#import scipy.io as sio

path_in = '../heidimod_data/'
path_out = '../heidimod_data/'
if not os.path.exists(path_out):
    os.makedirs(path_out)


#****************************************************************************
# data description
#****************************************************************************
NR = 50 # number of receivers
LEN = 10000 #8000# data length 
comp1 = 'Uz'


# read Heidimode data
dz_src = path_in  + '/arbzseis'
dz_data = np.fromfile(dz_src,dtype='>f') # big endian float (4 bytes)
dz_data = dz_data.reshape((NR,LEN), order="F")
dz_data = np.float32(dz_data)


# write SU-format binary
dz_dest = path_out  + '/' + comp1 + '_file_single.su'
stats = Stats()
stats.filename = dz_dest
stats.starttime = 0.0
stats.delta = 1.0e-4 # ObsPy DO NOT WORK WITH VERY SMALL TIMESTEP 10.0e-8
stats.npts = LEN
stream = Stream()
for i in range(NR):
    stream.append(Trace(data=dz_data[i,:], header=stats))
    #stream.append(Trace(data=dz_data[i,:]))
    print 'max of trace %i is %f' % (i, np.max(dz_data[i,:]))
print stream
#stream[0].plot()

stream.write(dz_dest,format='su')


