#!/usr/bin/env python

import numpy as np
import scipy.fftpack
import matplotlib.pyplot as plt
from src import readers, writers
from src.dispersion import get_dispersion


reader = getattr(readers, 'su')
writer = getattr(writers, 'su') 
def read(path, filename):
    traces = reader(path, filename)
    return traces

def read_ascii(path, NR, nt):
    from numpy import loadtxt
    from obspy.core import Stream, Stats, Trace
    dat_type = 'semd'
    comp1 = 'FXX'
    comp2 = 'FXY'
    stream = Stream()
    for rec_x in range(0,NR):
        file_name_in1 = path + 'P.R' + str(int(rec_x+1)) + '.' + comp1 + '.' + dat_type
        file_name_in2 = path + 'P.R' + str(int(rec_x+1)) + '.' + comp2 + '.' + dat_type
        xz1 = np.genfromtxt(file_name_in1)
        xz2 = np.genfromtxt(file_name_in2)
        deg = 0.0
        alpha = np.arctan(xz2[:nt,1]/(1.0e-40 + xz1[:nt,1])) # angle of projection
        direction = np.sign(np.cos(deg*np.pi/180.0)*xz1[:nt,1]*np.cos(alpha) + np.sin(deg*np.pi/180.0)*xz2[:nt,1]*np.cos(alpha))    
        data = direction*np.sqrt(xz1[:nt,1]**2 + xz2[:nt,1]**2)*np.cos(alpha) # scalar radial component

        stats = Stats()
        stats.filename = path + 'P.R' + str(int(rec_x+1))
        stats.starttime = xz1[0,0]
        stats.delta = xz1[1,0] - xz1[0,0]
        stats.npts = len(xz1[:nt,0])

        try:
            parts = filename.split('.')
            stats.network = parts[0]
            stats.station = parts[1]
            stats.channel = temp[2]
        except:
            pass

        stream.append(Trace(data=data[:], header=stats))

    return stream


if __name__ == '__main__':

    path = '../../data/guided_waves_pipe/data_L01/'
    NR = 24
    nt = 2000 # 6000
    u = read_ascii(path, NR, nt)
    dt = 3.5e-5
    dx = 10.0
    cmin = 50.
    cmax = 8000.0
    dc = 10.
    fmax = 600.0
    im, ax = plt.subplots(figsize=(7.0,5.0))

    ## phase shift method
    f,c,img,fmax_idx,U,t = get_dispersion(u,dx,cmin,cmax,dc,fmax)
    ax.imshow(img[:,:],aspect='auto',origin='lower',extent=(f[0],f[fmax_idx],c[0],c[-1]),interpolation='bilinear')
    ax.grid(linestyle='--',linewidth=2)

    ## slant stack method
    #p,tau,U_ptau = slant_stack(u,dx,cmax,fmax)
    #img,f = np.abs(get_fft(U_ptau,dt,nt))
    #df = f[1] - f[0]
    #fmax_idx = int(fmax//df)
    #ax.imshow(img[:,:fmax_idx],aspect='auto',origin='upper',extent=(f[0],f[fmax_idx],p[0],p[-1]),interpolation='bilinear')

    ax.set_xlabel('Frequency [kHz]', fontsize=14)
    ax.set_ylabel('Phase velocity [m/s]', fontsize=14)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 14)

    r = 10
    fig, axarr = plt.subplots(3,1)
    axarr[0].plot(t,u[r],label='time signal at receiver %i' % (r+1))
    axarr[0].legend()
    axarr[0].set_xlabel('Time [s]')
    axarr[1].plot(f,np.abs(U[r]),label='abs')
    axarr[1].plot(f,U[r].real,label='real')
    axarr[1].plot(f,U[r].imag,label='imag')
    axarr[1].legend()
    axarr[1].set_xlabel('Frequency [kHz]')
    axarr[1].set_xlim([0, 2.0*fmax])
    axarr[2].plot(f[:2*fmax_idx], np.angle(U[r,:2*fmax_idx]), label='phase')
    axarr[2].plot(f[:2*fmax_idx], np.unwrap(np.angle(U[r,:2*fmax_idx])), label='unwrapped phase')
    axarr[2].legend()
    axarr[2].set_xlabel('Frequency [kHz]')

    #plt.show()
    im.savefig('./pipe_syn_dispersion_curves_L01'  +  '.png',dpi=300)
