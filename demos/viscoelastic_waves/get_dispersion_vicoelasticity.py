#!/usr/bin/env python

from src import readers, writers
from src.dispersion import get_dispersion, get_fft, slant_stack_time
import numpy as np
import matplotlib.pyplot as plt


reader = getattr(readers, 'su')
writer = getattr(writers, 'su') 
def read(path, filename, endian):
    traces = reader(path, filename, endian)
    return traces


if __name__ == '__main__':

    path = '../../data/viscoelastic_waves/'
    comp = 'Uz'
    fname = comp + '_file_single.su'
    u50 = read(path, fname, 'big')
    u = u50[0:50] # select traces
    #print(u)

    nt = 10000 #30000 
    dt = 1.0e-4 #ms
    dx = 50.0 #mm
    cmin = 4100.0
    cmax = 5300.0
    dc = 5.0
    fmax = 200.0 #kHz
    fig1, ax = plt.subplots(figsize=(8.0,5.0))

    ## phase shift method
    f,c,img,fmax_idx,U,t = get_dispersion(u,dx,cmin,cmax,dc,fmax)
    cax = ax.imshow(img[:,:],aspect='auto',origin='lower',extent=(f[0],f[fmax_idx],c[0],c[-1]),interpolation='bilinear')
    ax.grid(linestyle='--',linewidth=1)
    #cbar = fig1.colorbar(cax, orientation='horizontal', label='wave energy')
    # extract maxima along frequency axis
    c_line = np.zeros(fmax_idx)
    for fi in range(fmax_idx):
        cmax_idx =np.argmax(img[:,fi])
        #print cmax_idx
        c_line[fi] = c[cmax_idx]
    ax.plot(range(fmax_idx), c_line, 'c--', linewidth=3.0)
    # labels and fontsize
    ax.set_xlabel('Frequency [kHz]', fontsize=14)
    ax.set_ylabel('Phase velocity [m/s]', fontsize=14)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 14)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 14)

    #fig,ax = plt.subplots()
    #ax.plot(range(fmax_idx), c_line, 'c--', linewidth=3.0)

    r = 0
    fig2, axarr = plt.subplots(3,1)
    axarr[0].plot(t,u[r],label='time signal at receiver %i' % (r+1))
    axarr[0].legend()
    axarr[0].set_xlabel('Time [s]')
    axarr[1].plot(f,np.abs(U[r]),label='abs')
    axarr[1].plot(f,U[r].real,label='real')
    axarr[1].plot(f,U[r].imag,label='imag')
    axarr[1].legend()
    axarr[1].set_xlabel('Frequency [kHz]')
    axarr[1].set_xlim([0, 2.0*fmax])
    axarr[2].plot(f,np.angle(U[r]), label='phase')
    axarr[2].plot(f[:2*fmax_idx],np.unwrap(np.angle(U[r,:2*fmax_idx])), label='unwrapped phase')
    axarr[2].set_xlim([0, 2.0*fmax])
    axarr[2].legend()
    axarr[2].set_xlabel('Frequency [kHz]')

    ## slant stack method
    #fig3, ax3 = plt.subplots(figsize=(8.0,5.0))
    #p,tau,U_ptau = slant_stack_time(u,dx,cmax)
    #img,f = np.abs(get_fft(U_ptau,dt,nt))
    #df = f[1] - f[0]
    #fmax_idx = int(fmax//df)
    #ax3.imshow(img[:,:fmax_idx],aspect='auto',origin='upper',extent=(f[0],f[fmax_idx],p[0],p[-1]),interpolation='bilinear')


    #plt.show()
    fig1.savefig('./viscous_dispersion_curves'  +  '.png',dpi=300)
