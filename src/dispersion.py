#!/usr/bin/env python

import numpy as np
import scipy.fftpack
import cmath
import matplotlib.pyplot as plt


def get_fft(traces, dt, nt):
    """ Get temporal Fourier transform for each of the traces
    """
    f = scipy.fftpack.fftfreq(nt,dt) #f = np.linspace(0.0, 1.0/(2.0*dt), nt//2)
    U = scipy.fftpack.fft(traces)
    if np.size(U.shape) > 1:
        return U[:,0:nt//2], f[0:nt//2]
    else:
        return U[0:nt//2], f[0:nt//2]

#def get_ifft(traces, df, nf):
def get_ifft(traces):
    """ Get temporal inverse Fourier transform for each of the traces
    """
    U = scipy.fftpack.ifft(traces)
    return U

def get_fk(traces, dt, nt, dx, nr):
    """ Apply 2D Fourier to convert time-offset data 
        into data in the frequency-wavenumber domain
    """
    f = scipy.fftpack.fftfreq(nt,dt) #f = np.linspace(0.0, 1.0/(2.0*dt), nt//2)
    x = scipy.fftpack.fftfreq(nr,dx) #f = np.linspace(0.0, 1.0/(2.0*dx), nr//2)
    k = 1.0/(x + dx)
    U = scipy.fftpack.fft2(traces)
    return U, k, f[0:nt//2]

def get_dispersion(traces,dx,cmin,cmax,dc,fmax):
    """ calculate dispersion curves after Park et al. 1998
    INPUTS
    traces: SU traces
    dx: distance between stations (m)
    cmax: upper velocity limit (m/s)
    fmax: upper frequency limit (Hz)
    OUTPUTS
    f: 1d array frequency vector
    c: 1d array phase velocity vector
    img: 2d array (c x f) dispersion image
    fmax_idx: integer index corresponding to the given fmax
    U: 2d array (nr x npts//2) Fourier transform of traces
    t: 1d array time vector
    """
    nr = len(traces) 
    dt = traces[0].stats.delta
    print('dt: ', dt)
    nt = traces[0].stats.npts
    print('nt: ', nt)
    t = np.linspace(0.0, nt*dt, nt)
    traces.detrend()
    traces.taper(0.05,type='hann')
    U, f = get_fft(traces, dt, nt)
    #dc = 10.0 # phase velocity increment
    c = np.arange(cmin,cmax,dc) # set phase velocity range
    df = f[1] - f[0]
    fmax_idx = int(fmax//df)
    print('Frequency resolution up to %5.2f kHz: %i bins' % (fmax, fmax_idx))
    print('Phase velocity resolution up to %5.2f m/s: %i bins' % (cmax, len(c)))
    img = np.zeros((len(c),fmax_idx))
    x = np.linspace(0.0, (nr-1)*dx, nr)
    for fi in range(fmax_idx): # loop over frequency range
        for ci in range(len(c)): # loop over phase velocity range
            k = 2.0*np.pi*f[fi]/(c[ci])
            img[ci,fi] = np.abs(np.dot(dx * np.exp(1.0j*k*x), U[:,fi]/np.abs(U[:,fi])))

    return f,c,img,fmax_idx,U,t

def slant_stack_frequency(traces, dx, cmax, fmax):
    """
    Slant stack method after McMechan and Yedlin 1981
    """
    nr = len(traces) 
    dt = traces[0].stats.delta
    nt = traces[0].stats.npts
    t = np.linspace(0.0, nt*dt, nt)
    traces.detrend()
    traces.taper(0.05,type='hann')
    U, f = get_fft(traces, dt, nt)  # u(r,t) -> u(r,f)
    df = f[1] - f[0]
    fmax_idx = int(fmax//df)

    dc = 50.0 # phase velocity increment
    c = np.arange(50.0,cmax,dc) # set phase velocity range
    p = 1.0/c # sorted in p_max,..,p_min order
    print('Phase slowness resolution up to %5.2f m/s: %i bins' % (cmax, len(c)))

    Upf = np.zeros((len(p), len(f)), dtype = complex)
    x = np.linspace(0.0, (nr-1)*dx, nr)
    for fi in range(fmax_idx): # loop over frequency range
        for pi in range(len(p)): # loop over slowness
            k = 2.0*np.pi*f[fi]*p[pi]
            Upf[pi,fi] = np.dot(dx * np.exp(1.0j*k*x), U[:,fi]/np.abs(U[:,fi]))

    Utaup = np.zeros((len(t)//2, len(p)), dtype = complex)
    for pi in range(len(p)):
        Utaup[:,pi] = get_ifft(Upf[pi,:])

    return f, p, fmax_idx, Upf[:,:fmax_idx], Utaup


#def inverse_slant_stack_frequency(Utaup, p):
#    """ 
#    Inverse slant stack
#    """
#    Upf = np.zeros(Utaup.transpose().shape, dtype = complex)
#    for pi in range(len(p)): # loop over slowness
#    Upf[pi,:] = get_fft(Utaup[:,pi]
#
#
#    #TO-DO
#
#    return Upf


def slant_stack_time(traces,dx,cmax):
    """
    Slant stack method after McMechan and Yedlin 1981
    """
    nr = len(traces) 
    dt = traces[0].stats.delta
    nt = traces[0].stats.npts
    t = np.linspace(0.0, nt*dt, nt)
    traces.detrend() 
    traces.taper(0.05,type='hann') 
    dc = 50.0 # phase velocity increment 
    c = np.arange(50.0,cmax,dc) # set phase velocity range 
    p = 1.0/c # sorted in p_max,..,p_min order 
    print('Phase slowness resolution up to %5.2f m/s: %i bins' % (cmax, len(c)))
    tau = np.linspace(0.0, np.max(t), len(t))
    U = np.zeros((len(p),len(tau)))
    u = np.array(traces) # obspy stream to 2d array
    x = np.linspace(0.0, (nr-1)*dx, nr) # array coordinates
    for pi in range(len(p)):
        print('slowness is %f: ' % p[pi])
        for taui in range(len(tau)):
        #   x = 0.0
        #   ux = 0.0
        #   for rcv_i in range(nr):
        #      ti = tau[taui] + p[pi]*x
        #      ti_idx = int(ti//dt)
        #      #print(ti, ti_idx)
        #      if ti_idx < nt:
        #           ux += traces[rcv_i].data[ti_idx]
        #      x += dx
            ti = tau[taui] + p[pi]*x
            ti_idx = (ti//dt).astype(int)
            for rcv_i in range(nr):
                #U[pi,taui] += (u[rcv_i,ti_idx[rcv_i]] if ti_idx[rcv_i] < nt else 0.0)
                if ((ti_idx[rcv_i] < nt) and (ti_idx[rcv_i] > 0)):
                    U[pi,taui] += u[rcv_i,ti_idx[rcv_i]]

    return p,tau,U

