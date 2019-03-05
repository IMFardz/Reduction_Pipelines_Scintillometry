import numpy as np
from astropy.time import Time

def read_and_bin_profile(fname,tbinned):
    gpu = np.load(fname)
    t = gpu['t']
    I = gpu['I']
    f = gpu['f']

    for i in range (len(t)):
        t[i] = np.float(t[i].value)
    t = np.array(t, dtype=float)
    t = Time(t,format='mjd')

    Ibinned = np.zeros((len(tbinned),I.shape[1], I.shape[2]))
    for i in range(I.shape[0]):
        idx = np.searchsorted(tbinned.mjd,t[i].mjd)
        Ibinned[idx-1,:]+=I[i,:,:]
    Ibinned = Ibinned/Ibinned.mean(axis=1,keepdims=True)
    return Ibinned,f, tbinned

def read_and_bin(fname,tbinned):
    gpu = np.load(fname)
    t = gpu['time']
    I = gpu['I']
    f = gpu['freq']

    t = Time(t,format='mjd')

    Ibinned = np.zeros((len(tbinned),I.shape[1]))
    for i in range(I.shape[0]):
        idx = np.searchsorted(tbinned.mjd,t[i].mjd)
        Ibinned[idx-1,:]+=I[i,:]
    Ibinned = Ibinned/Ibinned.mean(axis=1,keepdims=True)
    return Ibinned,f, tbinned
