"""This script is designed to make secondary spectra from the GMRT DataSet"""
import numpy as np
import astropy.units as u


def svd_model(arr, nmodes, phase_only=True):
    """
    Take time/freq visibilities SVD, zero out all but the largest
    mode, multiply original data by complex conjugate
    Parameters
    ----------
    arr : array_like Time/freq visiblity matrix
    Returns
    -------
    Original data array multiplied by the largest SVD mode conjugate
    """

    u,s,w = np.linalg.svd(arr)
    s[nmodes:] = 0.0
    S = np.zeros([len(u), len(w)], np.complex128)
    S[:len(s), :len(s)] = np.diag(s)

    model = np.dot(np.dot(u, S), w)

    #model = norm(model)  # This is the line I added

    return np.abs(model)


def make_secondary(dynamic, time, frequency, cropped=False):
    """
    Create a secondary spectrum from a given dynamic spectrum, time axis and
    frequency axis
    =====================
    Parameters:
    dynamic: 2D numpy array. Intensity. Axis are (frequency and time)
    time: array of times. Astropy quantities of times
    frequency: array of frequency. Astropy quantities of Hz
    =====================
    Return:
    2-D Secondary Spectrum, doppler frequency and delay axis.
    If cropped, only positive delays are outputted
    """
    secondary = np.fft.fftshift(np.fft.fft2(dynamic))
    doppmin = 50 / 1000 * u.Hz
    delaymin = 0.5 * frequency.size / (frequency.max()-frequency.min()).to(u.Hz)
    dopp = np.linspace(-1 * doppmin.value, doppmin.value, time.size) * u.Hz
    delay = np.linspace(-1 * delaymin.value, delaymin.value, frequency.size) * u.s
    if cropped:
        secondary = secondary[secondary.shape[0]:,:]
        delay = delay[delay.shape[0]:]
    return secondary, dopp, delay


def make_gate(I, on_bound, off_bound, svd=True):
    """
    Create a dynamic spectrum for given phase bounds defined by on.
    If svd is True, normalize via an SVD
    """
    # Create On and Off Gate
    on = I[:,:,on_bound[0]:on_bound[1]].mean(-1)
    off = np.delete(I, off_bound, axis=-1).mean(-1)

    # Divide out Off Gate
    gate = np.divide(on , off, where=off!=0)
    gate = gate.T

    if svd:
        mask = svd_model(gate, nmodes=1, phase_only=False)
        dynamic = np.divide(gate, mask, where=mask!=0) - 1
    else:
        dynamic = gate - 1

    return dynamic
