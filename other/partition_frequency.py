import numpy as np
from bin_time import read_and_bin
from astropy.time import Time
import glob

filename = "gk049c_gb_thr0_dynspec.npz"

print(filename)
file = np.load(filename)
time = file['time']
time = np.arange(time.min(), time.max(), step=(time.max()-time.min())/time.shape[0])
t = Time(time, format='mjd')
spec = file["I"]
newI, newF = read_and_bin(filename, t)
np.savez(filename, I=newI, time=np.array(t), freq = np.array(newF))
