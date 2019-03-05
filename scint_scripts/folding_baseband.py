import astropy.units as u
import astropy.constants as c
from baseband import mark4
from astropy.time import Time, TimeDelta
import numpy as np
import sys
import os
from glob import glob
import scintillometry_reduction as sr
from scintillometry.shaping import Reshape
import traceback

fdir = '/scratch/p/pen/simardda/B1133/gk049c/ar/'
fname = sys.argv[1]
output_name = fname[:-4]

# Load Data
frequency = np.array([[320.25], [320.25], [336.25], [336.25]]) * u.MHz
sideband = np.array([[-1, -1], [1, 1], [-1, -1], [1, 1]])
polarization = ['R', 'L']
dispersion_measure = 4.84066 * u.pc / u.cm**3
polyco_file = 'B1133_ar_polycopolyco_new.dat'
fullpol = False
print("Parameters set")

# Creating stream reader. For other formats, such as vdif, use vdif.open(...)
fh = mark4.open(fdir + fname, 'rs', decade=2010)
rh = Reshape(fh, (4, 2))
dt = TimeDelta(10, format='sec')

# Rounding Time
start = Time(rh.time)
start_time_str = start.iso.__str__()
new_time = Time(start_time_str, precision = -1)
new_time_str = new_time.iso.__str__()
start_time = Time(new_time_str) + dt

print("Opened stream reader with sample shape:", rh.sample_shape)
print("Starting at time:", start_time)

# Initial waterfall interpretor
WF = sr.Fold(rh, dispersion_measure, frequency, sideband, polyco_file, polarization, fullpol,start=start_time, nthreads=1)
print("Initialized waterfall interpretor with shape:", WF.integrator.shape)

# Determine how many samples to output at a time. I reccomend 1.
nsamples = WF.integrator.shape[0]
nsamples_per_output = 1
times = []

# Loop through integrator, creating one time bin at a time
try:
    output, t = WF.integrate_and_save(count=1)
    times.append(t)
    counter = 1
    while WF.integrator.tell() < nsamples - nsamples_per_output:
        sample, t = WF.integrate_and_save(count=1)
        print(t.yday)
        times.append(t)
        output = np.r_[output, sample]
        counter += 1

    np.savez(output_name, I=output, t=times, f=frequency)
except:
    print("Something went wrong. Likely, you inputted noise or too small of a sample set")
    print(traceback.format_exc())
    print(sys.exc_info()[0])
    np.savez(output_name, I=output, t=times, start=start_time)
