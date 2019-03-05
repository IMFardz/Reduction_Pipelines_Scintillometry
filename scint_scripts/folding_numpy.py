# Use scintillometry to fold a numpy array
import astropy.units as u
import astropy.constants as c
from baseband import mark4
from astropy.time import Time
import numpy as np
import sys
import os
from pulsar import predictor
from glob import glob
import traceback
from scintillometry import integration, generators

# User input
fname = sys.argv[1]
start_time = sys.argv[2]
freq_cent = sys.argv[3]

# Set constants
output_name = 'test'
date = '2018-02-28 '
dispersion_measure = 4.84066 * u.pc / u.cm**3
dt = (2.6/1000) * u.s
nphase = 1024
central_freq = np.int(freq_cent) * u.MHz
if freq_cent == '400':
    psr_polyco = predictor.Polyco('polycos/B1133_lower_gmrtpolyco_new.dat')
else:
    psr_polyco = predictor.Polyco('polycos/B1133_upper_gmrtpolyco_new.dat')

# Create stream generator
print("Loading data")
fh = np.load(fname)
I = fh['arr_0']
print("Shape of Data:{}".format(I.shape))
freq = np.linspace(-100, 100, I.shape[1]) * u.MHz + central_freq
samples_per_frame = np.int(np.round(10 / (2.6 / 1000)))
start = Time(date + start_time)

def my_source(sh):
    offset = sh.tell()
    return I[offset:offset+sh.samples_per_frame]

sh = generators.StreamGenerator(my_source, shape=I.shape, start_time=start, sample_rate=1/dt)

print("Creating Integrator")
integrator = integration.Fold(sh, n_phase=nphase, phase=psr_polyco, samples_per_frame=1, step=10*u.s, average=True)
print("Created Integrator with Shape:{}".format(integrator.shape))
# Read out one sample at  a time
# Determine how many samples to output at a time
nsamples = integrator.shape[0]
times = []

# Loop through integrator, creating one time bin at a time
print("Starting fold:")
try:
    output, t = integrator.read(count=1), integrator.time
    times.append(t)
    while integrator.tell() < nsamples - 1:
        sample, t = integrator.read(count=1), integrator.time
        print(t.yday)
        times.append(t)
        output = np.r_[output, sample]
    times = Time(times)
    np.savez(output_name, I=output, t=times, f=frequency)
except Exception as err:
    print("Something went wrong. Likely, you inputted noise or too small of a sample set")
    print(traceback.format_exc())
    print(sys.exc_info()[0])
    np.savez(output_name, I=output, t=times)
