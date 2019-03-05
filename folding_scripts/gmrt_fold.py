"""This script folds a portion of the GMRT Data"""
import numpy as np
import matplotlib.pyplot as plt
from pulsar import predictor
from astropy.time import Time, TimeDelta
import astropy.units as u
import sys

# User inputs
# Start time are taken from the header files
fname = sys.argv[1]
bw = sys.argv[2]
start_time = sys.argv[3]
plt_dir = "gmrt_plots/"

# Setting some parameters
n_phase = 1024
outname = fname.split("/")[-1][:-4]
date = '2018-02-28 '
print(date + start_time)
print("The outname for this file is:", outname)
if bw == 0:
    psr_polyco = predictor.Polyco('polycos/B1133_lower_gmrtpolyco_new.dat')
else:
    psr_polyco = predictor.Polyco('polycos/B1133_upper_gmrtpolyco_new.dat')

# Load data
file = np.load(fname)
fh = file['arr_0']
print("Shape of stream is: ", fh.shape)

# Set time step
orig_dt = 2.62144 / 1000
subint = 10 * u.s
dt = TimeDelta(orig_dt, format='sec')
print("Length of each time bin:", dt)

# Set the frequency axis
if bw == 0:
    freq = np.linspace(300, 500, fh.shape[1])
else:
    freq = np.linspace(550, 750, fh.shape[1])

# Make a reorganized array
# First column is time, second frequency and third phase
start = Time(date + start_time)
end = start + dt * fh.shape[0]
total_t = (end - start).to(u.s)
print("Length of time chunk:", (end - start).to(u.s))
total_tbins =  np.int((total_t / subint).value) + 1
folded = np.zeros(shape=(total_tbins, fh.shape[1], n_phase))
counts = np.zeros(shape=(total_tbins, fh.shape[1], n_phase))
print("Shape of output is:", folded.shape)

# Set values for current time, current time bin and time array
start = Time(date + start_time)
t = start.copy()
time = []
current_time_bin = Time(date + start_time, precision=-1)
current_time_bin = Time(current_time_bin.__str__())
time.append(current_time_bin)
t_bin = 0
p0 = psr_polyco(t).value

# Loop through each data point
for i in range(fh.shape[0]):
    p = psr_polyco(t).value
    # Calculate which index
    if t >= current_time_bin + subint:
        t_bin += 1
        current_time_bin += subint
        time.append(current_time_bin)
    phase = np.int(np.floor((p % 1) * n_phase))
    # Populate arrays
    folded[t_bin, :, phase] += fh[i]
    counts[t_bin, :, phase] += 1
    t += dt

# To compute the average, divide by number of samples per bins
time = Time(time)
folded = np.divide(folded, counts, where=counts!=0)

# Save the Data
np.savez(outname, I=folded, f=freq, t=time)


# Create a plot describing arrays
fig, sub = plt.subplots(1, 3, figsize=(16, 4))
# Plot the pulse profile
sub[0].plot(folded.mean(0).mean(0))
sub[0].set_title("Profile", outname)
# Plot the pulse stack
sub[1].imshow(folded.mean(1), aspect='auto', interpolation='none')
sub[1].set_title("Pulse Stack", outname)
# Plot the frequency vs. phase
a = folded.mean(0) / folded.mean(0).mean(-1, keepdims=True)
a = a.reshape(a.shape[0]//2, 2, a.shape[1]).mean(1)
a = a.reshape(a.shape[0], a.shape[1]//2, 2).mean(-1)
sub[2].imshow(a, aspect='auto', interpolation='none', origin='lower')
sub[2].set_title("Frequencs vs. Phase", outname)
plt.savefig(plt_dir + outname)
