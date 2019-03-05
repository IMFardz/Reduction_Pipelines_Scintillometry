"""Use this script to concatenate different sections of dataself.
Ensure that each segment has a separate time axis of array type time
"""
import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
from glob import glob
import sys

# Get files and determine outname
fdir = sys.argv[1]
band = sys.argv[2]
subint = 10 * u.s
outname = "B1133_GMRT_feb28_" + band
files = np.sort(glob(fdir + "*28feb2018.npz"))
print("Outname file name:{}".format(outname))

# Load time arries
print("Making time array")
times = []
for i in range(len(files)):
    print(files[i])
    fh = np.load(files[i])
    t = fh['t']
    times.append(t)

# Find the gap between each files
gaps = []
for i in range(1, len(times)):
    dt = times[i][0] - times[i - 1][-1]
    # Number of time bins between the two
    tbins = (dt.to(u.s) / subint).value
    print("Number of time bins between:{}".format(tbins))
    gaps.append(np.int(tbins))

# Load first time segment
print("Concatenating files")
sh = np.load(files[0])
output = sh["I"]

# Iteratively Adjoin Each Segment
for i in range(len(gaps)):
    ih = np.load(files[i + 1])
    I = ih["I"]
    gap = np.ones(shape=(gaps[i], I.shape[1], I.shape[2]))
    output = np.r_[output, gap]
    output = np.r_[output, I]

# Create final time array
print("Creating output time array")
pointer = times[0][0]
dt = TimeDelta(dt)
final_t = []
for i in range(output.shape[0]):
    final_t.append(pointer)
    pointer += dt

# Save
print("Saving Data")
np.savez(outname, I = output, t = Time(final_t))
