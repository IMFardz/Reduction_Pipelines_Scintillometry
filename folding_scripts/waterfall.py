import numpy as np
import baseband
import sys
import astropy.time as time
import astropy.units as u


def waterfall_by_sec(filename, seconds, time_step=0.002*u.s, Output_name="output", start_time=None):
   """Produce a waterfall given a from the baseband data
   === Parameters ===
   filename: Name of baseband data -> Code modelled after vdif format, hopefully should work fo
   mark5b too
   seconds: The number of seconds produced in the waterfall
   time_step: The length of each time bin in the waterfall: Defaults to 2 microseconds
   Output_name: The name of the Output waterfall: Defaults to output.npy
   Start: The start time for the waterfall. Defaults to the start of the scan
   """
   fh = baseband.open(filename)
   # Determine number of seconds in each time step
   dt = 0
   fh.seek(0)
   st = fh.time
   fn = fh.time
   while (time_step > fn - st):
      fh.read(512)
      fn = fh.time
      dt += 512
   print("Length of each time bin: ", fn - st)
   print("Number of bins in each time bin: ", dt)
   fh.seek(0)
   # Set start time:
   if start_time is None:
       start_time = fh.time
   fh.seek(start_time)
   print("Start time of waterfall; ", start_time)
   # Make one time bin at a time
   I = [] 
   while (seconds > fh.time - start_time):
        data = fh.read(dt) 
        data = data.reshape(-1, 512) # Not entirely confident about this line: Rebins nearby values 
        ft = np.fft.rfft(data, axis=-1)
        ft = np.fft.fftshift(ft)
        spec = ft.sum(axis=0)
        I.append(spec)
   return np.savez(Output_name, I = np.array(I))
    
