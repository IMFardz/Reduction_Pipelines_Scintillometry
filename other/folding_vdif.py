import numpy as np
import matplotlib.pyplot as plt
import baseband
import sys
import astropy.time as time
import astropy.units as u


def waterfall_by_sec(filename, seconds, dt=512, Output_name="output", start_time=None):
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
   fh.seek(0)
   st = fh.time
   fh.seek(dt)
   fn = fh.time
   sub_int = fn - st
   print("Length oft each time bin: ", sub_int)
   fh.seek(0)
   # Set start time:
   if start_time is None:
       start_time = fh.timec
   fh.seek(start_time)
   print("Start time of waterfall; ", start_tiime)
   # Make one time bin at a time
   I = [] 
   while (seconds > fh.time - start_time):
        data = fh.read(dt) 
        data = data.reshape(512, -1).mean(-1) # Not entirely confident about this line: Rebins nearby values 
        ft = np.fft.rfft(data, axis=-1)
        ft = np.fft.fftshift(ft)
        #spec = ft.sum(axis=0)
        I.append(ft)
   return np.savez(Output_name, I = np.array(I))


def fold_by_sec(filename, seconds, pulse_period, time_step=None, waterFall_name="WaterFall_Plot", fold_name="Fold", offset=0):
    fh = vdif.open(filename, "rs")
    I = []
    phase = np.arange(0, 512, 512/pulse_period)
    times_gates = np.arange(0, pulse_period, pulse_period/512)
    time = []
    fh.seek(offset)
    start_time = fh.time
    print("Start time of fold:", start_time)
    if time_step == None:
        dt = 79360
    else:
        dt = fk.seek(0)
        while (time_step > (fh.tell() - start_time).to(u.s).value):
            place_holder = fh.read(50)
            dt = fh.tell()
        fh.seek(0)
    fk.seek(offset)
    while (seconds > (fh.tell() - start_time).to(u.s).value):
        data = fh.read(dt)
        # Determine the phase bin to add to
        i = 0
        while (time_gates[i] < fh.time%pulse_period):
            i += 1;
        data = data.reshape(-1, 512)
        ft = np.fft.rfft(data, axis=-1)
        ft = np.fft.fftshift(ft)
        spec = ft.sum(axis=0)
        I.append(spec)
        if phase[i] == 0:
            phase[i] = [spec]
        else:
            phase[i].append(spec)
    final_time = fh.time
        
    phase = np.array(phase)
    time_cluster = fh.seek(0)
    while fh.time < time_gates[1]:
        fh.read(dt)
        time_cluster += dt

    phase = phase.reshape(phase.shape[0], -1, phase.shape[1]//int(time_cluster), phase.shape[2]).mean(axis=1)
    np.savez(waterFall_name, I=I)
    np.savez(fold_name, I=phase, t=np.arange(start_time, final_time, (final_time-start_time)/phase.shape[1]))
    
