import astropy.io.fits as pf
import numpy as np
from astropy.time import Time
import astropy.units as u, astropy.constants as c
import glob
import os

from psrfits_analysis import psrfits_data


path = sys.argv[1]
os.chdir(path)

files=np.sort(glob.glob('gk049c_gb_no00*.ar'))

time = []
I = []
freq = None
gates = {2: (), 5: , 8:() , 9:() , 12:(), 13:() , 14:() , 15:(), 16: (), 19: (), 20: (), 21: (), 22: (), 25: (), 26: (), 27: (), 31: (), 32: ()}
thr = {0:((400, 600), (720, 785)) ,1:((400, 600), (910, 960)) ,2:((400, 600),(720, 785)) ,3:((400, 600),(910, 960)) ,4:((400, 600), (720, 785)) ,5:((400, 600),(910, 960)) ,
       6:((400, 600), (720, 785)) ,7:((400, 600),(910, 960))}

for i in range(len(files)):
    thrnum = int(files[i][-4])
    scannum = int(files[i][12:16])
    print('File_Name: {}, File_Number: {}, Total_Files: {}'.format(files[i], i, len(files)))
    pfits = psrfits_data(files[i], remove_dispersion=True)
    dynspec = pfits.dynspec(offbins=throff[thrnum][0], onbin=gates[thrnum][0])
    for j in range(len(pfits.time)):
        time += [pfits.time[j].mjd]
        I += [dynspec[i,:]]
    freq = pfits.freq.to(u.MHz).value
print(np.array(I).shape)
print(np.array(freq).shape)
print(np.array(time).shape)
np.savez('gk049c_gb'+path[-3:]+'_dynspec.npz',
         I=np.array(I),time=np.array(time),freq=np.array(freq))
