import numpy as np
from baseband import vdif
import astropy.units as u
import sys

path = '/scratch/p/pen/simardda/B1133/gk049c/gb/'
filename = sys.argv[1]
thread_ids = np.arange(8)
outnames=[filename+'_pol_'+str(i) for i in range(2)]

# Characteristics of Data
nthread=8
nchan=1


try:
    fh = vdif.open(path+filename+'.vdif','rs',sample_rate=sample_rate,nthread=nthread,nchan=nchan)
    file_size = fh.fh_raw.seek(0, 2) // 2
    fh.fh_raw.seek(0)
    size = fh.seek(-1, 2)
    fh.seek(0)

    # Create stream reader for each polarization
    fos = [vdif.open(outname+'.vdif','ws',sample_rate=fh.sample_rate,nthread=4,nchan=1,
    samples_per_frame=fh.sample_rate, file_size=file_size) for outname in outnames]

    # Find each thread and make a vdif file for each polarization
    while fh.fh_raw.tell() < size:
        arr = fh.fh_raw.read_frame()
        arr.header.mutable = True
        old_thread = arr.header
        arr.header['thread_id'] = old_header // 2
        arr.header.mutable = False
        fos[old_thread % 2].write_frame(arr)
except KeyboardInterrupt:
    print('Closing files')
    fh.close()
    for i in range(len(fos)):
        fos[i].close()
    sys.exit()

fh.close()
for i in range(len(fos)):
    fos[i].close()
