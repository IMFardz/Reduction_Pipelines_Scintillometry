import numpy as np
from baseband import vdif
import astropy.units as u
import sys

# File path of basbeand data
path = '/Users/fardinsyed/Desktop/baseband/baseband/data/'
filename = sys.argv[1]
outnames = [filename + '_pol_' + str(i) for i in range(2)]

try:
    # Open stream reader and obtain file size
    fh = vdif.open(path+filename+'.vdif','rs')
    file_size = fh.fh_raw.seek(0, 2) // 2
    fh.fh_raw.seek(0)
    size = fh.seek(-1, 2)
    fh.seek(0)

    # Create stream writer for each polarization
    fos = [vdif.open(outname+'.vdif','ws',sample_rate=fh.sample_rate,nthread=4,
    nchan=1,samples_per_frame=fh.samples_per_frame, file_size=file_size) for outname in outnames]

    # Write into output streams frame by frame
    while fh.fh_raw.tell() < size:
        arr = fh.fh_raw.read_frame()
        arr.header.mutable = True
        old_thread = arr.header['thread_id']
        arr.header.update(thread_id = old_thread // 2)
        arr.header.mutable = False
        fos[old_thread % 2].fh_raw.write_frame(arr)

# Trying my best not to corrupt anything
except KeyboardInterrupt:
    print('Closing files')
    fh.close()
    for i in range(len(fos)):
        fos[i].close()
    sys.exit()

# Terminating
fh.close()
for i in range(len(fos)):
    fos[i].close()
