import astropy.io.fits as pf
import numpy as np
from astropy.time import Time
import astropy.units as u, astropy.constants as c
import glob

'''
A quick selection of code for reading in psrfits files and 
generating dynamic spectra.
'''

def __init__():
    return

class psrfits_data:
    '''
    A class for reading in psrfits files and generating a dynamic
    spectrum.

    Parameters:
    -----------
    fname: string, the full path to the psrfits file you would like to
           be read in.
    remove_dispersion: Boolean, optional, True if you would like the channel
           offsets due to dispersion to be removed.
    fref: astropy quantity, optional.  The frequency you would like to 
           dedispersion to.  If not passed, but remove_dispersion is True,
           the data will be dedispersed to a frequency of inf.
    '''

    def __init__(self,fname,remove_dispersion=False,fref=None):
        self.read(fname,remove_dispersion,fref)
            
    def read(self,fname,remove_dispersion,fref):
        '''
        Reads in the psrfits file.
        
        Parameters:
        -----------
        fname: string, the full path to the psrfits file you would like to
           be read in.
        remove_dispersion: Boolean, optional, True if you would like the channel
           offsets due to dispersion to be removed.
        fref: astropy quantity, optional.  The frequency you would like to 
           dedispersion to.  If not passed, but remove_dispersion is True,
           the data will be dedispersed to a frequency of inf.

        Sets:
        ----------
        self.data: numpy array, the flux table
        self.time: array of astropy times, the center time of each profile
        self.tint: astropy quantity, the integration time of each time bin
        self.freq: astropy quantity, the frequency of each channel
        '''
        with pf.open(fname) as hd:
            if hd['PRIMARY'].header['FITSTYPE'] != 'PSRFITS':
                raise ValueError('Not a psrfits file.')
            if hd['HISTORY'].data['NBIN'] == 0:
                raise ValueError('This data is in search mode. '
                              'Try the pypsrfits module by  '
                              'Paul Demorest, or use dspsr  '
                              'to create fold mode psrfits  '
                              'files.')
            self.freq = hd['SUBINT'].data['DAT_FREQ'][0]*u.MHz
            raw_data = hd['SUBINT'].data['DATA']
            shape    = raw_data.shape

            dat_offs = hd['SUBINT'].data['DAT_OFFS']
            dat_wts  = hd['SUBINT'].data['DAT_WTS']
            dat_scl  = hd['SUBINT'].data['DAT_SCL']
            self.data = (raw_data * np.repeat(dat_scl.reshape((shape[0],shape[1],shape[2],1)),
                                                  shape[3],axis=3) + 
                         np.repeat(dat_offs.reshape((shape[0],shape[1],shape[2],1)),
                                    shape[3],axis=3)) 
            self.data = self.data / dat_wts[:,np.newaxis,:,np.newaxis]
            self.data = self.data[:,:2,...].sum(axis=1)
            
            self.time = (hd['SUBINT'].data['OFFS_SUB']*u.s + 
                         Time(hd['PRIMARY'].header['STT_IMJD'],format='mjd') + 
                         hd['PRIMARY'].header['STT_SMJD']*u.s +
                         hd['PRIMARY'].header['STT_OFFS']*u.s)
            self.tint = (hd['SUBINT'].data['TSUBINT'])*u.s
            
            DM = hd['COHDDISP'].header['DM'] *u.pc *u.cm**(-3)
            f0 = float([h for h in hd['PSRPARAM'].data['PARAM'] if 'F0' in h][0].split()[1])
            dt_bin = 1/(f0*hd['SUBINT'].header['NBIN'])*u.s

            for i in np.arange(len(self.freq)):
                if fref is None:
                    dt = 4.18808*10.**(-3.)*u.s * DM.to(u.pc*u.cm**(-3)).value * (-1./self.freq[i].to(u.GHz).value**2.)
                else:
                    dt = 4.18808*10.**(-3.)*u.s * DM.to(u.pc*u.cm**(-3)).value * (1./fref.to(u.GHz).value**2.-1./self.freq[i].to(u.GHz).value**2.)
                nshift = np.int(np.rint(dt/dt_bin))
                self.data[:,i,:] = np.roll(self.data[:,i,:],nshift,axis=-1)
            
    def dynspec(self,offbins,onbins):
        '''
        Produces a dynamic spectrum from self.data
        
        Parameters:
        -----------
        onweights: array, float; the weight of each phase bin to generate 
           the on-pulse dynamic spectrum
        offweights: array, float; the weight of each phase bin to generate the
           off-pusle dynamic spectrum
           
        Returns:
        ---------
        the dynamic spectrum, un-normalized
        '''
        Ion  = self.data[:, :, onbins[0]:onbins[1]].mean(axis=-1)
        Ioff = self.data[:,:,offbins[0]:offbins[1]].mean(axis=-1)
        return (Ion/Ioff - 1)
