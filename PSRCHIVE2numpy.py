#!/usr/local/python
import psrchive as psr
import matplotlib.pyplot as plt
import numpy as np
import os,sys,math

#-- input should an psrchive archive file with full stokes --#
fname = sys.argv[1]

#-- output will be an numpy array with same input name--#
dynamicf = ".".join(fname.split(".")[:-1])+".data.npy"
freqf = ".".join(fname.split(".")[:-1])+".freq.npy"
timef = ".".join(fname.split(".")[:-1])+".time.npy"

#-- settings --#
#DM
dm = 563
#output nchan
NCHAN = 128


#-- load file and dedisperse --#
fpsr = psr.Archive_load(fname)
fpsr.dededisperse()
fpsr.set_dispersion_measure(dm)
fpsr.dedisperse()

#-- frequency scrunch and remove baseline --#
fpsr.fscrunch_to_nchan(NCHAN)
fpsr.remove_baseline()

#-- apply weights for RFI lines --#
ds = fpsr.get_data().squeeze()
w = fpsr.get_weights().flatten()
w = w/np.max(w)
idx = np.where(w==0)[0]
ds = np.multiply(ds, w[np.newaxis,:,np.newaxis])
ds[:,idx,:] = np.nan

#-- Get total intensity data (I) from the full stokes --#
data1 = ds[0,:,:]

#-- Get frequency axis values --#
freq = np.linspace(fpsr.get_centre_frequency()-abs(fpsr.get_bandwidth()/2),fpsr.get_centre_frequency()+abs(fpsr.get_bandwidth()/2),fpsr.get_nchan())
freq = freq[::-1]

#-- Get time axis --#
tbin = float(fpsr.integration_length()/fpsr.get_nbin()) 
taxis = np.arange(0,fpsr.integration_length(),tbin) 
# Convert to time to msec
taxis = taxis*1000

#-- Plot Dyanamic spectra --#
plt.ylabel("Frequency (MHz)")
plt.xlabel("Time (ms)")
plt.imshow(data1,aspect='auto',interpolation='none',extent=[0,max(taxis),min(freq),max(freq)])
plt.show()

#-- save numpy array--#
np.save(dynamicf,data1)
np.save(freqf,freq)
np.save(timef,taxis)
