#!/usr/bin/env python
"""
Take .fil file as an input and produces given number of sub-band .fil files using bldice
"""
import numpy as np
import sys,os
from sigpyproc.Readers import FilReader


fil_file = str(sys.argv[1]) 
nsub = int(sys.argv[2])

fhdr = FilReader(fil_file)

chans = fhdr.header['nchans']
fch1 = fhdr.header['fch1']
foff = fhdr.header['foff']

subsize = int(np.floor(chans/nsub))

for i in range(nsub): 
	print fch1+foff*i*subsize,fch1 + foff*(i+1)*subsize
		
