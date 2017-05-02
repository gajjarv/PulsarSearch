"""Check the quality of a .mask file based on how much of it is flagged for RFI."""
"""Command line execution should include this file name, 
a .mask file as the next argument, and chanfrac as the final argument.
ex: python rfi_quality_check.py -mask Lband_rfifind.mask"""
import rfifind_bandpass_on
import os
from functools import reduce
import numpy as np
from time import sleep
import sys

def rfi_check(base_name, mask_file, time, nchans, tsamp, chanfrac):
	#Create rfifind class and set the channels to zap.
	a = rfifind_bandpass_on.rfifind(mask_file)
	channels = float(nchans)
	intervals = float(tsamp) / float(time)
	a.get_bandpass()
	# Believed to be unnecessary --> a.set_zap_chans(100.0, 5.0, 2.0, False, 0.01, True, [])
	zapped = a.mask_zap_chans_per_int
	int_times = a.times
	os.system("mv pgplot.ps {0}.ps".format(base_name)) #Rename the .ps file that is outputted from the bandpass plotting function.

	total_zaps = 0
	
#	for item in zapped:
#		total_zaps += item.size
	zapcount = np.zeros(int(channels))
	data = list(np.concatenate(zapped,axis=0))
		
	zapcount,edges = np.histogram(data,int(channels))
	print zapcount 

	#kill_chans = list(reduce(np.intersect1d, zapped))
	avg_bad_chans = len(kill_chans) / channels #Give the average fraction of channels flagged across all intervals

	chan_killfile_name = base_name + "_chan.kill"
	chan_killfile = open(chan_killfile_name, "w")
	for i in range(int(channels)):
	        if i<10:
	                chan = "000" + str(i)
	        elif i<100:
	                chan = "00" + str(i)
	        elif i <1000:
	                chan = "0" + str(i)
	        else:
	                chan = str(i)
	        if i in kill_chans:
	                chan_killfile.write(chan + "          0\n")
	        else:
	                chan_killfile.write(chan + "          1\n")
	chan_killfile.close()

	time_killfile_name = base_name + "_time.kill"
	time_killfile = open(time_killfile_name, "w")
	bad_ints = 0
	for i in range(len(int_times)):
		item = int_times[i]
		interval = zapped[i]
		if interval.size / channels >= chanfrac:
			time_killfile.write("{0}    {1}\n".format(str(item), str(item + time)))
			bad_ints += 1
	time_killfile.close()

	percentage_bad_ints = float(bad_ints) / intervals

	#Calculate the percentage of channel-interval spaces that have been zapped.
	percent_flagged = (float(total_zaps) / (channels * intervals)) * 100
	print(str(percent_flagged) + "% of the data was flagged for RFI.")
	return percent_flagged, percentage_bad_ints


	
