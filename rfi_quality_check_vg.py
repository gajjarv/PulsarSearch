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
from collections import Counter
from itertools import groupby
from operator import itemgetter

def rfi_check(base_name, mask_file, rfitime, nchans, totime, chanfrac,intfrac):
	#Create rfifind class and set the channels to zap.
	a = rfifind_bandpass_on.rfifind(mask_file)
	channels = float(nchans)
	intervals = float(totime) / float(rfitime)
	a.get_bandpass()
	# Believed to be unnecessary --> a.set_zap_chans(100.0, 5.0, 2.0, False, 0.01, True, [])
	zapped = a.mask_zap_chans_per_int
	# Sometime there are repeted chan numbers in a given interval, so to get unique chan number per interval
	zapped = [np.unique(x) for x in zapped]	

	int_times = a.times
	if os.path.exists("pgplot.ps"):
		os.system("mv pgplot.ps {0}.ps".format(base_name)) #Rename the .ps file that is outputted from the bandpass plotting function.
	
	kill_chans = []
	data = list(np.concatenate(zapped,axis=0))

	#Get accetable bad intervals for each channel
	nbadtime = int(intfrac*len(int_times))	
	cdata = Counter(data)

	for i in range(int(nchans)): 
		#if(cdata[i] > nbadtime): kill_chans.append(i)
		# Output channel number is inverted to the original channel ordering, following corrects it (subtrating 1 if channel goes from 0 to nchan-1)
		if(cdata[i] > nbadtime): kill_chans.append(int(nchans)-i-1)
				
	avg_bad_chans = len(kill_chans) / channels #Give the average fraction of channels flagged across all intervals
	
	kill_chans_range = []

	for k,g in groupby(enumerate(kill_chans),lambda (i,x): i-x):
		temp = map(itemgetter(1),g)
		kill_chans_range.append(str(str(temp[0])+" "+str(temp[-1])))

	#print kill_chans_range

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
	kill_time_range = []
	for i in range(len(int_times)):
		#item = int_times[i]
		#interval = zapped[i]
		if (zapped[i].size/channels) >= chanfrac:
			#print i,zapped[i].size/channels,chanfrac
			temp = "{0}    {1}".format(str(int_times[i]), str(int_times[i] + rfitime))
			time_killfile.write(temp+"\n")
			kill_time_range.append((float(int_times[i]),float(int_times[i])+rfitime))
			bad_ints += 1

	#print kill_time_range

	time_killfile.close()

	percentage_bad_ints = float(bad_ints) / intervals

	total_zaps = bad_ints*len(kill_chans)

	#Calculate the percentage of channel-interval spaces that have been zapped.
	percent_flagged = (float(total_zaps) / (channels * intervals)) * 100
	print "%.1f%% of the data (%.1f%% chans and %.1f%% interval) was flagged for RFI." % (float(percent_flagged),float(avg_bad_chans*100),float(percentage_bad_ints*100))
	return percent_flagged, percentage_bad_ints,kill_chans,kill_chans_range,kill_time_range


	
