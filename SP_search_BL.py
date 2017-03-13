#!/usr/bin/env python
"""
Take in a Filterbank file (.fil) and check how much of it is flagged for RFI.
If the fraction of the data is less than the given threshold, 
then create a .bird file and .zap file.
"""
import sys
from argparse import ArgumentParser
import os
import rfi_quality_check
import bird_wrapper
import timeit
import optparse
from rfi_filter_vg import rfi_filter as rfi

def heimdall_run(fil_file,dmlo,dmhi):

	print "Running Heimdal with %f to %f DM range" % (lodm,hidm)
	
	os.system("heimdall -f %s -dm %f %f -v" % (fil_file,dmlo,dmhi));
	return

if __name__ == "__main__":

	parser = optparse.OptionParser()

	#parser = ArgumentParser(description = "Parser for inputs")
	parser.add_option("--fil", action='store', dest='fil_file', type=str,
                help="SIGPROC .fil file")
	parser.add_option("--lodm", action='store', dest='lodm', default=0.0, type=float,
                help="Heimdall: Low DM limit to search (Default: 0)")
	parser.add_option("--hidm", action='store', dest='hidm', default=1000.0, type=float,
                help="Heimdall: High DM limit to search (Default: 1000)")
	parser.add_option("--time", action='store', dest='time', default=2.0, type=float,
                help="RFIFIND: Seconds to integrate for stats and FFT calcs (Default: 2.0 seconds)")
	parser.add_option("--timesig", action='store', dest='timesig', default=10.0, type=float,
                help="RFIFIND: The +/-sigma cutoff to reject time-domain chunks (Default: 10)")
	parser.add_option("--freqsig", action='store', dest='freqsig', default=4.0, type=float,
                help="RFIFIND: The +/-sigma cutoff to reject freq-domain chunks (Default: 4)")
	parser.add_option("--chanfrac", action='store', dest='chanfrac', default=0.5, type=float,
                help="RFIFIND: The fraction of bad channels that will mask a full interval (Default: 0.5)")
	parser.add_option("--intfrac", action='store', dest='intfrac', default=0.3, type=float,
                help="RFIFIND: The fraction of bad intervals that will mask a full channel (Default: 0.3)")
	parser.add_option("--max_percent", action='store', dest='max_percent', default=20.0, type=float,
                help="Maximum percentage of flagged data allowed to pass through the filter. (Default: 20.0%)")
	parser.add_option("--mask", action='store_true', dest='mask',
                help='Use this flag to indicate whether a .mask file already exists for the given filterbank file.')
	parser.add_option("--sp", action='store_true', dest='sp',
                help='Use this flag for single-pulse searches instead of pulsar searches.')
	
		
	options,args = parser.parse_args()

	if (not options.fil_file):
                print 'Input file required.'
                print parser.print_help()
                sys.exit(1)

	fil_file = options.fil_file
	time = options.time
	timesig = options.timesig
	freqsig = options.freqsig
	chanfrac = options.chanfrac
	intfrac = options.intfrac
	max_percent = options.max_percent
	mask = options.mask
	sp = options.sp
	lodm = options.lodm
	hidm = options.hidm	
	
	rfi(fil_file, time, timesig, freqsig, chanfrac, intfrac, max_percent, mask, sp)
	heimdall_run(fil_file,lodm,hidm)

