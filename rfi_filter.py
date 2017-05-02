#!/usr/bin/env python
"""
Take in a Filterbank file (.fil) and check how much of it is flagged for RFI.
If the fraction of the data is less than the given threshold, 
then create a .bird file and .zap file.
"""
"""
Example: python rfi_filter.py -fil GBT_Lband_PSR.fil -time 30.0 -timesig 10.0 -freqsig 4.0 -chanfrac 0.7 -intfrac 0.3 -max_percent 20.0 -mask -sp
"""
from argparse import ArgumentParser
import os
import rfi_quality_check
import bird_wrapper
import timeit

def rfi_filter(fil_file, time, timesig, freqsig, chanfrac, intfrac, max_percent, mask, sp):
	start = timeit.default_timer()

	#Gather name and MDJ from the .fil file name.
	# """Name and MJD are hardcoded for now, until we involve SIGPROC."""
	# pulsar_name = "DIAG_PSR_J0034-0721"
	# pulsar_MJD = "57640_130983796298"
	file_name = fil_file[:-4] #file_name is filterbank file without .fil extension
	hdr_file_name = file_name+".hdr" #create name for .hdr file
	os.system("header {0} -source_name -tstart -tsamp -nchans > {1}.hdr".format(fil_file, file_name)) #create .hdr file with header information from .fil file
	hdr_file = open(hdr_file_name, "r")
	hdr_data = hdr_file.readlines()
	source_name = hdr_data[0].strip("\n")
	MJD = hdr_data[1].strip("\n").replace(".", "_")
	if (len(MJD) - MJD.index("_") - 1) > 4: #check to see if the MJD has more than 4 decimal places
		MJD = MJD[:MJD.index("_") + 5] #reduce the MJD to 4 decimal places
	base_name = source_name + "_" + MJD #create a base filename for files that will be created from the pipeline.
	tsamp = hdr_data[2].strip("\n")
	nchans = hdr_data[3].strip("\n")
	#Call PRESTO's rfifind command using the given inputs.
	if not mask:
		os.system("rfifind -time {0} -timesig {1} -freqsig {2} -chanfrac {3} -intfrac {4} -o {5} {6}".format(time, timesig, freqsig, chanfrac, intfrac, base_name, fil_file))
	#Create a string for the name of the .mask file, to be used later.
	mask_file = base_name + "_rfifind.mask"
	#Run the rfi_check command from the rfi_quality_check.py script to see what percentage of the data is flagged.
	percentage_flagged, percentage_bad_ints = rfi_quality_check.rfi_check(base_name, mask_file, time, nchans, tsamp, chanfrac)
	#See if the percentage of data flag exceeds the maximum allowed percentage input.
	if percentage_flagged > max_percent:
		print("File is bad. Too much data flagged.")
		return 
	else:
		#Create a .bird and .zap file because the .mask file passed the filter
		if not sp:
			dm_file = base_name + "_topo_DM0.00"
			bird_wrapper.create_bird_from_files(base_name, dm_file, fil_file, mask_file)
			#os.system("python bird_wrapper.py -name {0} -DM_name {1} -fil {2} -mask {3}".format(base_name, dm_file, fil_file, mask_file))

		tsamp = float(tsamp)
		channels = float(nchans)
		intervals = tsamp / time
		summary_file = open(base_name + "_summary.txt", "w")
		summary_file.write("{0}\n".format(fil_file))
		summary_file.write("Channels: {0}\n".format(str(channels)))
		summary_file.write("Intervals: {0}\n".format(str(intervals)))
		summary_file.write("Number of Samples: {0}\n".format(str(tsamp)))
		summary_file.write("Sample Time: {0} seconds\n".format(str(time)))
		summary_file.write("Percentage of Bad Data: {0}%\n".format(str(percentage_flagged)))
		summary_file.write("Percentage of Bad Intervals: {0}%\n".format(str(percentage_bad_ints)))
		summary_file.write("This filterbank file was analyzed using the rfifind program of PRESTO.\n")
		summary_file.write("The rfifind program produced a .mask file, which was analyzed to calculate the percentage of bad data.\n")
		summary_file.write("This file was created because the percentage of bad data was below {0}%.\n".format(max_percent))
		summary_file.write("There should be five other files (_chan.kill, _time.kill, .ps, .birds, .zaplist) to accompany this summary.txt file.\n")
		stop = timeit.default_timer()
		summary_file.write("Runtime: {0} seconds.".format(str(stop - start)))
		summary_file.close()
		return


if __name__ == "__main__":

	parser = ArgumentParser(description = "Parser for inputs")
	parser.add_argument("-fil", action='store', dest='fil_file', required=True, type=str,
	                help="The .fil filename")
	parser.add_argument("-time", action='store', dest='time', required=False, default=2.0, type=float,
	                help="Integration Time in seconds. Default is 2.0 seconds.")
	parser.add_argument("-timesig", action='store', dest='timesig', required=False, default=10.0, type=float,
	                help="Time of Signal in seconds. Default is 10.0 seconds.")
	parser.add_argument("-freqsig", action='store', dest='freqsig', required=False, default=4.0, type=float,
	                help="Frequency of Signal in Hz. Default is 4 Hz")
	parser.add_argument("-chanfrac", action='store', dest='chanfrac', required=False, default=0.5, type=float,
	                help="Threshold for Fraction of Channels Flagged. Default is 0.5.")
	parser.add_argument("-intfrac", action='store', dest='intfrac', required=False, default=0.3, type=float,
	                help="Threshold for Fraction of Intervals Flagged. Default is 0.3.")
	parser.add_argument("-max_percent", action='store', dest='max_percent', required=False, default=20.0, type=float,
	                help="Maximum percentage of flagged data allowed to pass through the filter. Default is 20.0 percent.")
	parser.add_argument("-mask", action='store_true', dest='mask',
	                help='Use this flag to indicate whether a .mask file already exists for the given filterbank file.')
	parser.add_argument("-sp", action='store_true', dest='sp',
	                help='Use this flag for single-pulse searches instead of pulsar searches.')

	args = parser.parse_args()

	fil_file = args.fil_file
	time = args.time
	timesig = args.timesig
	freqsig = args.freqsig
	chanfrac = args.chanfrac
	intfrac = args.intfrac
	max_percent = args.max_percent
	mask = args.mask
	sp = args.sp

	rfi_filter(fil_file, time, timesig, freqsig, chanfrac, intfrac, max_percent, mask, sp)










