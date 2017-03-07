from argparse import ArgumentParser
import os
import create_bird
"""
Take in a filterbank file, mask file, and DM file, and create a .bird file and a .zap file.
"""
"""
Example: python bird_wrapper.py -name Lband -DM_name Lband_topo_DM0.00 -fil GBT_Lband_PSR.fil -mask Lband_rfifind.mask 
"""

def prepdata(dm_name, mask_file, fil_file): #prepdata command in PRESTO
	dm = dm_name[-4:-1]
	os.system("prepdata -nobary -o {0} -dm {1} -mask {2} -numout 530000 {3}".format(dm_name, dm, mask_file, fil_file))

def real_fft(dat_file): #real_fft command in PRESTO
	os.system("realfft {0}".format(dat_file))

def accel_search(dat_file, num_harm = '4', zmax = '0'): #accel_search command in PRESTO
	os.system("accelsearch -numharm {0} -zmax {1} {2}".format(num_harm, zmax, dat_file))

def create_bird_from_files(name, dm_name, fil_file, mask_file):
	prepdata(dm_name, mask_file, fil_file)

	dat_file = dm_name + '.dat' #create string for name of .dat file, product of the prepdata command

	real_fft(dat_file)

	accel_search(dat_file)

	fft_file = dm_name + '.fft' #create string for name of .fft file, product of the accel_search command

	accel_file = dm_name + "_ACCEL_" + dm_name[-4] #create string for name of accel file. 

	"""Call create_bird.py, which constructs
	the .bird file and the .zap file."""
	create_bird.go(name, accel_file, fft_file)
	#os.system("python create_bird.py -name {0} -accel {1} -fft {2}".format(name, accel_file, fft_file))

if __name__ == '__main__':
	parser = ArgumentParser(description = "Parser for Name, .fil File, and .mask File")
	parser.add_argument("-name", action='store', dest='name', required=True, type=str,
	                help="A proper filename to use in creating the .bird file")
	parser.add_argument("-DM_name", action='store', dest='dm_name', required=True, type=str,
	                help="A proper filename for the 'prepdata' command, including the band and DM")
	parser.add_argument("-fil", action='store', dest='fil_file', required=True, type=str,
	                help="The .fil file name")
	parser.add_argument("-mask", action='store', dest='mask_file', required=True, type=str,
	                help="The .mask file name")
	args = parser.parse_args()

	name = args.name
	dm_name = args.dm_name
	fil_file = args.fil_file
	mask_file = args.mask_file

	create_bird_from_files(name, dm_name, fil_file, mask_file)



