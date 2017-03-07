"""Create a .bird file."""
"""
Make sure that the ACCEL file and .fft file are in the same directory in which you run this script.

Example: create_bird.py -name Lband -accel Lband_topo_DM0.00_ACCEL_0 -fft Lband_topo_DM0.00.fft
"""
import numpy as np 
import os, sys
from argparse import ArgumentParser

def read_accel_file(path):
	"""Read in an ACCEL file and obtain data from it."""

	with open(path) as f:
		raw_lines = [line.split() for line in f]

	num_harmonics = []
	frequencies = []
	widths = []

	done_reading = False
	i = 0
	while not done_reading:
		if raw_lines[i] == []:
			done_reading = True
		elif i >= 3:
			harm_count = raw_lines[i][4] 
			num_harmonics.append(harm_count)
			freq_and_width = raw_lines[i][6]
			l_index = freq_and_width.index('(')
			r_index = freq_and_width.index(')')
			frequency = freq_and_width[: l_index]
			width = freq_and_width[l_index + 1: r_index]
			frequencies.append(frequency)
			widths.append(width)	
		i += 1

	assert len(frequencies) == len(widths), "There are an unequal number of frequencies and widths."

	for i in range(len(frequencies)):
		places = len(frequencies[i]) - frequencies[i].index('.') - 1
		new_width = float(widths[i]) / (10 ** places)
		widths[i] = str(new_width)

	frequencies = ["#Freq"] + frequencies
	widths = ["Width"] + widths
	num_harmonics = ["#harm"] + num_harmonics
	grow = ["grow?"] + (['0'] * (len(frequencies)-1))	
	bary = ["bary?"] + (['0'] * (len(frequencies)-1))

	print("Gathering data from ACCEL file.")

	data = [frequencies, widths, num_harmonics, grow, bary]
	data = np.array(data)
	data = data.T
	#Transpose the data, so that it is in columns.
	return data

def write_birdie_file(data, datafile_path):
	"""Write the data into a .birds file."""
	bird_file = open(datafile_path, 'w+')
	#Open the ascii file.

	col_width = max(len(word) for row in data for word in row) + 2  #Padding in the .txt file
	for row in data:
		bird_file.write("".join(word.ljust(col_width) for word in row)+"\n")
		#Code is written line by line

	#Close the file.
	bird_file.close()

def zap(name, bird_file, fft_file):
	"""Execute the zaplist commands on the .birds and other files."""
	os.system("cp " + name + "_rfifind.inf " + name + ".inf")
	os.system("makezaplist.py " + bird_file)

	os.system("zapbirds -zap -zapfile " + name + ".zaplist " + fft_file)

def go(name, accel_file, fft_file):
	"""Run the entire script."""

	print("\n")
	print("Birdwriter:")

	assert 'ACCEL' in accel_file, "You did not pass an ACCEL file as your first argument."
	assert '.fft' in fft_file, " You did not pass an .fft file as your first argument."

	path = os.getcwd() + '/'
	accel_path = path + accel_file

	data = read_accel_file(accel_path)

	bird_file = name + ".birds"
	output_path = path + bird_file
	write_birdie_file(data, output_path)
	zap(name, bird_file, fft_file)

	print("Birdie & Zap File Created!")

if __name__ == "__main__":
	parser = ArgumentParser(description = "Parser for name, ACCEL File, and .fft File")
	parser.add_argument("-name", action='store', dest='name', required=True, type=str,
	                help="A proper filename to use to create the .bird file")
	parser.add_argument("-accel", action='store', dest='accel_file', required=True, type=str,
	                help="Name of the ACCEL file")
	parser.add_argument("-fft", action='store', dest='fft_file', required=True, type=str,
	                help="Name of the .fft_file file name")
	args = parser.parse_args()

	name = args.name
	accel_file = args.accel_file
	fft_file = args.fft_file
	go(name, accel_file, fft_file)









