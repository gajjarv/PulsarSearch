#!/usr/bin/env python

import math, numpy, sys, os, glob, time, datetime

if __name__ == "__main__":

	if len(sys.argv) != 6:
		print "function call:"
		print "  [1] path to files"
		print "  [2] repeating part in file names"
		print "  [3] starting time in s"
		print "  [4] stopping time in s"
		print "  [5] path to extract the data to"
		sys.exit()
	
	os.chdir(sys.argv[1])
	nTotalBlocs = 0
	nListBlocs = []
	nListBlocsCumul = []
	idx = 0

	flist = sorted(glob.glob(sys.argv[2] + "*.raw"))	# lists all files corresponding to entries
	nnumfiles =  len(flist)								# number of files to process
	print nnumfiles," files found in the dataset"		# print number of files

	for file in flist:				# .raw files starting with argv 2
		f = open(file,'rb')			# open as [b]inary and to [r]ead
		f.seek(0,0)					# 0=bof, 1=current, 2=eof
		currline = f.read(80)		# read 1st line
		nHeadLine = 0				# counter for number of lines in header
		# read header
		while (not 'END  ' in currline):
			# print currline
			if ('BLOCSIZE' in currline):
				subline = currline[10:]         # remove keyword from string
				subline.replace(' ', '')        # remove empty space
				nBlocsize = int(subline)        # convert string to integer
			if ('DIRECTIO' in currline):
				subline = currline[10:]         # remove keyword from string
				subline.replace(' ', '')        # remove empty space
				directio = int(subline)         # convert string to integer
			if ('CHAN_BW ' in currline):
				subline = currline[10:]         # remove keyword from string
				subline.replace(' ', '')        # remove empty space
				dchanbw = float(subline)		# convert string to float
			currline = f.read(80)
			nHeadLine = nHeadLine + 1			# count number of lines in header
		if directio == 1:
			nHeaderSize = 512*(math.floor((nHeadLine*80)/512.)+1)   # size of header
		if directio == 0:
			nHeaderSize = nHeadLine*80      # size of header
		f.close()
	
		nBlocs = numpy.ceil(os.path.getsize(file)/(nHeaderSize + nBlocsize))	# number of blocks in current file
		nTotalBlocs = nTotalBlocs + nBlocs	# total number of blocks
		nListBlocs.append(nBlocs)			# lists number of blocks per file
		nListBlocsCumul.append(nTotalBlocs)	# lists number of cumulative blocks

	# print "# of blocs = ",nTotalBlocs
	# print "total duration = ",nTotalBlocs*nBlocsize/64./4./abs(dchanbw)/1e6," seconds"

	if float(sys.argv[3]) < 0:	# verify starting time is >0
		print "starting time must be positive and < ",nTotalBlocs*nBlocsize/64./4./abs(dchanbw)/1e6," s"
		sys.exit()
	# verify stopping time is > starting time and < total duration
	if float(sys.argv[4]) <= float(sys.argv[3]) or float(sys.argv[4]) > (nTotalBlocs*nBlocsize/64./4./abs(dchanbw)/1e6):
		print "stopping time must be > ",sys.argv[3]," and < ",nTotalBlocs*nBlocsize/64./4./abs(dchanbw)/1e6," s"
		sys.exit()

	idx = 0
	while nListBlocsCumul[idx]*nBlocsize/64./4./abs(dchanbw)/1e6 <= float(sys.argv[3]):
		idx = idx + 1
	StartFile = idx
	if StartFile == 0:
		StartBlock = int(math.floor(float(sys.argv[3]) * nListBlocs[StartFile] / (nListBlocs[StartFile]*nBlocsize/64./4./abs(dchanbw)/1e6)))
	else:
		timefile = float(sys.argv[3]) - nListBlocsCumul[idx-1]*nBlocsize/64./4./abs(dchanbw)/1e6
		StartBlock = int(math.floor(timefile * nListBlocs[StartFile] / (nListBlocs[StartFile]*nBlocsize/64./4./abs(dchanbw)/1e6)))

	idx = 0
	while nListBlocsCumul[idx]*nBlocsize/64./4./abs(dchanbw)/1e6 <= float(sys.argv[4]):
		idx = idx + 1
	StopFile = idx
	if StopFile == 0:
		StopBlock = int(math.floor(float(sys.argv[4]) * nListBlocs[StopFile] / (nListBlocs[StopFile]*nBlocsize/64./4./abs(dchanbw)/1e6)))
	else:
		timefile = float(sys.argv[4]) - nListBlocsCumul[idx-1]*nBlocsize/64./4./abs(dchanbw)/1e6
		StopBlock = int(math.floor(timefile * nListBlocs[StopFile] / (nListBlocs[StopFile]*nBlocsize/64./4./abs(dchanbw)/1e6)))

	print "copy starts with file #",StartFile," - block #",StartBlock
	print "copy stops with file #",StopFile," - block #",StopBlock


	ts = time.time()
	st = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d_%H_%M_%S')

	newfilename = sys.argv[5]+st+".raw"
	print "writing data to "+newfilename
	nf = open(newfilename,'wb')	# open as [b]inary and to [w]rite

	for filenum in range(StartFile,StopFile+1):
		f = open(flist[filenum],'rb')
		if(StartFile == StopFile):
			for numblk in numpy.arange(StartBlock,StopBlock+1):
				f.seek(int(numblk*(nHeaderSize+nBlocsize)),0)
				# nf.write(f.read(int(nHeaderSize+nBlocsize)))
				nf.write(numpy.fromfile(f, dtype=numpy.uint8, count=int(nHeaderSize+nBlocsize)))
		elif(filenum == StartFile):
			for numblk in numpy.arange(StartBlock,nListBlocs[filenum]):
				f.seek(int(numblk*(nHeaderSize+nBlocsize)),0)
				# nf.write(f.read(int(nHeaderSize+nBlocsize)))
				nf.write(numpy.fromfile(f, dtype=numpy.uint8, count=int(nHeaderSize+nBlocsize)))
		elif(filenum == StopFile):
			for numblk in numpy.arange(0,StopBlock+1):
				f.seek(int(numblk*(nHeaderSize+nBlocsize)),0)
				# nf.write(f.read(int(nHeaderSize+nBlocsize)))
				nf.write(numpy.fromfile(f, dtype=numpy.uint8, count=int(nHeaderSize+nBlocsize)))
		
		else:
			for numblk in numpy.arange(0,nListBlocs[filenum]):
				f.seek(int(numblk*(nHeaderSize+nBlocsize)),0)
				# nf.write(f.read(int(nHeaderSize+nBlocsize)))
				nf.write(numpy.fromfile(f, dtype=numpy.uint8, count=int(nHeaderSize+nBlocsize)))
		
		f.close()

	nf.close()
		
