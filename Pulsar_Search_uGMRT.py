#!/usr/bin/env python
#
# Script to search pulsars in GMRT data
# Documentation can be found in Pulsar_seach_BL-GBT.txt file
#
# Author : Vishal Gajjar (1 June 2016) 
# Email : vishalg@berkeley.edu or vishal.gajjar2002@gmail.com
# v2 : digilab version
# v3 : RFIFIND and zapbirdie incorported
#-------------------------------------------------------------------------------------------#

import os, sys, glob
import numpy,time 
import sifting
import math
import optparse

#from ssps.group.groups import group_presto_data
#from ssps.presto.directories.base import PRESTODir

ddplanflag = 2
totime  = 0
pi = 3.1415926535897932384626433

class dedisp_plan:
    """
    class dedisp_plan(lodm, dmstep, dmsperpass, numpasses, numsub, downsamp)
       A class describing a de-dispersion plan for prepsubband in detail. 
       Imported from GBT_search pipe-line. 	 
    """
    def __init__(self, lodm, dmstep, dmsperpass, numpasses, numsub, downsamp):
        self.lodm = float(lodm)
        self.dmstep = float(dmstep)
        self.dmsperpass = int(dmsperpass)
        self.numpasses = int(numpasses)
        self.numsub = int(numsub)
        self.downsamp = int(downsamp)
        self.sub_dmstep = self.dmsperpass * self.dmstep
        self.dmlist = []  # These are strings for comparison with filenames
        self.subdmlist = []
        for ii in range(self.numpasses):
             self.subdmlist.append("%.2f"%(self.lodm + (ii+0.5)*self.sub_dmstep))
             lodm = self.lodm + ii * self.sub_dmstep
             dmlist = ["%.2f"%dm for dm in \
                       numpy.arange(self.dmsperpass)*self.dmstep + lodm]
             self.dmlist.append(dmlist)
#            print " %s " % dmlist 

def timed_execute(cmd, run_cmd=1):
    """
    timed_execute(cmd):
        Execute the command 'cmd' after logging the command
            to STDOUT.  Return the wall-clock amount of time
            the command took to execute.
    """
    sys.stdout.write("\n'"+cmd+"'\n")
    sys.stdout.flush()
    start = time.time()
    if run_cmd:  os.system(cmd) # for test all command are blocked
    end = time.time()
    return end - start


def radial_distance(RA1,DEC1,RA2,DEC2):
    """
     radial_distance: prog calculates radial distance between 
		given (RA1 DEC1) and (RA2 DEC2) position. It uses the formula 
	 	used in ATNF pulsar catalogue to calculate distances between 
		the pulsars. It is knowm as Haversine formula and it 
		can be found here. http://en.wikipedia.org/wiki/Great-circle_distance

    """
    RA1 = RA_to_rad(RA1)
    RA2 = RA_to_rad(RA2)
    DEC1 = DEC_to_rad(DEC1)
    DEC2 = DEC_to_rad(DEC2)
   
    dist = 2*math.asin(math.sqrt((math.sin((DEC1 - DEC2)/2))**2 + math.cos(DEC1)*math.cos(DEC2)*(math.sin((RA1-RA2)/2))**2))
 	   
    dist = rad_to_deg(dist)

    return dist

def RA_to_rad(RA):
    """
    RA_to_rad : converts RA which is in hh:mm:ss.ss formate into radian
    """
    RA = RA.split(":")	
    RA_H = float(RA[0])
    RA_M = float(RA[1])
    RA_S = float(RA[2])
  
    deg = (RA_H + RA_M/60 + RA_S/3600)*15
    rad = deg*pi/180.0
    return rad
  
def DEC_to_rad(DEC):
    """
    DEC_to_rad : converts radial dd:mm:ss.ss into radian
    """
    DEC = DEC.split(":")
    DEC_H = float(DEC[0])
    DEC_M = float(DEC[1])
    DEC_S = float(DEC[2])
   
    deg = (DEC_H + DEC_M/60 + DEC_S/3600)
    rad = deg*pi/180.0
    return rad


def rad_to_deg(rad):
    """
    deg_to_rad : converts radian to deg
    """
    deg = rad*180/pi
    return deg


if __name__ == "__main__":
	
	parser = optparse.OptionParser()

	parser.add_option('--ddplan',dest='ddplanflag', metavar='0/1/2',
          	help="DDplan flag,  0: To run DDplan script and determine Ddplan from it," 
		     "1: To use the hardwired DDplan default values, 2: For quick DM search, Default=2",
	        default=2, type='float')
	parser.add_option('--i',dest='inpath',metavar='INPATH',
       		help='Path and input .fil file',
	        type='string')
	parser.add_option('--o',dest='outpath',metavar='OUTPATH',
	        help='Output file Dir. If exists it will ovewrite it if not then it will create one.',
	        type='string')
	parser.add_option('--lodm',dest='lodm',metavar='LO_DM',
                help='low dm vlaue for quick DM search, only used when ddplan = 2, default = 0.0',
		default=0,type='float')
	parser.add_option('--dmstep',dest='dmstep',metavar='DM_STEP',
                help='dm step for quick DM search, only used when ddplan = 2, default = 1',
                default=1,type='float')
	parser.add_option('--ndm',dest='ndm',metavar='Num_DM',
                help='Number of DMs for quick DM search, only used when ddplan = 2, default = 150 Note: with mpi this number will be adjusted to be divisible by np',
                default=150,type='float')
	parser.add_option('--downsamp',dest='downsamp',metavar='DOWNSAMPLE',
                help='Down sample the input data for quick DM search, only used when ddplan = 2, default = 1',
                default=1,type='float')
	parser.add_option('--mpi',dest='mpiflag',metavar='MPI_FLAG',
                help='To use mpiprepsubband instead of normal prepsubband MPI_FLAG shouble be 1. default = 0 NOTE : mpiprepsubband will not work with --ddplan 0 or 1',
                default=0,type='int')
	parser.add_option('--np',dest='np',metavar='NUMBER_OF_PORT',
                help='Number of port to be used with mpiprepsubband. default = 8.',
                default=8,type='int')	

	parser.add_option('--mask',dest='maskfile',metavar='MASKFILE',
                help='.mask file to be used for the RFI flagging. No further RFI flagging will be done.',
                type='string',default='')
	
	parser.add_option('--zap',dest='zapfile',metavar='ZAPFILE',
                help='.zaplist file which has list of birdi found at zero DM. No further RFI flagging will be done.',
                type='string',default='')

	options, args = parser.parse_args()
	
	if (not options.inpath) or (not options.outpath):
	        print 'Input and output paths are required.'
        	print parser.print_help()
	        sys.exit(1)

	inpath = options.inpath
	outpath = options.outpath	
	ddplanflag = options.ddplanflag
	LODM = options.lodm
	DMSTEP = options.dmstep
	NDM = options.ndm
	DOWNSAMP = options.downsamp
	mpiflag = options.mpiflag	
	MASKFILE = options.maskfile
	zapfile = options.zapfile	

	print MASKFILE

	if(mpiflag==1):
		np = options.np - 1 

	        blk = 1000 
# 	As mpiprepsubband can not take more than 1000 dm at a time numdm will be divided in blks 
#	which should be divisible by np
	
	        blk = (int(blk/np))*np
        	tempdm = NDM % blk
	        if((tempdm % np) != 0):
        		 tempdm = (int(tempdm/np) + 1)*np
	        NDM = tempdm + int(NDM/blk)*blk

        	print "The Num DM will be adjusted to %d to be divisible by (%d - 1) nodes" % (NDM,np+1)

	HIDM =  LODM+NDM*DMSTEP

# 	Outpath can have / at the end or it can be without it. In both the cases it will process the script

	temp1 = inpath.split('/')
	temp2 = outpath.split('/')
	
	'''
	if(temp1[-1]==''): 
	# The inpath has / at the end
		indir = inpath 
	else: 
	# To add / at the end of path	
		indir = inpath + "/"
	'''	

	if(temp2[-1]==''):
                outdir = outpath
	else:
		outdir = outpath + "/"
	
	if(temp1[-1].endswith("fil")):
		filename = temp1[-1]
		indir = '/'.join(temp1[:-1]) + "/"
	else:
		print "Expected .fil file. Exit"
		sys.exit()
	
	#infile = indir + "*.sub[0-9]???"
        infile = indir + filename 
	print infile

	# Creating inf file 
	inffile = indir + '.'.join(filename.split('.')[:-1]) 
	cmd = "prepdata -numout 100 -o %s -dm 0.0 %s" % (inffile,inpath)
	prepdata_time = timed_execute(cmd,1)		

	inffile = inffile + ".inf"

	print inffile

	outfile = outdir + "Rfifind" 
	#outfile = "Rfifind"	

	if(not os.path.isdir(indir)):
		print "Input dir does not exist \n"
		sys.exit()
	if(not os.path.isdir(outdir)):
		print "Output dir does not exist. Created: " + outpath
		cmd="mkdir -p " + outpath
	        totime += timed_execute(cmd,1)
       
	fp = open(inffile,"r")
  	data = fp.read().split("\n")
	fp.close()
	lofq = float(data[15].split("=")[1])
	bandwidth  = float(data[16].split("=")[1])
	chan = int(data[17].split("=")[1])
	samptime = float(data[10].split("=")[1])
	RA = (data[4].split("=")[1])
	DEC = (data[5].split("=")[1])
	PSR = (data[3].split("=")[1])
	

#	print chan, bandwidth, lofq, samptime
#	print ddplanflag,inpath,outpath
#	print inffile[0],infile,outfile
	
	print "\nCHANNELS : %d" % (chan)
	print "Lowest Frequency : %f" % (lofq)
	print "BANDWIDTH : %f" % (bandwidth)
	print "SAMPLING TIME : %f" % (samptime)
	print "INPUT File/s :  %s" % (infile)
	print "OUTPUT Dir : %s" % (outdir)
	print "INF FILE : %s\n" % (inffile)
	
	report_file = outdir + "report.log"
        rfp = open(report_file,"w")
	rfp.write("This is report file\n")
	
	ddplan = []
        ddplans = {1:[]} # Here 1 is just the key of the variables 
 
	if(ddplanflag == 0):
		outddfile = outdir + "DDplan.out"

		cmd = "DDplan.py -f %s -b %s -n %s -t %f -r 0.2 -o %s -d 1000 > %s" % \
		       (lofq+(bandwidth/2), bandwidth, chan,samptime,outdir + "DDplan.eps",outddfile)

		timed_execute(cmd,1)
   
		fp = open(outddfile,"r")
		data = fp.read().split("\n")

#       	This is to take direct plan by runing DDplan.py command and using the output to do prepsubband
#               The order in which arg are passed ==> lodm  dmstep  dms/call  num_calls  num_sub  downsamp
#   		these lists will seperate the ddplan details and read the importain parts of it. It cant be used 
#		with mpiprepsubband		
		
		test1 = []
		test2 = []
		test3 = []
		for ii in range(0,len(data)-1):
		  for jj in data[ii].split(' '):
		   if(jj != '' ):
		       test1.append(jj)
		  test2.append(test1)
		  test1 = []
		for ii in range(0,len(test2)-1):
		    if(test2[ii] ==  ['Low', 'DM', 'High', 'DM', 'dDM', 'DownSamp', '#DMs', 'WorkFract']):
# 		       This reads all the raws with the corresponding DM values to make ddplan
		       for jj in range(ii+1,len(test2)-1):       
		          if(test2[jj] != []):
		            test3.append(test2[jj])

#		Here test3 array has all the important information which will be passed on 
#		to dedisp_plan to make suitable variables (i.e. ddplan.lodm etc...)
#					
		for ii in test3:
                    ddplans[1].append(dedisp_plan(float(ii[0]),float(ii[2]),float(ii[4]),1,chan,float(ii[3])))

 
	if(ddplanflag == 1): 
#        This is for hard coded DDplan for RSPA with HBA. Cant be used with mpiprepsubband  
		print "Hardwired values of the DDplan will be used\n"

                ddplans[1].append(dedisp_plan(0.00,0.02,7720,1,chan,1))
                ddplans[1].append(dedisp_plan(154.4,0.05,2655,1,chan,2))
                ddplans[1].append(dedisp_plan(287.15,0.1,2635,1,chan,4))
                ddplans[1].append(dedisp_plan(550.65,0.2,2247,1,chan,8))

	if(ddplanflag == 2):
		print "Quick DM search on the data with DM range of %.1f to %.1f with DM step of %.1f" % (LODM,HIDM,DMSTEP)	
		ddplans[1].append(dedisp_plan(LODM,DMSTEP,NDM,1,chan,DOWNSAMP))

	ddplans = ddplans[1]

	print "\nLODM DMSTEP DMPERPASS DOWNSAMP NUMPASS"
        for ddplan in ddplans:
                 print ddplan.lodm, ddplan.dmstep, ddplan.dmsperpass, ddplan.downsamp, ddplan.numpasses
	
	
        # rfifind (block ~2.6s)
        #cmd="rfifind -blocks 2048 -o %s %s > %s" % (outfile, infile, outfile+"_output.log")
	if(MASKFILE == ''):
		cmd="rfifind -time 2.0 -o %s %s > %s" % (outfile, infile, outfile+"_output.log")
        	print "Running RFI excision..."

        	rfi_time = timed_execute(cmd,1)
		rfi_time = 1
		rfp.write("RFI masking Time : %.2f\n" % (rfi_time))
        	totime += rfi_time
		maskfile=outfile+"_rfifind.mask"
	else:
		maskfile = MASKFILE		

	print "\nMASK FILE : %s\n" % (maskfile)
	print "\nZapfile : %s\n" % (zapfile)

	dmstrs = []
	totdm = 0
	tot_prep_time = 0

#       Prepsubband outfile name

        prepsubbandout = outdir + "prepsubband_out" 

	if(mpiflag != 1):
		for ddplan in ddplans:
			for passnum in range(ddplan.numpasses):
			     for dmstr in ddplan.dmlist[passnum]:
		                dmstrs.append(dmstr)
			print "Running prepsubband ...."
        	        totdm = ddplan.dmsperpass*ddplan.numpasses
 	
			for jj in range(0,((totdm/1000)+1)):
        	              if((totdm - jj*1000)>1000): # Because prepsubband can handle only 1000 dm values

                	          cmd = "prepsubband -nsub %d -mask %s -runavg -noclip -lodm %.2f -dmstep %.2f -numdms \
					    1000 -numout 835200 -downsamp %d -o %s %s >> %s " % (32*(jj+1),maskfile, \
                        	            ddplan.lodm+ddplan.dmstep*jj*1000 ,ddplan.dmstep,ddplan.downsamp, \
                                	    prepsubbandout,infile, prepsubbandout + ".log")

			      	  prepsubband_time = timed_execute(cmd,1)

				  tot_prep_time += prepsubband_time	

				  rfp.write("prepsubband time for DD range %.2f - %.2f w %.2f : %.2f\n" % \
				  (ddplan.lodm+ddplan.dmstep*jj*1000,(ddplan.lodm+ddplan.dmstep*jj*1000) + \
				  ddplan.dmstep*1000,ddplan.dmstep,prepsubband_time))

        	                  totime += prepsubband_time
				  prepsubband_time = 0

		 	      else:

        	                 # cmd = "prepsubband -mask %s -runavg -noclip -lodm %.2f -dmstep %.2f -numdms %d \
				 #	-downsamp %d -o %s %s >> %s " % (maskfile, \
                                 #           ddplan.lodm+ddplan.dmstep*jj*1000,ddplan.dmstep,totdm - jj*1000, \
                                 #           ddplan.downsamp,prepsubbandout,infile, prepsubbandout+ ".log")
				
				 # New V3 with numout otherwise it generates in realfft ERROR	
				  cmd = "prepsubband -nsub %d -mask %s -runavg -noclip -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -numout 835200 -o %s %s >> %s " % (32*(jj+1),maskfile, \
					ddplan.lodm+ddplan.dmstep*jj*1000,ddplan.dmstep,totdm - jj*1000, \
					ddplan.downsamp,prepsubbandout,infile, prepsubbandout+ ".log")

				 # cmd = "prepsubband -nsub %d -mask %s -runavg -noclip -lodm %.2f -dmstep %.2f -numdms %d \
				 #	    -downsamp %d -o %s %s >> %s " % (32*(jj+1),maskfile, \
                	         #           ddplan.lodm+ddplan.dmstep*jj*1000,ddplan.dmstep,totdm - jj*1000, \
				 # 	     ddplan.downsamp,prepsubbandout,infile, prepsubbandout+ ".log")
				  prepsubband_time = timed_execute(cmd,1)

				  tot_prep_time += prepsubband_time

                	          rfp.write("prepsubband time for DD range %.2f - %.2f w %.2f : %.2f\n" % \
                        	  (ddplan.lodm+ddplan.dmstep*jj*1000,(ddplan.lodm+ddplan.dmstep*jj*1000) + \
				  ddplan.dmstep*(totdm - jj*1000),ddplan.dmstep,prepsubband_time))

                     	          totime += prepsubband_time

                        	  prepsubband_time = 0

	if(mpiflag==1):

#	  For mpi prepsubband to work the number of DM should be divisible by number of nodes/cores - 1 that can be used.
#	  Here for given number of DMs, it will round off to the nearest value where it matches with the number divisible

        	for ddplan in ddplans:
                	for passnum in range(ddplan.numpasses):
	                     for dmstr in ddplan.dmlist[passnum]:
        	                dmstrs.append(dmstr)
                	print "Running mpiprepsubband ...."
	                totdm = ddplan.dmsperpass*ddplan.numpasses
        	        for jj in range(0,((totdm/blk)+1)):
                	      if((totdm - jj*blk)>blk): # Because prepsubband can handle only 1000 dm values

                        	  cmd = "mpirun -np %d mpiprepsubband -runavg -mask \
					     %s -noclip -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o %s %s >> %s " % \
                                	    (np + 1,maskfile, \
	                                    ddplan.lodm+ddplan.dmstep*jj*blk ,ddplan.dmstep,blk,ddplan.downsamp, \
        	                            prepsubbandout, infile, prepsubbandout + ".log")

                	          prepsubband_time = timed_execute(cmd,1)
                        	  tot_prep_time += prepsubband_time

               		          rfp.write("mpiprepsubband time for DD range %.2f - %.2f w %.2f : %.2f\n" % \
	                          (ddplan.lodm+ddplan.dmstep*jj*blk,(ddplan.lodm+ddplan.dmstep*jj*blk)+ ddplan.dmstep*blk \
				  ,ddplan.dmstep,prepsubband_time))
        	                  
				  totime += prepsubband_time

                	          prepsubband_time = 0

	                      else:

        	                  cmd = "mpirun -np %d mpiprepsubband -runavg -mask %s \
					    -noclip -lodm %.2f -dmstep %.2f -numdms %d -downsamp %d -o \
					    %s %s >> %s " % (np + 1, maskfile, \
                	                    ddplan.lodm+ddplan.dmstep*jj*blk,ddplan.dmstep,totdm - jj*blk, \
					    ddplan.downsamp,prepsubbandout, infile, prepsubbandout + ".log")

	                          prepsubband_time = timed_execute(cmd,1)
        	                  tot_prep_time += prepsubband_time

                	          rfp.write("mpiprepsubband time for DD range %.2f - %.2f w %.2f : %.2f\n" % \
                        	  (ddplan.lodm+ddplan.dmstep*jj*blk,(ddplan.lodm+ddplan.dmstep*jj*blk) + \
				  ddplan.dmstep*tempdm,ddplan.dmstep,prepsubband_time))

	                          totime += prepsubband_time
        	                  
				  prepsubband_time = 0

	base = prepsubbandout
	
	rfp.write("Total prepsubband time : %.2f\n" % (tot_prep_time))

	# Single pulse search
	cmd="single_pulse_search.py -t 10 %s > %s" % (base + "_DM*.dat", outdir + "singlepulse.log")
	print "Searching for single pulses..."
	sp_time = timed_execute(cmd,1)
	totime += sp_time
	rfp.write("Single_pulse_search time : %.2f\n" % (sp_time))
	# Periodicity search
	datfiles=glob.glob(outdir + "*_DM*.dat")
	test = outdir + "*_DM*.dat"
	print test  
	print "Running periodicity search..."

	# cleaning previous log files
	cmd="rm -f %s %s %s" % (base+"_realfft.log", base+"_rednoise.log", base+"_accelsearch.log")	
	totime += timed_execute(cmd,1)	

	counter=0
	numharm = 8
#	All the candidate above this level will be consider as hits
	sigma = 6

	zmax = 0

	flo = 1

	realfft_time = 0
	readnoise_time = 0
	accelsearch_time = 0	
	zaptime  = 0
	for dat in datfiles:
        	stem=dat.split(".dat")[0]
	        dm=float(stem.split("_DM")[1])
		print dat
        	cmd="realfft %s >> %s" % (dat, base+"_realfft.log")
		temp = timed_execute(cmd,1)	
		totime += temp
		realfft_time += temp 

		if(zapfile != ''):
			cmd = "zapbirds -zap -zapfile %s %s.fft" % (zapfile,stem)
			temp = timed_execute(cmd,1)
                	totime += temp		
			zaptime += temp			

	        cmd="rednoise %s.fft >> %s" % (stem, base+"_rednoise.log")
		temp = timed_execute(cmd,1)
		totime += temp
		readnoise_time += temp			
        	cmd="mv %s_red.fft %s.fft" % (stem, stem)
		totime += timed_execute(cmd,1)
	        cmd="accelsearch -numharm %d -sigma %f -zmax %d -flo %f %s.fft >> %s" % \
		(numharm, sigma, zmax, flo, stem, base+"_accelsearch.log")
		temp = timed_execute(cmd,1)
		accelsearch_time += temp
		totime += temp		
	
	rfp.write("Realfft total time : %.2f\n" % (realfft_time))
	rfp.write("Rednoise total time : %.2f\n" % (readnoise_time))
	rfp.write("Accelsearch total time : %.2f\n" % (accelsearch_time))
	
# Following will sort out the candidates 

        numhits_to_fold = 2 
        low_DM_cutoff = 1.0
	lo_accel_zmax = zmax  

#       read_candidate will generate collective information about found candidate in 
#	in all ACCEL files.    	
	lo_accel_cands = sifting.read_candidates(glob.glob(outdir + "*ACCEL_0"))
#	Remove candidates with same period and low significance. 
	if len(lo_accel_cands):
        	lo_accel_cands = sifting.remove_duplicate_candidates(lo_accel_cands) 
        if len(lo_accel_cands):
    	        lo_accel_cands = sifting.remove_DM_problems(lo_accel_cands, numhits_to_fold, dmstrs, low_DM_cutoff)
        if len(lo_accel_cands):
        	lo_accel_cands.sort(sifting.cmp_sigma)
	        sifting.write_candlist(lo_accel_cands,outdir+"/candidate_list_from_script")
	
	fp = open(outdir+"candidate_list_from_script","r")
	candfile = fp.read().split("\n")
	fp.close()

#    The actual path of the psr_cats should be specified here. 

	fp = open("/home/gajjar/SP_search_wrapper/PulsarSearch/psr_cats.txt","r")
	realpsr = fp.read().split("\n")
	fp.close()	

#       This part is to findout in all the candidates found how many of them are known strong pulsars. 
# 	In addition to that it also calculates the radial distances to these psr from the given data field
#       center. The Candidate DM and realpsr DM will be compared withing range of 5 and period will be 
#       consider in the range of 0.5 msec 
 
	test1 = []
	test2 = []
	
	for ii in range(0,len(realpsr)-1):
		for jj in realpsr[ii].split(' '):
		     if(jj != '' ):
                  	test1.append(jj)	
		test2.append(test1)
		test1 = []
	
	realpsr = test2
		
        realcand = []
        harmcand = []
	normalcand = []
	numcand = 0
	prepfold_time = 0

	for ii in range(len(candfile)):

#	only lines which has ACCEL word in it will be read from the candidate sorted file
		if 'ACCEL' in candfile[ii].split(" ")[0]:

		  temp = candfile[ii].split(" ")
		  temp1 = []

#	To remove unwanted spaces 	

		  for jj in range(len(temp)):
			  if (temp[jj] != ''):
                             temp1.append(temp[jj])	

#   	DM and period for the candidates sorted from the ACCEL search

		  DM = float(temp1[1])
		  Sigma = float(temp1[3])
		  period = float(temp1[7])/1000

#	Every candidate will be consider as possible candidate. If its period/harmonic and DM matches with any known pulsar
#       then flag will be removed. 

		  normalcandflag = 1	

#	Only the first 10 strong candidate will be folded

		  if(counter < 10): 
			  #cmd = "prepfold -mask %s -runavg -p %f -dm %f -noxwin -nosearch -o %s %s >> %s" % \
			  #			(maskfile,period,DM,outdir+"prepfold",infile,outdir+"prepfold.log")
			  cmd = "prepfold -mask %s -p %f -dm %f -noxwin -nosearch -o %s %s >> %s" % \
						(maskfile,period,DM,outdir+"prepfold",infile,outdir+"prepfold.log")
			  temp = timed_execute(cmd)
			  prepfold_time += temp
			  totime += temp
			  print "DM %f Period %f\n" % (DM, period)
 			  counter+=1
  
		  for kk in range(4,len(realpsr)-1):
		
			DM = float(DM)
			rDM = float(realpsr[kk][6]) 
		       	rPeriod = float(realpsr[kk][5])
			period = float(period)

#	If the DMs are in +/- 5 range and periods are in +/- 0.5 msec range

			if((DM > rDM - 5) and (DM < rDM + 5)):
			     if(abs(period - rPeriod) < 0.0005 ):
				   normalcandflag = 0
			  	   test1 = []
				   test1.append(realpsr[kk][1])
				   test1.append(realpsr[kk][5])
				   test1.append(realpsr[kk][6])
				   RA2 = realpsr[kk][3]
				   DEC2 = realpsr[kk][4]
				   dist = radial_distance(RA,DEC,RA2,DEC2)
				   test1.append(dist) 
				   
				   realcand.append(test1)
				   
				   print "period %f and rPeriod %f\n" % (period,rPeriod)			
				   print "%f" % (abs(period - rPeriod))	  
						     	
#			     For harmonic search actual candidate period will be compared upto 
#			     10th Harmonic of known pulsar period
 				
			     for jj in range(1,10):
				if(abs(period - rPeriod/jj) < 0.0005):
				  normalcandflag = 0	
				  test1 = []
				  test1.append(jj)
				  test1.append(period)	
				  test1.append(realpsr[kk][1])
			          test1.append(realpsr[kk][5])		 			
				  test1.append(realpsr[kk][6])
				  test1.append(DM)	
				  RA2 = realpsr[kk][3]
				  DEC2 = realpsr[kk][4]
				  dist = radial_distance(RA,DEC,RA2,DEC2)	
				  test1.append(dist)
				  test1.append(Sigma)	

				  harmcand.append(test1)

	
		  if(normalcandflag == 1):
			test1 = []
			numcand+=1		 
			test1.append(numcand)
			test1.append(DM)
			test1.append(period)
			test1.append(Sigma)

			normalcand.append(test1)	
	
	rfp.write("Total prepfold time : %.2f\n" % (prepfold_time))	
		
#	This will write the detail report file containing information about detection and 
# 	radial distances of the PSR detected. 

	str = "\n\nThe following are the results\n The Observered PSR for the given field was :" + PSR + "\n The detected PSRs are as follows \n"
	rfp.write(str)
	for ii in realcand:
	    str = ii[0] + " which has period " + ii[1] + " and DM of " + ii[2] + "\n" 
	    rfp.write(str)

	str = "\n Following are the Harmonics detected\n"
	rfp.write(str)
	str = "Harm_Num  PSR           PSR_period  DM        Detected_period   Detected_DM   Radial_Distance     Sigma\n"
	rfp.write(str)

	for ii in harmcand:
	    str = "%d         %s      %s    %s     %f           %.2f          %.2f             %.2f\n" % (ii[0],ii[2],ii[3],ii[4],ii[1],ii[5],ii[6],ii[7])
	    rfp.write(str)
	
	rfp.write("\n\nFollowing are the UPOs (Unidentified Pulsating Objects)\n")

	str = "Num  Detected_period   Detected_DM     Sigma\n"
	rfp.write(str)

	for ii in normalcand:
	    str = "%d         %f          %.2f          %.2f\n" % (ii[0],ii[2],ii[1],ii[3]) 
	    rfp.write(str)
	
        rfp.write("\n\nTotal time taken by the script : %f\n" % (totime))

	rfp.close()

	print realcand 
#	print harmcand	
#	print totime


