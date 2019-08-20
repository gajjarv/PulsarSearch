#!/usr/bin/env python
"""
Take .fil file as an input file and does single pulse search
"""
import matplotlib
matplotlib.use('pdf')
from argparse import ArgumentParser
import os,sys,math
import rfi_quality_check
import bird_wrapper
import timeit
import optparse
from rfi_filter_vg import rfi_filter as rfi
import numpy as np
import glob
from itertools import chain
sys.path.insert(0,'/home/vgajjar/SP_search_wrapper/PulsarSearch/robert_sp/')
import sp_cand_find as sp
import smtplib
from os.path import basename
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import COMMASPACE, formatdate
import subprocess as sb
from sigpyproc.Readers import FilReader
from PlotCand import extractPlotCand

def dm_delay(fl, fh, DM):
    kdm = 4148.808 # MHz^2 / (pc cm^-3)
    return abs(kdm * DM * (1.0 / (fl * fl) - 1 / (fh * fh)))

def emailsend(send_to,subject,msgtxt,files,candtxt):
	msg = MIMEMultipart()
	msg['From'] = 'blfrbcand@gmail.com'
	msg['To'] = send_to
	msg['Date'] = formatdate(localtime=True)
	msg['Subject'] = subject

	username = 'blfrbcand@gmail.com'
	password = 'breakthrough'

	msg.attach(MIMEText(msgtxt))
	if candtxt:
		with open(candtxt) as fp:
			msg.attach(MIMEText(fp.read()))

	for f in files or []:
		with open(f, "rb") as fil:
			part = MIMEApplication(fil.read(),Name=basename(f))
			part['Content-Disposition'] = 'attachment; filename="%s"' % basename(f)
			msg.attach(part)

	server = smtplib.SMTP('smtp.gmail.com:587')
	server.ehlo()
	server.starttls()
	server.login(username,password)
	server.sendmail(msg['From'], msg['To'], msg.as_string())
	server.quit()
	print("email sent")


def candplots(fil_file,source_name,snr_cut,filter_cut,maxCandSec,noplot,minMem,kill_chans,kill_time_range,nogpu,gcands,fl,fh,tint,Ttot,nchan):
	if(nogpu is not True):
		#os.chdir(basedir)
		#os.system("cd %s" % (basedir))
		#print "Inside : %s" % (basedir) 
		os.system("rm *_all.cand")
		os.system("rm *.ar *.norm")
		os.system("coincidencer *.cand")	
		os.system("trans_gen_overview.py -cands_file *_all.cand")
		os.system("mv overview_1024x768.tmp.png %s.overview.png" % (source_name))
		os.system("frb_detector_bl.py  -gdm 6 -cands_file *_all.cand -filter_cut %d -snr_cut %f -max_cands_per_sec %f -min_members_cut %f -verbose" % (filter_cut,snr_cut,maxCandSec,minMem))
		os.system("frb_detector_bl.py  -gdm 6 -cands_file *_all.cand -filter_cut %d -snr_cut %f -max_cands_per_sec %f -min_members_cut %f  > FRBcand" % (filter_cut,snr_cut,maxCandSec,minMem))
		if(os.stat("FRBcand").st_size is not 0):
			frb_cands = np.loadtxt("FRBcand",dtype={'names': ('snr','time','samp_idx','dm','filter','prim_beam'),'formats': ('f4', 'f4', 'i4','f4','i4','i4')})
		else:
			print "No candidate found"
			return
		#print frb_cands['time'],frb_cands['dm']
	else:
		if(gcands is not ""):
			os.system("rm *.ar *.norm")
			dt = np.dtype(dtype={'names': ('snr','time','samp_idx','dm','filter','prim_beam'),'formats': ('f4', 'f4', 'i4','f4','i4','i4')})
			frb_cands = np.zeros(len(gcands),dt)
			for i,dd in enumerate(gcands):
				frb_cands[i] = np.array([(dd.sigma,dd.time,dd.sample,dd.dm,dd.dfact,0)][0],dt)
		else:
			print "No candidate found"
			return

	extractPlotCand(fil_file,frb_cands,noplot,fl,fh,tint,Ttot,kill_time_range,kill_chans,source_name,nchan)

def heimdall_run(fil_file,dmlo,dmhi,base_name,snr_cut,dorfi,kill_chan_range):

	print "Running Heimdal with %f to %f DM range" % (lodm,hidm)
	#Test 
	#os.system("heimdall -zap_chans 1775 1942 -f %s -dm_tol 1.01 -dm %f %f -boxcar_max %f -output_dir %s/  -v" % (fil_file,dmlo,dmhi,boxcar_max,base_name));
	#Orig
	if dorfi is True:
		zapchan = ""
		#print kill_chan_range
		for r in kill_chan_range:
			zapchan = zapchan + " -zap_chans " + r 
		# After talking to AJ and SO and after much testing I found that 'rfi_no_narrow' works better. 
		cmd = "heimdall -f %s -rfi_tol 10 -dm_tol 1.15 -dm_pulse_width 100 -scrunching_tol 1.01 -rfi_no_narrow -dm_nbits 32 -dm %f %f -boxcar_max %f -output_dir %s -v %s" % (fil_file,dmlo,dmhi,boxcar_max,outdir,zapchan)		
		print cmd
		os.system(cmd)
	else:
		# After talking to AJ and SO
		cmd = "heimdall -f %s -dm_tol 1.15 -rfi_tol 10 -dm_pulse_width 100 -scrunching_tol 1.01 -rfi_no_narrow  -dm_nbits 32 -dm %f %f -boxcar_max %f -output_dir %s  -v" % (fil_file,dmlo,dmhi,boxcar_max,outdir)	
		print cmd
		os.system(cmd);
	return

def PRESTOsp(fil_file,dmlo,dmhi,outdir,snr_cut,zerodm,mask_file,base_name,nosearch):
	
	print "Running PRESTO with %f to %f DM range" % (lodm,hidm)	

	#prepsubband can only take 500 DMs to  dedisperse
	dmr = int(500)

	#DM Step
	dmstep = 1
	#Number of subbands
        nsub = 4096	

	if not nosearch:
		cmd = "rm *.dat prepsub*.inf *.singlepulse" 
		os.system(cmd)

	if zerodm:
		if mask_file:
			cmd = "prepsubband %s -lodm 0 -nsub %d -numdms 1 -dmstep 1 -mask %s -o prepsubband" % (fil_file,nsub,mask_file)	
		else:
			cmd = "prepsubband %s -lodm 0 -nsub %d  -numdms 1 -dmstep 1 -o prepsubband" % (fil_file,nsub)	

		#if not nosearch: os.system(cmd)
		os.system(cmd)

	for d in range(int(dmlo),int(dmhi),dmr):
		if(d+dmr>dmhi): dmhi1 = dmhi
		else: dmhi1 = d+dmr
		if(d>1000): dmstep = dmstep + 1
		if mask_file:
			cmd = "prepsubband %s -lodm %f -numdms %d -nsub %d -dmstep %d -mask %s -o prepsubband" % (fil_file,d,dmhi1-d,nsub,dmstep,mask_file)
		else:
			cmd = "prepsubband %s -lodm %f -numdms %d -nsub %d -dmstep %d -o prepsubband" % (fil_file,d,dmhi1-d,nsub,dmstep)
		if not nosearch: os.system(cmd)

	#Run single pulse search on the .dat files
	tint = 1 # Seconds. Time window to compare candidates 
	nhits_max = 12000
	dm_min = dmlo
	dm_max = dmhi

	fullfile = '%s_orig_cand.txt' %base_name
	allfile  = '%s_clean_cand.txt' %base_name
	tmp_file = '%s_FRBcand.txt' %base_name

	cmd = "single_pulse_search.py -m 0.01 -b -t %f prepsubband*.dat" % (snr_cut)

	if not nosearch: os.system(cmd)	

	cmd = "cat prepsubband*.singlepulse > All_cand.singlepulse"
	os.system(cmd)
	cand_file = "All_cand.singlepulse"
	sys.path.insert(0,'/home/vgajjar/SP_search_wrapper/PulsarSearch/robert_sp/')
	import sp_cand_find as sp	
	cands = sp.cands_from_file(cand_file, 0)
	print("%d total candidates" %len(cands))
	cands = sp.find_duplicates(cands, tint, 1000.0)

	if(len(cands)):
		if zerodm:
                	zdm = []
                	dupdms = np.array([ dd.dupes_dms for dd in cands ])
                	for i,d in enumerate(dupdms):
                        	if(0.0 not in d): zdm.append(i)
                	cands = [cands[ii] for ii in zdm]

	   	sp.write_cands( fullfile, cands )
    		ndupes = np.array([ dd.nhits for dd in cands ])
		#print 	dupes
    		yy = np.where( (ndupes > 0) & (ndupes <= nhits_max) )[0]
		#print len(yy)
	    	all_cands = [ cands[ii] for ii in yy ]
    		if len(all_cands): sp.write_cands( allfile, all_cands )
	
	    	dms = np.array([ dd.dm for dd in cands ])
		snrs = np.array([ dd.sigma for dd in cands ])
	    	xx = np.where( (ndupes > 0) & (ndupes <= nhits_max) & (dms >= dm_min) & (dms <= dm_max) & (snrs >= snr_cut))[0]
    		gcands = [ cands[ii] for ii in xx ]

		print("%d good candidates" %len(gcands))
		if(len(gcands)):
    			sp.write_cands(tmp_file, gcands)
	    		sp.make_nhits_plot(ndupes, nhits_max, base_name)	
		else: 	print "Nothing to plot"	
	else:
		gcands = []

	return gcands	

#def candplots_nogpu(fil_file,source_name,noplot,kill_chans,kill_time_range):	

def downsample(fil_file,inbits,inchans):

        #basename = ".".join(fil_file.split(".")[:-1])
	basename="downsampled"
	#sum_fil="/home/obs/sw/bl_sigproc/src/sum_fil"
	sum_fil="/home/vgajjar/bl_sigproc/src/sum_fil"
	if(inbits>8 and inchans < 8193):
		outf = basename + "_8bit.fil"
		cmd = sum_fil + " %s -o %s -obits 8 -qlen 10000" % (fil_file,outf)
	if(inbits>8 and inchans > 8192):
		outf = basename + "_8bit_2chan.fil"
		cmd = sum_fil + "  %s -o %s -obits 8 -fcollapse 2 -qlen 10000" % (fil_file,outf)
	if(inbits < 9 and inchans > 8192):
		outf = basename + "_2chan.fil"
		cmd = sum_fil + "  %s -o %s -obits 8  -fcollapse 2 -qlen 10000" % (fil_file,outf)
	if(inbits < 9 and inchans < 8193):
		outf = fil_file

	print cmd
	os.system(cmd)
        return outf

if __name__ == "__main__":

	parser = optparse.OptionParser()

	#parser = ArgumentParser(description = "Parser for inputs")
	parser.add_option("--fil", action='store', dest='fil_file', type=str,
                help="SIGPROC .fil file")
	parser.add_option("--outdir", action='store', dest='outdir', default="",type=str,
                help="Full Output directory where SOURCE_NAME_MJD folder will be created. (Default : .)")

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
	'''
	parser.add_option("--dozap", action='store_true', dest='sp',
                help='Use this flag to run zap (Default: Do not Run)')
	'''
	parser.add_option("--dorfi", action='store_true', dest='dorfi',
                help='Run RFIFIND (Default: Do not Run)')
	parser.add_option("--negDM", action='store_true', dest='negdm',
                help='Do all four types of negative search (Only works with SPANDAK; Default: Do not Run)')
	parser.add_option("--subBand", action='store_true', dest='negdm',
                help='Do sub-band search search (Only works with SPANDAK; Default: Do not Run)')

	parser.add_option("--lodm", action='store', dest='lodm', default=0.0, type=float,
                help="Heimdall: Low DM limit to search (Default: 0)")
        parser.add_option("--hidm", action='store', dest='hidm', default=1000.0, type=float,
                help="Heimdall: High DM limit to search (Default: 1000)")
        parser.add_option("--boxcar_max", action='store', dest='boxcar_max', default=16, type=float,
                help="Heimdall: Boxcar maximum window size to search (Default: 16)")
	parser.add_option("--nosearch", action='store_true', dest='nosearch',
                help='Do not run Heimdall (Default: Run)')

	parser.add_option("--snr_cut", action='store', dest='snr_cut', default=6.0, type=float,
                help="Post Heimdall: SNR cut for candidate selection (Default: 6.0)")	
	parser.add_option("--filter_cut", action='store', dest='filter_cut', default=16.0, type=int,
                help="Post Heimdall: Window size or filter cut for candidate selection (Default: 16.0)")
	parser.add_option("--maxCsec", action='store', dest='maxCandSec', default=2.0, type=float,
                help="Post Heimdall: Maximum allowed candidate per sec (Default: 2.0)")
	parser.add_option("--min_members_cut", action='store', dest='minMem', default=3.0, type=float,
                help="Post Heimdall: Number of required minimum memebers in a cluster for a real candidate (Default: 10.0)")
	parser.add_option("--noplot", action='store_true', dest='noplot',
                help='Do not run plot candidates (Default: Run)')
        parser.add_option("--nogpu", action='store_true', dest='nogpu',
                help='Run single_pulse_search.py (no heimdall)')
	parser.add_option("--zerodm", action='store_true', dest='zerodm',
                help='Remove zerodm candidates, (only work with nogpu)')

	parser.add_option("--email", action='store_true', dest='email',
                help='Send candidate file over email (no email)')	

	parser.add_option("--nodsamp", action='store_true', dest='nodsamp',
                help='Do not downsample; For the current 32-bit BL files,it is required (default: do downsample)')

	options,args = parser.parse_args()

	if (not options.fil_file):
                print 'Input file required.'
                print parser.print_help()
                sys.exit(1)

	nodsamp = options.nodsamp
        fil_file = os.path.abspath(options.fil_file)

	origf = FilReader(fil_file)
	inchans = origf.header['nchans']
        inbits = origf.header['nbits']

        # For 32-bit BL, downsample to 8-bit and create new file
        if(nodsamp is not True):
                if(inchans > 8192 or inbits > 8):
                        print "Running sum_fil"
                        fil_file = downsample(fil_file,inbits,inchans)

	fname = fil_file.split("/")[-1]
	print fil_file
	time = options.time
	timesig = options.timesig
	freqsig = options.freqsig
	chanfrac = options.chanfrac
	intfrac = options.intfrac
	max_percent = options.max_percent
	mask = options.mask
	#sp = options.sp
	sp = True
	lodm = options.lodm
	hidm = options.hidm	
	dorfi = options.dorfi
	snr_cut = options.snr_cut
	filter_cut = options.filter_cut
	boxcar_max = options.boxcar_max	
	maxCandSec = options.maxCandSec
	nosearch = options.nosearch
	noplot = options.noplot
	minMem = options.minMem
	nogpu = options.nogpu
	zerodm = options.zerodm
	email = options.email
	negdm = options.negdm
	#outdir = options.outdir

	if not options.outdir: outdir = os.getcwd()
	else: 
		outdir = options.outdir
		if(os.path.isdir(outdir) is not True):       
                	os.system("mkdir %s" % (outdir))
			
	#print "Output will go to %s" % (outdir)					
	
	os.chdir(outdir)

	f = FilReader(fil_file)
	nchan = f.header['nchans']
	fch1 = f.header['fch1']
	foff = f.header['foff']
	tint = f.header['tsamp']
        Ttot = f.header['tobs']

	print "\n Nchan : " + str(nchan) + \
              "\n High Freq (MHz) : " + str(fch1) + \
              "\n Chan. Bandwidth (MHz) : " + str(foff) + \
              "\n Integration time : " + str(tint) + \
              "\n Total time : " + str(Ttot) + "\n"
	
	fh = fch1
	fl = fch1 + (foff*nchan)
	
	fname = fname[:-4]
	hdr_file_name = fname+".hdr"
	#print fname,hdr_file_name	
	os.system("header {0} -source_name -tstart -tsamp -nchans > {1}.hdr".format(fil_file, fname))
	hdr_file = open(hdr_file_name, "r")
        hdr_data = hdr_file.readlines()
        source_name = hdr_data[0].strip("\n")
	if source_name == "": 
		source_name = fname.split(".")[0]
	source_name = source_name.replace(" ","_")
	#source_name = "fake"
	#print source_name

	MJD = hdr_data[1].strip("\n").replace(".", "_")
	MJDfloat = hdr_data[1].strip("\n")
	if (len(MJD) - MJD.index("_") - 1) > 4: #check to see if the MJD has more than 4 decimal places
                MJD = MJD[:MJD.index("_") + 5] #reduce the MJD to 4 decimal places
        base_name = source_name + "_" + MJD
	outdir = outdir + "/" + base_name 
	print "Output will go to %s"  % (outdir)
	if(os.path.isdir(base_name) is not True):	
		os.system("mkdir %s" % (base_name))
	basedir = os.getcwd() + "/" + base_name
	fil_file=os.path.abspath(fil_file)
	os.system("mv %s.hdr %s/" % (fname,basedir))
	os.chdir(basedir)
	if(dorfi is True):
		kill_chans,kill_chan_range,kill_time_range,mask_file = rfi(fil_file, time, timesig, freqsig, chanfrac, intfrac, max_percent, mask, sp)
	else : 
		kill_chans = []
		kill_chan_range = []
		kill_time_range = []
		mask_file = ""
	if(nogpu is not True):
		if(nosearch is not True):
			# IF running heimdall then remove old candidates 
                	os.system("rm %s/*.cand" % (outdir))
			heimdall_run(fil_file,lodm,hidm,outdir,boxcar_max,dorfi,kill_chan_range)
			#os.system("mv %s.* %s" % (base_name,base_name))
			#if(os.path.isfile("*.cand") is True):

		if filter(os.path.isfile,glob.glob("*.cand")):
			gcands = []
			candplots(fil_file,source_name,snr_cut,filter_cut,maxCandSec,noplot,minMem,kill_chans,kill_time_range,nogpu,gcands,fl,fh,tint,Ttot,nchan)
		else:	print "No heimdall candidate found"
	else:
		gcands = PRESTOsp(fil_file,lodm,hidm,outdir,snr_cut,zerodm,mask_file,base_name,nosearch)
		candplots(fil_file,source_name,snr_cut,filter_cut,maxCandSec,noplot,minMem,kill_chans,kill_time_range,nogpu,gcands,fl,fh,tint,Ttot,nchan)

	if(email is True):
		pdffile = source_name + "_frb_cand.pdf"
		pngfile = source_name + ".overview.png"
		#cmd = "convert *.ps %s" % (pdffile)
		pdffile = [pdffile]
		subject = "Broadband candidates from %s" % (base_name)
		cmd = "mjd2cal %s" % (MJDfloat)	
		p = sb.Popen([cmd],stdout=sb.PIPE,shell=True)
		mjdstr = [p.stdout.read()]
		msgtxt = "Broadband candidates found on %s" % (mjdstr[0])
		send_to = "vishal.gajjar2002@gmail.com"
		#print msgtxt,subject
		#candtxt = "DIAG_FRB130729_57913_0488_FRBcand.txt"
		candtxt = ""	
		emailsend(send_to,subject,msgtxt,pdffile,candtxt)
