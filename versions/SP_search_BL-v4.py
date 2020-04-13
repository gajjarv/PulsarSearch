#!/usr/bin/env python
"""
Take in a Filterbank file (.fil) and check how much of it is flagged for RFI.
If the fraction of the data is less than the given threshold, 
then create a .bird file and .zap file.
"""
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

def candplots(fil_file,source_name,snr_cut,filter_cut,maxCandSec,noplot,minMem,kill_chans,kill_time_range):
    #os.chdir(basedir)
    #os.system("cd %s" % (basedir))
    #print "Inside : %s" % (basedir) 
    os.system("rm *_all.cand")
    os.system("rm *.ar")
    os.system("coincidencer *.cand")	
    os.system("trans_gen_overview.py -cands_file *_all.cand")
    os.system("mv overview_1024x768.tmp.png %s.overview.png" % (source_name))
    os.system("frb_detector_bl.py -cands_file *_all.cand -filter_cut %d -snr_cut %f -max_cands_per_sec %f -min_members_cut %f -verbose" % (filter_cut,snr_cut,maxCandSec,minMem))
    os.system("frb_detector_bl.py -cands_file *_all.cand -filter_cut %d -snr_cut %f -max_cands_per_sec %f -min_members_cut %f  > FRBcand" % (filter_cut,snr_cut,maxCandSec,minMem))
    if(os.stat("FRBcand").st_size is not 0):
        frb_cands = np.loadtxt("FRBcand",dtype={'names': ('snr','time','samp_idx','dm','filter','prim_beam'),'formats': ('f4', 'f4', 'i4','f4','i4','i4')})
    else:
        print "No candidate found"
        return
    #print frb_cands['time'],frb_cands['dm']

    #Extract block to plot in seconds
    extime = 1.0

    if(noplot is not True):
        if(frb_cands.size > 1):
            frb_cands = np.sort(frb_cands)	
            frb_cands[:] = frb_cands[::-1]	
            for indx,frb in enumerate(frb_cands):
                time = frb['time']
                dm = frb['dm']
                stime = time-(extime/2)
                if(stime<0): stime = 0
                #if(any(l<=stime<=u for (l,u) in kill_time_ranges)):
                if(any(l<=time<=u for (l,u) in kill_time_range)):
                    print "Candidate inside bad-time range"
                else:
                    os.system("dspsrfil -S %f -c %f -T %f -D %f  -O %04d_%fsec_DM%f -e ar %s" % (stime,extime,extime,dm,indx,time,dm,fil_file))
        elif(frb_cands.size):
            time = float(frb_cands['time'])
            dm = float(frb_cands['dm'])
            stime = time-(extime/2)
            if(stime<0): stime = 0
            if(any(l<=time<=u for (l,u) in kill_time_range)):
                print "Candidate inside bad-time range"
            else:

                os.system("dspsrfil -cepoch=start -S %f -c %f -T %f -D %f  -O 0000_%fsec_DM%f -e ar %s" % (stime,extime,extime,dm,time,dm,fil_file))		
        else:
            print "No candidate found"
            return
        # If no kill_chans, do an automatic smoothing
        temp = ""
        #os.system("paz -r -b -L -m *.ar")
        if kill_chans: 	
            for k in kill_chans: temp = temp +" "+str(k)
            temp = "paz -z \"" + temp	+ "\" -m *.ar"
            print temp
            os.system(temp)	
        #os.system()
        #os.system("paz -r -b -L -m *.ar")
        #os.system("paz -Z '1775 1942' -m *.ar")
        os.system("psrplot -p F -j 'D, F 32, B 128' -D %s_frb_cand.ps/cps *.ar" % (source_name))


def heimdall_run(fil_file,dmlo,dmhi,base_name,boxcar_max,dorfi,kill_chan_range):

    print "Running Heimdal with %f to %f DM range" % (lodm,hidm)
    #Test 
    #os.system("heimdall -zap_chans 1775 1942 -f %s -dm_tol 1.01 -dm %f %f -boxcar_max %f -output_dir %s/  -v" % (fil_file,dmlo,dmhi,boxcar_max,base_name));
    #Orig
    if dorfi is True:
        zapchan = ""
        #print kill_chan_range
        for r in kill_chan_range:
            zapchan = zapchan + " -zap_chans " + r 
        # After talking to AJ and SO
        #cmd = "heimdall -f %s -scrunching 1 -scrunching_tol 1.05 -rfi_tol 5 -dm_nbits 32 -dm_pulse_width 1000 -dm_tol 1.05 -dm %f %f -boxcar_max %f -output_dir %s/  -v %s" % (fil_file,dmlo,dmhi,boxcar_max,base_name,zapchan)		
        cmd = "heimdall -f %s -scrunching 1 -rfi_tol 10 -dm_nbits 32 -dm %f %f -boxcar_max %f -output_dir %s  -v %s" % (fil_file,dmlo,dmhi,boxcar_max,outdir,zapchan)		
        print cmd
        os.system(cmd)
    else:
        # After talking to AJ and SO
        os.system("heimdall -f %s -scrunching 1 -rfi_tol 10 -dm_nbits 32 -dm %f %f -boxcar_max %f -output_dir %s -v" % (fil_file,dmlo,dmhi,boxcar_max,outdir));
        #os.system("heimdall -f %s -dm_tol 1.01 -dm %f %f -boxcar_max %f -output_dir %s/  -v" % (fil_file,dmlo,dmhi,boxcar_max,base_name));
    return

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
    parser.add_option("--min_members_cut", action='store', dest='minMem', default=10.0, type=float,
            help="Post Heimdall: Number of required minimum memebers in a cluster for a real candidate (Default: 10.0)")
    parser.add_option("--noplot", action='store_true', dest='noplot',
            help='Do not run plot candidates (Default: Run)')

    options,args = parser.parse_args()

    if (not options.fil_file):
        print 'Input file required.'
        print parser.print_help()
        sys.exit(1)

    fil_file = os.path.abspath(options.fil_file)
    fname = fil_file.split("/")[-1]
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
    #outdir = options.outdir

    if not options.outdir: outdir = os.getcwd()
    else: 
        outdir = options.outdir
        if(os.path.isdir(outdir) is not True):       
            os.system("mkdir %s" % (outdir))

    #print "Output will go to %s" % (outdir)					

    os.chdir(outdir)

    fname = fname[:-4]
    hdr_file_name = fname+".hdr"
    #print fname,hdr_file_name
    os.system("header {0} -source_name -tstart -tsamp -nchans > {1}.hdr".format(fil_file, fname))
    hdr_file = open(hdr_file_name, "r")
    hdr_data = hdr_file.readlines()
    source_name = hdr_data[0].strip("\n")
    source_name = source_name.replace(" ","_")
    #source_name = "fake"
    #print source_name

    MJD = hdr_data[1].strip("\n").replace(".", "_")
    if (len(MJD) - MJD.index("_") - 1) > 4: #check to see if the MJD has more than 4 decimal places
        MJD = MJD[:MJD.index("_") + 5] #reduce the MJD to 4 decimal places
    base_name = source_name + "_" + MJD
    outdir = outdir + "/" + base_name 
    print "Output will go to %s"  % (outdir)
    if(os.path.isdir(base_name) is not True):	
        os.system("mkdir %s" % (base_name))
    basedir = os.getcwd() + "/" + base_name
    os.system("mv %s.hdr %s/" % (fname,basedir))
    os.chdir(basedir)
    if(dorfi is True):
        kill_chans,kill_chan_range,kill_time_range = rfi(fil_file, time, timesig, freqsig, chanfrac, intfrac, max_percent, mask, sp)
    else : 
        kill_chans = []
        kill_chan_range = []
        kill_time_range = []
    if(nosearch is not True):
             # IF running heimdall then remove old candidates 
        os.system("rm %s/*.cand" % (outdir))
        heimdall_run(fil_file,lodm,hidm,outdir,boxcar_max,dorfi,kill_chan_range)
    #os.system("mv %s.* %s" % (base_name,base_name))
    #if(os.path.isfile("*.cand") is True):
    if filter(os.path.isfile,glob.glob("*.cand")):
        candplots(fil_file,source_name,snr_cut,filter_cut,maxCandSec,noplot,minMem,kill_chans,kill_time_range)
    else:	print "No heimdall candidate found"

