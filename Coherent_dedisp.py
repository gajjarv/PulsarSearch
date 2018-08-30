#!/usr/bin/env python

from argparse import ArgumentParser
import os,sys,math
import timeit
import optparse
import numpy as np
import glob
from itertools import chain
import smtplib
from os.path import basename
import subprocess as sb
import numpy as np
#import pypulse as pp
import matplotlib.pyplot as plt
from scipy.signal import detrend 
from peakutils.baseline import baseline
import psrchive as psr
import matplotlib.ticker as ticker

def run(gonogo,cmd):
	cmd = cmd.split(" ")
	if(gonogo): 
		#os.system(cmd)
		print cmd
		p = sb.Popen(cmd)
	else:	p = sb.Popen('echo')
	return p 

def simrun(execute,cmd):
	print cmd
	if(execute):
                os.system(cmd)

  	
def thresholding_algo(y, lag, threshold, influence):
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)):
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            if y[i] > avgFilter[i-1]:
                signals[i] = 1
            else:
                signals[i] = -1

            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            avgFilter[i] = np.mean(filteredY[(i-lag):i])
            stdFilter[i] = np.std(filteredY[(i-lag):i])
        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = np.mean(filteredY[(i-lag):i])
            stdFilter[i] = np.std(filteredY[(i-lag):i])

    return dict(signals = np.asarray(signals),
                avgFilter = np.asarray(avgFilter),
                stdFilter = np.asarray(stdFilter))

def find_nearidx(array,val):
	idx = (np.abs(array-val)).argmin()	
	return idx

def plotfig(filename,freql,freqh,snr,low,high,RM,polplot,ext,dm,PAfile,minlabx,minlaby):
    #Setting
    #polplot=1
    #minlabx = 1
    #minlaby = 1
    ticksize = 6
    fontsize = 6
    #ext = 10 # Extra phase bins around the pulse
    #-----------
    
    #data = pp.Archive(filename,prepare=False)
    #data1 = data.data[0][0]
    #time = np.sum(data1,axis=0)	
    #freq = data.freq[0]
   
    fpsr = psr.Archive_load(filename)
    fpsr.dededisperse()
    #fpsr.set_dispersion_measure(565)
    fpsr.set_dispersion_measure(dm)
    fpsr.dedisperse()
    if polplot:	
	    fpsr.set_rotation_measure(RM)
	    fpsr.defaraday() 
    fpsr.fscrunch_to_nchan(152)
    #fpsr.tscrunch(128)
    #fpsr.tscrunch_to_nsub(128)
    name = fpsr.get_source()
    fpsr.remove_baseline()
    #fpsr.set_rotation_measure(RM)
    #Profile before rotation
    ds1 = fpsr.get_data().squeeze()
    data2 = ds1[0,:,:]

    print low,high
#   If low and high are given do not center 
    if(low==0 and high==127):
	    fpsr.centre_max_bin()
 
    #fpsr.defaraday()
    freq = np.linspace(fpsr.get_centre_frequency()-abs(fpsr.get_bandwidth()/2),fpsr.get_centre_frequency()+abs(fpsr.get_bandwidth()/2),fpsr.get_nchan())
    freq = freq[::-1]
    ds = fpsr.get_data().squeeze()
    nbin = int(fpsr.get_nbin())

#   If low and high are given, use them. If not then define a window
    if(low==0 and high==127): 
	    low = int(nbin - nbin/2	- nbin*0.1)
	    high = int(nbin - nbin/2 + nbin*0.1)

    # Get weights 
    w = fpsr.get_weights().flatten()
    w = w/np.max(w) # Normalized it
    idx = np.where(w==0)[0]
    ds = np.multiply(ds, w[np.newaxis,:,np.newaxis]) # Apply it
    ds[:,idx,:] = np.nan
    data1 = ds[0,:,:]
    Qdata1 = ds[1,:,:]
    Udata1 = ds[2,:,:]
    Vdata1 = ds[3,:,:]
    Iprof = np.nanmean(data1[find_nearidx(freq,freqh):find_nearidx(freq,freql),:],axis=0)
    if polplot:
	Qdata = ds[1,find_nearidx(freq,freqh):find_nearidx(freq,freql),:]
	Qprof = np.nanmean(Qdata,axis=0)
	Udata = ds[2,find_nearidx(freq,freqh):find_nearidx(freq,freql),:]
	Uprof = np.nanmean(Udata,axis=0)
	Vdata = ds[3,find_nearidx(freq,freqh):find_nearidx(freq,freql),:]		
	Vprof = np.nanmean(Vdata,axis=0)
	Lprof = np.sqrt(pow(Qprof,2) + pow(Uprof,2))

	# Baseline removal from off-pulse
	if(nbin-high>50):
		Lprof = Lprof - np.nanmean(Lprof[high+ext:nbin-1])
		Vprof = Vprof - np.nanmean(Vprof[high+ext:nbin-1])
		Iprof = Iprof - np.nanmean(Iprof[high+ext:nbin-1]) 
	else:
		Lprof = Lprof - np.nanmean(Lprof[1:low-10])
		Vprof = Vprof - np.nanmean(Vprof[1:low-10])
		Iprof = Iprof - np.nanmean(Iprof[1:low-10])
		

    #print ds.shape 
    #sys.exit(1)

    #tbin = list(data.durations/data1[0].size)[0]
    tbin = float(fpsr.integration_length()/fpsr.get_nbin())
    #taxis = np.arange(0,data.durations,tbin)
    taxis = np.arange(0,fpsr.integration_length(),tbin)
    ptime = taxis[low-ext:high+ext]
    ptime=(ptime-np.mean(ptime))*1000 # msec
    y = Iprof
    #print len(y),len(ptime)
    spec = np.nanmean(data1[:,low:high],axis=1)
    Qspec = np.nanmean(Qdata1[:,low:high],axis=1)
    Uspec = np.nanmean(Udata1[:,low:high],axis=1)
    Vspec = np.nanmean(Vdata1[:,low:high],axis=1)
    #V Check
    plt.figure(2)
    plt.imshow(Vdata1,aspect='auto',interpolation='none')
    np.savetxt('Varray.txt',Vdata1)
    plt.show()
    sys.exit(0)	
    #Check end

    #offspec = np.nanmean(data1[:,high+40:high+40+(high-low)],axis=1)	  
 
    #oname = "".join(filename.split(".")[:-1]) + ".eps"
    specname = "".join(filename.split(".")[:-1]) + "_spectra.txt"	
    specoffname = "".join(filename.split(".")[:-1]) + "_spectra_off.txt"

    if(polplot):		
	oname = filename + ".withPA.eps"
    else:
	oname = filename + ".withoutPA.eps"
   
    #specname = filename + "_spectra.txt"
	 
    #np.savetxt(specname,spec,fmt="%.2f")	
    #np.savetxt(specname,np.c_[freq,spec],fmt="%.2f %.2f")	
    #np.savetxt(specoffname,np.c_[freq,offspec],fmt="%.2f %.2f") 
    np.savetxt(specname,np.c_[freq,spec,Vspec],fmt="%.2f %.2f %.2f")
	
    plt.rcParams['axes.linewidth'] = 0.5
    plt.subplots_adjust(hspace = .001)
    plt.subplots_adjust(wspace = .001)
    ax1 = plt.subplot2grid((7,1), (1,0), rowspan=2,colspan=1)

    #plt.subplot(211)
    plt.xlim(low-ext,high+ext)
    print low,high,ext,low-ext,high+ext
    plt.setp(ax1.get_xticklabels(), visible=False)
    #plt.setp(ax1.get_yticklabels(), visible=False
    if minlaby:	ax1.set_ylabel('Flux \n (mJy)',fontsize=fontsize, fontweight='bold')
    plt.tick_params(axis='both', which='major', labelsize=ticksize)
    #ax1.yaxis.set_ticks(np.arange(0,))

    if(max(y)>=1000): ax1.yaxis.set_major_locator(ticker.MultipleLocator(1000))
    elif(max(y)>500 and max(y)<1000): ax1.yaxis.set_major_locator(ticker.MultipleLocator(500))
    elif(max(y)>200 and max(y)<500): ax1.yaxis.set_major_locator(ticker.MultipleLocator(200))
    elif(max(y)>100 and max(y)<200): ax1.yaxis.set_major_locator(ticker.MultipleLocator(100))	
    else: ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))

    ax1.set_xticks([])
    #ax1.set_xlim(min(ptime),max(ptime))
    ax1.set_ylim([min(y)-min(y)/10,max(y)+max(y)/10])
    plt.yticks(rotation=90)
    ax1.tick_params(length=1, width=0.5)
    #plt.locator_params(axis='y', nticks=4)
    ax1.text(0.03,0.92,name,horizontalalignment='left',verticalalignment='top',transform=ax1.transAxes,fontsize=8,fontweight='bold')
  	
    '''	
    ax1.plot(np.arange(1, len(y)+1), y,linewidth=0.5,color='black')
    print min(y),max(y)
    if polplot:
	ax1.plot(np.arange(1, len(y)+1), Lprof,linewidth=0.5,color='red')
	ax1.plot(np.arange(1, len(y)+1), Vprof,linewidth=0.5,color='blue')
    '''

    #plt.subplot(212)
    #plt.step(np.arange(1, len(y)+1), result["signals"], color="red", lw=2)
    ax2 = plt.subplot2grid((7,1), (3,0), rowspan=4,colspan=1)
    #plt.xlim(low-20,high+20)
    pdata = data1[:,low-ext:high+ext]
    #ptime = taxis[low-ext:high+ext]
    lowedge = (taxis[low]-np.mean(ptime))*1000 # msec
    highedge = (taxis[high]-np.mean(ptime))*1000 # msec	
    #ptime=(ptime-np.mean(ptime))*1000 # msec 
    plt.yticks(rotation=90)
    ax2.tick_params(length=1, width=0.5)
    plt.tick_params(axis='both', which='major', labelsize=ticksize)
    if minlaby: ax2.set_ylabel('Frequency \n (MHz)',fontsize=fontsize, fontweight='bold')
    if minlabx: ax2.set_xlabel('Time (msec)',fontsize=fontsize, fontweight='bold')
    #plt.axvline(lowedge, color='b', linestyle='dashed', linewidth=2)
    #plt.axvline(highedge, color='b', linestyle='dashed', linewidth=2)
    plt.axhline(freqh,color='b', linestyle='dashed', linewidth=0.5)
    plt.axhline(freql, color='b', linestyle='dashed', linewidth=0.5)
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(1000))
    print min(freq),max(freq)	
    plt.imshow(pdata,aspect='auto',cmap='binary',extent=[min(ptime),max(ptime),min(freq),max(freq)],interpolation='none')

    if polplot:
    	ax3 = plt.subplot2grid((7,1), (0,0), rowspan=1,colspan=1) 	
	if(PAfile):
		d = np.loadtxt(PAfile,skiprows=1)
		Iprof3 = d[:,3]
		Qprof3 = d[:,4]
		Uprof3 = d[:,5]
		Lprof3 = np.sqrt(pow(Uprof3,2)+pow(Qprof3,2))
		Iprof2 = np.nanmean(data2[find_nearidx(freq,freqh):find_nearidx(freq,freql),:],axis=0)
		PA = d[:,7]
		PAerr = d[:,8]
		# To fix the PSRCHIVE truncation
		for i,pa in enumerate(PA): 
			if pa>0: PA[i] = PA[i]-180
		PA = np.roll(PA,nbin/2 - Iprof2.argmax())
		PAerr = np.roll(PAerr,nbin/2 - Iprof2.argmax())
		Lprof3 = np.roll(Lprof3,nbin/2 - Iprof2.argmax())
		
		Iprof3 = np.roll(Iprof3,nbin/2 - Iprof2.argmax())
		Iprof2 = np.roll(Iprof2,nbin/2 - Iprof2.argmax())
		Lprof = np.roll(Lprof,nbin/2 - y.argmax())
		y = np.roll(y,nbin/2 - y.argmax())
		SPEED_OF_LIGHT=300000000
		PA  = PA - RM*(SPEED_OF_LIGHT/((np.mean(freq))*1000000))*(SPEED_OF_LIGHT/((np.mean(freq))*1000000)) + 360
		#PA = PA - PAconst
		#print PA,PAerr,PAconst
	else:
		PA = np.arctan(Uprof/Qprof)
		PA = PA*0.5*180/math.pi + 60 #Not sure but this is offset I am getting from pav 
	#idx2 = np.where(Lprof>0)[0]
	idx2 = np.where(Lprof>snr*np.nanstd(Lprof))[0]
	#print idx1,PA[idx1]
	
	'''	
	# Diag	
	plt.figure(2)
	#plt.scatter(np.arange(len(PA)),PA)
	plt.plot(np.arange(1, len(y)+1), y,linewidth=0.5,color='black') # Original which gets plotted in final plot 
	plt.plot(np.arange(1, len(Iprof2)+1), Iprof2,linewidth=0.5,color='blue') # From input archive file but freq selection comes later (i.e. after dedispersion and RM correction)
	plt.plot(np.arange(1, len(Iprof3)+1), Iprof3,linewidth=0.5,color='red') # Iprof from PA file 
	plt.plot(np.arange(1, len(Iprof3)+1), Lprof,linewidth=0.5,color='green') # Linear pol which gets plotted
	plt.scatter(np.arange(1, len(y)+1), PA)
	plt.show()
	sys.exit(0)
	#Diag end
	'''

	#Print mean PA
	OnPA = PA[low-ext:high+ext]
        OnPAerr = PAerr[low-ext:high+ext]
	PAmean = OnPA[np.where((OnPA<120) & (OnPA>0))]
        PAmeanerr = OnPAerr[np.where((OnPA<120) & (OnPA>0))]
	#Working 
	#print "PA: " + str(np.nanmean(PAmean)) + " +/- " +str(np.nanstd(PAmean)) + " Mean Error " +str(np.nanmean(PAmeanerr))
	# Testing
	print filename + " PA: %.2f Std %.2f Uncertainty on avg: %.2f" % (np.average(PAmean,weights=1/(PAmeanerr*PAmeanerr)), np.nanstd(PAmean), 1/np.sqrt(np.sum(1/(PAmeanerr*PAmeanerr))))
	
	plt.yticks(rotation=90)
	plt.xlim(low-ext,high+ext)
	plt.ylim(-10,110)
	ax3.yaxis.set_major_locator(ticker.MultipleLocator(80))
	ax3.set_xticks([])
	ax3.tick_params(length=1, width=0.5)
	#plt.xticks([])
	plt.tick_params(axis='both', which='major', labelsize=ticksize)
	if minlaby:
		ax3.set_ylabel('PA \n (deg)',fontsize=fontsize, fontweight='bold')
		#plt.setp(ax3.get_yticklabels(), visible=False)
  	plt.setp(ax3.get_xticklabels(), visible=False) 
	#ax3.scatter(np.arange(1, len(y)+1)[idx2], PA[idx2], linewidth=0.5,s=0.5,color='black')      
	#ax3.errorbar(np.arange(1, len(y)+1)[idx2], PA[idx2], yerr=PAerr[idx2],color='black',markersize='0.001',capsize=0.1,elinewidth=0.01,fmt="o")
	ax3.errorbar(np.arange(1, len(y)+1), PA, yerr=PAerr,color='black',markersize='0.001',capsize=0.1,elinewidth=0.01,fmt="o")
	#ax3.plot(np.arange(1, len(y)+1), PA, linewidth=2,color='black')
    else:
        ax3 = plt.subplot2grid((7,1), (0,0), rowspan=1,colspan=1)
        #ax3.spines["top"].set_visible(False)       	
	#fig.patch.set_visible(False)
	ax3.axis('off')


    print min(y),max(y)
    if polplot:
        ax1.plot(np.arange(1, len(y)+1), Vprof,linewidth=0.5,color='blue')
        ax1.plot(np.arange(1, len(y)+1), Lprof,linewidth=0.5,color='red')
    ax1.plot(np.arange(1, len(y)+1), y,linewidth=0.5,color='black')

    '''
    ax3 = plt.subplot2grid((6,5), (2,4), rowspan=4,colspan=1)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.ylim(min(freq),max(freq))
    plt.plot(spec,freq,linewidth=2,color='black')
    '''
    #plt.show()
    fig = plt.gcf()
    fig.set_size_inches(1.4,3.5) 
    plt.savefig(oname, bbox_inches='tight',dpi=300)	


if __name__ == "__main__":

	parser = optparse.OptionParser()

	parser.add_option("-f", action='store', dest='infile', type=str, help="Spliced RAW file")
	
	parser.add_option("-o", action='store', dest='outdir', default="",type=str,help="Full Output directory (Default : .)")
	
	parser.add_option("-b", action='store', dest='nbin', default=16192, type=int,help="Number of time bins (DSPSR)")
   
	parser.add_option("-F", action='store', dest='fbin', default=1024, type=int,help="Number of freq bins (DSPSR)")

	parser.add_option("-D", action='store', dest='DM', default=565, type=float,help="DM")

	parser.add_option("-t", action='store', dest='totime', default=1.5, type=float,help="Total seconds of data in raw files")
	
	parser.add_option("--nodspsr", action='store_true', dest='nodspsr',help='Do not run DSPSR (Default: Run)')

	parser.add_option("--plot", action='store_true', dest='plot',help='Plot the final profile in eps file (Default: DO not plot)')

	parser.add_option("--fl", action='store', dest='freql', default=4000, type=float,help="Lower frequency for the pulse to add (this only affect the average pulse)")

	parser.add_option("--fh", action='store', dest='freqh', default=8000, type=float,help="Higher frequency for the pulse to add (this only affect the average pulse)")

	parser.add_option("-s", action='store', dest='snr', default=0.01, type=float,help="SNR threshold to detect pulse")

	options,args = parser.parse_args()

	infile = os.path.abspath(options.infile)
	nbin = options.nbin	
	fbin = options.fbin
	DM = options.DM
	totime = options.totime
	nodspsr = options.nodspsr
	plot = options.plot
	freql = options.freql
	freqh = options.freqh
	snr = options.snr

	gonogo=1
	if(nodspsr): gonogo=0
	
	if not options.outdir: outdir = os.getcwd()
        else:
                outdir = options.outdir
                if(os.path.isdir(outdir) is not True):
                        os.system("mkdir %s" % (outdir))	

	# Settings 
	execute = 0
	nprocess = 5 # Number of parallel dspsr # Not working
	threshold = 5 # Threshold to detect pulse 
	#	
	
	if (not options.infile or len(sys.argv)== 0):
                print 'Input file required.'
                print parser.print_help()
                sys.exit(1)

	cmd = 	"time dspsr -U 1024 -T %.1f -F 29696:D -K -d 4 -b 2048  -E /home/vgajjar/FRB121102.par -s -a psrfits -e fits %s" % (totime,infile)
	simrun(gonogo,cmd)
		
	print "DSPSR Finished"		
	
	os.chdir(outdir)
	
	cmd = "/home/vgajjar/Bandpass_correction/psrtools/normalize_rms spliced.fits"
	simrun(execute,cmd)
	
	cmd = "pam --setnchn 512 --setnbin 8196 -e fschTsch -D -p spliced.norm" # Total intensity dedisp file for plot
	#cmd = "pam --setnchn 512  -e fschTsch -D -p spliced.norm" # Total intensity dedidesp file for plot
	simrun(execute,cmd)
	
	if(plot): plotfig("spliced.fschTsch",freql,freqh,snr)		
						
