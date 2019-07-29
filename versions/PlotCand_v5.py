#!/usr/bin/env python
# A separate function to extract and plot
# heimdall candidate 
# This script is a modified version of the heimdall plotting scipt 'trans_freq_time.py' 
# 

import os,sys,math
import numpy as np
import glob
from itertools import chain
from os.path import basename
from itertools import tee, izip
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
#plt.ioff()
import psrchive as psr
from sigpyproc.Readers import FilReader

def find_nearidx(array,val):
        idx = (np.abs(array-val)).argmin()
        return idx

def dedispblock(ar,lodm,hidm):
    fpsr = psr.Archive_load(ar)
    toplot = []
    dmstep = 1
    dmrange = range(lodm,hidm,dmstep)
    for dm in dmrange:
        fpsr.remove_baseline()
        fpsr.set_dispersion_measure(dm)
        fpsr.dedisperse()
        ds = fpsr.get_data().squeeze()
        w = fpsr.get_weights().flatten()
        w = w/np.max(w) # Normalized it
        idx = np.where(w==0)[0]
        ds = np.multiply(ds, w[np.newaxis,:,np.newaxis]) # Apply it
        ds[:,idx,:] = np.nan
        data1 = ds[0,:,:]
        time = np.nanmean(data1[:,:],axis=0)
        toplot.append(time)

    tbin = float(fpsr.integration_length()/fpsr.get_nbin())        
    taxis = np.arange(0,fpsr.integration_length(),tbin)
    taxis=taxis*1000 #Get to msec
    toplot = np.array(toplot)
    toplot = [list(i) for i in zip(*toplot)]
    toplot = np.transpose(toplot)
    return toplot,taxis

def negDMplot(ar,FTdirection,nchan):
    fpsr = psr.Archive_load(ar)
    fpsr.remove_baseline()
    ds = fpsr.get_data().squeeze()
    w = fpsr.get_weights().flatten()
    w = w/np.max(w) # Normalized it
    idx = np.where(w==0)[0]
    ds = np.multiply(ds, w[np.newaxis,:,np.newaxis]) # Apply it
    ds[:,idx,:] = np.nan
    data = ds[0,:,:]
    if FTdirection == 'nT': 
    	ndata = data[...,::-1]
	print "Will be flipped in Time"
    elif FTdirection == 'nF': 
	ndata = data[::-1,...]
	print "Will be flipped in freq" 
    elif FTdirection == 'nTnF': 
	ndata = data[::-1,::-1]
	print "Will be flipped in time and freq"
    else: 
	ndata = data 
	print "No flip"
    return ndata

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def plotParaCalc(snr,filter,dm,fl,fh,tint,nchan):
        #Extract block factor plot in seconds
        extimefact = 1.0

        # Total extract time Calc
        # Extract according to the DM delay    
        cmd = 'dmsmear -f %f -b %f -n 2048 -d ' % (fl+(fh-fl)/2,fh-fl) + str(dm) + " -q 2>&1 "
        p = os.popen(cmd)
        cand_band_smear = p.readline().strip()
        p.close()
        extime = extimefact/2 + extimefact*float(cand_band_smear)
        if extime < 1.0: extime = 1.0

        # Tbin calc
        # For Filter widths startting from 2^0 to 2^12=4096
        #widths = [2048,2048,2048,1024,1024,512,512,256,256,128,128,64,32]
        #tbin = widths[filter]
        bin_width = tint * (2 ** filter)
	
	
	#So that we have at least 4 bins on pulse
	if filter <= 4 and snr > 20:
	        tbin = 4.0*int(extime / bin_width)
	else:
		tbin = 2.0*int(extime / bin_width)
	
        if tbin < 16:
            tbin = 16

	if tint > (extime/tbin):
	   tbin = int(extime/tint)	

        #Fbin Calc
        fbin = int(round(math.pow(float(snr)/4.0,2)))

	
	#if nchan is not power of 2, get fbin modulo of nchan
	i=0
	while nchan%(fbin+i): i+=1
	fbin+=i

	if fbin<16:
	    i=0
            while nchan%(16+i): i+=1
            fbin=i+16

        #fbin_base2 = int(round(math.log(fbin,2)))
        #fbin = pow(2,fbin_base2)
        if fbin > 512:
            #fbin = 512
	    i=0 
            while nchan%(512-i): i+=1
            fbin=512-i

	'''
	if(nchan%float(fbin)):
		# If the fbin is not modulo of number of channel, we select closed modulo nchan		
		fbin = nchan/2**np.argmin([abs(float(fbin)-nchan/2**i) for i in range(0,10)]) 
	'''

        # Fraction of extraction to plot each time calc
        if tbin>512:
            frac = np.linspace(0,1,np.ceil(tbin/512.0)) 
        else:
            frac = np.array([0,1])

        return tbin,fbin,extime,frac  

def extractPlotCand(fil_file,frb_cands,noplot,fl,fh,tint,Ttot,kill_time_range,kill_chans,source_name,nchan):
        
        # Half of this time will be subtracted from the Heimdall candidate time
        extimeplot = 1.0

        if(noplot is not True):
                if(frb_cands.size >= 1):
                        if(frb_cands.size>1):
                                frb_cands = np.sort(frb_cands)
                                frb_cands[:] = frb_cands[::-1]
                        if(frb_cands.size==1): frb_cands = [frb_cands]
                        for indx,frb in enumerate(frb_cands):
                                time = frb['time']
                                dm = frb['dm']
	    			filter = frb['filter']
				width = tint * (2 ** filter)*(10**3) # Width in msec
				snr = frb['snr']

                                tbin,fbin,extime,frac=plotParaCalc(snr,filter,dm,fl,fh,tint,nchan)
                                #print tbin,fbin,extime,frac
        
                                stime = time-(extimeplot*0.1) # Go 10% back data

                                if(stime<0): stime = 0
				if(stime+extime>=Ttot): extime=Ttot-stime
                                #if(any(l<=stime<=u for (l,u) in kill_time_ranges)):
                                if(any(l<=time<=u for (l,u) in kill_time_range) or extime < 1.0):
                                        print "Candidate inside bad-time range"
                                else:
                                        if(indx<1000):
                                                candname = '%04d' % (indx) + "_" + '%.3f' % (time) + "sec_DM" + '%.2f' % (dm) 
                                                cmd = "dspsr -cepoch=start -N %s" % (source_name) + \
                                                        " -b " + str(tbin) +   \
                                                        " -S " + str(stime) +  \
                                                        " -c " + str(extime) + \
                                                        " -T " + str(extime) + \
                                                        " -D " + str(dm) + \
                                                        " -O " + candname + " -e ar " + \
                                                        fil_file  
                                                print cmd        
                                                os.system(cmd)                 
                                                # If kill_chans, do first manual and then an automatic smoothing for remaining channels
                                                temp = ""
                                                if kill_chans:
                                                    for k in kill_chans: 
						    	if(k!=2048): temp = temp +" "+str(k)
                                                    cmd = "paz -z \"" + temp       + "\" -m %s.ar" % (candname)
                                                    print cmd
                                                    os.system(cmd)
						    cmd = "paz -r -b -L -m %s.ar" % (candname) 	
                                                    os.system(cmd)

                                                cmd = "pam --setnchn %d -m %s.ar" % (fbin,candname)
                                                print cmd
                                                os.system(cmd)

						# Correct the variable baseline, this script writes out .norm files 	
						cmd = "running_mean_sub %s.ar" % (candname)
						os.system(cmd)
                                                
                                                ar = candname + ".norm"

                                                lodm = int(dm-(dm*0.15))
                                                if lodm < 0: lodm = 0
                                                hidm = int(dm+(dm*0.15))
                                                print "Producing DM vs phase plot"
                                                toplot,taxis = dedispblock(ar,lodm,hidm)
                                                #plt.show()

                                                for i,j in pairwise(frac):
                                                    #New plotting technique
                                                    cmd = "psrplot -N 1x3 -p flux -p freq -p freq " + \
                                                          " -j ':1:dedisperse,F %d' -j ':2:F %d' " % (int(fbin),int(fbin)) + \
                                                          " -j :0:dedisperse -j :0:fscrunch " + \
                                                          " -c ':0:x:range=(%f,%f)' -c ':1:x:range=(%f,%f)'" % (i,j,i,j) + \
                                                          " -c ':2:x:range=(%f,%f)'" % (i,j) + \
                                                          " -c ':1:y:view=(0.1,1.13)' -c ':2:y:view=(0.13,1.08)'" + \
                                                          " -c ':0:set=pub,below:l=SNR: %.2f,ch=2,below:r=Wid: %.2f'" % (float(snr),float(width))  + \
                                                          " -c ':0:above:c=$file'" + \
                                                          " -c ':1:set=pub,above:c= ,ch=2,y:reverse=1'" + \
                                                          " -c ':2:set=pub,above:c= ,ch=2'" + \
                                                          " -c ':2:x:unit=ms,y:reverse=1' " + \
                                                          " -c ':1:cmap:map=heat' -c ':2:cmap:map=heat' -c ':1:crop=0.9' -c ':2:crop=0.9'" + \
                                                          " -D %s_%.2f.ps/cps %s.norm " % (candname,i,candname)

                                                    print cmd
                                                    os.system(cmd)

                                                    #DM vs time plot
                                                    plt.rcParams["figure.figsize"] = (20,8)
                                                    plt.rcParams.update({'font.size':22})
                                                    #plt.set_cmap('gray')
						    fig1 = plt.figure(1)	 
						    ax1 = fig1.add_subplot(1,1,1)
                                                    ax1.set_xlabel("Time (msec)")
                                                    ax1.set_ylabel("DM")
                                                    lti = int(np.floor(len(taxis)*i))
                                                    lhi = int(np.ceil(len(taxis)*j))
                                                    if(lhi>=len(taxis)): lhi = len(taxis)-1
                                                    lt = taxis[lti]
                                                    lh = taxis[lhi]
                                                    ax1.imshow(toplot[:,lti:lhi],extent=[lt,lh,lodm,hidm],origin='lower',aspect='auto')
						    #For testing
						    #ax1.axvline(((extime-0.5)/2 - 0.1)*1000,0,1,ls='--',lw=3,c='red')	
						    #ax1.axvline((time-stime+0.673478*0.39)*1000,0,1,ls='--',lw=3,c='red')
						    #ax1.set_xlim(lt,lh)
						    

                                                    pngfile = candname + ".dmspace_%.2f.png" % (i)
						    fig1.savefig(pngfile,format='png',bbox_inches='tight')	
						    #Negative DM plots if required      
                                                    FTdirection = source_name.split("_")[0]
                                                    if FTdirection in ['nT','nF','nTnF']:
							print "Will be plotting original axis direction" 
							fig2 = plt.figure(2)
							ax2=fig2.add_subplot(1,1,1)
							ax2.set_xlabel("Time (msec)")	
                                                        ndata=negDMplot(ar,FTdirection,fbin)
							ax2.set_ylabel("Frequency (channels)")
							ax2.imshow(ndata[:,lti:lhi],extent=[lt,lh,0,fbin],aspect='auto',interpolation='none')
							npngfile = candname + ".ndm_%.2f.png" % (i)
							fig2.savefig(npngfile,format='png',bbox_inches='tight')
							cmd = "convert -rotate 90 %s_%.2f.ps +append \( -trim -resize 560x700 %s +append \) -append \( -trim -resize 560x700 %s +append \) -append %s_%.2f.pdf" % (candname,i,npngfile,pngfile,candname,i)
                                                        os.system(cmd)
						    else:	
                                                    	cmd = "convert -rotate 90 %s_%.2f.ps +append \( -trim -resize 560x700 %s +append \) -append %s_%.2f.pdf" % (candname,i,pngfile,candname,i)
                                                    	os.system(cmd)

		        cmd = "gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=%s_frb_cand.pdf *sec*DM*.pdf" % (source_name)	  
                        print cmd
                        os.system(cmd)
                else:


                        print "No candidate found"
                        return

if __name__ == "__main__":

    fil_file = str(sys.argv[1]) # Filterbank file
    FinalList = str(sys.argv[2]) # Final list of candidate (output of FRB_detector_Fast.py)
     
    frb_cands = np.loadtxt(FinalList,dtype={'names': ('snr','time','samp_idx','dm','filter','prim_beam'),'formats': ('f4', 'f4', 'i4','f4','i4','i4')})

    #uGMRT
    #fl = 300
    #fh = 500
    #FAST
    fl = 1100
    fh = 1500
    noplot=0
    tint=0.000163
    Ttot = 80 # Total length of the file        
    kill_time_range=[]
    kill_chans=[]
    nchan = 2048
    source_name="Fake"
    
    f = FilReader(fil_file)
    nchan = f.header['nchans']
    fch1 = f.header['fch1']
    foff = f.header['foff']
    tint = f.header['tsamp']
    Ttot = f.header['tobs']
    fh = fch1
    fl = fch1 + (foff*nchan)
    source_name = f.header['source_name']
    extractPlotCand(fil_file,frb_cands,noplot,fl,fh,tint,Ttot,kill_time_range,kill_chans,source_name,nchan)

