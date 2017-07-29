from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
import os
import sys
import pickle
from collections import namedtuple
from astropy.io import fits
#from astropy.coordinates import SkyCoord
from astropy import units as u
import simple_dm as dedisp


def get_data(infile, tstart, tdur):
    hdulist = fits.open(infile)
    hdu1, hdu2 = hdulist

    tstop = tstart + tdur 
    
    dt = hdu2.header['tbin']
    freqs = hdu2.data[0]['dat_freq']
    df = freqs[1] - freqs[0]
    nchan = len(freqs)
    tsub = hdu2.data[0]['TSUBINT']
    nsblk = hdu2.header['nsblk']

    row_start = max(int(tstart / tsub), 0)
    row_stop  = max(int(np.ceil(tstop / tsub)), row_start + 1)
    trow_start = row_start * tsub

    dd = hdu2.data[row_start : row_stop]['data']
    dd = np.reshape(dd, ( (row_stop - row_start) * nsblk, nchan ))
    
    idx = max(int( (tstart - trow_start) / dt ), 0)
    ntx = int(tdur / dt)

    dd = dd[idx : idx + ntx, :]

    tt = np.arange(dd.shape[0]) * dt + trow_start + idx * dt

    hdulist.close()
    
    return tt, freqs, dd


def get_obs_info(infile):
    hdulist = fits.open(infile)
    hdu1, hdu2 = hdulist
    
    obs_info = {}
    hdr = hdu1.header
    keys = ['RA', 'DEC', 'STT_IMJD', 'STT_SMJD', 'STT_OFFS']
    for kk in keys:
        obs_info[kk] = hdr.get(kk)

    offs_sub = hdu2.data[0].field('OFFS_SUB')
    tsubint  = hdu2.data[0].field('TSUBINT')
    start_offs = offs_sub - 0.5 * tsubint
    obs_info['SEC_OFFS'] = start_offs

    stt_imjd = obs_info['STT_IMJD']
    stt_smjd = obs_info['STT_SMJD']
    stt_offs = obs_info['STT_OFFS']

    mjd_obs = np.float64(stt_imjd) + (stt_smjd + stt_offs) / (24. * 3600.0)
    mjd_dat = mjd_obs + start_offs / (24.0 * 3600.0)

    obs_info['MJD_OBS'] = mjd_obs
    obs_info['MJD_DAT'] = mjd_dat

    hdulist.close()

    return obs_info


def save_obs_info(fitsfile, basename):
    obs_info = get_obs_info(fitsfile)
    outfile = '%s_obs_info.p' %basename
    pickle.dump(obs_info, open(outfile, 'wb'))
    return


def avg_chan(dd, nchan=4):
    nsteps = int(dd.shape[1] / nchan)
    for nn in xrange(nsteps):
        dd[:,nn*nchan:(nn+1)*nchan] = \
            np.outer(np.mean(dd[:, nn*nchan:(nn+1)*nchan], axis=1), np.ones(nchan))
    return dd


def get_dedispersed_timeseries(tt, freqs, dat, dm):
    dt = tt[1] - tt[0]
    df = freqs[1] - freqs[0]
    ddm = dedisp.dedisperse_one(dat.T, dm, dt, df, freqs[0])
    return ddm


def get_ddm_spec_mod(tt, freqs, dat, dm, chan_weights=None, pad=5):
    dt = tt[1] - tt[0]
    df = freqs[1] - freqs[0]
    dd_dm = dedisp.dedisperse_dspec(dat.T, dm, dt, df, freqs[0]).T

    if chan_weights is None:
        ddm = np.mean(dd_dm, axis=1)
        xpk = np.argmax(ddm)
        spec = np.mean( dd_dm[xpk-1:xpk+2], axis=0 )
        tpk = tt[xpk]
        mod_idx = np.std(dd_dm, axis=1) / np.mean(dd_dm, axis=1)

    else:
        ddm = np.sum(dd_dm * chan_weights, axis=1) / np.sum(chan_weights)
        xpk = np.argmax(ddm)
        spec = np.mean(dd_dm[xpk-1:xpk+2], axis=0 ) * chan_weights
        tpk = tt[xpk]
        mom2 = np.sum(dd_dm**2.0 * chan_weights, axis=1) / np.sum(chan_weights)
        mom1 = np.sum(dd_dm * chan_weights, axis=1) / np.sum(chan_weights)
        mod_idx = np.sqrt((mom2 - mom1**2.0) / mom1**2.0)

    return tpk, ddm, spec, mod_idx


def make_plot(infile, tmid, tsearch, tshow, DM, beamnum=-999, 
              chan_weights=None, outfile=None, nchan_avg=1,
              snr_cut=4.0):
    tstart = tmid - tsearch * 0.5
    tt, freqs, dd = get_data(infile, tstart, tsearch)

    if nchan_avg > 1:
        dd = avg_chan(dd, nchan=nchan_avg)
    else: pass

    tpk, ddm, spec, midx = get_ddm_spec_mod(tt, freqs, dd, DM, chan_weights)
    tpk0, ddm0, spec0, midx0 = get_ddm_spec_mod(tt, freqs, dd, 0.0)

    if outfile is not None:
        plt.ioff()
    
    # Make Figure
    fig = plt.figure(figsize=(10,8))

    # Some params 
    dspec_snr = 25 * np.sqrt(nchan_avg)
    min_frac  = 0.25 #dspec_min = dspec_snr * min_frac
    cmap = plt.cm.magma_r
    #cmap = plt.cm.coolwarm

    # Dynamic Spectrum Axis
    bpad = 0.1
    wpad = hpad = 0.01
    wds = hds = 0.55
    xds = bpad 
    yds = bpad
    
    # Time Series Axis
    hts = 1.0 - 2 * bpad - hpad - hds
    wts = wds
    xts = xds
    yts = yds + hds + hpad 
    
    # Spectrum Axis
    hs = hds 
    ws = 1.0 - 2 * bpad - wpad - wds 
    xs = xds + wds + wpad 
    ys = yds 

    # Modulation Index / SNR Axis
    mod_pad = 0.035
    hm = hts - mod_pad
    wm = ws  - mod_pad 
    xm = xs  + mod_pad
    ym = yts + mod_pad

    ## Some Useful Limits ##
    df = freqs[1] - freqs[0]
    dt = tt[1] - tt[0]

    tlim = (- 0.2 * tshow,  0.8 * tshow)
    flim = (freqs[0] - 0.5 * df, freqs[-1] + 0.5 * df)

    ###  Make Dynamic Spectrum Plot ###
    ax_ds = fig.add_axes([xds, yds, wds, hds])

    if chan_weights is None:
        dd_sig = np.std(dd)
    else:
        xx = np.where(chan_weights)[0]
        dd_sig = np.std(dd[:, xx])

    #im_ext = [tt[0] - 0.5 * dt, tt[-1] + 0.5 * dt,
    #          freqs[0] - 0.5 * df, freqs[-1] + 0.5 * df]

    im_ext = [tt[0] - tpk, tt[-1] - tpk, freqs[0], freqs[-1]]

    
    ax_ds.imshow(dd.T / dd_sig, interpolation='nearest', origin='lower', 
                 aspect='auto', extent = im_ext, cmap=cmap, 
                 vmin= -1 * min_frac * dspec_snr * dd_sig, 
                 vmax = dspec_snr * dd_sig)

    # Add Dispersion Sweep
    offset = dt * 2
    dm_curve = dedisp.dm_delay(freqs, freqs[-1], DM)
    ax_ds.plot(dm_curve - offset, freqs, lw=1, c='k')
    ax_ds.plot(dm_curve + offset, freqs, lw=1, c='k')

    # Set limits
    ax_ds.set_xlim(tlim)
    #ax_ds.set_xlim(-2, 2)
    ax_ds.set_ylim(flim)

    # Add labels
    ax_ds.set_ylabel('Frequency (MHz)', fontsize=14, labelpad=10)
    ax_ds.set_xlabel('Time Offset (s)', fontsize=14, labelpad=10)

    ### Make Time Series Plot ###
    ax_ts = fig.add_axes([xts, yts, wts, hts])
    #ddm_sig = np.std(ddm)
    ddm_sig = np.std(ddm[tt > tpk + 0.01])
    snrs = ddm / ddm_sig
    xpk = np.where(tt == tpk)[0]

    ax_ts.plot(tt - tpk, snrs, c='k', ls='steps', lw=1)
    #ax_ts.plot(tt - tpk, ddm0 / ddm_sig, c='LimeGreen') 
    ax_ts.axvline(0, lw=3, c='r', alpha=0.2)
    ax_ts.axhline(y=0, lw=3, c='k', alpha=0.2)

    ax_ts.set_xlim(tlim)
    ax_ts.set_ylim(None, (int(1.1 * np.max(snrs)/5)+1) * 5)
    ax_ts.set_xticklabels([])
    ax_ts.set_ylabel('SNR', fontsize=14)
    ax_ts.set_title('De-dispersed Time Series', fontsize=14)


    ### Make Spectrum Plot ###
    ax_s = fig.add_axes([xs, ys, ws, hs])    
    ax_s.plot(spec / dd_sig, freqs, c='k', ls='steps', lw=1)
    #ax_s.plot(spec0 / dd_sig, freqs, c='LimeGreen')
    ax_s.axvline(x=0, lw=3, c='k', alpha=0.2)
    
    ax_s.set_ylim(flim)
    spec_max = (int(1.2 * np.max(spec/dd_sig)/5)+1) * 5
    spec_min = -(int(0.1 * np.max(spec/dd_sig)/5)+1) * 5
    ax_s.set_xlim(spec_min, spec_max)
    ax_s.set_yticklabels([])
    ax_s.set_xlabel('SNR', fontsize=14)
    ax_s.set_ylabel('De-dispersed Spectrum', fontsize=14, 
                    rotation=-90, labelpad=20)
    ax_s.yaxis.set_label_position("right")

    ### Make Mod Index Plot ###
    ax_m = fig.add_axes([xm, ym, wm, hm])

    ax_m.plot(midx, snrs, ls='', marker='o', c='k', alpha=0.5)
    #ax_m.plot(midx0, ddm0/ddm_sig, ls='', marker='o', c='LimeGreen', alpha=0.5)
    ax_m.plot(midx[xpk], snrs[xpk], ls='', marker='o', c='r')

    if chan_weights is None:
        midx_cut = np.sqrt(len(freqs) / nchan_avg) / snr_cut
    else:
        midx_cut = np.sqrt(np.sum(chan_weights) / nchan_avg) / snr_cut

    ax_m.text(0.95, 0.95, r"$m_{\rm I} = %.2f$" %midx[xpk],
              ha='right', va='top', transform=ax_m.transAxes)

    ax_m.axvline(x=midx_cut, c='b', lw=3, alpha=0.2)
    ax_m.axhline(y=snr_cut, c='b', lw=3, alpha=0.2)

    ax_m.set_xlim(0, 3 * midx_cut)
    ax_m.set_ylim(0, (int(1.2 * np.max(snrs)/5)+1) * 5)

    ax_m.set_ylabel('SNR', rotation=-90, labelpad=20, fontsize=14)
    ax_m.set_xlabel("Mod Index", fontsize=14)

    ax_m.yaxis.set_label_position("right")
    ax_m.xaxis.set_label_position("top")

    # Set the figure title
    title_str = r"$t_{\rm pk} = %.3f\, {\rm s}$" %(tpk) + "     " +\
                r"${\rm DM} = %.1f\, {\rm pc\, cm}^{-3}$" %(DM) + "     " +\
                r"${\rm Beam} = %d$" %(beamnum)
    fig.suptitle(title_str, fontsize=16)
    
    # Add a footer with the file name
    fig.text(0.02, 0.02, infile.split('/')[-1], 
             va='bottom', ha='left', fontsize=10)
    
    if outfile is not None:
        plt.savefig(outfile, dpi=100, bbox_inches='tight')
        plt.close()
        plt.ion()

    else:
        plt.show()
    return    


def dspec_stats(infile, tmid, tsearch, tshow, DM, chan_weights=None):
    tstart = tmid - tsearch * 0.5
    tt, freqs, dd = get_data(infile, tstart, tsearch)
    tpk, ddm, spec, midx = get_ddm_spec_mod(tt, freqs, dd, DM, chan_weights)
    tpk0, ddm0, spec0, midx0 = get_ddm_spec_mod(tt, freqs, dd, 0.0)

    dd_sig = np.std(dd)
    ddm_sig = np.std(ddm)
    snrs = ddm / ddm_sig 
    xpk = np.where(tt == tpk)[0][0]
    
    return tpk, snrs[xpk], midx[xpk]
    

def get_file_name(indir, basename, beamnum, nzeros=4):
    return "{}/{}_beam{:0{}}.fits".format(indir, basename, beamnum, nzeros)


def plots_from_candfile(candfile, fitsdir, basename, tsearch, tshow, 
                        nzeros=4, chan_weights=None, fixdm=None, nchan_avg=1,
                        snr_cut=4.0):
    dat = np.loadtxt(candfile)
    times = dat[:, 0]
    beams = dat[:, 1].astype('int')
    dms   = dat[:, 2]
    snrs  = dat[:, 3]
    N = len(dat)

    if fixdm is not None:
        dms = dms * 0 + fixdm

    print("\nMaking %d Dynamic Spectrum Plots" %N)
    for ii in xrange(N):
        print("%d%%" %(int(100.0 * (ii+1) / float(N))), end="...")
        sys.stdout.flush()
        infile = get_file_name(fitsdir, basename, beams[ii], nzeros=nzeros)
        outfile = "%s_beam%04d_T%08.3f_dspec.png" %( \
            basename, beams[ii], times[ii])
        if nchan_avg > 1:
            outfile = outfile.rstrip('.png') + "_avgchan%d.png" %nchan_avg
        make_plot(infile, times[ii], tsearch, tshow, dms[ii], beamnum=beams[ii],
                  chan_weights=chan_weights, outfile=outfile, nchan_avg=nchan_avg,
                  snr_cut=snr_cut)   
        print("\n")
    return


def plots_from_candlist(candlist, fitsdir, basename, tsearch, tshow, 
                        nzeros=4, chan_weights=None, fixdm=None, nchan_avg=1,
                        snr_cut=4.0):
    times = array_attr(candlist, 'tt')
    dms   = array_attr(candlist, 'dm')
    N = len(candlist)

    if fixdm is not None:
        dms = dms * 0 + fixdm

    print("\nMaking %d Dynamic Spectrum Plots" %N)
    for ii in xrange(N):
        print("%d%%" %(int(100.0 * (ii+1) / float(N))), end="...")
        sys.stdout.flush()
        infile = "%s/%s.fits" %(fitsdir, basename)
        outfile = "%s_T%08.3f_dspec.png" %(basename, times[ii])
        if nchan_avg > 1:
            outfile = outfile.rstrip('.png') + "_avgchan%d.png" %nchan_avg
        make_plot(infile, times[ii], tsearch, tshow, dms[ii],
                  chan_weights=chan_weights, outfile=outfile, 
                  nchan_avg=nchan_avg, snr_cut=snr_cut)    
    print("\n")
    return


def stats_from_candfile(candfile, fitsdir, basename, tsearch, tshow, 
                        nzeros=4, chan_weights=None, fixdm=None):
    dat = np.loadtxt(candfile)
    times = dat[:, 0]
    dms_det   = dat[:, 2]
    snrs  = dat[:, 3]
    N = len(dat)

    Cand = namedtuple('Cand', ['tt', 'beam', 'dm', 'dm_det', 'snr', 'fit_tt', 
                               'fit_snr', 'fit_midx'])
    candlist = []

    if fixdm is None:
        dms = dms_det
    else:
        dms = np.ones(len(dms_det)) * fixdm

    for ii in xrange(N):
        #infile = get_file_name(fitsdir, basename, 0, nzeros=nzeros)
        infile = "%s/%s.fits" %(fitsdir, basename)
        ftt, fss, fmm = dspec_stats(infile, times[ii], tsearch, tshow, dms[ii],
                                    chan_weights=chan_weights)
        cc = Cand(tt=times[ii], beam=0, dm=dms[ii], dm_det=dms_det[ii], 
                  snr=snrs[ii], fit_tt = ftt, fit_snr=fss, fit_midx=fmm)
        candlist.append(cc)
    return candlist


def cands_from_candfile(candfile, fixdm=None):
    dat = np.loadtxt(candfile)
    times = dat[:, 0]
    beams = dat[:, 1].astype('int')
    dms_det   = dat[:, 2]
    snrs  = dat[:, 3]
    N = len(dat)

    if dat.shape[-1] == 6:
        fit_snr = dat[:, 4]
        fit_midx = dat[:, 5]
    else:
        fit_snr = np.zeros(len(snrs))
        fit_midx = np.zeros(len(snrs))

    Cand = namedtuple('Cand', ['tt', 'beam', 'dm', 'dm_det', 'snr',
                               'fit_snr', 'fit_midx'])
    candlist = []

    if fixdm is None:
        dms = dms_det
    else:
        dms = np.ones(len(dms_det)) * fixdm


    for ii in xrange(N):
        cc = Cand(tt=times[ii], beam=beams[ii], dm=dms[ii], dm_det=dms_det[ii], 
                  snr=snrs[ii], fit_snr=fit_snr[ii], fit_midx=fit_midx[ii])
        candlist.append(cc)
    return candlist


def array_attr(clist, attr):
    if hasattr(clist[0], attr):
        return np.array([ getattr(cc, attr) for cc in clist ])
    else:
        print("Object has no attribute: %s" %attr)
        return


def write_select_cands(cands, outfile):
    fout = open(outfile, 'w')
    hdr = "#{:<12}{:<10}{:<10}{:<10}{:<10}{:<10}".format(\
        "Time", "Beam", "DM", "SNR", "New SNR", "Mod Index")
    fout.write(hdr + "\n")

    fit_snrs = array_attr(cands, 'fit_snr')
    idx = np.argsort(fit_snrs)
    sorted_cands = [ cands[ii] for ii in idx[::-1] ]
    
    for cc in sorted_cands:
        outstr = "{:<12.3f}{:<10d}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.2f}".format(\
            cc.tt, cc.beam, cc.dm_det, cc.snr, cc.fit_snr, cc.fit_midx)
        fout.write(outstr + "\n")
    fout.close()
    return


if __name__ == '__main__':
    pass
