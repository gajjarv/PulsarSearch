import numpy as np
import os
import sys
from copy import deepcopy as dcopy
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
import matplotlib as mpl
import shutil
from astropy import units as u
#from astropy.coordinates import SkyCoord
#from astropy.coordinates import Angle

from itertools import combinations

import get_dspec as sp_plt


class Pulse(object):
    def __init__(self, beam, line):
        dm, sigma, time, sample, dfact = line.split()
        dfact.rstrip('\n')
        self.beam   = beam
        self.dm     = np.float(dm)
        self.sigma  = np.float(sigma)
        self.time   = np.float(time)
        self.sample = np.int(sample)
        self.dfact  = np.int(dfact)
        self.line   = line
        self.ra     = 0
        self.dec    = 0

        self.dupes = []
        self.dupes_beams  = [self.beam]
        self.dupes_sigmas = [self.sigma]
        self.dupes_times  = [self.time]
        self.dupes_dms    = [self.dm]
        self.dupes_dfacts = [self.dfact]
        self.dupes_maxsnr = self.sigma
        self.dupes_ras    = []
        self.dupes_decs   = []

        self.dupes_nbeams = 0
        self.dupes_dist_mean = -1
        self.dupes_dist_std =  -1
        self.nhits = 0

        self.tbin  = -1
        self.dmbin = -1

        self.dm_list = []
        self.ndms    = 0

    def add_dupes(self, newpulse):
        self.dupes_beams.append(newpulse.beam)
        self.dupes_sigmas.append(newpulse.sigma)
        self.dupes_times.append(newpulse.time)
        self.dupes_dms.append(newpulse.dm)
        self.dupes_dfacts.append(newpulse.dfact)
        if newpulse.sigma > self.dupes_maxsnr:
            self.dupes_maxsnr = newpulse.sigma

    def beam_to_radec(self, coords):
        self.ra  = coords[self.beam][0]
        self.dec = coords[self.beam][1]
        dbeams = self.dupes_beams
        nbeams = len(dbeams)

        if self.nhits:
            for dd in dbeams:
                self.dupes_ras.append( coords[dd][0] )
                self.dupes_decs.append( coords[dd][1] )

    def calc_beams(self):
        self.dupes_nbeams = len(set(self.dupes_beams))
        if self.dupes_nbeams == 0:
            print self.dupes_beams

    def get_pdist_stats(self):
        if self.nhits > 1:
            idx = np.arange(self.nhits, dtype='int')
            pidx = np.array([ii for ii in combinations(idx,2)])
            ras  = np.array(self.dupes_ras)
            decs = np.array(self.dupes_decs)
            dra = (ras[pidx[:, 0]] - ras[pidx[:, 1]]) * 3600.
            ddec = (decs[pidx[:, 0]] - decs[pidx[:, 1]]) * 3600.
            dd = np.sqrt( dra**2.0 + ddec**2.0 )
            self.dupes_dist_mean = np.mean(dd)
            self.dupes_dist_std  = np.std(dd)
        else:
            self.dupes_dist_mean = 0.0
            self.dupes_dist_std = 0.0


    def bin_time(self, dt):
        self.tbin = int(self.time / dt)

    def bin_dm(self, ddm):
        self.dmbin = int(self.dm / ddm)

    def bin_time_dm(self, dt, ddm):
        self.bin_time(dt)
        self.bin_dm(ddm)

    def reset_dms(self):
        self.dm_list = []
        self.ndms = 0

    def reset_bin(self):
        self.tbin  = -1
        self.dmbin = -1
        self.dupes = []
        self.nhits = 0

    def reset_dupes(self):
        self.dupes = []
        self.nhits = 0

    def reset_all(self):
        self.reset_dms()
        self.reset_bin()
        self.reset_dupes()

    def __str__(self):
        return self.line

    def __repr__(self):
        out_str = "Pulse(T=%.3f, "    %(self.time) +\
                        "DM=%.2f, "   %(self.dm)   +\
                        "beam=%d, "   %(self.beam) +\
                        "width=%d, "  %(self.dfact) +\
                        "sigma=%.2f) " %(self.sigma)
        return out_str


def cmp_pulse(p1, p2):
    cmp_val = cmp(p1.sample, p2.sample)
    if cmp_val == 0:
        cmp_val = cmp(p1.dm, p2.dm)
        if cmp_val == 0:
            cmp_val = -cmp(p1.sigma, p2.sigma)
        else: pass
    else: pass
    return cmp_val


def cmp_pulse2(p1, p2):
    cmp_val = cmp(p1.sample, p2.sample)
    if cmp_val == 0:
        cmp_val = -cmp(p1.sigma, p2.sigma)
        if cmp_val == 0:
            cmp_val = cmp(p1.beam, p2.beam)
        else: pass
    else: pass
    return cmp_val


def cmp_bins(p1, p2):
    """
    Sort by time bin and DM bin
    """
    cmp_val = cmp(p1.tbin, p2.tbin)
    if cmp_val == 0:
        cmp_val = cmp(p1.dmbin, p2.dmbin)
        if cmp_val == 0:
            cmp_val = -cmp(p1.sigma, p2.sigma)
        else: pass
    else: pass
    return cmp_val


def cmp_tbin_sigma(p1, p2):
    """
    Sort by time bin then by sigma
    """
    cmp_val = cmp(p1.tbin, p2.tbin)
    if cmp_val == 0:
        cmp_val = -cmp(p1.sigma, p2.sigma)
    else: pass
    return cmp_val


def cmp_beam_tbin_sigma(p1, p2):
    """
    Sort by beam, time bin, then sigma
    """
    cmp_val = cmp(p1.beam, p2.beam)
    if cmp_val == 0:
        cmp_val = cmp(p1.tbin, p2.tbin)
        if cmp_val == 0:
            cmp_val = -cmp(p1.sigma, p2.sigma)
        else: pass
    else: pass
    return cmp_val


def cands_from_file(infile, beam):
    """
    Return list of Pulse objects from the file "infile".
    Assumes this file is in the same format as a PRESTO
    *.singlepulse file.
    """
    plist = []
    for line in file(infile):
        if line[0] == "#":
            continue
        else: pass
        plist.append(Pulse(beam, line))
    return plist


def cands_from_many_files(beam_nums, cands_dir):
    """
    Assumes all the files have the name format "beam%03d.cands"
    and are located in cands_dir
    """
    plist = []
    for bnum in beam_nums:
        sp_file = "%s/beam%04d.cands" %(cands_dir, bnum)
        # check that file exists
        if not os.path.isfile(sp_file):
            sp_file = "%s/beam%03d.cands" %(cands_dir, bnum)
            if not os.path.isfile(sp_file):
                print("File not found: %s" %sp_file)
                continue
            else: pass
        else: pass
        plist += cands_from_file(sp_file, bnum)
    return plist


def find_duplicates(plist, dt, ddm):
    # Bin the candidates
    for pp in plist:
        pp.bin_time_dm(dt, ddm)
    # Sort the candidates
    pl = dcopy(sorted(plist, cmp=cmp_bins))
    # Loop through and add duplicates
    tbin  = -1
    dmbin = -1
    top_ii = -1
    for ii, pulse in enumerate(pl):
        if (pulse.tbin == tbin and pulse.dmbin == dmbin):
            #pl[top_ii].dupes.append(pulse.beam)
            pl[top_ii].add_dupes(pulse)
        else:
            top_ii = ii
            tbin   = pulse.tbin
            dmbin  = pulse.dmbin
        # Update hits (dupes will be 0)
        pl[top_ii].nhits += 1
    return pl


def reset_bins(plist):
    for pp in plist:
        pp.reset_bin()
    return plist


def reset_dupes(plist):
    for pp in plist:
        pp.reset_dupes()
    return plist


def reset_all(plist):
    for pp in plist:
        pp.reset_all()
    return plist


def get_best_DM(plist, dt):
    """
    For pulses that are in the same time bin and beam,
    select the dm that maximizes the source sigma
    """
    # set time bin
    for pp in plist:
        pp.bin_time(dt)

    # sort by beam, time bin, sigma
    ppi = dcopy(sorted(plist, cmp=cmp_beam_tbin_sigma))
    beam = -1
    tbin = -1
    top_ii = -1
    for ii, pulse in enumerate(ppi):
        if pulse.beam != beam:
            beam = pulse.beam
            tbin = pulse.tbin
            top_ii = ii
        else:
            if pulse.tbin == tbin:
                ppi[top_ii].dm_list.append(pulse.dm)
            else:
                top_ii = ii
                tbin   = pulse.tbin
        ppi[top_ii].ndms += 1
    return ppi


def filter_plist(plist, hit_min=1, hit_max=10, ndms_min=1):
    ppi = get_best_DM(plist)
    ppi = [pp for pp in ppi if pp.nhits >= hit_min ]
    ppi = [pp for pp in ppi if pp.nhits <= hit_max ]
    ppi = [pp for pp in ppi if pp.ndms >= ndms_min ]
    return ppi


def add_time_allcands(cands, dt):
    for cc in cands:
        cc.time += dt
    return cands


def dm_delay(freq_lo, freq_hi, DM):
    return 4.15e3 * (freq_lo**-2.0 - freq_hi**-2.0) * DM


def attrarr(obj_list, attr):
    if hasattr(obj_list[0], attr):
        out_arr = np.array([ getattr(bb, attr) for bb in obj_list ])
        return out_arr
    else:
        print("List has no attribute \"%s\" " %attr)
        return


def write_cands_save(outfile, cands):
    """
    Output cand info sorted by snr
    """
    sigmas = attrarr(cands, 'sigma')
    idx = np.argsort(sigmas)
    idx = idx[::-1]

    fout = open(outfile, 'w')
    hdr = '#{:<12}{:<10}{:<10}{:<10}'.format( \
        'Time', 'Beam', 'DM', 'SNR')
    fout.write(hdr + '\n')
    for xx in idx:
        out_str = '{:<12.3f}{:<10d}{:<10.2f}{:<10.2f}'.format( \
            cands[xx].time, cands[xx].beam, cands[xx].dm, cands[xx].sigma)
        fout.write(out_str + '\n')
    fout.close()
    return


def write_cands(outfile, cands):
    """
    Output cand info sorted by snr
    """
    sigmas = attrarr(cands, 'sigma')
    idx = np.argsort(sigmas)
    idx = idx[::-1]

    fout = open(outfile, 'w')
    hdr = '#{:<12}{:<10}{:<10}{:<10}{:<10}'.format( \
        'Time', 'Beam', 'DM', 'SNR', 'Nhits')
    fout.write(hdr + '\n')
    for xx in idx:
        out_str = '{:<12.3f}{:<10d}{:<10.2f}{:<10.2f}{:<10d}'.format( \
            cands[xx].time, cands[xx].beam, cands[xx].dm, cands[xx].sigma,
            cands[xx].nhits)
        fout.write(out_str + '\n')
    fout.close()
    return



def make_midx_plot(midxs, snrs, midx_cut, snr_cut, basename):
    plt.ioff()
    xx = np.where( (midxs < midx_cut) & (snrs > snr_cut) )[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(midxs, snrs, ls='', marker='o', c='b')
    if len(xx):
        ax.plot(midxs[xx], snrs[xx], ls='', marker='o', c='r')
    ax.axvline(x=midx_cut, ls='--', lw=2, c='k')
    ax.axhline(y=snr_cut, ls='--', lw=2, c='k')
    ax.set_xlabel("Spectral Modulation Index")
    ax.set_ylabel("Single Pulse SNR")
    ax.set_title("Candidate Selection for %s" %basename)
    plt.savefig("%s_modindex.png" %basename, dpi=100, bbox_inches='tight')
    plt.close()
    plt.ion()
    return


def make_nhits_plot(ndupes, nhitmax, basename):
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    xx = np.where(ndupes)[0]
    nmax = np.max(ndupes[xx])
    nmin = np.min(ndupes[xx])
    logbinmax = np.ceil(np.log10(nmax))
    bb = 10.0**np.linspace(0, logbinmax, 50)
    ax.hist(ndupes[xx], bins=bb, log=True)
    ax.set_xscale('log')
    ax.set_ylim(0.3)
    ax.axvline(x=nhitmax, ls='--', lw=3, c='r')

    ax.set_xlabel("Nhits (Nbeams x NDMs)")
    ax.set_ylabel("Counts")

    ax.set_title("%s: %d cands, %d groups with Nhit < %d" \
                 %(basename, len(ndupes), len(xx), nhitmax),
                 fontsize=12)
    outname = "%s_hits.png" %basename
    plt.savefig(outname, dpi=100, bbox_inches='tight')
    plt.close()
    plt.ion()
    return


def get_beamlist(candfile):
    beams = []
    with open(candfile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", '\n', ' ']:
                continue
            else:
                pass
            beams.append(int(line.split()[1]))
    beams = np.array(beams)
    return beams


def copy_beams(beamlist, basename, indir, outdir):
    for bb in beamlist:
        fitsfile = "%s_beam%04d.fits" %(basename, bb)
        inpath = "%s/%s" %(indir, fitsfile)
        if os.path.isfile(inpath):
            print("FOUND: %s" %inpath)
            outpath = "%s/%s" %(outdir, fitsfile)
            shutil.copy2(inpath, outpath)
        else:
            print("MISSING: %s" %inpath)
    return


def copy_beams_from_candfile(candfile, basename, indir, outdir):
    beams = get_beamlist(candfile)
    bb, cc = np.unique(beams, return_counts=True)
    print("%d Beams, %d Unique Beams" %(len(beams), len(cc)))
    xx = np.where(cc > 1)[0]
    if len(xx):
        for idx in xx:
            print("BEAM %d -- %d hits" %(bb[idx], cc[idx]))
    copy_beams(bb, basename, indir, outdir)
    return

if __name__ == "__main__":
    
     print "Single pulse searching...."	     

