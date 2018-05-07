#!/usr/bin/env python 

import sys
import numpy as np

MAX_DM = 4000

class Classifier(object):
    def __init__(self):
        self.nbeams      = 13
        self.snr_cut     = 6.5
        self.members_cut = 3
        self.nbeams_cut  = 3
        self.dm_cut      = 1.5
        self.filter_cut  = 10
        self.beam_mask   = (1<<13) - 1
        self.filter_max  = 12
        
    def is_masked(self, beam):
        return ((1<<beam) & self.beam_mask) == 0
    
    def is_hidden(self, cand):
        # self.is_noise(cand) | 
        return ( (cand['snr'] < self.snr_cut) |
                 (cand['filter'] > self.filter_cut) |
                 self.is_masked(cand['beam']) |
                 ((self.is_masked(cand['beam']) != True) &
                  (cand['beam'] != cand['prim_beam'])) )
    
    def is_noise(self, cand):
        return cand['members'] < self.members_cut

    def is_fat(self, cand):
        return cand['filter'] >= self.filter_max
    
    def count_nbeams(self, mask):
        n = 0
        for i in range(self.nbeams):
            n += (mask & (1<<i)) > 0
        return n
            
    def is_coinc_rfi(self, cand):
        nbeams = self.count_nbeams(cand['beam_mask'] & self.beam_mask)
        return nbeams > self.nbeams_cut
    
    def is_lowdm_rfi(self, cand):
        return cand['dm'] < self.dm_cut

class TimeDMPlot(object):
    def __init__(self, g, multiplot):
        self.g = g
        self.dm_base = 1.0
        self.snr_min = 6.0
        self.multiplot = multiplot
        
    def plot(self, data):
        self.g.reset()
        if self.multiplot:
          self.g('set size 1.0,0.5')
        else:
          self.g('set size 1.0,1.0')
        self.g('set origin 0.0,0.0')
        self.g('unset key')
        self.g('set autoscale x')
        self.g('set logscale y')
        self.g('set logscale y2')
        self.g('set yrange[1.0:'+str(MAX_DM * 1.1)+']')
        self.g('set y2range[1.0:'+str(MAX_DM * 1.1)+']')
        self.g('set cbrange[-0.5:12.5]')
        self.g('set palette positive nops_allcF maxcolors 13 gamma 1.5 color model RGB')
        self.g("set palette defined ( 0 'green', 1 'cyan', 2 'magenta', 3 'orange' )")
        self.g('unset colorbox')
        self.g('set grid noxtics nomxtics ytics mytics lt 9 lw 0.2')
        self.g('set ytics 10')
        self.g('set mytics 10')
        self.g('set y2tics 10 out mirror format ""')
        self.g('set my2tics 10')
        self.g('set xtics auto')
        self.g('set x2tics auto out mirror format ""')
        self.g('set mxtics 4')
        self.g('set mx2tics 4')
        self.g('set xlabel "Time [s]"')
        self.g('set ylabel "DM + 1 [pc cm^{-3}]"')
        self.g('min(x,y) = x<=y?x:y')
        self.g('max(x,y) = x>=y?x:y')

        to_plot = []

        if (len(data['noise']['snr']) > 0):        
            to_plot.append(Gnuplot.Data(data['noise']['snr'],
                                 data['noise']['time'],
                                 data['noise']['dm'],
                                 using="2:($3+%f):(($1-%f)/2.0+0.5)" \
                                     % (self.dm_base,self.snr_min),
                                 with_="p pt 2 lt 9 lw 0.5 ps variable", inline=True))

        if (len(data['fat']['snr']) > 0):        
            to_plot.append(Gnuplot.Data(data['fat']['snr'],
                                 data['fat']['time'],
                                 data['fat']['filter'],
                                 data['fat']['dm'],
                                 using="2:($4+%f):(min(($1-%f)/2.0+0.9,5)):3" \
                                     % (self.dm_base,self.snr_min),
                                 with_="p pt 6 lw 0.5 lt palette ps variable", inline=True))
            to_plot.append(Gnuplot.Data(data['fat']['beam'],
                                       data['fat']['time'],
                                       data['fat']['dm'],
                                       using='2:($3+%f):(sprintf("%%d",$1+1))' \
                                           % (self.dm_base),
                                       with_='labels center font ",7" offset 0,0.05 textcolor rgbcolor "gray"', inline=True))


        if (len(data['coinc']['snr']) > 0):
            to_plot.append(Gnuplot.Data(data['coinc']['snr'],
                                 data['coinc']['time'],
                                 data['coinc']['filter'],
                                 data['coinc']['dm'],
                                 using="2:($4+%f):(min(($1-%f)/2.0+0.9,5)):3" \
                                     % (self.dm_base,self.snr_min),
                                 with_="p pt 3 lw 0.25 lt palette ps variable", inline=True))
            
        if (len(data['lowdm']['snr']) > 0):        
            to_plot.append(Gnuplot.Data(data['lowdm']['snr'],
                                 data['lowdm']['time'],
                                 data['lowdm']['filter'],
                                 data['lowdm']['dm'],
                                 using="2:($4+%f):(min(($1-%f)/2.0+0.9,5)):3" \
                                     % (self.dm_base,self.snr_min),
                                 with_="p pt 6 lt palette ps variable", inline=True))
            to_plot.append(Gnuplot.Data(data['lowdm']['beam'],
                                       data['lowdm']['time'],
                                       data['lowdm']['dm'],
                                       using='2:($3+%f):(sprintf("%%d",$1+1))' \
                                           % (self.dm_base),
                                       with_='labels center font ",7" offset 0,0.05 textcolor rgbcolor "black"', inline=True))

        if (len(data['valid']['snr']) > 0):
            to_plot.append(Gnuplot.Data(data['valid']['snr'],
                                 data['valid']['time'],
                                 data['valid']['filter'],
                                 data['valid']['dm'],
                                 using="2:($4+%f):(min(($1-%f)/2.0+0.9,5)):3" \
                                     % (self.dm_base,self.snr_min),
                                 with_="p pt 7 lt palette ps variable", inline=True))
            to_plot.append(Gnuplot.Data(data['valid']['beam'],
                                       data['valid']['time'],
                                       data['valid']['dm'],
                                       using='2:($3+%f):(sprintf("%%d",$1+1))' \
                                           % (self.dm_base),
                                       with_='labels center font ",7" offset 0,0.05 textcolor rgbcolor "black"', inline=True))
        self.g.plot(*to_plot)


class DMSNRPlot(object):
    def __init__(self, g):
        self.g = g
        self.dm_base = 1.0
        self.snr_base = 5.9
        self.max_filter = 12
	self.dt = 163.84e-6
    def plot(self, data):
        self.g.reset()
        self.g('set size 0.5,0.5')
        self.g('set origin 0.5,0.5')
        self.g('unset key')
        self.g('set xrange[1.0:'+str(MAX_DM)+']')
        self.g('set x2range[1.0:'+str(MAX_DM)+']')
        self.g('set logscale x')
        self.g('set logscale x2')
        self.g('set xtics 10')
        self.g('set mxtics 10')
        self.g('set x2tics 10 out mirror format ""')
        self.g('set mx2tics 10')
        self.g('set xtics mirror')
        self.g('set grid noytics nomytics xtics mxtics lt 9 lw 0.2')
        self.g('set logscale y')
        self.g('set logscale y2')
        self.g('set yrange[0.1:100]')
        self.g('set y2range[0.1:100]')
        self.g('set cbrange[-0.5:12.5]')
        self.g('set palette positive nops_allcF maxcolors 13 gamma 1.5 color model RGB')
        self.g("set palette defined ( 0 'green', 1 'cyan', 2 'magenta', 3 'orange' )")
        self.g('set colorbox')
        self.g('snr_min = %f' % self.snr_base)
        self.g('unset mytics')
        self.g('unset my2tics')
        self.g('set ytics ("6 " 6-snr_min, "6.1 " 6.1-snr_min, "6.4 " 6.4-snr_min, "7 " 7-snr_min, "8 " 8-snr_min, "10 " 10-snr_min, "13 " 13-snr_min, "20 " 20-snr_min, "40 " 40-snr_min, "100 " 100-snr_min)')
        self.g('set y2tics ("6" 6-snr_min, "6.1 " 6.1-snr_min, "6.4 " 6.4-snr_min, "7 " 7-snr_min, "8 " 8-snr_min, "10 " 10-snr_min, "13 " 13-snr_min, "20 " 20-snr_min, "40 " 40-snr_min, "100 " 100-snr_min) out mirror format ""')
        #self.g('set cbtics 1 format "2^%g"')
        self.g('set cbtics 1 format ""')
        filter_tics = [1000*self.dt * 2**i for i in range(self.max_filter+1)]
        #self.g('set cbtics add ("64 us" 0, "128 us" 1, "256 us" 2, "512 us" 3, "1 ms" 4, "2 ms" 5, "4 ms" 6, "8 ms" 7, "16 ms" 8, "32 ms" 9, "64 ms" 10, "128 ms" 11, "256 ms" 12)')
        self.g('set cbtics add ('+', '.join(['"%.4g" %i'%(x,i) for i,x in enumerate(filter_tics)])+')')
        self.g('set xlabel "DM+1 [pc cm^{-3}]"')
        self.g('set ylabel "SNR"')
        #self.g('set cblabel "log_{2} boxcar width"')
        self.g('set cblabel "Boxcar width [ms]"')
  
        to_plot = []
            
        if (len(data['noise']['snr']) > 0):
            to_plot.append(Gnuplot.Data(data['noise']['snr'],
                                     data['noise']['dm'],
                                     using="($2+%f):($1-%f)" \
                                        % (self.dm_base,self.snr_base),
                                     with_="p pt 2 ps 0.5 lt 9", inline=True))
        if (len(data['fat']['snr']) > 0):
            to_plot.append(Gnuplot.Data(data['fat']['snr'],
                                     data['fat']['dm'],
                                     data['fat']['filter'],
                                     using="($2+%f):($1-%f):3" \
                                         % (self.dm_base,self.snr_base),
                                     with_="p pt 6 ps 0.8 lw 0.3 lt palette", inline=True))
        if (len(data['lowdm']['snr']) > 0):
            to_plot.append(Gnuplot.Data(data['lowdm']['snr'],
                                     data['lowdm']['dm'],
                                     data['lowdm']['filter'],
                                     using="($2+%f):($1-%f):3" \
                                         % (self.dm_base,self.snr_base),
                                     with_="p pt 6 ps 0.8 lt palette", inline=True))
        if (len(data['coinc']['snr']) > 0):
            to_plot.append(Gnuplot.Data(data['coinc']['snr'],
                                     data['coinc']['dm'],
                                     data['coinc']['filter'],
                                     using="($2+%f):($1-%f):3" \
                                         % (self.dm_base,self.snr_base),
                                     with_="p pt 3 ps 0.8 lw 0.3 lt palette", inline=True))
        if (len(data['valid']['snr']) > 0):
            to_plot.append(Gnuplot.Data(data['valid']['snr'],
                                 data['valid']['dm'],
                                 data['valid']['filter'],
                                 using="($2+%f):($1-%f):3" \
                                     % (self.dm_base,self.snr_base),
                                 with_="p pt 7 ps 0.8 lt palette", inline=True))
        self.g.plot(*to_plot)

class DMHistogram(object):
    def __init__(self, cands=None):
        self.dm_min   = 0.10
        self.dm_max   = MAX_DM * 1.01
        self.min_bins = 30
        self.hist     = None
        if cands is not None:
            self.build(cands)
            
    def build(self, cands):
        # cands = cands[cands["filter"] < 11]
        import math
        N = len(cands)
        log_dm_min = math.log10(self.dm_min)
        log_dm_max = math.log10(self.dm_max)
        nbins    = max(self.min_bins, 2*int(math.sqrt(N)))
        binwidth = (log_dm_max - log_dm_min) / nbins
        bins_    = 10**(log_dm_min + (np.arange(nbins)+0.5)*binwidth)
        dms      = np.maximum(cands['dm'], self.dm_min)
        log_dms  = np.log10(dms)
        vals, edges = np.histogram(log_dms, bins=nbins,
                                   range=(log_dm_min,log_dm_max))
        self.hist = np.rec.fromrecords(np.column_stack((bins_, vals)),
                                       names=('bins', 'vals'))

class DMHistPlot(object):
    def __init__(self, g):
        self.g = g
        self.dm_base = 1.0
        self.snr_base = 5.9
        self.max_filter = 12
        self.dt = 163.84e-6
    def plot(self, data):
        self.g.reset()
        self.g('set size 0.5,0.5')
        self.g('set origin 0.0,0.5')
        self.g('unset key')
        self.g('set logscale x')
        self.g('set xrange[1:'+str(MAX_DM)+']')
        self.g('set logscale x2')
        self.g('set x2range[1:'+str(MAX_DM)+']')
        self.g('set logscale y')
        self.g('set logscale y2')
        self.g('set yrange [1:2000]')
        self.g('set y2range [1:2000]')
        self.g('set xtics 10')
        self.g('set mxtics 10')
        self.g('set x2tics 10 out mirror format ""')
        self.g('set mx2tics 10')
        self.g('set ytics 10')
        self.g('set y2tics 10 out mirror format ""')
        self.g('set mytics 10')
        self.g('set my2tics 10')
        self.g('set grid noxtics nomytics xtics mxtics lt 9 lw 0.2')
        #self.g('set key inside top center horizontal samplen 2 maxcols 6')
        self.g('set key inside top center horizontal samplen 2')
        self.g('set xlabel "DM+1 [pc cm^{-3}]')
        self.g('set ylabel "Candidate count"')
        
        beams = []
        for b,beam_hist in enumerate(data):
            beams.append( Gnuplot.Data(beam_hist['bins'],
                                       beam_hist['vals'],
                                       using="($1+%f):2" \
                                           % (self.dm_base),
                                       with_='histeps lw %i lt 1 lc %i' \
                                           % (1+(b+1<8),b+1),
                                       title=str(b+1),inline=True))
        self.g.plot(*beams)

class TextOutput(object):
    def __init__(self):
        self.dm_base = 1.0
        self.snr_min = 6.0
    def print_text(self, data, max_cands):
        sys.stdout.write("SNR\tTime\tDM\n")
        for (i, item) in enumerate(data['valid']['snr']):
            sys.stdout.write (str(data['valid']['snr'][i]) + "\t" + str(data['valid']['time'][i]) + "\t" + str(data['valid']['dm'][i]) + "\n")
    def print_xml(self, data, max_cands):
        # get indicie list for sorting via snr
        snr_sorted_indices = [i[0] for i in sorted(enumerate(data['valid']['snr']), key=lambda x:x[1],reverse=True)]

        cand_i = 0
        for i in snr_sorted_indices:
            if cand_i >= max_cands:
                return
            else:
                cand_i += 1
            sys.stdout.write ("<candidate snr='" + str(data['valid']['snr'][i]) + \
                                       "' time='" + str(data['valid']['time'][i]) + \
                                       "' dm='" + str(data['valid']['dm'][i]) + \
                                       "' samp_idx='" + str(data['valid']['samp_idx'][i]) + \
                                       "' filter='" + str(data['valid']['filter'][i]) + \
                                       "' prim_beam='" + str(data['valid']['prim_beam'][i] + 1) + "'/>\n")




if __name__ == "__main__":
    import argparse
    import Gnuplot
    
    parser = argparse.ArgumentParser(description="Generates data for Heimdall overview plots.")
    parser.add_argument('-cands_file', default="all_candidates.dat")
    parser.add_argument('-nbeams', type=int, default=13)
    parser.add_argument('-snr_cut', type=float)
    parser.add_argument('-beam_mask', type=int, default=(1<<13)-1)
    parser.add_argument('-nbeams_cut', type=int, default=3)
    parser.add_argument('-members_cut', type=int, default=3)
    parser.add_argument('-dm_cut', type=float, default=1.5)
    parser.add_argument('-filter_cut', type=int, default=99)
    parser.add_argument('-filter_max', type=int, default=12)
    parser.add_argument('-min_bins', type=int, default=30)
    parser.add_argument('-resolution', default="1024x768")
    parser.add_argument('-std_out', action="store_true")
    parser.add_argument('-skip_rows', type=int, default=0)
    parser.add_argument('-just_time_dm', action="store_true")
    parser.add_argument('-cand_list_xml', action="store_true")
    parser.add_argument('-max_cands', type=int, default=20)
    parser.add_argument('-no_plot', action="store_true")
    parser.add_argument('-interactive', action="store_true")
    parser.add_argument('-verbose', action="store_true")
    args = parser.parse_args()
    
    filename = args.cands_file
    nbeams = args.nbeams
    interactive = args.interactive
    std_out = args.std_out
    skip_rows = args.skip_rows
    just_time_dm = args.just_time_dm
    verbose = args.verbose
    cand_list_xml = args.cand_list_xml
    max_cands = args.max_cands
    no_plot = args.no_plot
    resolution = args.resolution
    res_parts = resolution.split("x")
    if (len(res_parts) != 2):
      sys.stderr.write("ERROR: resolution must be of form 1024x768")
      sys.exit(1)

    res_x = res_parts[0]
    res_y = res_parts[1]
    
    # Load candidates from all_candidates file
    all_cands = \
        np.loadtxt(filename,
                   dtype={'names': ('snr','samp_idx','time','filter',
                                    'dm_trial','dm','members','begin','end',
                                    'nbeams','beam_mask','prim_beam',
                                    'max_snr','beam'),
                          'formats': ('f4', 'i4', 'f4', 'i4',
                                      'i4', 'f4', 'i4', 'i4', 'i4',
                                      'i4', 'i4', 'i4',
                                      'f4', 'i4')},
                   skiprows=skip_rows)

    # Adjust for 0-based indexing
    all_cands['prim_beam'] -= 1
    all_cands['beam'] -= 1

    if verbose:
      sys.stderr.write ("Loaded %i candidates\n" % len(all_cands))
    
    classifier = Classifier()
    classifier.nbeams = args.nbeams
    classifier.snr_cut = args.snr_cut
    classifier.beam_mask = args.beam_mask
    classifier.nbeams_cut = args.nbeams_cut
    classifier.members_cut = args.members_cut
    classifier.dm_cut = args.dm_cut
    classifier.filter_cut = args.filter_cut
    classifier.filter_max = args.filter_max
    
    # Filter candidates based on classifications
    if verbose:
      sys.stderr.write ("Classifying candidates...\n")
    categories = {}
    is_hidden = classifier.is_hidden(all_cands)
    is_noise  = (is_hidden==False) & classifier.is_noise(all_cands)
    is_coinc  = (is_hidden==False) & (is_noise ==False) & classifier.is_coinc_rfi(all_cands)
    is_fat    = (is_hidden==False) & (is_noise ==False) & (is_coinc ==False) & classifier.is_fat(all_cands)
    is_lowdm  = (is_hidden==False) & (is_noise ==False) & (is_fat   ==False) & (is_coinc ==False) & classifier.is_lowdm_rfi(all_cands)
    is_valid  = (is_hidden==False) & (is_noise ==False) & (is_fat   ==False) & (is_coinc ==False) & (is_lowdm ==False)
    categories["hidden"] = all_cands[is_hidden]
    categories["noise"]  = all_cands[is_noise]
    categories["coinc"]  = all_cands[is_coinc]
    categories["fat"]    = all_cands[is_fat]
    categories["lowdm"]  = all_cands[is_lowdm]
    categories["valid"]  = all_cands[is_valid]
    
    if verbose:
      sys.stderr.write ( "Classified %i as hidden\n" % len(categories["hidden"]))
      sys.stderr.write ( "           %i as noise spikes\n" % len(categories["noise"]))
      sys.stderr.write ( "           %i as coincident RFI\n" % len(categories["coinc"]))
      sys.stderr.write ( "           %i as fat RFI\n" % len(categories["fat"]))
      sys.stderr.write ( "           %i as low-DM RFI\n" % len(categories["lowdm"]))
      sys.stderr.write ( "           %i as valid candidates\n" % len(categories["valid"]))
    
    if verbose:
      sys.stderr.write ( "Building histograms...\n")
    beam_hists = []
    for beam in range(nbeams):
        cands = all_cands[all_cands['beam'] == beam]
        beam_hists.append(DMHistogram(cands).hist)

    if cand_list_xml:
      if verbose:
        sys.stderr.write ( "Generating text only listing on stdout:\n")
      text_output = TextOutput()
      text_output.print_xml(categories, max_cands)

    if not no_plot:    
      # Generate plots
      if verbose:
        sys.stderr.write ( "Generating plots...\n")
      g = Gnuplot.Gnuplot(debug=0)
      if not interactive:
        g('set terminal pngcairo enhanced font "arial,10" size ' + res_x + ', ' + res_y)
        if std_out:
          g('set output')
          if verbose:
            sys.stderr.write ( "Writing binary image data to STDOUT\n")
        else:
          g('set output "overview_' + resolution + '.tmp.png"')
          if verbose:
            sys.stderr.write ( "Writing plots to overview_" + resolution + ".tmp.png\n")
      else:
        g('set terminal x11 size '+ res_x + ', ' + res_y)

      if just_time_dm:
        timedm_plot = TimeDMPlot(g, False)
        timedm_plot.plot(categories)
      else:
        # g('set terminal x11 size '+ res_x + ', ' + res_y)
        g('set multiplot')
        if verbose:
          sys.stderr.write ( "Gen TimeDM\n")
        timedm_plot = TimeDMPlot(g, True)
        if verbose:
          sys.stderr.write ( "Gen DMSNR\n")
        dmsnr_plot  = DMSNRPlot(g)
        if verbose:
          sys.stderr.write ( "Gen DMHist\n")
        dmhist_plot = DMHistPlot(g)
        if verbose:
          sys.stderr.write ( "Plot TimeDM\n")
        timedm_plot.plot(categories)
        if verbose:
          sys.stderr.write ( "Plot DMSNR\n")
        dmsnr_plot.plot(categories)
        if verbose:
          sys.stderr.write ( "Plot BeamHists\n")
        dmhist_plot.plot(beam_hists)
        g('unset multiplot')
      
      if interactive:
          raw_input('Please press return to close...\n')
        
    if verbose:
      sys.stderr.write ( "Done\n")
