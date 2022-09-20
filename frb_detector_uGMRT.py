#!/usr/bin/env python 

#
# Modify the orignal script to suit single beam FRB search from GBT
#

import sys, math
import numpy as np
from math import sin, pi

class Classifier(object):
    def __init__(self, gdm):
        self.nbeams      = 13
        self.snr_cut     = 10.0
        self.members_cut = 3
        self.nbeams_cut  = 4
        self.dm_cut      = gdm;
        self.filter_cut  = 8
        self.beam_mask   = (1<<13) - 1  # 1111111111111 allow all beams
	# For 13 beam mask. These number are in decimal but if we convert them to binnary they generate valid beam masks. 
        self.valid_masks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 15, 16, 17, 24, 25, 29, 32, 33, 48, 49, 57, 64, 65, 66, 67, 71, 96, 97, 99, 113, 128, 130, 192, 194, 195, 256, 258, 260, 262, 263, 386, 512, 516, 520, 524, 525, 772, 1024, 1032, 1040, 1048, 1049, 1544, 2048, 2064, 2080, 2096, 2097, 3088, 4096, 4128, 4160, 4192, 4193, 4288, 6176]

    def is_hidden(self, cand):
        return ( (cand['snr'] < self.snr_cut) |
                 (cand['filter'] > self.filter_cut) | (cand['filter'] <= 2) )

    def is_noise(self, cand):
        return cand['members'] < self.members_cut

    # test if candidate is galactic or low dm in BL case
    def is_galactic(self, cand):
      #Orig
      return cand['dm'] <= self.dm_cut
      #return (cand['dm'] <= 70) | (cand['dm'] >= 72)

    def is_not_adjacent(self, cand, epoch, half_time): 	 
      return (cand['time'] > epoch + half_time) | (cand['time'] < epoch - half_time)

    # count the maximum time
    def min_time(self, cand):
        return np.amin(cand['time'])

    # count the maximum time
    def max_time(self, cand):
        return np.amax(cand['time'])


class TextOutput(object):
    def __init__(self):
        self.dm_base = 1.0
        self.snr_min = 6.0

    def print_html(self, data):
        if len(data['valid']) > 0:
            sys.stdout.write("<table width='100%' border=1 cellpadding=4px cellspacing=4px>\n")
            sys.stdout.write("<tr><th align=left>SNR</th><th align=left>Time</th><th align=left>DM</th><th align=left>Filter [ms]</th><th align=left>Beam</th></tr>\n")
            for (i, item) in enumerate(data['valid']['snr']):
                sys.stdout.write ("<tr>" + \
                                  "<td>" + str(data['valid']['snr'][i]) + "</td>" + \
                                  "<td>" + str(data['valid']['time'][i]) + "</td>" + \
                                  "<td>" + str(data['valid']['dm'][i]) + "</td>" + \
                                  "<td>" + str(0.064 * (2 **data['valid']['filter'][i])) + "</td>" + \
                                  "<td>" + str(data['valid']['prim_beam'][i]+1) + "</td>" + \
                                  "</tr>\n")
            sys.stdout.write("</table>\n")

    def print_text(self, data):
        
        cand_type = 'valid'
    
        if len(data[cand_type]) > 0:

          # get indicies list for sorting via time
          sorted_indices = [i[0] for i in sorted(enumerate(data[cand_type]['time']), key=lambda x:x[1])]

          #sys.stdout.write ( "SNR\tTIME\tSAMP\tDM\tFILTER\tBEAM\n")
          for i in sorted_indices:
                sys.stdout.write (str(data[cand_type]['snr'][i]) + "\t" + \
                                  str(data[cand_type]['time'][i]) + "\t" + \
                                  str(data[cand_type]['samp_idx'][i]) + "\t" + \
                                  str(data[cand_type]['dm'][i]) + "\t" + \
                                  str(data[cand_type]['filter'][i]) + "\t" + \
                                  str(data[cand_type]['prim_beam'][i]+1) + \
                                  "\n")

    def print_xml(self, data):
        # get indicie list for sorting via snr
        snr_sorted_indices = [i[0] for i in sorted(enumerate(data['valid']['snr']), key=lambda x:x[1],reverse=True)]

        cand_i = 0
        for i in snr_sorted_indices:
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
    
    parser = argparse.ArgumentParser(description="Detects FRB's in candidates file")
    parser.add_argument('-gdm',default=1, type=float)
    parser.add_argument('-cands_file', default="all_candidates.dat")

    parser.add_argument('-snr_cut', type=float, default=10)
    parser.add_argument('-filter_cut', type=int, default=8)
    parser.add_argument('-nbeams_cut', type=int, default=4)
    parser.add_argument('-beam_mask', type=int, default=(1<<13)-1)

    parser.add_argument('-max_cands_per_sec', type=float, default=2)
    parser.add_argument('-cand_list_xml', action="store_true")
    parser.add_argument('-cand_list_html', action="store_true")
    parser.add_argument('-verbose', action="store_true")
    parser.add_argument('-min_members_cut',type=float,default=3)
    args = parser.parse_args()
	 
    max_cands_per_second = args.max_cands_per_sec
    filename = args.cands_file
    verbose = args.verbose
    cand_list_xml = args.cand_list_xml
    cand_list_html = args.cand_list_html

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
                                      'f4', 'i4')})

    # Adjust for 0-based indexing
    all_cands['prim_beam'] -= 1
    all_cands['beam'] -= 1

    # to clear the 17th bit (RFI tag)
    clear_rfi_mask = 0b10001111111111111;

    all_cands['beam_mask'] &= clear_rfi_mask

    if verbose:
      sys.stderr.write ("Loaded %i candidates\n" % len(all_cands))
    
    classifier = Classifier(math.fabs(args.gdm))
    classifier.snr_cut = args.snr_cut
    classifier.filter_cut = args.filter_cut
    classifier.nbeams_cut = args.nbeams_cut
    classifier.beam_mask = args.beam_mask
    classifier.members_cut = args.min_members_cut    
    # Filter candidates based on classifications
    if verbose:
      sys.stderr.write ("Classifying candidates...\n")

    categories = {}

    is_hidden      = classifier.is_hidden(all_cands)
    is_noise       = classifier.is_noise(all_cands)
    #is_coinc_dumb  = classifier.is_coinc_rfi_dumb(all_cands)
    #is_coinc_smart = classifier.is_coinc_rfi_smart(all_cands)
    is_galactic    = classifier.is_galactic(all_cands)
    is_valid       = (is_hidden == False) & (is_noise == False) & (is_galactic == False)

    categories["hidden"]      = all_cands[is_hidden]
    categories["noise"]       = all_cands[(is_hidden == False) & is_noise]
    #categories["coinc_dumb"]  = all_cands[(is_hidden == False) & (is_noise == False) & is_coinc_dumb]
    #categories["coinc_smart"] = all_cands[(is_hidden == False) & (is_noise == False) & (is_coinc_dumb == False) & is_coinc_smart]
    categories["galactic"]    = all_cands[(is_hidden == False) & (is_noise == False) & is_galactic]
    categories["valid"]       = all_cands[is_valid]

    pre_valid = len(categories["valid"])

    # for valid events, check the event rate around the time of the event
    if len(categories['valid']) > 0 & False:
      min_time = classifier.min_time(all_cands)
      max_time = classifier.max_time(all_cands)

      # look in a 8 second window around the event for excessive RFI
      event_time = 8

      # Here we check if the event rate excedes set maximum number of events per sec
	
      for (i, item) in reversed(list(enumerate(categories['valid']))):
        #cands_per_second = classifier.events_per_sec (all_cands, item['time'], min_time, max_time)
        #print "cands_per_second(orig)="+str(cands_per_second)

        epoch = item['time']
        half_time = event_time / 2.0
        new_time = 0

        if (epoch - half_time) < min_time:
          new_time += (epoch - min_time)
        else:
          new_time += half_time

        if (epoch + half_time) > max_time:
          new_time += (max_time - epoch)
        else:
          new_time += half_time

        half_time = new_time / 2.0

        is_not_adjacent = (is_noise == False) & classifier.is_not_adjacent(all_cands, epoch, half_time) 
        is_valid = (is_noise == False) & (is_not_adjacent == False)

        event_sum = float(np.count_nonzero(is_valid))
        cands_per_second = event_sum / new_time

        if cands_per_second < 0:
          sys.stderr.write ( "cands_per_second = %f \n" % ( cands_per_second )) 
          cands_per_second = max_cands_per_second
        if verbose:
          sys.stderr.write ( "cands_per_second around %f was %f [max = %f]\n" % ( item['time'], cands_per_second, max_cands_per_second )) 
        #I commented this for testing 
        #if cands_per_second >= max_cands_per_second:
        #  categories['valid'] = np.delete(categories['valid'], i, axis=0)
	#if(categories['valid']['time'] < item['time']+0.0015/2.0 and categories['valid']['time'] > item['time'] - 0.0015/2.0 ):
	#	print categories['valid']['snr']
	#I think rather than rejecting all candidates in a block where cand_per_sec exceed limit. 
	period=1.0/max_cands_per_second
	in_same_grp =(categories['valid']['time'] < item['time'] + period) & (categories['valid']['time'] > item['time'] - period)
	if item['snr'] < np.max(categories['valid'][in_same_grp]['snr']):
		categories['valid'] = np.delete(categories['valid'], i, axis=0)


    rfi_storm = pre_valid - len(categories["valid"])

    if verbose:
      sys.stderr.write ( "Classified %i as hidden \n" % len(categories["hidden"]))
      sys.stderr.write ( "           %i as noise spikes\n" % len(categories["noise"]))
      #sys.stderr.write ( "           %i as coinc RFI [nbeam > %i]\n" % (len(categories["coinc_dumb"]), classifier.nbeams_cut))
      #sys.stderr.write ( "           %i as coinc RFI [mask_rule]\n" % len(categories["coinc_smart"]))
      sys.stderr.write ( "           %i as low DM spikes\n" % len(categories["galactic"]))
      sys.stderr.write ( "           %i as RFI storm\n" % rfi_storm)
      sys.stderr.write ( "           %i as valid FRB candidates\n" % len(categories["valid"]))

    text_output = TextOutput()

    if cand_list_xml:
      text_output.print_xml(categories)
    elif cand_list_html:
      text_output.print_html(categories)
    else:
      text_output.print_text(categories)

    if verbose:
      sys.stderr.write ( "Done\n")


