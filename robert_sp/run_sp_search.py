import numpy as np
import subprocess as sub
import time
import os
import glob
import sys
import shutil


def time_check_call(cmd_str, shell=True, debug=False):
    tstart = time.time()
    retval = 0
    try:
        if debug:
            print cmd_str
        else:
            retval = sub.check_call(cmd_str, shell=shell)
    except sub.CalledProcessError:
        retval = 1
        print("cmd failed: %s" %cmd_str)
    dt = tstart - time.time()
    return dt, retval


def param_string_from_dict(pdict):
    pstr = ""
    for kk, vv in pdict.items():
        pstr += "-%s %s " %(kk, str(vv))
    return pstr


def dedisperse_beam(fitsfile, out_dir, out_name, rfi_mask=None, pdict={}, debug=False):
    """
    This version writes the dat files to a local scratch directory, 
    then moves them to the output directory.  This is just a work-around 
    to avoid overly long names (ie, the whole path) in the inf files
    """
    err_val = 0

    param_str = param_string_from_dict(pdict)
    cmd_str = "prepsubband -o %s " %(out_name)
    cmd_str += param_str

    # if rfi_mask exists, add it to command string
    print rfi_mask
    if rfi_mask is not None and os.path.isfile(rfi_mask):
        cmd_str += "-mask %s " %rfi_mask
        print "Using RFI mask: %s" %rfi_mask
        sys.stdout.flush()
    else:
        print "No RFI mask!!!"
        sys.stdout.flush()
    
    # if fits file exists, add to command string
    if os.path.isfile(fitsfile):
        cmd_str += fitsfile
        print cmd_str
    else:
        err_val += 1
        
    # if no errors, try running function
    if err_val == 0:
        dt, retval = time_check_call(cmd_str, shell=True, debug=debug)
        err_val += retval
    else:
        dt = 0
    
    # Now gather the files to move to the output directory
    dat_files = glob.glob("%s*dat" %out_name)
    inf_files = glob.glob("%s*inf" %out_name)

    print(dat_files)
    print(inf_files)

    # Check to see if output dir exists, if not make it
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else: pass

    for ii in xrange(len(dat_files)):
        shutil.move(dat_files[ii], out_dir)
        shutil.move(inf_files[ii], out_dir)

    return dt, err_val


def run_single_pulse_search(dat_files, pdict={}, debug=False):
    err_val = 0
    
    param_str = param_string_from_dict(pdict)
    
    # check that data files have been created
    file_err = [ not (os.path.isfile(dat)) for dat in dat_files ]
    if np.any(file_err):
        err_val += 1
        print("Not all files exist!")
        sys.stdout.flush()
    else: pass

    # if files are good, build cmd string and run 
    cmd_str = "single_pulse_search.py %s %s" %(param_str, ' '.join(dat_files))
    print(cmd_str)

    if err_val == 0:
        dt, retval = time_check_call(cmd_str, shell=True, debug=debug)
        err_val += retval
    else:
        dt = 0
    
    return dt, err_val




def combine_sp_files(sp_files, outfile):
    sp_files.sort()
    fout = open(outfile, 'w')
    for ii, spfile in enumerate(sp_files):
        jj = 0
        for line in file(spfile):
            if (ii==0 and jj==0) or (jj>0):
                fout.write(line)
            else:
                pass
            jj += 1
    fout.close()
    return


def combine_sp_files_fromdir(out_dir, out_base):
    # Make sure beams_dir exists
    if not os.path.exists(out_dir):
        print("%s does not exist!  Exiting..." %out_dir)
        return
    else: pass

    sp_files = glob.glob("%s/*singlepulse" %(out_dir))
    sp_files.sort()
    
    outfile = "%s/%s.cands" %(out_dir, out_base)
    combine_sp_files(sp_files, outfile)
    return



if __name__ == "__main__":

    do_dm = 1
    do_sp = 1

    fits_base = 'test4'
    fits_file = "%s.fits" %(fits_base)
    
    top_dir = '/home/vishal/SP_pipeline/sp_scripts/'

    out_dir = '%s/search' %top_dir

    mask_dir = '%s/rfi_mask' %top_dir
    mask_file = '%s/%s_rfifind.mask' %(mask_dir, fits_base)
    
    ddm_pars_fine_full = {'numdms' : 10,
                          'dmstep' : 50,
                          'lodm'   : 1500.0, 
                          'nsub'   : 32}

    
    ddm_pars = ddm_pars_fine_full

    #sp_pars = {'f' : '', 
    #           'm' : 0.02}


    sp_pars = {'m' : 0.02, 
               'b' : ''  , 
               't' : 5 }

    # De-disperse
    tdm_start = time.time()
    if do_dm:
        dt_dm, err_dm = dedisperse_beam(fits_file, out_dir, fits_base, 
                                        rfi_mask=mask_file, pdict=ddm_pars, 
                                        debug=False)
    tdm = time.time() - tdm_start

    # Single Pulse Search
    tsp_start = time.time()
    if do_sp:
        tsp_start = time.time()
        dat_base = fits_base
        dat_files = glob.glob("%s/%s*dat" %(out_dir, dat_base))
        dat_files.sort()
        dt_sp, err_sp = run_single_pulse_search(dat_files, pdict=sp_pars, debug=False)
        combine_sp_files_fromdir(out_dir, fits_base)
        
    tsp = time.time() - tsp_start
    

    print("Dedispersion took %.1f min" %(tdm / 60.0))
    print("SP searching took %.1f min" %(tsp / 60.0))
        
