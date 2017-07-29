import numpy as np
import matplotlib.pyplot as plt
import sys

kdm = 4148.808 # MHz^2 / (pc cm^-3)

def dm_delay(f1, f2, DM, kdm=kdm):
    return kdm * DM * (1.0 / (f1 * f1) - 1 / (f2 * f2))


def deltaT(ichan, dt, df, f0, kdm=kdm):
    return (kdm / dt) * ((f0 + ichan * df)**-2.0 - f0**-2)


def dmdt(DM, ichan, dt, df, f0, kdm=kdm):
    return int(np.round( DM * deltaT(ichan, dt, df, f0, kdm=kdm)))


def fake_dspec(tt, t0, W, f0, df, nchan, DM):
    pp = np.exp(-(tt - t0)**2.0 / (2 * W**2.0))
    dt = tt[1] - tt[0]
    nt = len(tt)
    dsamp = np.array([ dmdt(DM, ichan, dt, df, f0) for ichan in xrange(nchan) ])
    offs = dsamp - np.min(dsamp)
    dspec = np.zeros(shape=(nchan, len(tt)))
    for jj in xrange(nchan):
        dspec[jj, offs[jj]:] = pp[:nt-offs[jj]]
    return dspec
    

def dedisperse_one(dspec, dm, dt, df, f0, kdm=kdm):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = np.array([ dmdt(dm, ichan, dt, df, f0, kdm=kdm) for ichan in xrange(nchans) ])
    dsamps -= np.min(dsamps)
    tpad = np.max(dsamps)
    outarr = np.zeros( nt + tpad )
    for ii in xrange(nchans):
        osl = slice(tpad - dsamps[ii], nt + tpad - dsamps[ii])
        outarr[osl] += dspec[ii]
    return outarr[tpad:nt + tpad] / float(nchans)


def dedisperse_naive(dspec, dms, dt, df, f0, kdm=kdm):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    outarr = np.zeros( shape=(len(dms), nt) )
    for ii in xrange(len(dms)):
        outarr[ii] = dedisperse_one(dspec, dms[ii], dt, df, f0, kdm=kdm)
    return outarr


def dedisperse_dspec(dspec, dm, dt, df, f0, kdm=kdm):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = np.array([ dmdt(dm, ichan, dt, df, f0, kdm=kdm) for ichan in xrange(nchans) ])
    dsamps -= np.min(dsamps)
    tpad = np.max(dsamps)
    outarr = np.zeros( (nchans, nt + tpad) )
    for ii in xrange(nchans):
        osl = slice(tpad - dsamps[ii], nt + tpad - dsamps[ii])
        outarr[ii, osl] = dspec[ii]
    return outarr[:, tpad:nt + tpad]
