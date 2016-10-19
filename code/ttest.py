#!/usr/bin/env python

import sys
import itertools
import operator
import numpy as np
import matplotlib.pyplot as plt

from obspy.core import UTCDateTime, read
from obspy.signal.invsim import paz_to_freq_resp
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from matplotlib.mlab import csd
import matplotlib as mpl

mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font', size=12)


def computeresp(resp, delta, lenfft):
    respval = paz_to_freq_resp(resp['poles'],
                               resp['zeros'],
                               resp['sensitivity']*resp['gain'],
                               t_samp=delta,
                               nfft=lenfft, freq=False)
    respval = np.absolute(respval*np.conjugate(respval))
    respval = respval[1:]
    return respval


def cp(tr1, tr2, lenfft, lenol, delta):
    sr = 1./delta
    cpval, fre = csd(tr1.data, tr2.data,
                     NFFT=lenfft, Fs=sr, noverlap=lenol,
                     scale_by_freq=True)
    fre = fre[1:]
    cpval = cpval[1:]
    return cpval, fre


def getFile(sta, loc, chan, stime):
    #fileStr = '/tr1/telemetry_days/' + sta + '/'
    #fileStr += str(stime.year) + '/' + str(stime.year) + '_'
    # fileStr = str(stime.julday).zfill(3) # + '/'
    fileStr = sta + '_'
    fileStr += loc + '_' + chan + '.512.seed'
    return fileStr


def compNoise(idx, length, overlap, delta, instresp):
    cpFix = lambda i1, i2: cp(st[i1], st[i2], length,
                              overlap, delta)

    # We could do these as permutations but that gets confusing
    # Instead we will just hard code the indices
    pp, f = cpFix(idx, idx)
    if idx == 0:
        noisetemp = pp - \
                    cpFix(1, 0)[0]*cpFix(0, 2)[0]/cpFix(1, 2)[0]
    elif idx == 1:
        noisetemp = pp - \
                    cpFix(2, 1)[0]*cpFix(1, 0)[0]/cpFix(2, 0)[0]
    elif idx == 2:
        noisetemp = pp - \
                    cpFix(1, 2)[0]*cpFix(2, 0)[0]/cpFix(1, 0)[0]
    else:
        print 'Bad index, crash landing.'
        sys.exit()
    # Convert to acceleration
    noisetemp *= (2.*np.pi*f)**2
    pp *= (2.*np.pi*f)**2
    # Remove the response
    noisetemp *= 1./instresp
    pp *= 1./instresp
    # Convert to dB
    noisetemp = 10.*np.log10(np.absolute(noisetemp))
    pp = 10.*np.log10(np.absolute(pp))
    return noisetemp, pp, f


def mkFigLabel(st):
    figlabel = str(st[0].stats.starttime.year)
    figlabel += str(st[0].stats.starttime.julday)
    figlabel += str(st[0].stats.starttime.hour)
    figlabel += str(st[0].stats.starttime.minute)
    for idx in range(0, 3):
        figlabel += st[idx].stats.station + st[idx].stats.location
        figlabel += st[idx].stats.channel
    return figlabel


def pltTseries(st):
    # Make a time series plot
    # Grab a legend
    titlelegend = 'Time Series: ' + str(st[0].stats.starttime.year)
    titlelegend += ' ' + str(st[0].stats.starttime.julday) + ' '
    titlelegend += str(st[0].stats.starttime.hour) + ':'
    titlelegend += str(st[0].stats.starttime.minute) + ':'
    titlelegend += str(st[0].stats.starttime.second)
    titlelegend += ' ' + str(st[0].stats.npts*delta) + ' seconds'

    tval = np.arange(0, st[0].stats.npts) / st[0].stats.sampling_rate
    plt.figure(1)

    colors = ['r', 'b', 'g']
    for idx in range(0, 3):

        # Grab the data and detrend
        stplt = (st[idx].detrend()).data
        plt.subplot(3, 1, idx+1)
        plt.plot(tval, stplt, colors[idx], label='TSeries ' +
                 st[idx].stats.station + ' ' + st[idx].stats.location +
                 ' ' + st[idx].stats.channel)
        plt.legend(prop={'size': 12}, loc=2)
        plt.xlim((0, np.amax(tval)))
        if idx == 2:
            plt.ylabel('Counts')
    plt.xlabel('Time (s)')
    plt.suptitle(titlelegend, fontsize=12)
    figlabel = mkFigLabel(st)
    plt.savefig('TSERIES' + figlabel + '.jpg', format='jpeg', dpi=400)
    plt.clf()

    # Do not bother returning the figure handle just save it
    return


def pltNoise(noises, st):
    # Make a title
    titlelegend = st[0].stats.channel + ' Self-Noise Start Time: '
    titlelegend += str(st[0].stats.starttime.year) + ' '
    titlelegend += str(st[0].stats.starttime.julday) + ' '
    titlelegend += str(st[0].stats.starttime.hour) + ':'
    titlelegend += str(st[0].stats.starttime.minute) + ':'
    titlelegend += str(st[0].stats.starttime.second) + ' Duration: '
    titlelegend += str(int(st[0].stats.npts*delta/(60*60))) + ' Hours'

    colors = ['r', 'b', 'g']
    plt.figure(1)
    plt.subplot(1, 1, 1)
    plt.title(titlelegend, fontsize=12)
    for idx, noise in enumerate(noises):
        lbstr = st[idx].stats.station + ' ' + st[idx].stats.location
        plt.plot(1./noise[2], noise[1], colors[idx],
                 label='PSD ' + lbstr, linewidth=1.5)
        plt.plot(1./noise[2], noise[0], colors[idx] + ':',
                 label='Noise ' + lbstr, linewidth=1.5)

    # Grab the noise model
    NLNMper, NLNMpower = get_nlnm()
    NHNMper, NHNMpower = get_nhnm()
    plt.plot(NLNMper, NLNMpower, 'k', label='NHNM/NLNM', linewidth=2.0)
    plt.plot(NHNMper, NHNMpower, 'k', linewidth=2.0)
    plt.legend(frameon=True)
    plt.xlabel('Period (s)')
    plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz)$')
    plt.xscale('log')
    # We might want to limit this to the sample rate
    plt.xlim((1/20., 200.))
    figlabel = mkFigLabel(st)
    plt.savefig('NOISE' + figlabel + '.jpg', format='jpeg', dpi=400)

    # Don't bother keeping the figure handle
    plt.clf()
    return


def writeResults(noises, st, chan):
    pers = [(.1, 30.), (30., 100.)]
    resultslabel = mkFigLabel(st)
    f = open(resultslabel + 'RESULTS', 'w+')
    for per in pers:
        for idx, noise in enumerate(noises):
            noidx = (per[0] <= 1./noise[2]) & (per[1] <= 1./noise[2])
            ps = noise[1][noidx]
            ns = noise[0][noidx]
            valstr = str(per[0]) + ' to ' + str(per[1]) + \
                ' s sensor ' + str(idx+1) + ': ' + chan + \
                ' ' + str(np.average(ns.real)) + '+/-' + \
                str(np.std(ns.real)) + '\n'
            f.write('Self-Noise ' + valstr)
            valstr = str(per[0]) + ' to ' + str(per[1]) + \
                ' s sensor ' + str(idx+1) + ': ' + chan + \
                ' ' + str(np.average(ps.real)) + '+/-' + \
                str(np.std(ps.real)) + '\n'
            f.write('PSD Noise ' + valstr)
    f.close()
    return

def grabPZ(sensor):
    # Do we want to maintain this or pull this in.
    if sensor == 'T-Compact':
        paz = {'zeros': [0.j, 0.j, -392. + 0.j, -1960. + 0.j,
                     -1490. + 1740.j, -1490. - 1740.j],
           'poles': [-0.3691 + 0.03702j, -0.3691 - 0.03702j,
                     -343. + 0.j, -370. + 467.j, -370. - 467.j,
                     -836. + 1522.j, -836. - 1522.j, -4900. + 4700.j,
                     -4900. - 4700.j, -6900. + 0.j, -15000. + 0.j],
           'gain': 4.344928*10**17, 'sensitivity': 754.3*2.**26/40.}
    else:
        'sensor type is not implemented'
        sys.exit()
    return paz


class SNCLSen:
    # Here is the SNCL class
    def __init__(self, net, sta, loc, chans, sensor):
        self.network = net
        self.station = sta
        self.location = loc
        self.channels = chans
        self.sensor = sensor


if __name__ == '__main__':

    # Here is a list of SNCLSen objcts for now we can use getFile
    # a different approach would be to grab the data and have these
    # fields get populated

    # No one seems to happily use argparse so that might be too much
    # work
    sncls =[]
    sncls.append(SNCLSen('XX', 'TST5', '00',
                   ['BH0', 'BH1', 'BH2'], 'T-Compact'))
    sncls.append(SNCLSen('XX', 'TST5', '10',
                   ['BH0', 'BH1', 'BH2'], 'T-Compact'))
    sncls.append(SNCLSen('XX', 'TST6', '00',
                   ['BH0', 'BH1', 'BH2'], 'T-Compact'))

    # Here is the start time of the test
    stime = UTCDateTime('2016-196T01:00:00')
    etime = stime + 6.*60.*60.

    # Here are the channels, the stations, and locations
    # We can always put the sensor under test as the first station

    # This piece we will want to think about moving forward.  For now
    # we can run the test with each channel and assume they are the same
    # we might want to add rotation later
    chans = sncls[0].channels

    #stas = [sncl.network + '_' + sncl.station for sncl in sncls]
    stas = [sncl.station for sncl in sncls]

    locs = [sncl.location for sncl in sncls]

    # Here are the parameters
    length = 4*4096
    overlap = 2*1024

    # Grab the data and do the test by channel
    for chan in chans:

        fileStrs = itertools.starmap(getFile,
                                     zip(stas, locs,
                                         [chan]*3, [stime]*3))

        print '\n ----------------------------- \n'


        print '\n ----------------------------- \n'

        st = reduce(operator.add, map(read, fileStrs))
        st.trim(starttime=stime, endtime=etime)
        st.merge()

        # Time to compute some PSD and estimate the noise
        delta = st[0].stats.delta

        paz = grabPZ(sncls[0].sensor)
        # Need delta from first t series
        instresp = computeresp(paz, delta, length)

        compNoiseRed = lambda x: compNoise(x, length,
                                           overlap, delta, instresp)

        # Compute the noise
        noises = map(compNoiseRed, range(0, 3))

        # Plot the noise
        pltNoise(noises, st)
        # Make a time series plot
        pltTseries(st)

        writeResults(noises, st, chan)
