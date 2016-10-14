#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import matplotlib as mpl

#TODO encapsulate this, no scripting in middle of definitions!
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font', size=12)

def mkFigLabel(st):
    figlabel = str(st[0].stats.starttime.year)
    figlabel += str(st[0].stats.starttime.julday)
    figlabel += str(st[0].stats.starttime.hour)
    figlabel += str(st[0].stats.starttime.minute)
    for idx in range(0, 3):
        figlabel += st[idx].stats.station + st[idx].stats.location
        figlabel += st[idx].stats.channel
    return figlabel


def plot_time_series(st, delta):
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
    plt.savefig('plots/TSERIES' + figlabel + '.jpg', format='jpeg', dpi=400)
    plt.clf()

    # Do not bother returning the figure handle just save it
    return


def plot_noise(noises, st, delta):
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
    plt.savefig('plots/NOISE' + figlabel + '.jpg', format='jpeg', dpi=400)

    # Don't bother keeping the figure handle
    plt.clf()
    return
