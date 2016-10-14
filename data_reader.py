#!/usr/bin/env python

from sensor import *
import itertools
import operator
from obspy.core import UTCDateTime, read

def getFile(sta, loc, chan, stime):
    # TODO get more options to find files
    #fileStr = '/tr1/telemetry_days/' + sta + '/'
    #fileStr += str(stime.year) + '/' + str(stime.year) + '_'
    # fileStr = str(stime.julday).zfill(3) # + '/'
    fileStr = ''
    fileStr += 'data/'
    fileStr += sta + '_'
    fileStr += loc + '_' + chan + '.512.seed'
    return fileStr

def get_seismic_data_stream(sncls,
                            start_time,
                            end_time,
                            channel_index = None):
    """Takes mseed files and retuns Obspy seismic data stream"""

    # Assume that all station have the same channels
    # Assume 3 channels
    # TODO: check

    channels_len = len(sncls[0].channels)
    channels = []

    if channel_index == None:
        channels = sncls[0].channels
    else:
        channels = [ sncls[0].channels[channel_index] ] * channels_len

    stas = [sncl.station for sncl in sncls]
    locs = [sncl.location for sncl in sncls]

    fileStrs = itertools.starmap(getFile, zip(stas,
                                              locs,
                                              channels,
                                              [start_time] * channels_len))

    st = reduce(operator.add, map(read, fileStrs))
    st.trim(starttime=start_time, endtime=end_time)
    # removes redudant data or traces in the stream
    st.merge()

    return st
