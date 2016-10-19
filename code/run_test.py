#!/usr/bin/env python

import numpy as np

from obspy.core import UTCDateTime

from sensor import Sensor, SNCL
from data_reader import get_seismic_data_stream
from grab_sensors import grab_sensor
from results_file_writer import write_results
from self_noise_test import SelfNoiseTest
from seismic_data_plotter import *

if __name__ == '__main__':

    # B broadband
    # H high gain

    network = 'XX'
    channels = ['BH0', 'BH1', 'BH2']

    sncls =[]
    sncls.append( SNCL('TST5', network, channels, '00') )
    sncls.append( SNCL('TST5', network, channels, '10') )
    sncls.append( SNCL('TST6', network, channels, '00') )

    # TODO: encapsulate test
    # start time of the test
    start_time = UTCDateTime('2016-196T01:00:00')
    # end time of the test
    end_time = start_time + 6.*60.*60.

    # TODO: grab sensor from database, file, etc, or from metadata using sncl
    sensor = grab_sensor()

    # what is the right place to put this global? maybe here
    # the other thing is to encapsulate the fft itself with our fixed
    # version with less options
    fft_length = 4*4096
    fft_overlap = 2*1024

    self_noise_test = SelfNoiseTest(fft_length, fft_overlap)

    # TODO: encapsulate everything about the test inside the class
    # tests = [self_noise_test]

    # Grab the data and do the test by channel
    for channel_index, channel in enumerate(channels):

        print '\n ----------------------------- \n'
        print '\n ' + channel + '\n'

        st = get_seismic_data_stream(sncls,
                                     start_time,
                                     end_time,
                                     channel_index)

        # delta from first t series
        delta = st[0].stats.delta

        # Compute PSD
        # TODO: this only needs to be computed once per sensor
        # is property of the sensor, not of the channel
        # needs delta, which right now comes from st
        instresp = sensor.compute_response(delta, fft_length)

        # Compute the noise
        compNoiseReduced = lambda x: self_noise_test.compute_noise(x,
                                                                   instresp,
                                                                   st,
                                                                   delta)

        noises = map(compNoiseReduced, range(0, 3))

        # Plot the noise
        plot_noise(noises, st, delta)
        # Make a time series plot
        plot_time_series(st, delta)

        # Write results to file
        results_label = mkFigLabel(st)
        write_results(noises, results_label, channel)
