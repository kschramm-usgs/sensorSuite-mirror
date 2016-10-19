#!/usr/bin/env python

import sys
from obspy.signal.invsim import paz_to_freq_resp
import numpy as np

class SNCL:
    """Station Network Channels Location."""

    def __init__(self, sta, net, channels, loc):
        self.station = sta
        self.network = net
        self.channels = channels
        self.location = loc


class Sensor:
    """Seismic sensor."""

    def __init__(self, name, zeros, poles, gain, sensitivity):

        self.zeros = zeros
        self.poles = poles
        self.gain = gain
        self.sensitivity = sensitivity

    def compute_response(self, delta, fft_length):
        """Convert Poles and Zeros (PAZ) to frequency response.
           and ?

        Parameters
        ----------
        paz : type
            Description of parameter `x`.
        """

        respval = paz_to_freq_resp(self.poles,
                                   self.zeros,
                                   self.sensitivity * self.gain,
                                   t_samp=delta,
                                   nfft=fft_length, freq=False)

        respval = np.absolute(respval*np.conjugate(respval))
        respval = respval[1:]
        return respval
