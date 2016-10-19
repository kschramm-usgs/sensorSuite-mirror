#!/usr/bin/env python

import sys
import numpy as np

from matplotlib.mlab import csd

class SensorTest():

    # TODO move delta here, when delta is available before stream
    def __init__(self, fft_length, fft_overlap):
        # Fast Fourier Transform parameters
        self.fft_length = fft_length
        self.fft_overlap = fft_overlap

    def cross_power(self, tr1, tr2, delta):
        """Calculate cross power."""
        sr = 1./delta
        cpval, fre = csd(tr1.data, tr2.data,
                         NFFT=self.fft_length, Fs=sr, noverlap=self.fft_overlap,
                         scale_by_freq=True)
        fre = fre[1:]
        cpval = cpval[1:]
        return cpval, fre

class SelfNoiseTest(SensorTest):

    def compute_noise(self, idx, instresp, st, delta):
        """Compute noise."""
        cpFix = lambda i1, i2: self.cross_power(st[i1], st[i2], delta)

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
