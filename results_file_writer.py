#!/usr/bin/env python

import numpy as np

def write_results(noises, results_label, chan):
    """Write results to file"""
    pers = [(.1, 30.), (30., 100.)]
    f = open('results/' + results_label + 'RESULTS', 'w+')
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
