#!/usr/bin/env python

from sensor import Sensor

def grab_sensor():

    # TODO: take this values from mseed,
    # database or GUI application

    # T-Compact
    zeros = [0.j,
             0.j,
             -392. + 0.j,
             -1960. + 0.j,
             -1490. + 1740.j,
             -1490. - 1740.j]

    poles = [-0.3691 + 0.03702j,
             -0.3691 - 0.03702j,
             -343. + 0.j,
             -370. + 467.j,
             -370. - 467.j,
             -836. + 1522.j,
             -836. - 1522.j,
             -4900. + 4700.j,
             -4900. - 4700.j,
             -6900. + 0.j,
             -15000. + 0.j]

    gain = 4.344928*10**17

    sensitivity = 754.3*2.**26/40.

    sensor_name = 'T-Compact'

    sensor = Sensor(sensor_name,
                    zeros,
                    poles,
                    gain,
                    sensitivity)

    return sensor
