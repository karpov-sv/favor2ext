#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import posixpath, glob, datetime, os, sys

from astropy.io import fits

from favor2 import fix_remote_path

class Calibrator:
    def __init__(self, shutter=0, channel=None):
        if not channel:
            self.channel = self.get_channel_id()
        else:
            self.channel = channel

        self.shutter = shutter
        self.header = None

        self.load_calibrations()

    def get_channel_id():
        import socket
        import re

        m = re.match("^mmt(\d+)", socket.gethostname())

        if m:
            id = int(m.group(1))

            if id:
                #print "Channel ID=%d" % id
                channel = id
            else:
                channel = 0

        return channel

    def load_calibrations(self):
        base = 'calibrations/shutter_%d_channel_%d' % (self.shutter, self.channel)
        self.bias = self.load_file(base + '_bias.fits')
        self.mask = self.load_file(base + '_mask.fits')
        self.thresh = self.load_file(base + '_thresh.fits')

        self.scale1_0 = self.load_file(base + '_scale1_0.fits')
        self.scale1_1 = self.load_file(base + '_scale1_1.fits')
        self.scale1_2 = self.load_file(base + '_scale1_2.fits')
        self.scale1_3 = self.load_file(base + '_scale1_3.fits')
        self.scale1_4 = self.load_file(base + '_scale1_4.fits')

        self.scale2_0 = self.load_file(base + '_scale2_0.fits')
        self.scale2_1 = self.load_file(base + '_scale2_1.fits')
        self.scale2_2 = self.load_file(base + '_scale2_2.fits')
        self.scale2_3 = self.load_file(base + '_scale2_3.fits')
        self.scale2_4 = self.load_file(base + '_scale2_4.fits')

    def load_file(self, filename):
        fullname = fix_remote_path(filename, self.channel)

        #print filename, fullname, posixpath.exists(fullname)

        if self.header is None:
            self.header = fits.getheader(fullname, -1)

        return fits.getdata(fullname, -1).astype(np.double)

    def calibrate(self, image, header=None, mask_nans=False, mask_negative=True):
        saturation = self.header['SATURATE']
        bias0 = self.header['BIAS0']

        if header and header.get('SIMPLEPREAMP', 2) == 1: # 11-bit mode, GLOBAL shutter
            saturation = 1500

        v = image.astype(np.float) - bias0
        idx0 = v <= 0.1
        v[idx0] = 0.1

        lv = np.log10(v)
        lv2 = lv*lv
        lv3 = lv2*lv
        lv4 = lv2*lv2

        v1 = self.scale1_0*lv4 + self.scale1_1*lv3 + self.scale1_2*lv2 + self.scale1_3*lv + self.scale1_4
        idx = v1 > 0
        v1[idx] = v[idx]/v1[idx]

        if mask_negative:
            v1[v1 <= 0] = saturation

        v2 = self.scale2_0*lv4 + self.scale2_1*lv3 + self.scale2_2*lv2 + self.scale2_3*lv + self.scale2_4
        idx = v2 > 0
        v2[idx] = v[idx]/v2[idx]

        if mask_negative:
            v2[v2 <= 0] = saturation

        idx = v > self.thresh
        result = v1
        result[idx] = v2[idx]

        result -= self.bias

        # Saturation
        result[~np.isfinite(result)] = saturation

        if mask_negative:
            result[idx0] = saturation

        if mask_nans:
            result[self.mask > 0] = saturation
        else:
            result[self.mask > 1] = saturation

        result[result > saturation] = saturation

        if header is not None:
            header = header.copy()
            header['SATURATE'] = saturation

            return result,header

        return result

_calibrators = {}

def calibrate(image, header=None, channel=None, shutter=None, dark=None, calibrator=None):
    if channel is None:
        channel = header['CHANNEL ID']
    if shutter is None:
        shutter = header['SHUTTER']

    key = '%d_%d' % (channel, shutter)

    if calibrator is None:
        if _calibrators and key in _calibrators:
            calibrator = _calibrators[key]
        else:
            calibrator = Calibrator(shutter=shutter, channel=channel)

    if header and header.get('SERIAL') in ['SCC-01798']:
    # Temporarily disable 2nd channel with 10th channel camera
        calibrator = None

    if calibrator:
        if not key in _calibrators:
            _calibrators[key] = calibrator

        if header is not None:
            image,header = calibrator.calibrate(image, header)
        else:
            image = calibrator.calibrate(image)

        if dark is not None:
            dark = calibrator.calibrate(dark)
            image -= dark
    else:
        if dark is not None:
            image -= dark

    if header is not None:
        return image,header
    else:
        return image
