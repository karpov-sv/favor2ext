from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

import tempfile, datetime, posixpath, shutil, re, os

import ephem
import socket
import re

from db import DB

def get_time_from_string(filename):
    m = re.match("^(\d\d)\.(\d\d)\.(\d\d\d\d) (\d\d):(\d\d):(\d\d).(\d\d\d)", filename)
    if m:
        t = [int(m.group(i)) for i in range(1,8)]
        return datetime.datetime(t[2], t[1], t[0], t[3], t[4], t[5], t[6]*1000)
    else:
        return None

def get_time_from_filename(filename):
    filename = posixpath.split(filename)[-1]
    m = re.match("^(\d\d\d\d)(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)(\d\d\d)", filename)
    if m:
        t = [int(m.group(i)) for i in range(1,8)]
        return datetime.datetime(t[0], t[1], t[2], t[3], t[4], t[5], t[6]*1000)
    else:
        return None

def split_path(filename):
    paths = []

    while posixpath.basename(filename):
        paths.append(posixpath.basename(filename))
        filename = posixpath.dirname(filename)

    paths.reverse()

    return paths

def fix_remote_path(filename, channel_id=0):
    return filename

def get_time_from_night(night):
    m = re.match("^(\d\d\d\d)_(\d\d)_(\d\d)", night)
    if m:
        t = [int(m.group(i)) for i in range(1,4)]
        return datetime.datetime(t[0], t[1], t[2])

def parse_time(string):
    return datetime.datetime.strptime(string, '%d.%m.%Y %H:%M:%S.%f')

def get_night(t):
    et = t - datetime.timedelta(hours=12)
    return "%04d_%02d_%02d" % (et.year, et.month, et.day)

class Favor2(DB):
    def __init__(self, base='.', readonly=False, latitude=43.649861, longitude=41.4314722, elevation=2030, **kwargs):
        DB.__init__(self, **kwargs)

        self.base = base

        self.obs = ephem.Observer()
        self.obs.lat = np.deg2rad(latitude)
        self.obs.lon = np.deg2rad(longitude)
        self.obs.elevation = elevation

        self.moon = ephem.Moon()
        self.sun = ephem.Sun()

    def tempdir(self, prefix='favor2'):
        return tempfile.mkdtemp(prefix=prefix)

    def find_image(self, type='survey', time=None, channel=None, shutter=None, filter=None, pos0=None, pos1=None, header=None, before=True, debug=False):
        if channel is None and header is not None:
            channel = header.get('CHANNEL ID', 0)

        if shutter is None and header is not None:
            shutter = header.get('SHUTTER', 0)

        if type in ['masterflat']:
            if pos0 is None and header is not None:
                pos0 = header.get('MIRROR_POS0', 0)

            if pos1 is None and header is not None:
                pos1 = header.get('MIRROR_POS1', 0)

        if type not in ['dark', 'masterdark']:
            if filter is None and header is not None:
                filter = header.get('FILTER', 'Clear')

        if time is None and header is not None:
            time = header['TIME']

        if isinstance(time, basestring) or isinstance(time, str):
            time = parse_time(time)

        conditions,args = [],[]

        if type is not None:
            conditions.append('type=%s')
            args.append(type)

        if channel is not None:
            conditions.append('channel=%s')
            args.append(channel)

        if shutter is not None:
            conditions.append('shutter=%s')
            args.append(shutter)

        if filter is not None:
            conditions.append('filter=%s')
            args.append(filter)

        if pos0 is not None:
            conditions.append('pos0=%s')
            args.append(pos0)

        if pos1 is not None:
            conditions.append('pos1=%s')
            args.append(pos1)

        res = None

        if before:
            res = self.query('SELECT filename FROM IMAGES WHERE ' + ' AND '.join(conditions + ['time <= %s']) + ' ORDER BY time DESC LIMIT 1;', args + [time], debug=debug)

        if not res or not before:
            res = self.query('SELECT filename FROM IMAGES WHERE ' + ' AND '.join(conditions + ['time > %s']) + ' ORDER BY time DESC LIMIT 1;', args + [time], debug=debug)

        return res

    def get_night(self, t):
        et = t - datetime.timedelta(hours=12)
        return "%04d_%02d_%02d" % (et.year, et.month, et.day)

    # Normalize image using dark and flat frames from the archive
    def prepare_image(self, filename, outname):
        print("Normalizing %s -> %s" % (filename, outname))

        img = fits.open(filename)
        time = get_time_from_string(img[-1].header['TIME'])
        channel_id = int(img[-1].header['CHANNEL_ID']) or -1

        dark_filename = self.find_image(time, 'dark', channel_id)
        if dark_filename:
            print("Subtracting DARK: %s" % (dark_filename))
            dark = fits.open(posixpath.join(self.base, dark_filename))

            if dark[-1].shape == img[-1].shape:
                img[-1].data -= dark[-1].data
                img[-1].header['comment'] = 'Dark frame subtracted: %s' % (dark_filename)

                flat_filename = self.find_image(time, 'flat', channel_id)
                if flat_filename:
                    print("Normalizing to FLAT: %s" % (flat_filename))
                    flat = pyfits.open(posixpath.join(self.base, flat_filename))

                    if flat[-1].shape == img[-1].shape:
                        img[-1].data *= np.mean(flat[-1].data - dark[-1].data)/(flat[-1].data - dark[-1].data)
                        img[-1].header['comment'] = 'Flatfield corrected: %s' % (flat_filename)

        img.writeto(outname, overwrite=True)

    # Register the image in the database
    def register_image(self, filename, imtype='avg', time=None, channel_id=0, filter=None, ra0=None, dec0=None):
        print("Inserting %s (type=%s) into database" % (filename, imtype))

        img = fits.open(filename)
        h = img[-1].header

        if not time:
            time = get_time_from_string(h['TIME'])
        elif type(time) == str:
            time = get_time_from_string(time)

        if not channel_id:
            channel_id = h['CHANNEL ID']

        if not filter:
            filter = h['FILTER']

        night = self.get_night(time)
        width = img[-1].shape[1]
        height = img[-1].shape[0]

        if not ra0 or not dec0:
            w = WCS(header=h)
            if w.wcs.has_cd():
                [ra0],[dec0] = w.all_pix2world([0.5*width], [0.5*height], 0)
            else:
                ra0,dec0 = 0,0

        # header to dict
        keywords = dict(h)
        for kw in ["SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "EXTEND", "COMMENT", "BZERO", "BSCALE", "HISTORY"]:
            while keywords.pop(kw, None):
                pass

        # Fix the filename
        if '/MMT/' in filename:
            filename = filename[7:]
        elif '/data/' in filename:
            filename = filename[6:]

        # FIXME: escape quotes!
        keywords = ",".join(['"%s"=>"%s"'%(k,keywords[k]) for k in keywords.keys()])
        #print(keywords)

        self.query('DELETE FROM images WHERE type=%s AND filename=%s', (imtype, filename))

        self.query("INSERT INTO images (filename, night, channel_id, time, type, ra0, dec0, width, height, keywords, filter) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s::hstore, get_filter_id(%s))",
                   (filename, night, channel_id, time, imtype, ra0, dec0, width, height, keywords, filter))
