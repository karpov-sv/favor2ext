#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import posixpath, glob, datetime, os, sys, shutil, random

from astropy.io import fits

import sep

from favor2 import Favor2, parse_time
from calibrate import calibrate

def process_cfg(cfg, dbname=None, dbhost=None, replace=False):
    print("Processing channel %d / shutter %d / pos %d %d / filter %s / %s - %s: %d frames" % (cfg['channel'], cfg['shutter'], cfg['pos0'], cfg['pos1'], cfg['filter'], cfg['night1'], cfg['night2'], cfg['count']))

    favor2 = Favor2(dbname=dbname, dbhost=dbhost)

    res = favor2.query('select time,night,filename from images where channel=%s and shutter=%s and pos0=%s and pos1=%s and filter=%s and (type=\'survey\' or type=\'widefield\') and night>=%s and night<=%s order by time', (cfg['channel'], cfg['shutter'], cfg['pos0'], cfg['pos1'], cfg['filter'], cfg['night1'], cfg['night2']))

    night0 = res[0]['night']

    header0 = fits.getheader(res[0]['filename'], -1)
    header0['type'] = 'masterflat'
    header0['NFRAMES'] = cfg['count']
    header0['NIGHT1'] = cfg['night1']
    header0['NIGHT2'] = cfg['night2']

    darkname = favor2.find_image(res[0]['time'], type='masterdark', shutter=cfg['shutter'], channel=cfg['channel'])
    dark = fits.getdata(darkname, -1)

    basename = 'calibrations/masterflats/superflat_channel_%d_shutter_%d_pos_%d_%d_%s_%s.fits' % (cfg['channel'], cfg['shutter'], cfg['pos0'], cfg['pos1'], cfg['filter'], night0)
    if posixpath.exists(basename) and not replace:
        print('Skipping', basename)
        return

    random.shuffle(res)

    images = []
    coadd,counts = None,0

    for i,r in enumerate(res):
        filename = r['filename']

        image,header = fits.getdata(filename,-1), fits.getheader(filename, -1)
        image,header = calibrate(image, header, dark=dark)

        # Mask stars in the image
        bg = sep.Background(image)
        mask = image > bg.back() + 2.0*bg.rms()
        image[mask] = np.nan
        image /= np.nanmedian(image)

        images.append(image)

        if len(images) == 5:
            flat = np.nanmedian(images, axis=0)
            flat /= np.nanmedian(flat)
            idx = np.isfinite(flat)
            images = []

            if coadd is None:
                coadd = np.zeros_like(image)
                counts = np.zeros_like(image, dtype=np.int)

            coadd[idx] += flat[idx]
            counts[idx] += 1

            print(i, '/', len(res), basename)

    # Make masterflat and interpolate the gaps
    idx = counts > 5
    masterflat = np.ones_like(coadd)*np.nan
    masterflat[idx] = coadd[idx]/counts[idx]

    idx |= masterflat < 0.1

    bg = sep.Background(masterflat, mask=~idx)
    masterflat[~idx] = bg.back()[~idx]

    fits.writeto(basename, masterflat, header0, overwrite=True)

    print(basename)

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-d', '--db', help='Database name', action='store', dest='db', type='str', default='favor2')
    parser.add_option('-H', '--host', help='Database host', action='store', dest='dbhost', type='str', default=None)
    parser.add_option('-n', '--nthreads', help='Number of threads to use', action='store', dest='nthreads', type='int', default=0)
    # parser.add_option('-t', '--type', help='Type', action='store', dest='type', default='survey')
    parser.add_option('-r', '--replace', help='Replace', action='store_true', dest='replace', default=False)

    (options, args) = parser.parse_args()

    favor2 = Favor2(dbname=options.db, dbhost=options.dbhost)
    cfgs0 = favor2.query('select * from (select channel,shutter,pos0,pos1,filter,count(*) from images where (type=\'survey\' or type=\'widefield\') and filter=\'Clear\' group by channel,shutter,pos0,pos1,filter order by channel,shutter,pos0,pos1) t where count>100 order by count')

    cfgs = []

    for cfg in cfgs0:
        res0 = favor2.query('select night,count(*) from images where channel=%s and shutter=%s and pos0=%s and pos1=%s and filter=%s and (type=\'survey\' or type=\'widefield\') group by night order by night', (cfg['channel'], cfg['shutter'], cfg['pos0'], cfg['pos1'], cfg['filter']))

        is_first = True

        night1,night2,N = None,None,0

        for r0 in res0:
            if night1 is None:
                night1 = r0['night']
            night2 = r0['night']
            N += r0['count']

            if N > 3000:
                cfgs.append({'night1':night1, 'night2':night2, 'count':N, 'channel':cfg['channel'], 'shutter':cfg['shutter'], 'pos0':cfg['pos0'], 'pos1':cfg['pos1'], 'filter':cfg['filter']})
                night1,night2,N = None,None,0
                is_first = False

        if is_first:
            cfgs.append({'night1':night1, 'night2':night2, 'count':N, 'channel':cfg['channel'], 'shutter':cfg['shutter'], 'pos0':cfg['pos0'], 'pos1':cfg['pos1'], 'filter':cfg['filter']})

    print(len(cfgs), 'configurations loaded')

    if options.nthreads > 1:
        import multiprocessing
        from functools import partial

        pool = multiprocessing.Pool(options.nthreads)
        # Make wrapper function to pass our arguments inside worker processes
        fn = partial(process_cfg, dbname=options.db, dbhost=options.dbhost, replace=options.replace)
        pool.map(fn, cfgs, 1)

        pool.close()
        pool.join()

    else:
        for cfg in cfgs:
            process_cfg(cfg, dbname=options.db, dbhost=options.dbhost, replace=options.replace)
