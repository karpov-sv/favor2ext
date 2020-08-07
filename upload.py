#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import posixpath, glob, sys

from astropy.wcs import WCS
from astropy.io import fits

import warnings
from astropy.wcs import FITSFixedWarning
warnings.simplefilter('ignore', FITSFixedWarning)

from esutil import coords

# from calibrate import crop_overscans

from favor2 import Favor2, get_night, parse_time

def process_file(filename, night=None, favor2=None, verbose=False):
    if favor2 is None:
        favor2 = Favor2()

    header = fits.getheader(filename, -1)

    if night is None:
        time = parse_time(header['TIME'])
        night = get_night(time)

    if verbose:
        print(night,header['TYPE'])

    # Skip old master calibrations
    # if header['TYPE'] in ['mdark', 'mflat']:
    #     return None

    image = fits.getdata(filename, -1).astype(np.double)

    # Frame dimensions
    width,height = header['NAXIS1'],header['NAXIS2']

    # image,header = crop_overscans(image, header, subtract=False)
    image -= header.get('BASELINE', 100.0)

    # Clean up the header a bit
    header.remove('HISTORY', remove_all=True, ignore_missing=True)
    header.remove('COMMENT', remove_all=True, ignore_missing=True)
    header.remove('', remove_all=True, ignore_missing=True)
    for _ in header.keys():
        if _ and _[0] == '_':
            header.remove(_, remove_all=True, ignore_missing=True)

    type = header.get('TYPE', 'unknown')

    if type not in ['dark', 'masterdark', 'flat', 'skyflat', 'superflat', 'masterflat'] and header.get('CTYPE1'):
        wcs = WCS(header)
        if not wcs:
            print('No WCS information in', filename)
            return None

        ra,dec = wcs.all_pix2world([0, image.shape[1], 0.5*image.shape[1]], [0, image.shape[0], 0.5*image.shape[0]], 0)
        radius = 0.5*coords.sphdist(ra[0], dec[0], ra[1], dec[1])[0]
        ra0,dec0 = ra[2],dec[2]

        # Frame footprint
        ra,dec = wcs.all_pix2world([0, 0, image.shape[1], image.shape[1]], [0, image.shape[0], image.shape[0], 0], 0)
        footprint = "(" + ",".join(["(%g,%g)" % (_,__) for _,__ in zip(ra, dec)]) + ")"

        # Frame footprint at +10 pixels from the edge
        ra,dec = wcs.all_pix2world([10, 10, image.shape[1]-10, image.shape[1]-10], [10, image.shape[0]-10, image.shape[0]-10, 10], 0)
        footprint10 = "(" + ",".join(["(%g,%g)" % (_,__) for _,__ in zip(ra, dec)]) + ")"

    else:
        # Should we really discard WCS for non-object frames?
        ra0,dec0,radius = 0,0,0
        footprint,footprint10 = None,None

    filter = header.get('FILTER', 'unknown')
    time = parse_time(header['TIME'])

    exposure = header.get('EXPOSURE')
    shutter = header.get('SHUTTER')
    channel = header.get('CHANNEL ID')

    mean = np.mean(image)

    keywords = dict(header)

    favor2.query('INSERT INTO images (filename,night,time,channel,type,filter,exposure,shutter,ra,dec,radius,width,height,mean,footprint,footprint10,keywords) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s) ON CONFLICT (filename) DO NOTHING', (filename,night,time,channel,type,filter,exposure,shutter,ra0,dec0,radius,width,height,mean,footprint,footprint10,keywords))

    return {'filename':filename, 'night':night, 'time':time, 'channel':channel, 'type':type, 'filter':filter, 'shutter':shutter, 'ra0':ra0, 'dec0':dec0, 'radius':radius, 'exposure':exposure, 'width':width, 'height':height, 'mean':mean, 'keywords':keywords}

def process_dir(dir, dbname='favor2'):
    favor2 = Favor2()
    favor2.conn.autocommit = False

    # Night
    night = posixpath.split(dir)[-1]

    print(night, '/', dir)
    files = glob.glob('%s/20*.fits' % dir)
    files.sort()

    res = favor2.query('SELECT filename FROM images WHERE night=%s', (night,), simplify=False)
    filenames = [_['filename'] for _ in res]

    j = 0

    for j,filename in enumerate(files):
        if filename in filenames:
            continue
        try:
            result = process_file(filename, night=night, favor2=favor2)

            sys.stdout.write('\r  %d / %d - %s' % (j, len(files), filename))
            sys.stdout.flush()

        except KeyboardInterrupt:
                raise

        except:
            import traceback
            print("Exception while processing", filename)
            traceback.print_exc()
            pass

        #break

    favor2.conn.commit()

    if j:
        print

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage="usage: %prog [options] arg")
    parser.add_option('-n', '--nthreads', help='Number of threads to use', action='store', dest='nthreads', type='int', default=1)
    parser.add_option('-d', '--db', help='Database name', action='store', dest='db', type='str', default='favor2')
    parser.add_option('-f', '--files', help='Process files instead of directories', action='store_true', dest='process_files', default=False)
    parser.add_option('-r', '--replace', help='Replace already existing records in database', action='store_true', dest='replace', default=False)

    (options,args) = parser.parse_args()

    dirs = args

    if not dirs:
        dirs = glob.glob('AVG/*')

    dirs.sort(reverse=False)

    print(len(dirs), "dirs")

    if options.process_files:
        # Process single files
        favor2 = Favor2()
        for i,filename in enumerate(args):
            try:
                print(i, '/', len(args), filename)
                if options.replace:
                    favor2.query('DELETE FROM images WHERE filename=%s', (filename,))
                process_file(filename, favor2=favor2, verbose=True)
            except KeyboardInterrupt:
                raise
            except:
                import traceback
                print("Exception while processing", filename)
                traceback.print_exc()

    else:
        # Process directories
        if options.nthreads > 1:
            import multiprocessing
            from functools import partial

            pool = multiprocessing.Pool(options.nthreads)
            # Make wrapper function to pass our arguments inside worker processes
            fn = partial(process_dir, dbname=options.db)
            pool.map(fn, dirs, 1)

            pool.close()
            pool.join()

        else:
            for dirname in dirs:
                process_dir(dirname)
