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

from calibrate import crop_overscans

from fram import Fram, get_night, parse_iso_time

def process_file(filename, night=None, site=None, fram=None, verbose=False):
    if fram is None:
        fram = Fram()

    if site is None:
        # Simple heuristics to derive the site name
        for _ in ['auger', 'cta-n', 'cta-s0', 'cta-s1']:
            if _ in filename:
                site = _
                break

    header = fits.getheader(filename)

    if night is None:
        time = parse_iso_time(header['DATE-OBS'])
        if header.get('LONGITUD') is not None:
            night = get_night(time, lon=header['LONGITUD'])
        else:
            night = get_night(time, site=site)

    if verbose:
        print(night,site,header['IMAGETYP'])

    # Skip old master calibrations
    if header['IMAGETYP'] in ['mdark', 'mflat']:
        return None

    # Skip frames acquired by rts2-scriptexec
    if header.get('TARGET') is None or header.get('CCD_NAME') is None:
        return None

    image = fits.getdata(filename)

    # Original (uncropped) dimensions
    width,height = header['NAXIS1'],header['NAXIS2']

    image,header = crop_overscans(image, header, subtract=False)

    cropped_width,cropped_height = image.shape[1],image.shape[0]

    # Clean up the header a bit
    header.remove('HISTORY', remove_all=True, ignore_missing=True)
    header.remove('COMMENT', remove_all=True, ignore_missing=True)
    header.remove('', remove_all=True, ignore_missing=True)
    for _ in header.keys():
        if _ and _[0] == '_':
            header.remove(_, remove_all=True, ignore_missing=True)

    type = header.get('IMAGETYP', 'unknown')

    if type == 'object' and header.get('CTYPE1'):
        wcs = WCS(header)
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

    target = header.get('TARGET')

    ccd = header.get('CCD_NAME')
    serial = header.get('product_id')
    filter = header.get('FILTER', 'unknown')
    time = parse_iso_time(header['DATE-OBS'])

    exposure = header.get('EXPOSURE')
    binning = header.get('BINNING')

    mean = np.mean(image)
    median = np.median(image)

    keywords = dict(header)

    fram.query('INSERT INTO images (filename,night,time,target,type,filter,ccd,serial,site,ra,dec,radius,exposure,width,height,cropped_width,cropped_height,footprint,footprint10,binning,mean,median,keywords) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s) ON CONFLICT (filename) DO NOTHING', (filename,night,time,target,type,filter,ccd,serial,site,ra0,dec0,radius,exposure,width,height,cropped_width,cropped_height,footprint,footprint10,binning,mean,median,keywords))

    return {'filename':filename, 'night':night, 'time':time, 'target':target, 'type':type, 'filter':filter, 'ccd':ccd, 'serial':serial, 'site':site, 'ra0':ra0, 'dec0':dec0, 'radius':radius, 'exposure':exposure, 'width':width, 'height':height, 'binning':binning, 'mean':mean, 'median':median}

def process_dir(dir, dbname='fram'):
    fram = Fram()
    fram.conn.autocommit = False

    site = None
    # Simple heuristics to derive the site name
    for _ in ['auger', 'cta-n', 'cta-s0', 'cta-s1']:
        if _ in dir:
            site = _
            break

    # Night
    night = posixpath.split(dir)[-1]

    print(night, '/', dir)
    files = glob.glob('%s/[0-9]*/*/*.fits' % dir)
    files += glob.glob('%s/darks/*/*.fits' % dir)
    files += glob.glob('%s/skyflats/*/*.fits' % dir)
    files.sort()

    res = fram.query('SELECT filename FROM images WHERE night=%s AND site=%s', (night,site), simplify=False)
    filenames = [_['filename'] for _ in res]

    j = 0

    for j,filename in enumerate(files):
        if filename in filenames:
            continue
        if 'focusing' in filename:
            continue
        # if 'bad' in filename:
        #     continue
        try:
            result = process_file(filename, night=night, site=site, fram=fram)

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

    fram.conn.commit()

    if j:
        print

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage="usage: %prog [options] arg")
    parser.add_option('-n', '--nthreads', help='Number of threads to use', action='store', dest='nthreads', type='int', default=1)
    parser.add_option('-d', '--db', help='Database name', action='store', dest='db', type='str', default='fram')
    parser.add_option('-f', '--files', help='Process files instead of directories', action='store_true', dest='process_files', default=False)
    parser.add_option('-r', '--replace', help='Replace already existing records in database', action='store_true', dest='replace', default=False)

    (options,args) = parser.parse_args()

    dirs = args

    if not dirs:
        dirs = glob.glob('/mnt/data0/auger/2*/*')
        dirs += glob.glob('/mnt/data0/cta-n/2*/*')
        dirs += glob.glob('/mnt/data2/cta-s0/2*/*')
        dirs += glob.glob('/mnt/data2/cta-s1/2*/*')

    dirs.sort(reverse=True)

    print(len(dirs), "dirs")

    if options.process_files:
        # Process single files
        fram = Fram()
        for i,filename in enumerate(args):
            try:
                print(i, '/', len(args), filename)
                if options.replace:
                    fram.query('DELETE FROM images WHERE filename=%s', (filename,))
                process_file(filename, fram=fram, verbose=True)
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
