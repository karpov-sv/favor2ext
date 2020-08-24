#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals

import os, sys, glob

import numpy as np
import posixpath, glob, sys, tempfile, shutil

from astropy.wcs import WCS
from astropy.io import fits

from astropy.stats import mad_std

import warnings
from astropy.wcs import FITSFixedWarning
warnings.simplefilter(action='ignore', category=FITSFixedWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

warnings.catch_warnings()

import statsmodels.api as sm
from scipy.spatial import cKDTree

import survey, calibrate
from match import Match

from favor2 import Favor2, get_night, parse_time

import sep
from esutil import htm

def process_file(filename, favor2=None, fram=None, verbose=False, replace=False, dbname=None, dbhost=None):
    #### Some parameters
    aper = 2.0
    bkgann = None
    order = 4
    bg_order = 4
    color_order = 0
    sn = 5

    #### Preparation
    header = fits.getheader(filename, -1)

    if header['TYPE'] not in ['survey', 'imaging', 'widefield', 'Swift', 'Fermi', 'test']:
        return

    channel = header.get('CHANNEL ID')
    fname = header.get('FILTER', 'unknown')
    time = parse_time(header['TIME'])
    shutter = header.get('SHUTTER', -1)

    if fname not in ['Clear']:
        return

    if fname == 'Clear':
        effective_fname = 'V'
    else:
        effective_fname = fname

    night = get_night(time)

    dirname = 'photometry/%s' % (night)
    basename = posixpath.splitext(posixpath.split(filename)[-1])[0]
    basename = dirname + '/' + basename
    catname = basename + '.cat'

    if not replace and posixpath.exists(catname):
        return

    if verbose:
        print(filename, channel, night, fname, effective_fname)

    image = fits.getdata(filename, -1).astype(np.double)

    if favor2 is None:
        favor2 = Favor2(dbname=options.db, dbhost=options.dbhost)

    if fram is None:
        fram = Favor2(dbname='fram', dbhost=options.dbhost)

    #### Basic calibration
    darkname = favor2.find_image('masterdark', header=header, debug=False)
    flatname = favor2.find_image('masterflat', header=header, debug=False)

    if darkname:
        dark = fits.getdata(darkname)
    else:
        dark = None

    if flatname:
        flat = fits.getdata(flatname)
    else:
        flat = None

    if dark is None or flat is None:
        survey.save_objects(catname, None)
        return

    image,header = calibrate.calibrate(image, header, dark=dark)

    #### Basic masking
    mask = image > 0.9*header['SATURATE']
    fmask = ~np.isfinite(flat) | (flat < 0.5)
    dmask = dark > 10.0*mad_std(dark) + np.nanmedian(dark)

    image[~fmask] *= np.median(flat[~fmask])/flat[~fmask]

    #### WCS
    wcs = WCS(header)
    pixscale = np.hypot(wcs.pixel_scale_matrix[0,0], wcs.pixel_scale_matrix[0,1])
    gain = header.get('GAIN', 1.0)
    if gain > 100:
        gain /= 1000

    #### Background mask
    mask_bg = np.zeros_like(mask)
    mask_segm = np.zeros_like(mask)

    bg2 = sep.Background(image, mask=mask|mask_bg, bw=64, bh=64)

    for _ in xrange(3):
        bg1 = sep.Background(image, mask=mask|mask_bg, bw=256, bh=256)

        ibg = bg2.back() - bg1.back()

        tmp = np.abs(ibg - np.median(ibg)) > 5.0*mad_std(ibg)
        mask_bg |= survey.dilate(tmp, np.ones([50, 50]))

    mask_bg = survey.dilate(tmp, np.ones([50, 50]))

    if np.sum(mask_bg|mask|fmask|dmask)/mask_bg.shape[0]/mask_bg.shape[1] > 0.4:
        print(100*np.sum(mask_bg|mask|fmask|dmask)/mask_bg.shape[0]/mask_bg.shape[1], '% of image masked, skipping', filename)
        survey.save_objects(catname, None)
        return

    # Frame footprint at +10 pixels from the edge
    ra,dec = wcs.all_pix2world([10, 10, image.shape[1]-10, image.shape[1]-10], [10, image.shape[0]-10, image.shape[0]-10, 10], 0)
    footprint = "(" + ",".join(["(%g,%g)" % (_,__) for _,__ in zip(ra, dec)]) + ")"

    #### Catalogue
    ra0,dec0,sr0 = survey.get_frame_center(header=header)
    cat = fram.get_stars(ra0, dec0, sr0, catalog='gaia', extra=['g<14', 'q3c_poly_query(ra, dec, \'%s\'::polygon)' % footprint], limit=1000000)

    if verbose:
        print(len(cat['ra']), 'star positions from Gaia down to g=%.1f mag' % np.max(cat['g']))

    ## Detection of blended and not really needed stars in the catalogue
    h = htm.HTM(10)
    m = h.match(cat['ra'], cat['dec'], cat['ra'], cat['dec'], 2.0*aper*pixscale, maxmatch=0)
    m = [_[m[2]>1e-5] for _ in m]

    blended = np.zeros_like(cat['ra'], dtype=np.bool)
    notneeded = np.zeros_like(cat['ra'], dtype=np.bool)

    for i1,i2,dist in zip(*m):
        if dist*3600 > 0.2*aper*pixscale:
            if cat['g'][i1] - cat['g'][i2] < 3:
                blended[i1] = True
                blended[i2] = True
            else:
                # i1 is fainter by more than 3 mag
                notneeded[i1] = True

        if dist*3600 < 0.2*aper*pixscale:
            if cat['g'][i1] > cat['g'][i2]:
                notneeded[i1] = True

    cat,blended = [_[~notneeded] for _ in cat,blended]

    #### Background subtraction
    bg = sep.Background(image, mask=mask|dmask|fmask|mask_bg, bw=64, bh=64)
    image1 = image - bg.back()

    #### Detection of all objects on the frame
    obj0 = sep.extract(image1, err=bg.rms(), thresh=2, minarea=3, mask=mask|dmask|fmask|mask_bg, filter_kernel=None, clean=False)
    obj0 = obj0[(obj0['x'] > 10) & (obj0['y'] > 10) & (obj0['x'] < image.shape[1]-10) & (obj0['y'] < image.shape[0]-10)]
    obj0 = obj0[obj0['flag'] <= 1] # We keep only normal and blended oblects

    fields = ['ra', 'dec', 'fluxerr', 'mag', 'magerr', 'flags']
    obj0 = np.lib.recfunctions.append_fields(obj0, fields, [np.zeros_like(obj0['x'], dtype=np.int if _ in ['flags'] else np.double) for _ in fields], usemask=False)
    obj0['ra'],obj0['dec'] = wcs.all_pix2world(obj0['x'], obj0['y'], 0)

    if verbose:
        print(len(obj0['x']), 'objects detected on the frame')

    ## Filter out objects not coincident with catalogue positions
    h = htm.HTM(10)
    m = h.match(obj0['ra'], obj0['dec'], cat['ra'], cat['dec'], aper*pixscale)

    nidx = np.isin(np.arange(len(obj0['ra'])), m[0], invert=True)
    obj0 = obj0[nidx]

    if verbose:
        print(len(obj0['x']), 'are outside catalogue apertures')

    ## Filter out not too significant detections
    obj0['flux'],obj0['fluxerr'],obj0['flags'] = sep.sum_circle(image1, obj0['x'], obj0['y'], aper, err=bg.rms(), gain=gain, mask=mask|dmask|fmask|mask_bg, bkgann=bkgann)
    obj0['flags'] |= obj0['flag']

    obj0 = obj0[obj0['flux'] > 0]
    obj0 = obj0[obj0['fluxerr'] > 0]

    obj0['mag'] = -2.5*np.log10(obj0['flux'])
    obj0['magerr'] = 2.5/np.log(10)*obj0['fluxerr']/obj0['flux']

    obj0 = obj0[1.0/obj0['magerr'] > sn]

    if verbose:
        print(len(obj0['x']), 'are with S/N >', sn)

    #### Forced catalogue photometry
    xc,yc = wcs.all_world2pix(cat['ra'], cat['dec'], 0)

    flux,fluxerr,flag = sep.sum_circle(image1, xc, yc, aper, err=bg.rms(), gain=gain, mask=mask|dmask|fmask|mask_bg, bkgann=bkgann)
    flag[blended] |= sep.OBJ_MERGED # Blending flag like in SEP

    mag,magerr = np.zeros_like(flux), np.zeros_like(flux)
    fidx = np.isfinite(flux) & np.isfinite(fluxerr)
    fidx[fidx] &= (flux[fidx] > 0)

    mag[fidx] = -2.5*np.log10(flux[fidx])
    magerr[fidx] = 2.5/np.log(10)*fluxerr[fidx]/flux[fidx]

    fidx[fidx] &= (magerr[fidx] > 0)
    fidx[fidx] &= 1/magerr[fidx] > sn

    obj = {'x':xc[fidx], 'y':yc[fidx], 'ra':cat['ra'][fidx], 'dec':cat['dec'][fidx],
           'mag':mag[fidx], 'magerr':magerr[fidx], 'flags':flag[fidx],
           'aper':aper}

    ## Merge detected objects into the mix
    for _ in ['x', 'y', 'ra', 'dec', 'mag', 'magerr', 'flags']:
        obj[_] = np.concatenate((obj[_], obj0[_]))

    if len(obj['x']) < 1000:
        print('Only', len(obj['x']), 'objects on the frame, skipping', filename)
        survey.save_objects(catname, None)
        return

    bgflux,_,__ = sep.sum_circle(bg.back(), obj['x'], obj['y'], aper, err=bg.rms(), gain=gain, mask=mask|dmask|fmask|mask_bg, bkgann=bkgann)
    obj['bg'] = bgflux/np.pi/aper**2
    obj['fwhm'] = 2.0*sep.flux_radius(image1, obj['x'], obj['y'], 2.0*aper*np.ones_like(obj['x']), 0.5, mask=mask|dmask|fmask|mask_bg)[0]

    #### Catalogue matching
    cidx = (cat['good'] == 1) & (cat['var'] == 0)
    if np.sum(cidx & fidx & (cat['multi_70'] == 0)) > 2000:
        cidx &= (cat['multi_70'] == 0)
        obj['cat_multi'] = 70
    elif np.sum(cidx & fidx & (cat['multi_45'] == 0)) > 1000:
        cidx &= (cat['multi_45'] == 0)
        obj['cat_multi'] = 45
    else:
        cidx &= (cat['multi_30'] == 0)
        obj['cat_multi'] = 30

    if verbose:
            print('Using %d arcsec avoidance radius' % obj['cat_multi'])

    # We match with very small SR to only account for manually placed apertures
    if verbose:
        print('Trying full fit:', np.sum(fidx), 'objects,', np.sum(cidx), 'stars')

    match = Match(width=image.shape[1], height=image.shape[0])

    if not match.match(obj=obj, cat=cat[cidx], sr=1e-5, filter_name='V', order=order, bg_order=bg_order, color_order=color_order, verbose=False) or match.ngoodstars < 500:
        if verbose:
            print(match.ngoodstars, 'good matches, retrying without additive background term')

        if not match.match(obj=obj, cat=cat[cidx], sr=1e-5, filter_name='V', order=order, bg_order=None, color_order=color_order, verbose=False) or match.ngoodstars < 100:
            if verbose:
                print('Matching failed for', filename)
            survey.save_objects(catname, None)
            return

    if verbose:
        print(match.ngoodstars, 'good matches, std =', match.std)

    #### Store objects to file
    try:
        os.makedirs(dirname)
    except:
        pass

    obj['mag_limit'] = match.mag_limit
    obj['color_term'] = match.color_term

    obj['filename'] = filename
    obj['night'] = night
    obj['channel'] = channel
    obj['filter'] = fname
    obj['cat_filter'] = match.cat_filter_name
    obj['time'] = time

    obj['mag_id'] = match.mag_id

    obj['good_idx'] = match.good_idx
    obj['calib_mag'] = match.mag
    obj['calib_magerr'] = match.magerr

    obj['std'] = match.std
    obj['nstars'] = match.ngoodstars

    survey.save_objects(catname, obj, header=header)

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage="usage: %prog [options] arg")
    parser.add_option('-d', '--db', help='Database name', action='store', dest='db', type='str', default='favor2')
    parser.add_option('-H', '--host', help='Database host', action='store', dest='dbhost', type='str', default=None)
    parser.add_option('-r', '--replace', help='Replace already existing records in database', action='store_true', dest='replace', default=False)
    parser.add_option('-v', '--verbose', help='Verbose', action='store_true', dest='verbose', default=False)

    (options,files) = parser.parse_args()

    favor2 = None # Favor2(dbname=options.db, dbhost=options.dbhost)
    fram = None # Favor2(dbname='fram', dbhost=options.dbhost)

    for i,filename in enumerate(files):
        if len(files) > 1:
            print(i, '/', len(files), filename)
        try:
            process_file(filename, favor2=favor2, fram=fram, verbose=options.verbose, replace=options.replace, dbname=options.db, dbhost=options.dbhost)
        except:
            print('\nException while processing:', filename)
            raise
