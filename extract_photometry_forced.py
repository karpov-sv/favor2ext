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

# from scipy.interpolate import griddata
# from loess.loess_2d import loess_2d

# Flags, mostly analogous to SEP ones
FLAG_NORMAL = 0x0
FLAG_BLENDED = 0x1
FLAG_TRUNCATED = 0x2
FLAG_MASKED = 0x20
FLAG_NO_BACKGROUND = 0x100
FLAG_BAD_CALIBRATION = 0x200

def process_file(filename, favor2=None, verbose=False, replace=False, dbname=None, dbhost=None, photodir='photometry'):
    #### Some parameters
    aper = 2.0
    bkgann = None
    order = 4
    bg_order = 4
    color_order = 2
    sn = 5

    if not posixpath.exists(filename):
        return None

    # Rough but fast checking of whether the file is already processed
    if not replace and posixpath.exists(photodir + '/' + filename.split('/')[-2] + '/' + posixpath.splitext(posixpath.split(filename)[-1])[0] + '.cat'):
        return

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

    dirname = '%s/%s' % (photodir, night)
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

    # Check whether the calibration failed
    if 'SATURATE' not in header:
        print('Calibration failed for', filename)
        survey.save_objects(catname, None)
        return

    #### Basic masking
    mask = image > 0.9*header['SATURATE']
    fmask = ~np.isfinite(flat) | (flat < 0.5)
    dmask = dark > 10.0*mad_std(dark) + np.nanmedian(dark)

    if header.get('BLEMISHCORRECTION', 1):
        # We have to mask blemished pixels
        blemish = fits.getdata('calibrations/blemish_shutter_%d_channel_%d.fits' % (header['SHUTTER'], header['CHANNEL ID']))
        if verbose:
            print(100*np.sum(blemish>0)/blemish.shape[0]/blemish.shape[1], '% pixels blemished')
        dmask |= blemish > 0

    image[~fmask] *= np.median(flat[~fmask])/flat[~fmask]

    #### WCS
    wcs = WCS(header)
    pixscale = np.hypot(wcs.pixel_scale_matrix[0,0], wcs.pixel_scale_matrix[0,1])
    gain = 0.67 if header.get('SHUTTER') == 0 else 1.9

    #### Background mask
    mask_bg = np.zeros_like(mask)
    mask_segm = np.zeros_like(mask)

    # bg2 = sep.Background(image, mask=mask|mask_bg, bw=64, bh=64)

    # for _ in xrange(3):
    #     bg1 = sep.Background(image, mask=mask|mask_bg, bw=256, bh=256)

    #     ibg = bg2.back() - bg1.back()

    #     tmp = np.abs(ibg - np.median(ibg)) > 5.0*mad_std(ibg)
    #     mask_bg |= survey.dilate(tmp, np.ones([50, 50]))

    # mask_bg = survey.dilate(tmp, np.ones([50, 50]))

    # Large objects?..
    bg = sep.Background(image, mask=mask|dmask|fmask|mask_bg|mask_segm, bw=128, bh=128)
    image1 = image - bg.back()
    obj0,segm = sep.extract(image1, err=bg.rms(), thresh=10, minarea=10, mask=mask|dmask|fmask|mask_bg, filter_kernel=None, clean=False, segmentation_map=True)

    mask_segm = np.isin(segm, [_+1 for _,npix in enumerate(obj0['npix']) if npix > 500])
    mask_segm = survey.dilate(mask_segm, np.ones([20, 20]))

    if np.sum(mask_bg|mask_segm|mask|fmask|dmask)/mask_bg.shape[0]/mask_bg.shape[1] > 0.4:
        print(100*np.sum(mask_bg|mask_segm|mask|fmask|dmask)/mask_bg.shape[0]/mask_bg.shape[1], '% of image masked, skipping', filename)
        survey.save_objects(catname, None)
        return
    elif verbose:
        print(100*np.sum(mask_bg|mask_segm|mask|fmask|dmask)/mask_bg.shape[0]/mask_bg.shape[1], '% of image masked')

    # Frame footprint at +10 pixels from the edge
    ra,dec = wcs.all_pix2world([10, 10, image.shape[1]-10, image.shape[1]-10], [10, image.shape[0]-10, image.shape[0]-10, 10], 0)
    footprint = "(" + ",".join(["(%g,%g)" % (_,__) for _,__ in zip(ra, dec)]) + ")"

    #### Catalogue
    ra0,dec0,sr0 = survey.get_frame_center(header=header)
    cat = favor2.get_stars(ra0, dec0, sr0, catalog='gaia', extra=['g<14', 'q3c_poly_query(ra, dec, \'%s\'::polygon)' % footprint], limit=1000000)

    if verbose:
        print(len(cat['ra']), 'star positions from Gaia down to g=%.1f mag' % np.max(cat['g']))

    ## Detection of blended and not really needed stars in the catalogue
    h = htm.HTM(10)
    m = h.match(cat['ra'], cat['dec'], cat['ra'], cat['dec'], 2.0*aper*pixscale, maxmatch=0)
    m = [_[m[2]>1e-5] for _ in m]

    blended = np.zeros_like(cat['ra'], dtype=np.bool)
    notneeded = np.zeros_like(cat['ra'], dtype=np.bool)

    for i1,i2,dist in zip(*m):
        if dist*3600 > 0.5*aper*pixscale:
            if cat['g'][i1] - cat['g'][i2] < 3:
                blended[i1] = True
                blended[i2] = True
            else:
                # i1 is fainter by more than 3 mag
                notneeded[i1] = True

        if dist*3600 < 0.5*aper*pixscale:
            if cat['g'][i1] > cat['g'][i2]:
                notneeded[i1] = True

    cat,blended = [_[~notneeded] for _ in cat,blended]

    #### Background subtraction
    bg = sep.Background(image, mask=mask|dmask|fmask|mask_bg|mask_segm, bw=128, bh=128)
    # bg = sep.Background(image, mask=mask|dmask|fmask|mask_bg|mask_segm, bw=32, bh=32)
    image1 = image - bg.back()

    #### Detection of all objects on the frame
    obj0,segm = sep.extract(image1, err=bg.rms(), thresh=2, minarea=3, mask=mask|dmask|fmask|mask_bg|mask_segm, filter_kernel=None, clean=False, segmentation_map=True)
    obj0 = obj0[(obj0['x'] > 10) & (obj0['y'] > 10) & (obj0['x'] < image.shape[1]-10) & (obj0['y'] < image.shape[0]-10)]
    obj0 = obj0[obj0['flag'] <= 1] # We keep only normal and blended oblects

    fields = ['ra', 'dec', 'fluxerr', 'mag', 'magerr', 'flags', 'cat']
    obj0 = np.lib.recfunctions.append_fields(obj0, fields, [np.zeros_like(obj0['x'], dtype=np.int if _ in ['flags', 'cat'] else np.double) for _ in fields], usemask=False)
    obj0['ra'],obj0['dec'] = wcs.all_pix2world(obj0['x'], obj0['y'], 0)
    obj0['flags'] = obj0['flag']

    if verbose:
        print(len(obj0['x']), 'objects detected on the frame')

    ## Filter out objects not coincident with catalogue positions
    h = htm.HTM(10)
    m = h.match(obj0['ra'], obj0['dec'], cat['ra'], cat['dec'], aper*pixscale)

    nidx = np.isin(np.arange(len(obj0['ra'])), m[0], invert=True)
    obj0 = obj0[nidx]

    if verbose:
        print(len(obj0['x']), 'are outside catalogue apertures')

    # Catalogue stars
    xc,yc = wcs.all_world2pix(cat['ra'], cat['dec'], 0)
    obj = {'x':xc, 'y':yc, 'ra':cat['ra'], 'dec':cat['dec']}

    obj['flags'] = np.zeros_like(xc, dtype=np.int)
    obj['flags'][blended] |= FLAG_BLENDED

    obj['cat'] = np.ones_like(xc, dtype=np.int)

    for _ in ['mag', 'magerr', 'flux', 'fluxerr']:
        obj[_] = np.zeros_like(xc)

    # Merge detected objects
    for _ in ['x', 'y', 'ra', 'dec', 'flags', 'mag', 'magerr', 'flux', 'fluxerr', 'cat']:
        obj[_] = np.concatenate((obj[_], obj0[_]))

    if verbose:
        print(len(obj['x']), 'objects for photometry')

    # Simple aperture photometry
    obj['flux'],obj['fluxerr'],flag = sep.sum_circle(image1, obj['x'], obj['y'], aper, err=bg.rms(), gain=gain, mask=mask|dmask|fmask|mask_bg|mask_segm, bkgann=bkgann)
    obj['flags'] |= flag
    # Normalize flags
    obj['flags'][obj['flags'] & sep.APER_TRUNC] |= FLAG_TRUNCATED
    obj['flags'][obj['flags'] & sep.APER_ALLMASKED] |= FLAG_MASKED
    obj['flags'] &= FLAG_NORMAL | FLAG_BLENDED | FLAG_TRUNCATED | FLAG_MASKED | FLAG_NO_BACKGROUND | FLAG_BAD_CALIBRATION

    area,_,_ = sep.sum_circle(np.ones_like(image1), obj['x'], obj['y'], aper, err=bg.rms(), gain=gain, mask=mask|dmask|fmask|mask_bg|mask_segm, bkgann=bkgann)

    # Simple local background estimation
    bgflux,bgfluxerr,bgflag = sep.sum_circann(image1, obj['x'], obj['y'], 10, 15, err=bg.rms(), gain=gain, mask=mask|dmask|fmask|mask_bg|mask_segm|(segm>0))
    bgarea,_,_ = sep.sum_circann(np.ones_like(image1), obj['x'], obj['y'], 10, 15, err=bg.rms(), gain=gain, mask=mask|dmask|fmask|mask_bg|mask_segm|(segm>0))

    bgidx = np.isfinite(bgarea) & np.isfinite(area)
    bgidx[bgidx] &= (bgarea[bgidx] > 10) & (area[bgidx] > 1)

    obj['flux'][bgidx] -= bgflux[bgidx]*area[bgidx]/bgarea[bgidx]
    obj['flags'][~bgidx] |= FLAG_NO_BACKGROUND # No local background

    obj['deltabgflux'] = np.zeros_like(obj['x'])
    obj['deltabgflux'][bgidx] = bgflux[bgidx]*area[bgidx]/bgarea[bgidx]

    fidx = np.isfinite(obj['flux']) & np.isfinite(obj['fluxerr'])
    fidx[fidx] &= (obj['flux'][fidx] > 0)

    obj['mag'][fidx] = -2.5*np.log10(obj['flux'][fidx])
    obj['magerr'][fidx] = 2.5/np.log(10)*obj['fluxerr'][fidx]/obj['flux'][fidx]

    fidx[fidx] &= (obj['magerr'][fidx] > 0)
    fidx[fidx] &= 1/obj['magerr'][fidx] > sn

    for _ in obj.keys():
        if hasattr(obj[_], '__len__'):
            obj[_] = obj[_][fidx]

    obj['aper'] = aper

    if verbose:
        print(len(obj['x']), 'objects with S/N >', sn)

    if len(obj['x']) < 1000:
        print('Only', len(obj['x']), 'objects on the frame, skipping', filename)
        survey.save_objects(catname, None)
        return

    obj['fwhm'] = 2.0*sep.flux_radius(image1, obj['x'], obj['y'], 2.0*aper*np.ones_like(obj['x']), 0.5, mask=mask|dmask|fmask|mask_bg|mask_segm)[0]

    #### Check FWHM of all objects and select only 'good' ones
    idx = obj['flags'] == 0
    idx &= obj['magerr'] < 1/20

    fwhm0 = survey.fit_2d(obj['x'][idx], obj['y'][idx], obj['fwhm'][idx], obj['x'], obj['y'], weights=1/obj['magerr'][idx])

    fwhm_idx = np.abs(obj['fwhm'] - fwhm0 - np.median((obj['fwhm'] - fwhm0)[idx])) < 3.0*mad_std((obj['fwhm'] - fwhm0)[idx])
    obj['flags'][~fwhm_idx] |= FLAG_BLENDED

    #### Catalogue matching
    idx = obj['flags']
    m = htm.HTM(10).match(obj['ra'], obj['dec'], cat['ra'], cat['dec'], 1e-5)
    fidx = np.in1d(np.arange(len(cat['ra'])), m[1]) # Stars that got successfully measured and not blended

    cidx = (cat['good'] == 1) & (cat['var'] == 0)
    cidx &= np.isfinite(cat['B']) & np.isfinite(cat['V']) # & np.isfinite(cat['lum'])
    cidx[cidx] &= ((cat['B'] - cat['V'])[cidx] > -0.5) & ((cat['B'] - cat['V'])[cidx] < 2.0)
    # cidx[cidx] &= (cat['lum'][cidx] > 0.3) & (cat['lum'][cidx] < 30)

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
            print(np.sum(obj['flags'] == 0), 'objects without flags')
            print('Amount of good stars:',
                  np.sum(cidx & fidx & (cat['multi_70'] == 0)),
                  np.sum(cidx & fidx & (cat['multi_45'] == 0)),
                  np.sum(cidx & fidx & (cat['multi_30'] == 0)))
            print('Using %d arcsec avoidance radius' % obj['cat_multi'])

    # We match with very small SR to only account for manually placed apertures
    if verbose:
        print('Trying full fit:', len(obj['x']), 'objects,', np.sum(cidx), 'stars')

    match = Match(width=image.shape[1], height=image.shape[0])

    prev_ngoodstars = len(obj['x'])

    for iter in range(10):
        if not match.match(obj=obj, cat=cat[cidx], sr=1e-5, filter_name='V', order=order, bg_order=bg_order, color_order=color_order, verbose=False) or match.ngoodstars < 500:
            if verbose:
                print(match.ngoodstars, 'good matches, matching failed for', filename)
            survey.save_objects(catname, None)
            return

        if verbose:
            print(match.ngoodstars, 'good matches, std =', match.std)

        if match.ngoodstars == prev_ngoodstars:
            if verbose:
                print('Converged on iteration', iter)
            break
        prev_ngoodstars = match.ngoodstars

        # Match good objects with stars
        oidx = obj['flags'] == 0
        oidx1,cidx1,dist1 = htm.HTM(10).match(obj['ra'][oidx], obj['dec'][oidx], cat['ra'][cidx], cat['dec'][cidx], 1e-5)

        x = obj['x'][oidx][oidx1]
        y = obj['y'][oidx][oidx1]
        cbv = match.color_term[oidx][oidx1]
        cbv2 = match.color_term2[oidx][oidx1]
        cbv3 = match.color_term3[oidx][oidx1]
        bv = (cat['B'] - cat['V'])[cidx][cidx1]
        cmag = cat[match.cat_filter_name][cidx][cidx1]
        mag = match.mag[oidx][oidx1] + bv*cbv + bv**2*cbv2 + bv**3*cbv3
        magerr = np.hypot(obj['magerr'][oidx][oidx1], 0.02)

        dmag = mag-cmag
        ndmag = ((mag-cmag)/magerr)

        idx = cmag < match.mag_limit[oidx][oidx1]

        x,y,cbv,cbv2,cbv3,bv,cmag,mag,magerr,dmag,ndmag = [_[idx] for _ in [x,y,cbv,cbv2,cbv3,bv,cmag,mag,magerr,dmag,ndmag]]

        # Match all objects with good objects
        xy = np.array([x,y]).T
        xy0 = np.array([obj['x'], obj['y']]).T

        kd = cKDTree(xy)

        dist,m = kd.query(xy0, 101)
        dist = dist[:,1:]
        m = m[:,1:]

        vchi2 = mad_std(ndmag[m]**2, axis=1)

        # Mark regions of too sparse or too noisy matches as bad
        obj['flags'][vchi2 > 5] |= FLAG_BAD_CALIBRATION
        # obj['flags'][vchi2 > np.median(vchi2) + 5.0*mad_std(vchi2)] |= FLAG_BAD_CALIBRATION
        obj['flags'][dist[:,10] > np.median(dist[:,10]) + 10.0*mad_std(dist[:,10])] |= FLAG_BAD_CALIBRATION

    match.good_idx = (obj['flags'] & FLAG_BAD_CALIBRATION) == 0
    if verbose:
        print(np.sum(match.good_idx), 'of', len(match.good_idx), 'stars are good')

    #### Store objects to file
    try:
        os.makedirs(dirname)
    except:
        pass

    obj['mag_limit'] = match.mag_limit
    obj['color_term'] = match.color_term
    obj['color_term2'] = match.color_term2
    obj['color_term3'] = match.color_term3

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
    parser.add_option('-p', '--photodir', help='Base dir to store photometry', action='store', dest='photodir', type='str', default='photometry')
    parser.add_option('-r', '--replace', help='Replace already existing records in database', action='store_true', dest='replace', default=False)
    parser.add_option('-v', '--verbose', help='Verbose', action='store_true', dest='verbose', default=False)

    (options,files) = parser.parse_args()

    if len(files) > 1:
        favor2 = Favor2(dbname=options.db, dbhost=options.dbhost)
    else:
        favor2 = None

    for i,filename in enumerate(files):
        if len(files) > 1:
            print(i, '/', len(files), filename)
        try:
            process_file(filename, favor2=favor2, verbose=options.verbose, replace=options.replace, dbname=options.db, dbhost=options.dbhost, photodir=options.photodir)
        except KeyboardInterrupt:
            raise
        except:
            print('\nException while processing:', filename, file=sys.stderr)
            import traceback
            traceback.print_exc()
