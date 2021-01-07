from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import posixpath, glob, datetime, os, sys, tempfile, shutil

from astropy.wcs import WCS
from astropy.io import fits
from astropy.stats import mad_std

import warnings
from astropy.wcs import FITSFixedWarning
warnings.simplefilter(action='ignore', category=FITSFixedWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

from esutil import coords, htm
import statsmodels.api as sm
from scipy.spatial import cKDTree
from scipy.stats import binned_statistic
from scipy.interpolate import interp1d
from scipy.optimize import bisect
from scipy.special import legendre as leg

import sep
try:
    import cPickle as pickle
except:
    import pickle
import json

try:
    import cv2
    # Much faster dilation
    dilate = lambda image,mask: cv2.dilate(image.astype(np.uint8), mask).astype(np.bool)
except:
    from scipy.signal import fftconvolve
    dilate = lambda image,mask: fftconvolve(image, mask, mode='same') > 0.9

def get_frame_center(filename=None, header=None, wcs=None, width=None, height=None):
    if not wcs:
        if header:
            wcs = WCS(header=header)
        elif filename:
            header = fits.getheader(filename, -1)
            wcs = WCS(header=header)

    if (not width or not height) and header:
        width = header['NAXIS1']
        height = header['NAXIS2']

    [ra1],[dec1] = wcs.all_pix2world([0], [0], 1)
    [ra0],[dec0] = wcs.all_pix2world([width/2], [height/2], 1)

    sr = coords.sphdist(ra0, dec0, ra1, dec1)[0]

    return ra0, dec0, sr

def blind_match_objects(obj, order=4, extra="", verbose=False, fix=True, sn=20):
    dir = tempfile.mkdtemp(prefix='astrometry')
    wcs = None
    binname = None
    ext = 0

    for path in ['.', '/usr/local', '/opt/local']:
        if os.path.isfile(posixpath.join(path, 'astrometry', 'bin', 'solve-field')):
            binname = posixpath.join(path, 'astrometry', 'bin', 'solve-field')
            break

    if binname:
        idx = obj['magerr']<1/sn
        columns = [fits.Column(name='XIMAGE', format='1D', array=obj['x'][idx]+1),
                   fits.Column(name='YIMAGE', format='1D', array=obj['y'][idx]+1),
                   fits.Column(name='FLUX', format='1D', array=obj['flux'][idx])]
        tbhdu = fits.BinTableHDU.from_columns(columns)
        filename = posixpath.join(dir, 'list.fits')
        tbhdu.writeto(filename, overwrite=True)
        extra += " --x-column XIMAGE --y-column YIMAGE --sort-column FLUX --width %d --height %d" % (np.ceil(max(obj['x']+1)), np.ceil(max(obj['y']+1)))

        wcsname = posixpath.split(filename)[-1]
        tmpname = posixpath.join(dir, posixpath.splitext(wcsname)[0] + '.tmp')
        wcsname = posixpath.join(dir, posixpath.splitext(wcsname)[0] + '.wcs')

        if verbose:
            print("%s -D %s --no-verify --overwrite --no-plots -T %s %s" % (binname, dir, extra, filename))

        os.system("%s -D %s --no-verify --overwrite --no-plots -T %s %s" % (binname, dir, extra, filename))

        if order:
            order_str = "-t %d" % order
        else:
            order_str = "-T"

        if os.path.isfile(wcsname):
            shutil.move(wcsname, tmpname)
            os.system("%s -D %s --overwrite --no-plots %s %s --verify %s %s" % (binname, dir, order_str, extra, tmpname, filename))

            if os.path.isfile(wcsname):
                header = fits.getheader(wcsname)
                wcs = WCS(header)

                if fix and wcs:
                    obj['ra'],obj['dec'] = wcs.all_pix2world(obj['x'], obj['y'], 0)
    else:
        print("Astrometry.Net binary not found")

    #print order
    shutil.rmtree(dir)

    return wcs

def radectoxyz(ra_deg, dec_deg):
    ra  = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    xyz = np.array((np.cos(dec)*np.cos(ra),
                    np.cos(dec)*np.sin(ra),
                    np.sin(dec)))

    return xyz

def xyztoradec(xyz):
    ra = np.arctan2(xyz[1], xyz[0])
    ra += 2*np.pi * (ra < 0)
    dec = np.arcsin(xyz[2] / np.linalg.norm(xyz, axis=0))

    return (np.rad2deg(ra), np.rad2deg(dec))

def make_series(mul=1.0, x=1.0, y=1.0, order=1, sum=False, legendre=False, zero=True):
    if zero:
        res = [np.ones_like(x)*mul]
    else:
        res = []

    for i in xrange(1,order+1):
        maxr = i+1
        if legendre:
            maxr = order+1

        for j in xrange(maxr):
            #print i, '-', i - j, j
            if legendre:
                res.append(mul * leg(i)(x) * leg(j)(y))
            else:
                res.append(mul * x**(i-j) * y**j)
    if sum:
        return np.sum(res, axis=0)
    else:
        return res

def make_kernel(r0=1.0, ext=1.0):
    x,y = np.mgrid[np.floor(-ext*r0):np.ceil(ext*r0+1), np.floor(-ext*r0):np.ceil(ext*r0+1)]
    r = np.hypot(x,y)
    image = np.exp(-r**2/2/r0**2)

    return image

def fit_2d(x, y, z, x0, y0, weights=None, order=4, niter=3):
    # Fit
    sx = (x - np.mean(x))/np.max(x)
    sy = (y - np.mean(y))/np.max(y)

    X = make_series(1.0, sx, sy, order=order)
    X = np.vstack(X).T
    Y = z

    idx = np.isfinite(Y)

    for i in range(niter):
        if weights is not None:
            C = sm.WLS(Y[idx], X[idx], weights=weights[idx]).fit()
        else:
            C = sm.RLM(Y[idx], X[idx]).fit()

        YY = np.sum(X*C.params, axis=1)
        idx = np.abs(Y-YY - np.median((Y-YY)[idx])) < 3.0*mad_std((Y-YY)[idx])

    # Predict
    sx = (x0 - np.mean(x))/np.max(x)
    sy = (y0 - np.mean(y))/np.max(y)

    X = make_series(1.0, sx, sy, order=order)
    X = np.vstack(X).T

    return np.sum(X*C.params, axis=1)

def get_objects_sep(image, header=None, mask=None, thresh=4.0, aper=3.0, bkgann=None, r0=0.5, gain=1, edge=0, minnthresh=2, minarea=5, relfluxradius=2.0, wcs=None, use_fwhm=False, use_mask_bg=False, use_mask_large=False, subtract_bg=True, npix_large=300, sn=10.0, verbose=True, get_fwhm=False, **kwargs):
    if r0 > 0.0:
        kernel = make_kernel(r0)
    else:
        kernel = None

    if verbose:
        print("Preparing background mask")

    if mask is None:
        mask = np.zeros_like(image, dtype=np.bool)

    mask_bg = np.zeros_like(mask)
    mask_segm = np.zeros_like(mask)

    if use_mask_bg:
        # Simple heuristics to mask regions with rapidly varying background
        if verbose:
            print("Masking rapidly changing background")

        bg2 = sep.Background(image, mask=mask|mask_bg, bw=64, bh=64)

        for _ in xrange(3):
            bg1 = sep.Background(image, mask=mask|mask_bg, bw=256, bh=256)

            ibg = bg2.back() - bg1.back()

            tmp = np.abs(ibg - np.median(ibg)) > 5.0*mad_std(ibg)
            mask_bg |= dilate(tmp, np.ones([100, 100]))

        mask_bg = dilate(tmp, np.ones([100, 100]))

    if verbose:
        print("Building background map")

    bg = sep.Background(image, mask=mask|mask_bg, bw=64, bh=64)
    if subtract_bg:
        image1 = image - bg.back()
    else:
        image1 = image.copy()

    sep.set_extract_pixstack(image.shape[0]*image.shape[1])

    if use_mask_large:
        # Mask regions around huge objects as they are most probably corrupted by saturation and blooming
        if verbose:
            print("Extracting initial objects")

        obj0,segm = sep.extract(image1, err=bg.rms(), thresh=thresh, minarea=minarea, mask=mask|mask_bg, filter_kernel=kernel, segmentation_map=True, **kwargs)

        if verbose:
            print("Dilating large objects")

        mask_segm = np.isin(segm, [_+1 for _,npix in enumerate(obj0['npix']) if npix > npix_large])
        mask_segm = dilate(mask_segm, np.ones([50, 50]))

    if verbose:
        print("Extracting final objects")

    obj0 = sep.extract(image1, err=bg.rms(), thresh=thresh, minarea=minarea, mask=mask|mask_bg|mask_segm, filter_kernel=kernel, **kwargs)

    if use_fwhm:
        # Estimate FHWM and use it to get optimal aperture size
        idx = obj0['flag'] == 0
        fwhm = 2.0*np.sqrt(np.hypot(obj0['a'][idx], obj0['b'][idx])*np.log(2))
        fwhm = 2.0*sep.flux_radius(image1, obj0['x'][idx], obj0['y'][idx], relfluxradius*fwhm*np.ones_like(obj0['x'][idx]), 0.5, mask=mask)[0]
        fwhm = np.median(fwhm)

        aper = max(1.5*fwhm, aper)

        if verbose:
            print("FWHM = %.2g, aperture = %.2g" % (fwhm, aper))

    # Windowed positional parameters are often biased in crowded fields, let's avoid them for now
    # xwin,ywin,flag = sep.winpos(image1, obj0['x'], obj0['y'], 0.5, mask=mask)
    xwin,ywin = obj0['x'], obj0['y']

    # Filter out objects too close to frame edges
    idx = (np.round(xwin) > edge) & (np.round(ywin) > edge) & (np.round(xwin) < image.shape[1]-edge) & (np.round(ywin) < image.shape[0]-edge) # & (obj0['flag'] == 0)

    if minnthresh:
        idx &= (obj0['tnpix'] >= minnthresh)

    if verbose:
        print("Measuring final objects")

    flux,fluxerr,flag = sep.sum_circle(image1, xwin[idx], ywin[idx], aper, err=bg.rms(), gain=gain, mask=mask|mask_bg|mask_segm, bkgann=bkgann)
    # For debug purposes, let's make also the same aperture photometry on the background map
    bgflux,bgfluxerr,bgflag = sep.sum_circle(bg.back(), xwin[idx], ywin[idx], aper, err=bg.rms(), gain=gain, mask=mask|mask_bg|mask_segm)

    bg = bgflux/np.pi/aper**2

    # Fluxes to magnitudes
    mag,magerr = np.zeros_like(flux), np.zeros_like(flux)
    mag[flux>0] = -2.5*np.log10(flux[flux>0])
    # magerr[flux>0] = 2.5*np.log10(1.0 + fluxerr[flux>0]/flux[flux>0])
    magerr[flux>0] = 2.5/np.log(10)*fluxerr[flux>0]/flux[flux>0]

    # better FWHM estimation - FWHM=HFD for Gaussian
    if get_fwhm:
        fwhm = 2.0*sep.flux_radius(image1, xwin[idx], ywin[idx], relfluxradius*aper*np.ones_like(xwin[idx]), 0.5, mask=mask)[0]
    else:
        fwhm = np.zeros_like(xwin[idx])

    flag |= obj0['flag'][idx]

    # Quality cuts
    fidx = (flux > 0) & (magerr < 1.0/sn)

    if wcs is None and header is not None:
        # If header is provided, we may build WCS from it
        wcs = WCS(header)

    if wcs is not None:
        # If WCS is provided we may convert x,y to ra,dec
        ra,dec = wcs.all_pix2world(obj0['x'][idx], obj0['y'][idx], 0)
    else:
        ra,dec = np.zeros_like(obj0['x'][idx]),np.zeros_like(obj0['y'][idx])

    if verbose:
        print("All done")

    return {'x':xwin[idx][fidx], 'y':ywin[idx][fidx], 'flux':flux[fidx], 'fluxerr':fluxerr[fidx], 'mag':mag[fidx], 'magerr':magerr[fidx], 'flags':obj0['flag'][idx][fidx]|flag[fidx], 'ra':ra[fidx], 'dec':dec[fidx], 'bg':bg[fidx], 'fwhm':fwhm[fidx], 'aper':aper, 'bkgann':bkgann, 'a':obj0['a'][idx][fidx], 'b':obj0['b'][idx][fidx], 'theta':obj0['theta'][idx][fidx]}

def get_objects_sextractor(image, header=None, mask=None, thresh=2.0, aper=3.0, r0=0.5, bkgann=None, gain=1, edge=0, minarea=5, wcs=None, sn=3.0, verbose=False, checkimages=[], extra_params=[], extra_opts={}, _workdir=None, _tmpdir=None):
    # Find the binary
    binname = None
    for path in ['.', '/usr/bin', '/usr/local/bin', '/opt/local/bin']:
        for exe in ['sex', 'sextractor', 'source-extractor']:
            if os.path.isfile(posixpath.join(path, exe)):
                binname = posixpath.join(path, exe)
                break

    if binname is None:
        if verbose:
            print("Can't find SExtractor binary")
        return None

    workdir = _workdir if _workdir is not None else tempfile.mkdtemp(prefix='sex', dir=_tmpdir)
    obj = None

    # Prepare
    imagename = posixpath.join(workdir, 'image.fits')
    fits.writeto(imagename, image, header, overwrite=True)

    opts = {
        'VERBOSE_TYPE': 'QUIET',
        'DETECT_MINAREA': minarea,
        'GAIN': gain,
        'DETECT_THRESH': thresh,
        'WEIGHT_TYPE': 'BACKGROUND',
        'MASK_TYPE': 'NONE', # both 'CORRECT' and 'BLANK' seem to cause systematics
        # 'MASK_TYPE': 'CORRECT', # both 'CORRECT' and 'BLANK' seem to cause systematics
    }

    if mask is None:
        mask = np.zeros_like(image, dtype=np.bool)

    flagsname = posixpath.join(workdir, 'flags.fits')
    fits.writeto(flagsname, mask.astype(np.int16), overwrite=True)
    opts['FLAG_IMAGE'] = flagsname

    if np.isscalar(aper):
        opts['PHOT_APERTURES'] = aper*2 # SExtractor expects diameters, not radii
        size = ''
    else:
        opts['PHOT_APERTURES'] = ','.join([str(_*2) for _ in aper])
        size = '[%d]' % len(aper)

    checknames = [posixpath.join(workdir, _.replace('-', 'M_') + '.fits') for _ in checkimages]
    if checkimages:
        opts['CHECKIMAGE_TYPE'] = ','.join(checkimages)
        opts['CHECKIMAGE_NAME'] = ','.join(checknames)

    params = ['MAG_APER'+size, 'MAGERR_APER'+size, 'FLUX_APER'+size, 'FLUXERR_APER'+size, 'X_IMAGE', 'Y_IMAGE', 'ERRX2_IMAGE', 'ERRY2_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'FLUX_RADIUS', 'FWHM_IMAGE', 'FLAGS', 'IMAFLAGS_ISO', 'BACKGROUND']
    params += extra_params
    paramname = posixpath.join(workdir, 'cfg.param')
    open(paramname, 'w').write("\n".join(params))
    opts['PARAMETERS_NAME'] = paramname

    catname = posixpath.join(workdir, 'out.cat')
    opts['CATALOG_NAME'] = catname
    opts['CATALOG_TYPE'] = 'FITS_LDAC'

    if not r0:
        opts['FILTER'] = 'N'
    else:
        kernel = make_kernel(r0, ext=1.0)
        kernelname = posixpath.join(workdir, 'kernel.txt')
        np.savetxt(kernelname, kernel/np.sum(kernel), fmt=b'%.6f', header='CONV NORM', comments='')
        opts['FILTER'] = 'Y'
        opts['FILTER_NAME'] = kernelname

    opts.update(extra_opts)

    # Build the command line
    # FIXME: quote strings!
    cmd = binname + ' ' + imagename + ' ' + ' '.join(['-%s %s' % (_,opts[_]) for _ in opts.keys()])
    if not verbose:
        cmd += ' > /dev/null 2>/dev/null'
    if verbose:
        print(cmd)

    # Run the command!

    res = os.system(cmd)

    if res == 0:
        data = fits.getdata(catname, -1)

        idx = (data['X_IMAGE'] > edge) & (data['X_IMAGE'] < image.shape[1] - edge)
        idx &= (data['Y_IMAGE'] > edge) & (data['Y_IMAGE'] < image.shape[0] - edge)

        if np.isscalar(aper):
            idx &= data['MAGERR_APER'] < 1.0/sn
            idx &= data['FLUX_APER'] > 0
        else:
            idx &= np.all(data['MAGERR_APER'] < 1.0/sn, axis=1)
            idx &= np.all(data['FLUX_APER'] > 0, axis=1)

        data = data[idx]

        if wcs is None and header is not None:
            wcs = WCS(header)

        if wcs is not None:
            ra,dec = wcs.all_pix2world(data['X_IMAGE'], data['Y_IMAGE'], 1)
        else:
            ra,dec = np.zeros_like(data['X_IMAGE']),np.zeros_like(data['Y_IMAGE'])

        data['FLAGS'][data['IMAFLAGS_ISO'] > 0] |= 256

        obj = {
            'x': data['X_IMAGE']-1, 'y': data['Y_IMAGE']-1,
            'xerr': np.sqrt(data['ERRX2_IMAGE']), 'yerr': np.sqrt(data['ERRY2_IMAGE']),
            'flux': data['FLUX_APER'], 'fluxerr': data['FLUXERR_APER'],
            'mag': data['MAG_APER'], 'magerr': data['MAGERR_APER'],
            'flags': data['FLAGS'], 'ra':ra, 'dec': dec,
            'bg': data['BACKGROUND'], 'fwhm': data['FWHM_IMAGE'],
            'aper': aper, 'a': data['A_IMAGE'], 'b': data['B_IMAGE'], 'theta': data['THETA_IMAGE'],
        }

        for _ in extra_params:
            obj[_] = data[_]
    else:
        if verbose:
            print("Error", res, "running SExtractor")

    result = obj

    if checkimages:
        result = [result]

        for name in checknames:
            result.append(fits.getdata(name))

    if _workdir is None:
        shutil.rmtree(workdir)

    return result

def match_objects(obj, cat, sr, fname='V', order=4, bg_order=None, cmag=None, color_order=None, thresh=5.0, clim=None, sn=10, mag_idx=0):
    x0,y0,width,height = np.mean(obj['x']), np.mean(obj['y']), np.max(obj['x'])-np.min(obj['x']), np.max(obj['y'])-np.min(obj['y'])

    # Match stars
    h = htm.HTM(10)
    m = h.match(obj['ra'],obj['dec'], cat['ra'],cat['dec'], sr, maxmatch=0)
    oidx = m[0]
    cidx = m[1]
    dist = m[2]

    if fname == 'B':
        cmag,cmagerr = cat['B'], cat['Berr']
    elif fname == 'V':
        cmag,cmagerr = cat['V'], cat['Verr']
    elif fname == 'R' or fname == 'N':
        cmag,cmagerr = cat['R'], cat['Rerr']
    elif fname == 'I':
        cmag,cmagerr = cat['I'], cat['Ierr']
    elif fname == 'g':
        cmag,cmagerr = cat['g'], cat['gerr']
    elif fname == 'r':
        cmag,cmagerr = cat['r'], cat['rerr']
    elif fname == 'i':
        cmag,cmagerr = cat['i'], cat['ierr']
    elif fname == 'z':
        cmag,cmagerr = cat['z'], cat['zerr']

    color = cat['B'] - cat['V']

    if cmag is not None:
        cmag = cmag[cidx]
        cmagerr = cmagerr[cidx]
        color = color[cidx]
    else:
        print('Unsupported filter:', fname)
        return None

    # Optimal fit limit, as mean mag where S/N = 10
    idx = (1.0/obj['magerr'][oidx] > 5) & (1.0/obj['magerr'][oidx] < 15)
    if np.sum(idx) > 10:
        X = make_series(1.0, obj['x'][oidx], obj['y'][oidx], order=order)
        X = np.vstack(X).T
        Y = cmag
        weights = 1/obj['magerr'][oidx]**2
        Cmag_lim = sm.WLS(Y[idx], X[idx], weights=weights[idx]).fit()
        cmag_lim = np.sum(X*Cmag_lim.params, axis=1)
    else:
        # print('Not enough matches with SN~10:', np.sum(idx))
        Cmag_lim = None
        cmag_lim = 99.0*np.ones_like(cmag)

    # User provided upper fit limit
    if clim is not None:
        idx = cmag < clim

        oidx,cidx,dist,cmag,cmagerr,cmag_lim,color = [_[idx] for _ in [oidx,cidx,dist,cmag,cmagerr,cmag_lim,color]]

    x,y = obj['x'][oidx],obj['y'][oidx]
    oflags = obj['flags'][oidx]
    omag,omagerr = obj['mag'][oidx],obj['magerr'][oidx]
    if len(obj['mag'].shape) > 1:
        omag,omagerr = omag[:,mag_idx],omagerr[:,mag_idx]

    x = (x - x0)*2/width
    y = (y - y0)*2/height

    tmagerr = np.hypot(cmagerr, omagerr)
    tmagerr = np.hypot(tmagerr, 0.02)

    delta_mag = cmag - omag
    weights = 1.0/tmagerr**2

    X = make_series(1.0, x, y, order=order)
    if bg_order is not None and bg_order >= 0:
        X += make_series(-2.5/np.log(10)/10**(-0.4*omag), x, y, order=bg_order)

    if color_order is not None and color_order >= 0:
        X += make_series(color, x, y, order=color_order)

    X = np.vstack(X).T
    Y = delta_mag

    # X /= tmagerr[:, None]
    # Y /= tmagerr

    idx0 = (oflags == 0) & (cmag < cmag_lim)
    if sn and sn>0:
        idx0 &= (omagerr < 1/sn)

    idx = idx0.copy()

    # Do the fitting
    for iter in range(3):
        if len(X[idx]) < 3:
            print("Fit failed - %d objects" % len(X[idx]))
            return None

        C = sm.WLS(Y[idx], X[idx], weights=weights[idx]).fit()
        # C = sm.RLM(Y[idx], X[idx]).fit()

        YY = np.sum(X*C.params, axis=1)
        idx = idx0.copy()
        if thresh and thresh > 0:
            idx &= (np.abs((Y-YY)/tmagerr) < thresh)

    # Apply the fit for all objects
    x = (obj['x'] - x0)*2/width
    y = (obj['y'] - y0)*2/height

    if len(obj['mag'].shape) > 1:
        mag = obj['mag'][:,mag_idx]
        magerr = obj['magerr'][:,mag_idx]
    else:
        mag = obj['mag']
        magerr = obj['magerr']

    # Magnitude estimation for every position
    X = make_series(1.0, x, y, order=order)
    if bg_order is not None and bg_order >= 0:
        X += make_series(-2.5/np.log(10)/10**(-0.4*obj['mag']), x, y, order=bg_order)

    X = np.vstack(X).T
    YY1 = np.sum(X*C.params[0:X.shape[1]], axis=1)

    if len(obj['mag'].shape) > 1:
        mag = obj['mag'][:,mag_idx] + YY1
    else:
        mag = obj['mag'] + YY1

    # Additive flux component for every position
    if bg_order is not None and bg_order >= 0:
        Xbg = make_series(1, x, y, order=bg_order)
        Xbg = np.vstack(Xbg).T
        delta_flux = np.sum(Xbg*C.params[X.shape[1]-Xbg.shape[1]:X.shape[1]], axis=1)
    else:
        delta_flux = np.zeros_like(mag)

    # Color term for every position
    if color_order is not None and color_order >= 0:
        Xc = make_series(1, x, y, order=color_order)
        Xc = np.vstack(Xc).T

        Cbv = np.sum(Xc*C.params[X.shape[1]:], axis=1)
    else:
        Cbv = np.zeros_like(mag)

    # Approx magnitude limit (S/N=10) at every position
    if Cmag_lim is not None:
        Xlim = make_series(1.0, x, y, order=order)
        Xlim = np.vstack(Xlim).T
        cmag_lim = np.sum(Xlim*Cmag_lim.params, axis=1)
    else:
        cmag_lim = 99.0*np.ones_like(x)

    # Simple analysis of proximity to "good" points
    mx,my = obj['x'][oidx][idx], obj['y'][oidx][idx]
    kdo = cKDTree(np.array([obj['x'], obj['y']]).T)
    kdm = cKDTree(np.array([mx, my]).T)

    mr0 = np.sqrt(width*height/np.sum(idx))
    m = kdm.query_ball_tree(kdm, 5.0*mr0)

    dists = []
    for i,ii in enumerate(m):
        if len(ii) > 1:
            d1 = [np.hypot(mx[i]-mx[_], my[i]-my[_]) for _ in ii]
            d1 = np.sort(d1)

            dists.append(d1[1])
    mr1 = np.median(dists)

    m = kdo.query_ball_tree(kdm, 5.0*mr1)
    midx = np.array([len(_)>1 for _ in m])

    return {
        # Matching indices and a distance in degrees
        'cidx':cidx, 'oidx':oidx, 'dist':dist,
        # Pixel and sky coordinates of matched stars, as well as their flags
        'x':obj['x'][oidx], 'y':obj['y'][oidx], 'flags':oflags,
        'ra':obj['ra'][oidx], 'dec':obj['dec'][oidx],
        # All catalogue magnitudes of matched stars
        'cB':cat['B'][cidx], 'cV':cat['V'][cidx], 'cR':cat['R'][cidx], 'cI':cat['I'][cidx],
        'cBerr':cat['Berr'][cidx], 'cVerr':cat['Verr'][cidx], 'cRerr':cat['Rerr'][cidx], 'cIerr':cat['Ierr'][cidx],
        # Catalogue magnitudes of matched stars in proper filter
        'cmag':cmag, 'cmagerr':cmagerr, 'tmagerr':tmagerr,
        # Model zero point for all objects, and their corrected magnitudes
        'mag0':YY1, 'mag':mag,
        'Cbv':Cbv, 'delta_flux': delta_flux,
        'mag_idx': mag_idx,
        'Y':Y, 'YY':YY,
        'cmag_lim':cmag_lim,
        # Subset of matched stars used in the fits
        'idx':idx,
        # Subset of all objects at 'good' distances from matched ones
        'midx':midx, 'mr0':mr0, 'mr1':mr1}

def fix_wcs(obj, cat, sr, header=None, use_header_wcs=False, maxmatch=1, order=6, fix=True):
    '''Get a refined WCS solution based on cross-matching of objects with catalogue on the sphere.
    Uses external 'fit-wcs' binary from Astrometry.Net suite'''

    if header is not None:
        width,height = header['NAXIS1'],header['NAXIS2']

        if use_header_wcs:
            wcs = WCS(header)

            if wcs:
                obj['ra'],obj['dec'] = wcs.all_pix2world(obj['x'], obj['y'], 0)
    else:
        width,height = int(np.max(obj['x'])),int(np.max(obj['y']))

    h = htm.HTM(10)
    oidx,cidx,dist = h.match(obj['ra'], obj['dec'], cat['ra'], cat['dec'], sr, maxmatch=1)

    dir = tempfile.mkdtemp(prefix='astrometry')
    wcs = None
    binname = None

    for path in ['.', '/usr/local', '/opt/local']:
        if os.path.isfile(posixpath.join(path, 'astrometry', 'bin', 'fit-wcs')):
            binname = posixpath.join(path, 'astrometry', 'bin', 'fit-wcs')
            break

    if binname:
        columns = [fits.Column(name='FIELD_X', format='1D', array=obj['x'][oidx] + 1),
                   fits.Column(name='FIELD_Y', format='1D', array=obj['y'][oidx] + 1),
                   fits.Column(name='INDEX_RA', format='1D', array=cat['ra'][cidx]),
                   fits.Column(name='INDEX_DEC', format='1D', array=cat['dec'][cidx])]
        tbhdu = fits.BinTableHDU.from_columns(columns)
        filename = posixpath.join(dir, 'list.fits')
        wcsname = posixpath.join(dir, 'list.wcs')

        tbhdu.writeto(filename, overwrite=True)

        os.system("%s -c %s -o %s -W %d -H %d -C -s %d" % (binname, filename, wcsname, width, height, order))

        if os.path.isfile(wcsname):
            header = fits.getheader(wcsname)
            wcs = WCS(header)

            if fix and wcs:
                obj['ra'],obj['dec'] = wcs.all_pix2world(obj['x'], obj['y'], 0)

    else:
        print("Astrometry.Net binary not found")

    #print order
    shutil.rmtree(dir)

    return wcs

def fix_distortion(obj, cat, header=None, wcs=None, width=None, height=None, dr=3.0):
    '''Refines object sky coordinates based on cross-match with the catalogue in pixel space.
    Does not use any external tools but can't provide a fixed WCS solution'''
    if header:
        wcs = WCS(header)
        width = header['NAXIS1']
        height = header['NAXIS2']

    if not wcs or not width or not height:
        print("Nothing to fix")
        return

    kdo = cKDTree(np.array([obj['x'], obj['y']]).T)
    xc,yc = wcs.all_world2pix(cat['ra'], cat['dec'], 0)
    kdc = cKDTree(np.array([xc,yc]).T)

    m = kdo.query_ball_tree(kdc, dr)
    nm = np.array([len(_) for _ in m])

    # Distortions
    dx = np.array([obj['x'][_] - xc[m[_][0]] if nm[_] == 1 else 0 for _ in xrange(len(m))])
    dy = np.array([obj['y'][_] - yc[m[_][0]] if nm[_] == 1 else 0 for _ in xrange(len(m))])

    # Normalized coordinates
    x = (obj['x'] - width/2)*2.0/width
    y = (obj['y'] - height/2)*2.0/height

    X = make_series(1.0, x, y, order=6)

    X = np.vstack(X).T
    Yx = dx
    Yy = dy

    # Use only unique matches
    idx = (nm == 1) # & (distm < 2.0*np.std(distm))
    Cx = sm.WLS(Yx[idx], X[idx]).fit()
    Cy = sm.WLS(Yy[idx], X[idx]).fit()

    YYx = np.sum(X*Cx.params, axis=1)
    YYy = np.sum(X*Cy.params, axis=1)

    obj['ra'],obj['dec'] = wcs.all_pix2world(obj['x']-YYx, obj['y']-YYy, 0)

def save_objects(filename, obj, header=None, wcs=None):
    # Create directory hierarchy as necessary
    dirname = posixpath.split(filename)[0]

    try:
        os.makedirs(dirname)
    except:
        pass

    if obj is None:
        # Create dummy empty file
        with open(filename, 'w') as f:
            pass
        return

    header1 = fits.Header() if header is None else header.copy()
    # Clean header from unnecessary stuff
    header1.remove('HISTORY', ignore_missing=True, remove_all=True)
    header1.remove('COMMENT', ignore_missing=True, remove_all=True)
    for _ in [_ for _ in header1.keys() if _[0] == '_']:
        header1.remove(_, ignore_missing=True, remove_all=True)

    if wcs is not None:
        header1.update(wcs.to_header(relax=True))

    columns = []

    meta_len = 0

    for name in obj.keys():
        # Scalar type or array?
        is_scalar = np.isscalar(obj[name])
        if type(obj[name]) is datetime.datetime:
            is_scalar = True

        if not np.isscalar(obj[name]) and hasattr(obj[name], '__len__') and len(obj[name]) == len(obj['x']):
            # Add new data column
            size = 1 if len(obj[name].shape) == 1 else obj[name].shape[1]
            fmt = 'L' if obj[name].dtype == np.bool else 'D'
            col = fits.Column(name=name, array=obj[name], format='%d%s' % (size,fmt))
            columns.append(col)

        else:
            if np.isscalar(obj[name]):
                value = obj[name]
                vtype = 'scalar'
            elif type(obj[name]) is datetime.datetime:
                value = obj[name].strftime('%Y-%m-%dT%H:%M:%S.%f')
                vtype = 'datetime'
            else:
                value = json.dumps(obj[name])
                vtype = 'json'

            header1['METAN%d' % meta_len] = name
            header1['METAV%d' % meta_len] = value
            header1['METAT%d' % meta_len] = vtype

            meta_len += 1

    if meta_len:
        header1['METALEN'] = meta_len
        header1.comments['METALEN'] = 'Number of metadata entries'

    hdu1 = fits.BinTableHDU.from_columns(columns, header=header1)
    hdu1.writeto(filename, overwrite=True)

def load_objects(filename, get_header=False):
    obj,header = None,None

    try:
        data,header = fits.getdata(filename, -1), fits.getheader(filename, -1)

        obj = {}

        for col in data.columns:
            obj[col.name] = data[col.name]

        for i in range(header.get('METALEN', 0)):
            name = header.get('METAN%d' % i)
            value = header.get('METAV%d' % i)
            vtype = header.get('METAT%d' % i)

            if vtype == 'datetime':
                value = datetime.datetime.strptime(value, '%Y-%m-%dT%H:%M:%S.%f')
            elif vtype == 'json':
                value = json.loads(value)

            if name:
                obj[name] = value

    except IOError:
        pass
    except:
        import traceback
        traceback.print_exc()

    if get_header:
        return obj, header
    else:
        return obj
