from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from esutil import coords, htm
import statsmodels.api as sm
from scipy.spatial import cKDTree

from survey import make_series

class Match:
    def __init__(self, width=None, height=None):
        self.width,self.height = width, height
        self.x0,self.y0 = width/2 if width else None, height/2 if height else None
        self.C = None

    def match(self, obj=None, cat=None, sr=5./3600, verbose=False, predict=True,
              ra=None, dec=None, x=None, y=None, mag=None, magerr=None, flags=None,
              filter_name='V', order=4, bg_order=None, color_order=None,
              hard_mag_limit=99, mag_id=0, magerr0=0.02, sn=None, thresh=5.0):

        """Match a set of points with catalogue"""

        self.success = False
        self.ngoodstars = 0

        self.order = order
        self.bg_order = bg_order
        self.color_order = color_order
        self.mag_id = mag_id

        self.filter_name = filter_name
        if filter_name in ['B', 'V', 'R', 'I', 'g', 'r', 'i', 'z']:
            # Generic names
            cmag,cmagerr = cat[filter_name], cat[filter_name + 'err']
            self.cat_filter_name = filter_name
        elif filter_name == 'Clear':
            # Mini-MegaTORTORA
            cmag,cmagerr = cat['V'], cat['Verr']
            self.cat_filter_name = 'V'
        elif filter_name == 'N':
            # FRAMs
            cmag,cmagerr = cat['R'], cat['Rerr']
            self.cat_filter_name = 'R'
        else:
            if verbose:
                print('Unsupported filter name: %s' % filter_name)
            return False

        # TODO: make it configurable?..
        color = cat['B'] - cat['V']
        self.cat_color_name = 'B - V'

        # Objects to match
        if obj is not None:
            ra = obj['ra']
            dec = obj['dec']
            x = obj['x']
            y = obj['y']
            mag = obj['mag']
            magerr = obj['magerr']
            flags = obj['flags']

        else:
            if ra is None or dec is None or x is None or y is None or mag is None:
                raise ValueError('Data for matching are missing')

            if magerr is None:
                magerr = np.ones_like(mag)*np.std(mag)

            if flags is None:
                flags = np.zeros_like(ra, dtype=np.int)

        if self.width is None or self.height is None:
            self.x0,self.y0,self.width,self.height = np.mean(x), np.mean(y), np.max(x)-np.min(x), np.max(y) - np.min(y)

        # Match stars
        h = htm.HTM(10)
        oidx,cidx,dist = h.match(ra, dec, cat['ra'],cat['dec'], sr, maxmatch=0)

        if verbose:
            print(len(oidx), 'matches between', len(ra), 'objects and', len(cat['ra']), 'stars, sr = %.1f arcsec' % (3600.0*sr))

        self.oidx,self.cidx,self.dist = oidx, cidx, dist

        self.cmag = cmag[cidx]
        self.cmagerr = cmagerr[cidx]
        self.color = color[cidx]

        self.ox,self.oy = x[oidx], y[oidx]
        self.oflags = flags[oidx]
        self.omag,self.omagerr = mag[oidx], magerr[oidx]
        if len(self.omag.shape) > 1:
            # If we are given a multi-aperture magnitude column
            self.omag,self.omagerr = self.omag[:, mag_id], self.omagerr[:, mag_id]

        # Scaled spatial coordinates for fitting
        sx = (self.ox - self.x0)*2/self.width
        sy = (self.oy - self.y0)*2/self.height

        # Optimal magnitude cutoff for fitting, as a mean mag where S/N = 10
        idx = (1.0/self.omagerr > 5) & (1.0/self.omagerr < 15)
        if np.sum(idx) > 10:
            X = make_series(1.0, sx, sy, order=order)
            X = np.vstack(X).T
            Y = self.cmag

            self.C_mag_limit = sm.RLM(Y[idx], X[idx]).fit()
            mag_limit = np.sum(X*self.C_mag_limit.params, axis=1)
        else:
            if verbose:
                print('Not enough matches with SN~10:', np.sum(idx))
            self.C_mag_limit = None
            mag_limit = 99.0*np.ones_like(cmag)

        self.zero = self.cmag - self.omag # We will build a model for this variable

        self.zeroerr = np.hypot(self.omagerr, self.cmagerr)
        self.zeroerr = np.hypot(self.zeroerr, magerr0)

        self.weights = 1.0/self.zeroerr**2

        X = make_series(1.0, sx, sy, order=self.order)
        if self.bg_order is not None:
            X += make_series(-2.5/np.log(10)/10**(-0.4*self.omag), sx, sy, order=self.bg_order)

        if self.color_order is not None:
            X += make_series(self.color, sx, sy, order=self.color_order)

        X = np.vstack(X).T

        self.idx0 = (self.oflags == 0) & (self.cmag < hard_mag_limit) & (self.cmag < mag_limit)

        if sn is not None:
            self.idx0 &= (self.omagerr < 1.0/sn)

        # Actual fitting
        self.idx = self.idx0.copy()

        for iter in range(3):
            if np.sum(self.idx) < 3:
                if verbose:
                    print("Fit failed - %d objects" % np.sum(self.idx))
                return False

            self.C = sm.WLS(self.zero[self.idx], X[self.idx], weights=self.weights[self.idx]).fit()

            self.zero_model = np.sum(X*self.C.params, axis=1)

            self.idx = self.idx0.copy()
            if thresh and thresh > 0:
                self.idx &= (np.abs((self.zero - self.zero_model)/self.zeroerr) < thresh)

        self.std = np.std((self.zero - self.zero_model)[self.idx])
        self.ngoodstars = np.sum(self.idx)
        self.success = True

        if verbose:
            print('Fit finished:', self.ngoodstars, 'stars, rms', self.std)

        if predict:
            self.predict(obj=obj, x=x, y=y, mag=mag, magerr=magerr, mag_id=mag_id, verbose=verbose)

        return True

    def predict(self, obj=None, x=None, y=None, mag=None, magerr=None, mag_id=0, verbose=False):
        """Compute the magnitudes in catalogue system based on a model already built"""

        if self.C is None:
            return

        if obj is not None:
            x = obj['x']
            y = obj['y']
            mag = obj['mag']
            magerr = obj['magerr']
        elif x is None or y is None or mag is None or magerr is None:
            raise ValueError('Missing data for objects')

        self.x,self.y = x, y

        # Scaled spatial coordinates
        sx = (x - self.x0)*2/self.width
        sy = (y - self.y0)*2/self.height

        if len(mag.shape) > 1:
            mag = mag[:,mag_id]
            magerr = magerr[:,mag_id]

        # Magnitude estimation for every position
        X = make_series(1.0, sx, sy, order=self.order)
        if self.bg_order is not None:
            X += make_series(-2.5/np.log(10)/10**(-0.4*mag), sx, sy, order=self.bg_order)

        X = np.vstack(X).T
        self.mag0 = np.sum(X*self.C.params[0:X.shape[1]], axis=1)

        self.mag = mag + self.mag0
        self.magerr = magerr

        # Additive flux component for every position
        if self.bg_order is not None:
            Xbg = make_series(1.0, sx, sy, order=self.bg_order)
            Xbg = np.vstack(Xbg).T
            self.delta_flux = np.sum(Xbg*self.C.params[X.shape[1]-Xbg.shape[1]:X.shape[1]], axis=1)
        else:
            self.delta_flux = np.zeros_like(mag)

        # Color term for every position
        if self.color_order is not None:
            Xc = make_series(1.0, sx, sy, order=self.color_order)
            Xc = np.vstack(Xc).T

            self.color_term = np.sum(Xc*self.C.params[X.shape[1]:], axis=1)
        else:
            self.color_term = np.zeros_like(mag)

        # Approx magnitude limit (S/N=10) at every position
        if self.C_mag_limit is not None:
            Xlim = make_series(1.0, sx, sy, order=self.order)
            Xlim = np.vstack(Xlim).T
            self.mag_limit = np.sum(Xlim*self.C_mag_limit.params, axis=1)
        else:
            self.mag_limit = 99.0*np.ones_like(mag)

        # Simple analysis of proximity to "good" points used for model building
        mx,my = self.ox[self.idx], self.oy[self.idx]
        kdo = cKDTree(np.array([x, y]).T)
        kdm = cKDTree(np.array([mx, my]).T)

        # Mean distance between "good" points
        mr0 = np.sqrt(self.width*self.height/np.sum(self.idx))

        # Closest "good" points
        m = kdm.query_ball_tree(kdm, 5.0*mr0)

        dists = []
        for i,ii in enumerate(m):
            if len(ii) > 1:
                d1 = [np.hypot(mx[i] - mx[_], my[i] - my[_]) for _ in ii]
                d1 = np.sort(d1)

                dists.append(d1[1])
        mr1 = np.median(dists)

        # Closest "good" points to objects
        m = kdo.query_ball_tree(kdm, 5.0*mr1)

        self.good_idx = np.array([len(_) > 1 for _ in m])

        if verbose:
            print(np.sum(self.good_idx), 'of', len(mag), 'objects are at good distances from model points')
