#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals

import psycopg2, psycopg2.extras
import datetime

import numpy as np

class DB:
    """Class encapsulating the connection to PostgreSQL database"""
    def __init__(self, dbname='favor2', dbhost='', dbport=0, dbuser='', dbpassword='', readonly=False):
        connstring = "dbname=" + dbname
        if dbhost:
            connstring += " host="+dbhost
        if dbport:
            connstring += " port=%d" % dbport
        if dbuser:
            connstring += " user="+dbuser
        if dbpassword:
            connstring += " password='%s'" % dbpassword

        self.connect(connstring, readonly)

    def connect(self, connstring, readonly=False):
        self.conn = psycopg2.connect(connstring)
        self.conn.autocommit = True
        self.conn.set_session(readonly=readonly)
        psycopg2.extras.register_default_jsonb(self.conn)
        # FIXME: the following adapter is registered globally!
        psycopg2.extensions.register_adapter(dict, psycopg2.extras.Json)

        self.connstring = connstring
        self.readonly = readonly

    def query(self, string="", data=(), simplify=True, debug=False, array=False):
        if self.conn.closed:
            print("Re-connecting to DB")
            self.connect(self.connstring, self.readonly)

        cur = self.conn.cursor(cursor_factory = psycopg2.extras.DictCursor)

        if debug:
            print(cur.mogrify(string, data))

        if data:
            cur.execute(string, data)
        else:
            cur.execute(string)

        try:
            result = cur.fetchall()
            # Simplify the result if it is simple
            if array:
                # Code from astrolibpy, https://code.google.com/p/astrolibpy
                strLength = 10
                __pgTypeHash = {
                    16:bool,18:str,20:'i8',21:'i2',23:'i4',25:'|S%d'%strLength,700:'f4',701:'f8',
                    1042:'|S%d'%strLength,#character()
                    1043:'|S%d'%strLength,#varchar
                    1114:'|O',#datetime
                    1700:'f8' #numeric
                }

                desc = cur.description
                names = [d.name for d in desc]
                formats = [__pgTypeHash.get(d.type_code, '|O') for d in desc]

                table = np.recarray(shape=(cur.rowcount,), formats=formats, names=names)

                for i,v in enumerate(result):
                    table[i] = tuple(v)

                return table
            elif simplify and len(result) == 1:
                if len(result[0]) == 1:
                    return result[0][0]
                else:
                    return result[0]
            else:
                return result
        except:
            # Nothing returned from the query
            #import traceback
            #traceback.print_exc()
            return None

    def get_stars(self, ra0=0, dec0=0, sr0=0, limit=10000, catalog='pickles', extra=[], extrafields=None, debug=False):
        # Code from astrolibpy, https://code.google.com/p/astrolibpy
        strLength = 10
        __pgTypeHash = {
            16:bool,18:str,20:'i8',21:'i2',23:'i4',25:'|S%d'%strLength,700:'f4',701:'f8',
            1042:'|S%d'%strLength,#character()
            1043:'|S%d'%strLength,#varchar
            1114:'|O',#datetime
            1700:'f8' #numeric
        }

        substr = ""

        if catalog == 'tycho2':
            substr = "0.76*bt+0.24*vt as b , 1.09*vt-0.09*bt as v, 0 as r"

        # TODO: Do we really need brightess ordered output?..
        order = ""
        if False:
            if catalog in ['tycho2', 'apass', 'pickles']:
                order = "ORDER BY v"
            elif catalog in ['twomass']:
                order = "ORDER BY j"
            elif catalog in ['atlas']:
                order = "ORDER BY g"

        # if catalog == 'gaia':
        #     if extra and type(extra) == str:
        #         extra = [extra]
        #     elif not extra:
        #         extra = []
        #     extra += ["lum > 0.3 AND lum < 30 AND bp_rp_excess > 1.0 + 0.015*bp_rp*bp_rp AND bp_rp_excess < 1.3 + 0.06*bp_rp*bp_rp "]

        if extra and type(extra) == list:
            extra_str = " AND " + " AND ".join(extra)
        elif extra and type(extra) == str:
            extra_str = " AND " + extra
        else:
            extra_str = ""

        if extrafields:
            if substr:
                substr = substr + "," + extrafields
            else:
                substr = extrafields

        if substr:
            substr = "," + substr

        cur = self.conn.cursor(cursor_factory = psycopg2.extras.DictCursor)
        string = "SELECT * " + substr + " FROM " + catalog + " cat WHERE q3c_radial_query(ra, dec, %s, %s, %s) " + extra_str + " " + order + " LIMIT %s;"
        data = (ra0, dec0, sr0, limit)

        if debug:
            print(cur.mogrify(string, data))

        cur.execute(string, data)

        desc = cur.description
        names = [d.name for d in desc]
        formats = [__pgTypeHash.get(d.type_code, '|O') for d in desc]

        table = np.recarray(shape=(cur.rowcount,), formats=formats, names=names)

        for i,v in enumerate(cur.fetchall()):
            table[i] = tuple(v)

        # Add some computed fields to the table
        if catalog == 'pickles':
            # Pickles - has computed B, V, R and measured J, H, K, Bt, Vt
            bt = table['bt']
            vt = table['vt']

            # http://www.aerith.net/astro/color_conversion.html
            v = vt + 0.00097 - 0.1334 * (bt-vt) + 0.05486 * (bt-vt)**2 - 0.01998 * (bt-vt)**3
            b = v + (bt-vt) - 0.007813 * (bt-vt) - 0.1489 * (bt-vt)**2 + 0.03384 * (bt-vt)**3

            # quick and dirty calibration using Landolt (2009,2013)
            # b = bt + 0.02927929 - 0.16965093*(bt-vt) - 0.07067606*(bt-vt)**2
            # v = vt + 0.01082083 - 0.09096735*(bt-vt)

            # My cross-calibration using Gaia DR2 (see below)
            Cb = [-0.0245864 , -0.06181789, -0.18266344,  0.0063554 ]
            Cv = [-0.02894483,  0.05962909, -0.22627449,  0.10120628]
            Cr = [ 0.01460373, -0.70958196,  0.10229897,  0.02253856]

            bt_vt = bt - vt

            b = bt + Cb[0] + Cb[1]*bt_vt + Cb[2]*bt_vt**2 + Cb[3]*bt_vt**3
            v = vt + Cv[0] + Cv[1]*bt_vt + Cv[2]*bt_vt**2 + Cv[3]*bt_vt**3
            r = vt + Cr[0] + Cr[1]*bt_vt + Cr[2]*bt_vt**2 + Cr[3]*bt_vt**3

            table = np.lib.recfunctions.append_fields(table,
                        ['B', 'V', 'R', 'I', 'Berr', 'Verr', 'Rerr', 'Ierr'],
                        [
                            # B
                            b,
                            # 0.760*table['bt'] + 0.240*table['vt'], # B = 0.76*BT+0.24*VT
                            # V
                            v,
                            # table['vt'] - 0.090*(table['bt'] - table['vt']), # V = VT -0.090*(BT-VT)
                            # R
                            r,
                            # I
                            table['vt'] - 0.090*(table['bt'] - table['vt']) - 1.6069*(table['j'] - table['k']) + 0.0503, # V - Ic = 1.6069 * (J - Ks) + 0.0503
                            # Berr
                            np.hypot(table['ebt'], table['evt']),
                            # Verr
                            np.hypot(table['ebt'], table['evt']),
                            # Rerr
                            np.hypot(table['ebt'], table['evt']),
                            # Ierr
                            np.hypot(table['ebt'], table['evt']),
                        ],
                        [np.double, np.double, np.double, np.double, np.double, np.double, np.double, np.double])

        elif catalog == 'apass':
            # APASS - has measured B, V, g', r', i'
            # g'r'i' to gri - http://classic.sdss.org/dr7/algorithms/jeg_photometric_eq_dr1.html
            g = table['g'] + 0.060*(table['g']-table['r'] - 0.53)
            r = table['r'] + 0.035*(table['r']-table['i'] - 0.21)
            i = table['i'] + 0.041*(table['r']-table['i'] - 0.21)

            table = np.lib.recfunctions.append_fields(table,
                        ['B', 'V', 'R', 'I', 'Berr', 'Verr', 'Rerr', 'Ierr'],
                        [
                            # B
                            table['b'],
                            # V
                            table['v'],
                            # R
                            table['r'] - 0.257*(r-i) - 0.152, # R-r = (-0.257 +/- 0.004)*(r-i) + (0.152 +/- 0.002)
                            # I
                            table['v'] - 0.671*(g-i) - 0.359, # V-I = (0.671 +/- 0.002)*(g-i) + (0.359 +/- 0.002) if g-i <= 2.1
                            # Berr
                            table['berr'],
                            # Verr
                            table['verr'],
                            # Rerr
                            table['rerr'],
                            # Ierr
                            table['ierr'],
                        ],
                        [np.double, np.double, np.double, np.double, np.double, np.double, np.double, np.double])

        elif catalog == 'atlas':
            # ATLAS-refcat2 - has measured Gaia, GaiaBP, GaiaRP, g, r, i, z (PanSTARRS ones!)
            # Official transformations from https://arxiv.org/pdf/1203.0297.pdf
            B = table['g'] + 0.212 + 0.556*(table['g'] - table['r']) + 0.034*(table['g'] - table['r'])**2
            V = table['r'] + 0.005 + 0.462*(table['g'] - table['r']) + 0.013*(table['g'] - table['r'])**2
            R = table['r'] - 0.137 - 0.108*(table['g'] - table['r']) - 0.029*(table['g'] - table['r'])**2
            I = table['i'] - 0.366 - 0.136*(table['g'] - table['r']) - 0.018*(table['g'] - table['r'])**2

            # Alternative transfromation from https://arxiv.org/pdf/1706.06147.pdf, Stetson, seems better with Landolt
            B = table['g'] + 0.199 + 0.540*(table['g'] - table['r']) + 0.016*(table['g'] - table['r'])**2
            V = table['g'] - 0.020 - 0.498*(table['g'] - table['r']) - 0.008*(table['g'] - table['r'])**2
            R = table['r'] - 0.163 - 0.086*(table['g'] - table['r']) - 0.061*(table['g'] - table['r'])**2
            I = table['i'] - 0.387 - 0.123*(table['g'] - table['r']) - 0.034*(table['g'] - table['r'])**2

            table = np.lib.recfunctions.append_fields(table,
                        ['B', 'V', 'R', 'I', 'Berr', 'Verr', 'Rerr', 'Ierr', 'zerr'],
                        [
                            # corrected by me for ~0.02 mag median difference with Landolt (1992) standards
                            # B
                            B,
                            # V
                            V,
                            # R
                            R,
                            # I
                            I,
                            # Berr
                            table['dg'],
                            # Verr
                            table['dg'],
                            # Rerr
                            table['dr'],
                            # Ierr
                            table['di'],
                            # zerr
                            table['dz'],
                        ],
                        [np.double, np.double, np.double, np.double, np.double, np.double, np.double, np.double, np.double])

        elif catalog == 'gaia':
            # Gaia DR2 data
            # My simple fit using Landolt standards
            # pB = [ 0.00827462, -0.07061174,  0.37224837,  0.66911232,  0.01612951]
            # pV = [-0.05053524,  0.27906535, -0.34027179,  0.3705785 , -0.03900042]
            # pR = [-0.03439028,  0.18250704, -0.17043058, -0.22162674,  0.0165021 ]
            # pI = [ 0.02478455, -0.16225798,  0.46062943, -1.07182209,  0.09876639]
            pB,pCB = [-0.05927724559795761, 0.4224326324292696, 0.626219707920836, -0.011211539139725953], [876.4047401692277, 5.114021693079334, -2.7332873314449326, 0]
            pV,pCV = [0.0017624722901609662, 0.15671377090187089, 0.03123927839356175, 0.041448557506784556],[98.03049528983964, 20.582521666713028, 0.8690079603974803, 0]
            pR,pCR = [0.02045449129406191, 0.054005149296716175, -0.3135475489352255, 0.020545083667168156], [347.42190542330945, 39.42482430363565, 0.8626828845232541, 0]
            pI,pCI = [0.005092289380850884, 0.07027022935721515, -0.7025553064161775, -0.02747532184796779], [79.4028706486939, 9.176899238787003, -0.7826315256072135, 0]
            g = table['g']
            bp_rp = table['bp'] - table['rp']

            # https://www.cosmos.esa.int/web/gaia/dr2-known-issues#PhotometrySystematicEffectsAndResponseCurves
            gcorr = g.copy()
            gcorr[(g>2)&(g<6)] = -0.047344 + 1.16405*g[(g>2)&(g<6)] - 0.046799*g[(g>2)&(g<6)]**2 + 0.0035015*g[(g>2)&(g<6)]**3
            gcorr[(g>6)&(g<16)] = g[(g>6)&(g<16)] - 0.0032*(g[(g>6)&(g<16)] - 6)
            gcorr[g>16] = g[g>16] - 0.032
            g = gcorr

            Cstar = table['bp_rp_excess'] - np.polyval([-0.00445024,  0.0570293,  -0.02810592,  1.20477819], bp_rp)

            B = g + np.polyval(pB, bp_rp) + np.polyval(pCB, Cstar)
            V = g + np.polyval(pV, bp_rp) + np.polyval(pCV, Cstar)
            R = g + np.polyval(pR, bp_rp) + np.polyval(pCR, Cstar)
            I = g + np.polyval(pI, bp_rp) + np.polyval(pCI, Cstar)

            err = np.sqrt(3)*table['dg']

            table = np.lib.recfunctions.append_fields(table,
                        ['B', 'V', 'R', 'I', 'Berr', 'Verr', 'Rerr', 'Ierr'],
                        [
                            B,
                            V,
                            R,
                            I,
                            err,
                            err,
                            err,
                            err,
                        ],
                        [np.double, np.double, np.double, np.double, np.double, np.double, np.double, np.double])

        elif catalog == 'gaiaedr3':
            # Gaia EDR3 data
            # My simple fit using Landolt standards
            pB,pCB = [-1.10797403e-01,  5.05107200e-01,  6.08461421e-01, -2.29323596e-03], -7.10179154
            pV,pCV = [-0.00815956,  0.18593159,  0.02125959,  0.01791784],  0.25652809
            pR,pCR = [0.02395484,  0.05862712, -0.32681951,  0.0100306],   2.00476047
            pI,pCI = [-5.58513407e-04,  9.78202588e-02, -7.53451026e-01, -1.12096873e-02], -1.22532996

            g = table['g']
            bp_rp = table['bp'] - table['rp']
            bp_rp_excess = table['bp_rp_excess']

            B = g + np.polyval(pB, bp_rp) + pCB*bp_rp_excess
            V = g + np.polyval(pV, bp_rp) + pCV*bp_rp_excess
            R = g + np.polyval(pR, bp_rp) + pCR*bp_rp_excess
            I = g + np.polyval(pI, bp_rp) + pCI*bp_rp_excess

            err = np.sqrt(3)*table['dg']

            table = np.lib.recfunctions.append_fields(table,
                        ['B', 'V', 'R', 'I', 'Berr', 'Verr', 'Rerr', 'Ierr'],
                        [
                            B,
                            V,
                            R,
                            I,
                            err,
                            err,
                            err,
                            err,
                        ],
                        [np.double, np.double, np.double, np.double, np.double, np.double, np.double, np.double])
        return table
