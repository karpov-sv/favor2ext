from __future__ import absolute_import, division, print_function, unicode_literals

from django.http import HttpResponse, FileResponse
from django.template.response import TemplateResponse
from django.shortcuts import redirect
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_protect
from django.db.models import Q

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import json
import datetime

from astropy.time import Time
from astropy.stats import mad_std

from .models import Photometry

def radectoxieta(ra, dec, ra0=0, dec0=0):
    ra,dec = [np.asarray(_) for _ in ra,dec]
    delta_ra = np.asarray(ra - ra0)

    delta_ra[(ra < 10) & (ra0 > 350)] += 360
    delta_ra[(ra > 350) & (ra0 < 10)] -= 360

    xx = np.cos(dec*np.pi/180)*np.sin(delta_ra*np.pi/180)
    yy = np.sin(dec0*np.pi/180)*np.sin(dec*np.pi/180) + np.cos(dec0*np.pi/180)*np.cos(dec*np.pi/180)*np.cos(delta_ra*np.pi/180)
    xi = (xx/yy)

    xx = np.cos(dec0*np.pi/180)*np.sin(dec*np.pi/180) - np.sin(dec0*np.pi/180)*np.cos(dec*np.pi/180)*np.cos(delta_ra*np.pi/180)
    eta = (xx/yy)

    xi *= 180./np.pi
    eta *= 180./np.pi

    return xi,eta

def get_lc(request):
    lc = Photometry.objects.order_by('time')

    night = request.GET.get('night')
    if night and night != 'all':
        lc = lc.filter(night=night)

    night1 = request.GET.get('night1')
    if night1:
        lc = lc.filter(night__gte=night1)

    night2 = request.GET.get('night2')
    if night2:
        lc = lc.filter(night__lte=night2)

    # Filter out bad data
    # lc = lc.filter(Q(night__lt='20190216') | Q(night__gt='20190222'))

    channel = request.GET.get('channel')
    if channel and channel != 'all':
        lc = lc.filter(channel=channek)

    fname = request.GET.get('filter')
    if fname and fname != 'all':
        lc = lc.filter(filter=fname)

    magerr = request.GET.get('magerr')
    if magerr:
        magerr = float(magerr)
        lc = lc.filter(magerr__lt=magerr)

    nstars = request.GET.get('nstars')
    if nstars:
        nstars = int(nstars)
        lc = lc.filter(nstars__gte=nstars)

    ra = float(request.GET.get('ra'))
    dec = float(request.GET.get('dec'))
    sr = float(request.GET.get('sr', 0.01))

    # Lc with centers within given search radius
    lc = lc.extra(where=["q3c_radial_query(ra, dec, %s, %s, %s)"], params=(ra, dec, sr))

    return lc

# TODO: move to lcs.py
def filter_lc(mag, cbv, bv=None):
    if bv is None:
        bv = get_bv(mag, cbv)

    return mag + bv*cbv

def get_bv(mag, cbv):
    X = np.vstack([np.ones_like(mag), cbv]).T
    Y = mag

    idx = np.ones_like(Y, dtype=np.bool)
    for i in xrange(3):
        C = np.linalg.lstsq(X[idx], Y[idx])[0]
        YY = np.sum(X*C, axis=1)
        idx = np.abs(Y-YY) < 3.0*mad_std(Y-YY)

    return -C[1]

def lc(request, mode="jpg", size=800):
    if True:
        lc = get_lc(request)

        times = np.array([_.time for _ in lc])
        channels = np.array([_.channel for _ in lc])
        # filters = np.array([_.filter for _ in lc])
        filters = np.array(['V' for _ in times])
        ras = np.array([_.ra for _ in lc])
        decs = np.array([_.dec for _ in lc])
        mags = np.array([_.mag for _ in lc])
        magerrs = np.array([_.magerr for _ in lc])

        flags = np.array([_.flags for _ in lc])
        fwhms = np.array([_.fwhm for _ in lc])
        color_terms = np.array([_.color_term for _ in lc])

        stds = np.array([_.std for _ in lc])
        nstars = np.array([_.nstars for _ in lc])
    else:
        time0 = datetime.datetime(2018, 1, 1)

        times = np.array([time0 + datetime.timedelta(hours=8)*_ for _ in xrange(1000)])
        channels = np.repeat(1, len(times))
        filters = np.repeat('V', len(times))
        ras = float(request.GET.get('ra')) + np.random.normal(0, 0.001, size=len(times))
        decs = float(request.GET.get('dec')) + np.random.normal(0, 0.001, size=len(times))

        mags = np.random.normal(0, 0.03, size=len(times)) + 10.0 + 1.0*np.sin(2*np.pi*np.array([(_-time0).total_seconds() for _ in times])/3600/24/3/np.e)
        magerrs = np.repeat(0.03, len(times))

        flags = np.repeat(0, len(times))
        fwhms = np.random.normal(2.0, size=len(times))
        color_terms = np.repeat(0.0, len(times))

        stds = np.random.normal(0.8, 0.05, size=len(times))
        nstars = np.repeat(1000, len(times))

    mjds = Time(times).mjd if len(times) else []

    cols = np.array([{'B':'blue', 'V':'green', 'R':'red', 'I':'orange', 'z':'magenta'}.get(_, 'black') for _ in filters])

    ra = float(request.GET.get('ra'))
    dec = float(request.GET.get('dec'))
    sr = float(request.GET.get('sr', 0.01))
    name = request.GET.get('name')

    if request.GET.get('nofiltering'):
        filtering = False
    else:
        filtering = True
    bv_forced = request.GET.get('bv_forced', None)

    # Quality cuts
    idx0 = np.ones_like(mags, dtype=np.bool)
    if filtering:
        idx0 &= flags < 2

        for _ in range(3):
            idx0 &= stds < np.median(stds[idx0]) + 3.0*mad_std(stds[idx0])

        for _ in range(3):
            idx0 &= fwhms < np.median(fwhms[idx0]) + 3.0*mad_std(fwhms[idx0])

        for _ in range(3):
            idx0 &= color_terms < np.median(color_terms[idx0]) + 3.0*mad_std(color_terms[idx0])

    if bv_forced:
        bv = float(bv_forced)
    else:
        bv = get_bv(mags[idx0], color_terms[idx0])
    mags = filter_lc(mags, color_terms, bv)

    context = {}

    context['ra'] = ra
    context['dec'] = dec
    context['sr'] = sr
    context['filtering'] = filtering

    title = '%s - %.4f %.3f %.3f - %d pts - mean=%.2f std=%.2f B-V=%.2f' % (request.GET.get('name'), ra, dec, sr, len(mags), np.mean(mags), np.std(mags), bv)

    xi,eta = radectoxieta(ras, decs, ra, dec)
    xi *= 3600
    eta *= 3600

    if mode == 'jpeg':
        # Plot lc
        fig = Figure(facecolor='white', dpi=72, figsize=(size/72,0.5*size/72), tight_layout=True)
        ax = fig.add_subplot(111)
        ax.grid(True, alpha=0.1, color='gray')

        for fn in np.unique(filters):
            idx = idx0 & (filters == fn)

            if len(mags[idx]) < 2:
                continue

            ax.errorbar(times[idx], mags[idx], magerrs[idx], fmt='.', color=cols[idx][0], capsize=0, alpha=0.3)
            ax.scatter(times[idx], mags[idx], marker='.', c=cols[idx][0])
            ax.invert_yaxis()

        ax.invert_yaxis()

        ax.set_title(title)

        canvas = FigureCanvas(fig)

        response = HttpResponse(content_type='image/jpeg')
        canvas.print_jpg(response)

        return response

    elif mode == 'json':
        lcs = []

        for fn in np.unique(filters):
            idx = idx0 & (filters == fn)

            if len(mags[idx]) < 2:
                continue

            times_idx = [_.isoformat() for _ in times[idx]]

            lcs.append({'filter': fn, 'color': cols[idx][0], 'bv': bv,
                        'times': times_idx, 'mjds': list(mjds[idx]), 'xi': list(xi[idx]), 'eta': list(eta[idx]),
                        'mags': list(mags[idx]), 'magerrs': list(magerrs[idx]), 'flags': list(flags[idx]),
                        'fwhms': list(fwhms[idx]), 'channels': list(channels[idx]), 'color_terms': list(color_terms[idx]),
                        'stds': list(stds[idx]), 'nstars': list(nstars[idx])})

        data = {'name': name, 'title': title, 'ra': ra, 'dec': dec, 'sr': sr, 'lcs': lcs}

        return HttpResponse(json.dumps(data), content_type="application/json")

    elif mode == 'text':
        response = HttpResponse(request, content_type='text/plain')

        response['Content-Disposition'] = 'attachment; filename=lc_full_%s_%s_%s.txt' % (ra, dec, sr)

        print('# Date Time MJD Channel Filter Mag Magerr Flags Std Nstars FWHM ColorTerm', file=response)

        for _ in xrange(len(times)):
            print(times[_], mjds[_], channels[_], filters[_], mags[_], magerrs[_], flags[_], stds[_], nstars[_], fwhms[_], color_terms[_], file=response)

        return response

    elif mode == 'mjd':
        response = HttpResponse(request, content_type='text/plain')

        response['Content-Disposition'] = 'attachment; filename=lc_mjd_%s_%s_%s.txt' % (ra, dec, sr)

        print('# MJD Mag Magerr', file=response)

        idx = idx0

        for _ in xrange(len(times[idx])):
            print(mjds[idx][_], mags[idx][_], magerrs[idx][_], file=response)

        return response
