from __future__ import absolute_import, division, print_function, unicode_literals

from django.http import HttpResponse
from django.template.response import TemplateResponse
from django.shortcuts import redirect
from django.urls import reverse
from django.views.decorators.cache import cache_page
from django.views.decorators.csrf import csrf_protect
from django.db.models import Count, Min, Sum, Avg

import datetime, re, urllib

from .models import Images
from .utils import permission_required_or_403, redirect_get, db_query

# Favor2ext modules
from .parent.resolve import resolve

# @cache_page(3600)
def index(request):
    context = {}

    context['nights'] = db_query('select min(night),max(night),count(*) from images', ())

    return TemplateResponse(request, 'index.html', context=context)

# @cache_page(3600)
#@csrf_protect
def search(request, mode='images'):
    context = {}

    message,message_cutout = None,None

    if request.method == 'GET' and 'coords' in request.GET:# and (request.GET.get('coords') or request.GET.get('filename')):
        # Form submission handling

        params = {}

        for _ in ['type', 'shutter', 'filter', 'night1', 'night2', 'channel', 'maxdist', 'filename', 'coords', 'magerr', 'nstars', 'bv_forced', 'nofiltering', 'detailed']:
            if request.GET.get(_) and request.GET.get(_) != 'all':
                params[_] = request.GET.get(_)

        coords = request.GET.get('coords')
        if request.GET.get('sr_value'):
            sr = float(request.GET.get('sr_value', 0.1))*{'arcsec':1/3600, 'arcmin':1/60, 'deg':1}.get(request.GET.get('sr_units', 'deg'), 1)
            params['sr'] = sr
            params['sr_value'] = request.GET.get('sr_value')
            params['sr_units'] = request.GET.get('sr_units')
        else:
            sr = 0

        if coords:
            name,ra,dec = resolve(coords)
        else:
            name,ra,dec = None,None,None

        if name:
            params['name'] = name
            params['ra'] = ra
            params['dec'] = dec

        if name or mode == 'images':
            if mode == 'cutouts':
                # Search cutouts only
                if sr > 1:
                    params['sr'] = 1

                return redirect_get('images_cutouts',  get=params)

            elif mode == 'photometry':
                # Search photometry database
                if sr > 5/60:
                    params['sr'] = 5/60
                    params['sr_value'] = 1
                    params['sr_units'] = 'arcmin'

                context['lc'] = reverse('photometry_lc') + '?' + urllib.urlencode(params)
                context['lc_json'] = reverse('photometry_json') + '?' + urllib.urlencode(params)
                context['lc_text'] = reverse('photometry_text') + '?' + urllib.urlencode(params)
                context['lc_mjd'] = reverse('photometry_mjd') + '?' + urllib.urlencode(params)

            elif mode == 'images':
                # Search full images
                if name and not sr:
                    context['message'] = "Search radius not set"
                else:
                    return redirect_get('images',  get=params)

        else:
            context['message'] = "Can't resolve query position: " + coords

        context.update(params)

    # Possible values for fields
    types = Images.objects.distinct('type').values('type')
    context['types'] = types

    shutters = Images.objects.distinct('shutter').values('shutter')
    context['shutters'] = shutters

    channels = Images.objects.distinct('channel').values('channel')
    context['channels'] = channels

    filters = Images.objects.distinct('filter').values('filter')
    context['filters'] = filters

    if mode == 'cutouts':
        return TemplateResponse(request, 'cutouts.html', context=context)
    elif mode == 'photometry':
        return TemplateResponse(request, 'photometry.html', context=context)
    else:
        return TemplateResponse(request, 'search.html', context=context)
