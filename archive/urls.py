from __future__ import absolute_import, division, print_function, unicode_literals

from django.http import HttpResponse
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.urls import include, re_path

from . import settings

from . import views
from . import views_images
from . import views_photometry

urlpatterns = [
    # Index
    re_path(r'index^$', views.index, name="index"),
    re_path(r'^$', views.index, name="index"),

    # Images
    re_path(r'^images/?$', views_images.images_list, name='images'),

    # Nights
    re_path(r'^nights/?$', views_images.images_nights, name='nights'),

    # Detailed image view
    re_path(r'^images/(?P<id>\d+)/?$', views_images.image_details, name='image_details'),
    re_path(r'^images/(?P<id>\d+)/download$', views_images.image_download, name='image_download'),
    re_path(r'^images/(?P<id>\d+)/download/processed$', views_images.image_download, {'raw':False}, name='image_download_processed'),
    re_path(r'^images/(?P<id>\d+)/full$', views_images.image_preview, name='image_full'),
    re_path(r'^images/(?P<id>\d+)/view$', views_images.image_preview, {'size':800}, name='image_view'),
    re_path(r'^images/(?P<id>\d+)/preview$', views_images.image_preview, {'size':128}, name='image_preview'),
    # Image analysis
    re_path(r'^images/(?P<id>\d+)/bg$', views_images.image_analysis, {'mode':'bg'}, name='image_bg'),
    re_path(r'^images/(?P<id>\d+)/fwhm$', views_images.image_analysis, {'mode':'fwhm'}, name='image_fwhm'),
    re_path(r'^images/(?P<id>\d+)/wcs$', views_images.image_analysis, {'mode':'wcs'}, name='image_wcs'),
    re_path(r'^images/(?P<id>\d+)/filters$', views_images.image_analysis, {'mode':'filters'}, name='image_filters'),
    re_path(r'^images/(?P<id>\d+)/zero$', views_images.image_analysis, {'mode':'zero'}, name='image_zero'),

    # Cutouts
    re_path(r'^images/cutouts/?$', views_images.images_cutouts, name='images_cutouts'),
    re_path(r'^images/(?P<id>\d+)/cutout$', views_images.image_cutout, name='image_cutout'),
    re_path(r'^images/(?P<id>\d+)/cutout/preview$', views_images.image_cutout, {'size':300}, name='image_cutout_preview'),
    re_path(r'^images/(?P<id>\d+)/cutout/download$', views_images.image_cutout, {'mode':'download'}, name='image_cutout_download'),

    # Photometry
    # re_path(r'^photometry/?$', views_photometry.photometry, name='photometry'),
    re_path(r'^photometry/lc$', views_photometry.lc, {'mode': 'jpeg'}, name='photometry_lc'),
    re_path(r'^photometry/json$', views_photometry.lc, {'mode': 'json'}, name='photometry_json'),
    re_path(r'^photometry/text$', views_photometry.lc, {'mode': 'text'}, name='photometry_text'),
    re_path(r'^photometry/mjd$', views_photometry.lc, {'mode': 'mjd'}, name='photometry_mjd'),

    # Search
    re_path(r'^search/?$', views.search, name='search'),
    re_path(r'^search/cutouts/?$', views.search, {'mode':'cutouts'}, name='search_cutouts'),
    re_path(r'^search/photometry/?$', views.search, {'mode':'photometry'}, name='search_photometry'),

    # Robots
    re_path(r'^robots.txt$', lambda r: HttpResponse("User-agent: *\nDisallow: /\n", content_type="text/plain")),

    # Markdown
    #re_path(r'^about/(?P<path>.*)$', views_markdown.markdown_page, {'base':'about'}, name="markdown"),
]

urlpatterns += staticfiles_urlpatterns()

if settings.DEBUG:
    import debug_toolbar
    urlpatterns = [
        re_path(r'^__debug__/', include(debug_toolbar.urls)),
    ] + urlpatterns
