from __future__ import absolute_import, division, print_function, unicode_literals

from django.http import HttpResponse
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf.urls import include, url

from . import settings

from . import views
from . import views_images
# from . import views_photometry

urlpatterns = [
    # Index
    url(r'index^$', views.index, name="index"),
    url(r'^$', views.index, name="index"),

    # Images
    url(r'^images/?$', views_images.images_list, name='images'),

    # Nights
    url(r'^nights/?$', views_images.images_nights, name='nights'),

    # Detailed image view
    url(r'^images/(?P<id>\d+)/?$', views_images.image_details, name='image_details'),
    url(r'^images/(?P<id>\d+)/download$', views_images.image_download, name='image_download'),
    url(r'^images/(?P<id>\d+)/download/processed$', views_images.image_download, {'raw':False}, name='image_download_processed'),
    url(r'^images/(?P<id>\d+)/full$', views_images.image_preview, name='image_full'),
    url(r'^images/(?P<id>\d+)/view$', views_images.image_preview, {'size':800}, name='image_view'),
    url(r'^images/(?P<id>\d+)/preview$', views_images.image_preview, {'size':128}, name='image_preview'),
    # Image analysis
    url(r'^images/(?P<id>\d+)/bg$', views_images.image_analysis, {'mode':'bg'}, name='image_bg'),
    url(r'^images/(?P<id>\d+)/fwhm$', views_images.image_analysis, {'mode':'fwhm'}, name='image_fwhm'),
    url(r'^images/(?P<id>\d+)/wcs$', views_images.image_analysis, {'mode':'wcs'}, name='image_wcs'),
    url(r'^images/(?P<id>\d+)/filters$', views_images.image_analysis, {'mode':'filters'}, name='image_filters'),
    url(r'^images/(?P<id>\d+)/zero$', views_images.image_analysis, {'mode':'zero'}, name='image_zero'),

    # Cutouts
    url(r'^images/cutouts/?$', views_images.images_cutouts, name='images_cutouts'),
    url(r'^images/(?P<id>\d+)/cutout$', views_images.image_cutout, name='image_cutout'),
    url(r'^images/(?P<id>\d+)/cutout/preview$', views_images.image_cutout, {'size':300}, name='image_cutout_preview'),
    url(r'^images/(?P<id>\d+)/cutout/download$', views_images.image_cutout, {'mode':'download'}, name='image_cutout_download'),

    # Photometry
    # url(r'^photometry/?$', views_photometry.photometry, name='photometry'),
    # url(r'^photometry/lc$', views_photometry.lc, {'mode': 'jpeg'}, name='photometry_lc'),
    # url(r'^photometry/json$', views_photometry.lc, {'mode': 'json'}, name='photometry_json'),
    # url(r'^photometry/text$', views_photometry.lc, {'mode': 'text'}, name='photometry_text'),
    # url(r'^photometry/mjd$', views_photometry.lc, {'mode': 'mjd'}, name='photometry_mjd'),

    # Search
    url(r'^search/?$', views.search, name='search'),
    url(r'^search/cutouts/?$', views.search, {'mode':'cutouts'}, name='search_cutouts'),
    # url(r'^search/photometry/?$', views.search, {'mode':'photometry'}, name='search_photometry'),

    # Robots
    url(r'^robots.txt$', lambda r: HttpResponse("User-agent: *\nDisallow: /\n", content_type="text/plain")),

    # Markdown
    #url(r'^about/(?P<path>.*)$', views_markdown.markdown_page, {'base':'about'}, name="markdown"),
]

urlpatterns += staticfiles_urlpatterns()

if settings.DEBUG:
    import debug_toolbar
    urlpatterns = [
        url(r'^__debug__/', include(debug_toolbar.urls)),
    ] + urlpatterns
