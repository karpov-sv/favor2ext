# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey has `on_delete` set to the desired behavior.
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from __future__ import unicode_literals

from django.db import models
from django.contrib.postgres.fields import JSONField

class Images(models.Model):
    id = models.IntegerField(primary_key=True)
    filename = models.TextField(unique=True, blank=True, null=True)
    night = models.TextField(blank=True, null=True)
    time = models.DateTimeField(blank=True, null=True)
    channel = models.IntegerField(blank=True, null=True)
    type = models.TextField(blank=True, null=True)
    filter = models.TextField(blank=True, null=True)
    exposure = models.FloatField(blank=True, null=True)
    shutter = models.IntegerField(blank=True, null=True)
    pos0 = models.IntegerField(blank=True, null=True)
    pos1 = models.IntegerField(blank=True, null=True)
    ra = models.FloatField(blank=True, null=True)
    dec = models.FloatField(blank=True, null=True)
    radius = models.FloatField(blank=True, null=True)
    width = models.IntegerField(blank=True, null=True)
    height = models.IntegerField(blank=True, null=True)
    # footprints skipped
    mean = models.FloatField(blank=True, null=True)
    keywords = JSONField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'images'
        app_label = 'favor2'

class Photometry(models.Model):
    time = models.DateTimeField(blank=True, null=True, primary_key=True)
    channel = models.IntegerField(blank=True, null=True)
    # night = models.TextField(blank=True, null=True)
    # filter = models.TextField(blank=True, null=True)

    ra = models.FloatField(blank=True, null=True)
    dec = models.FloatField(blank=True, null=True)
    mag = models.FloatField(blank=True, null=True)
    magerr = models.FloatField(blank=True, null=True)

    x = models.FloatField(blank=True, null=True)
    y = models.FloatField(blank=True, null=True)
    fwhm = models.FloatField(blank=True, null=True)

    flags = models.FloatField(blank=True, null=True)
    color_term = models.FloatField(blank=True, null=True)
    color_term2 = models.FloatField(blank=True, null=True)
    color_term3 = models.FloatField(blank=True, null=True)

    std = models.FloatField(blank=True, null=True)
    nstars = models.FloatField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'photometry_all'
        app_label = 'favor2'
