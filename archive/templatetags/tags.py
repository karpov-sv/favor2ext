from __future__ import absolute_import, division, print_function, unicode_literals

from django import template
from django.urls import get_script_prefix

register = template.Library()

@register.simple_tag(takes_context=True)
def get_root(context):
    #return context['request'].path
    root = get_script_prefix()

    if root != '/':
        return root
    else:
        return ''
