#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import posixpath, glob, sys

from astropy.io import fits

from resolve import resolve
from favor2 import Favor2

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage="usage: %prog [options] arg")
    # 'Object in the field' search
    parser.add_option('-o', '--object', help='Object name, to be visible on all frames', action='store', dest='object', type='str', default=None)

    # 'Frame center in the cone' search
    parser.add_option('--ra', help='Center RA', action='store', dest='ra', type='float', default=None)
    parser.add_option('--dec', help='Center Dec', action='store', dest='dec', type='float', default=None)
    parser.add_option('--sr', help='Search radius', action='store', dest='sr', type='float', default=None)

    # Refinement
    parser.add_option('-c', '--channel', help='Channel', action='store', dest='channel', type='int', default=None)
    parser.add_option('-s', '--shutter', help='Shutter (0=Rolling, 1=Global)', action='store', dest='shutter', type='int', default=None)
    parser.add_option('-t', '--type', help='Image type', action='store', dest='type', type='str', default=None)
    parser.add_option('-f', '--filter', help='Filter', action='store', dest='filter', type='str', default=None)

    parser.add_option('-n', '--night', help='Night of observations', action='store', dest='night', type='str', default=None)
    parser.add_option('--night1', help='First night of observations', action='store', dest='night1', type='str', default=None)
    parser.add_option('--night2', help='Last night of observations', action='store', dest='night2', type='str', default=None)

    parser.add_option('--latest', help='Show latest images first', action='store_true', dest='latest', default=False)

    # Connection
    parser.add_option('-d', '--db', help='Database name', action='store', dest='db', type='str', default='favor2')
    parser.add_option('-H', '--host', help='Database host', action='store', dest='dbhost', type='str', default=None)

    (options,args) = parser.parse_args()

    wheres,wargs = [],[]

    if options.object:
        target_name,target_ra,target_dec = resolve(options.object)

        if target_name:
            print('Object resolved to:', target_name, 'at', target_ra, target_dec, file=sys.stderr)
            wheres += ['q3c_radial_query(ra, dec, %s, %s, radius)']
            wargs += [target_ra, target_dec]
            wheres += ['q3c_poly_query(%s, %s, footprint10)']
            wargs += [target_ra, target_dec]
        else:
            print('Can\'t resolve:', object, file=sys.stderr)
            sys.exit(1)

    elif options.ra is not None and options.dec is not None and options.sr is not None:
        print('Searching for images with centers within', options.sr, 'deg around ', options.ra, options.dec, file=sys.stderr)
        wheres += ['q3c_radial_query(ra, dec, %s, %s, %s)']
        wargs += [options.ra, options.dec, options.sr]

    if options.channel is not None:
        print('Searching for images from channel', options.channel, file=sys.stderr)
        wheres += ['channel=%s']
        wargs += [options.channel]

    if options.shutter is not None:
        print('Searching for images with shutter', options.shutter, file=sys.stderr)
        wheres += ['shutter=%s']
        wargs += [options.shutter]

    if options.type is not None:
        print('Searching for images with type', options.type, file=sys.stderr)
        wheres += ['type=%s']
        wargs += [options.type]

    if options.filter is not None:
        print('Searching for images with filter', options.filter, file=sys.stderr)
        wheres += ['filter=%s']
        wargs += [options.filter]

    if options.night is not None:
        print('Searching for images from night', options.night, file=sys.stderr)
        wheres += ['night=%s']
        wargs += [options.night]

    if options.night1 is not None:
        print('Searching for images night >=', options.night1, file=sys.stderr)
        wheres += ['night>=%s']
        wargs += [options.night1]

    if options.night2 is not None:
        print('Searching for images night <=', options.night2, file=sys.stderr)
        wheres += ['night<=%s']
        wargs += [options.night2]

    favor2 = Favor2(dbname=options.db, dbhost=options.dbhost)

    if not favor2:
        print('Can\'t connect to the database', file=sys.stderr)
        sys.exit(1)

    res = favor2.query('SELECT filename FROM images WHERE ' + ' AND '.join(wheres) + ' ORDER BY time ' + ('DESC' if options.latest else 'ASC'), wargs)
    print(len(res), 'images found', file=sys.stderr)

    for r in res:
        print(r['filename'])
