#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import posixpath, glob, sys, os

import cPickle as pickle

from StringIO import StringIO

from favor2 import Favor2
import survey

def touch(filename):
    with open(filename, 'a'):
        pass

    os.utime(filename, None)

if __name__ == '__main__':
    from optparse import OptionParser

    parser = OptionParser(usage="usage: %prog [options] arg")
    parser.add_option('-d', '--db', help='Database name', action='store', dest='db', type='str', default='favor2')
    parser.add_option('-H', '--host', help='Database host', action='store', dest='dbhost', type='str', default=None)
    parser.add_option('-t', '--table', help='Database table', action='store', dest='table', type='str', default='photometry_staging')
    parser.add_option('-r', '--replace', help='Replace already existing records in database', action='store_true', dest='replace', default=False)
    parser.add_option('-i', '--ignore', help='Ignore upload status file', action='store_true', dest='ignore', default=False)
    parser.add_option('-v', '--verbose', help='Verbose', action='store_true', dest='verbose', default=False)

    (options,files) = parser.parse_args()

    print('Uploading to table', options.table)

    favor2 = Favor2(dbname=options.db, dbhost=options.dbhost)
    favor2.conn.autocommit = False
    cur = favor2.conn.cursor()

    N = 0

    s = StringIO()
    filenames = []

    for i,filename in enumerate(files):
        if not options.replace and not options.ignore and posixpath.exists(filename + '.upload'):
            # print('Skipping', filename, 'as upload file exists')
            continue
            # pass

        obj = survey.load_objects(filename)

        if obj is None:
            if not options.replace:
                touch(filename + '.upload')
            continue

        if len(files) > 1:
            print(i, '/', len(files), filename, len(obj['ra']))
            sys.stdout.flush()

        if False and favor2.query('SELECT EXISTS(SELECT time FROM ' + options.table + ' WHERE time=%s AND channel=%s);', (obj['time'], obj['channel']), simplify=True):
            if not options.replace:
                touch(filename + '.upload')

            # print('Skipping', filename, 'as already in DB')
            continue

        for i in xrange(len(obj['ra'])):
            s.write("\t".join(str(_) for _ in (obj['time'], obj['channel'],
                                               obj['ra'][i], obj['dec'][i],
                                               obj['calib_mag'][i], obj['calib_magerr'][i],
                                               obj['x'][i], obj['y'][i],
                                               obj['fwhm'][i], int(obj['flags'][i]),
                                               obj['color_term'][i], obj['color_term2'][i], obj['color_term3'][i],
                                               obj['std'], obj['nstars'])) + '\n')

        filenames.append(filename)

        N += 1

        if N % 100 == 0:
            s.seek(0)
            cur.copy_from(s, options.table, sep='\t', columns=['time', 'channel',
                                                               'ra', 'dec', 'mag', 'magerr',
                                                               'x', 'y', 'fwhm', 'flags',
                                                               'color_term', 'color_term2', 'color_term3',
                                                               'std', 'nstars'], size=65535000)

            print('committing...')

            favor2.conn.commit()

            for fn in filenames:
                if not posixpath.exists(fn + '.upload'):
                    touch(fn + '.upload')

            s = StringIO()
            filenames = []

    if N > 0:
        print('committing...')
        favor2.conn.commit()
        for fn in filenames:
            if not posixpath.exists(fn + '.upload'):
                touch(fn + '.upload')

    print(N, '/', len(files), 'files uploaded')
