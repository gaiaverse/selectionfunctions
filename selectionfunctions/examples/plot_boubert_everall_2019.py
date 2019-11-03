#!/usr/bin/env python
#
# plot_boubert_everall_2019.py
# An example of how to query the Gaia DR2 selection function of
# Boubert & Everall (2019, submitted).
#
# Copyright (C) 2019  Douglas Boubert
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from __future__ import print_function

import numpy as np
import os.path

try:
    import PIL.Image
except ImportError as error:
    print('This example requires Pillow or PIL.\n'
          'See <http://pillow.readthedocs.io/en/stable/installation.html>.')
    raise error

from selectionfunctions.source_base import Source
import astropy.units as u

from selectionfunctions.boubert_everall_2019 import BoubertEverall2019Query


def numpy2pil(a, vmin, vmax):
    a = np.clip((a - vmin) / (vmax - vmin), 0., 1.)
    a = (254.99 * a).astype('u1')
    return PIL.Image.fromarray(a)


def main():
    w,h = (2056,1024)
    l_0 = 130.

    # Set up Bayestar query object
    print('Loading Boubert & Everall 2019 selection function...')
    boubert_everall_2019 = BoubertEverall2019Query()

    # Create a grid of coordinates
    print('Creating grid of coordinates...')
    l = np.linspace(-180.+l_0, 180.+l_0, 2*w)
    b = np.linspace(-90., 90., 2*h+2)
    b = b[1:-1]
    l,b = np.meshgrid(l, b)

    l += (np.random.random(l.shape) - 0.5) * 360./(2.*w)
    b += (np.random.random(l.shape) - 0.5) * 180./(2.*h)

    sf = np.empty(l.shape+(3,), dtype='f8')

    for k,G in enumerate([21.0, 21.2, 21.4]):
        # d = 5.    # We'll query integrated reddening to a distance of 5 kpc
        sources = Source(l*u.deg, b*u.deg, photometry={'gaia_g':G*np.ones(l.shape)}, frame='galactic')

        # Get the dust median reddening at each coordinate
        print('Querying map...')
        sf[:,:,k] = boubert_everall_2019.query(sources)

    # Convert the output array to a PIL image and save
    print('Saving image...')
    img = numpy2pil(sf[::-1,::-1,:], 0., 1)
    img = img.resize((w,h), resample=PIL.Image.LANCZOS)
    fname = 'boubert_everall_2019.png'
    img.save(fname)

    return 0


if __name__ == '__main__':
    main()
