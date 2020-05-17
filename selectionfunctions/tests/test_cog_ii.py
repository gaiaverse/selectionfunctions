#!/usr/bin/env python
#
# test_cog_ii.py
# Test query code for the Boubert & Everall (2020) selection function.
#
# Copyright (C) 2020  Douglas Boubert & Andrew Everall.
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

from __future__ import print_function, division

import unittest

import numpy as np
import astropy.coordinates as coords
import astropy.units as units
import os
import re
import time

from .. import cog_ii
from ..std_paths import *
from ..source import Source

class TestCoGII(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        print('Loading cog_ii.dr2_sf query object ...')
        t0 = time.time()

        # Set up query object
        self._cog_ii_dr2_sf = cog_ii.dr2_sf(crowding=True)

        t1 = time.time()
        print('Loaded cog_ii.dr2_sf test function in {:.5f} s.'.format(t1-t0))

    def test_plot_selection_function(self):
        # Draw random coordinates, both above and below dec = -30 degree line
        n_pix = 100000
        ra = -180. + 360.*np.random.random(n_pix)
        dec = -45. + 90.*np.random.random(n_pix)    # 45 degrees above/below
        G = 23.5*np.random.random(n_pix)
        c = Source(ra, dec, photometry={'gaia_g':G}, frame='icrs', unit='deg')

        sf_calc = self._cog_ii_dr2_sf(c)

        import matplotlib.pyplot as plt
        plt.hexbin(G,sf_calc,mincnt=1)
        #plt.savefig('./tests/test.png',dpi=500)
        plt.show()

    def test_bounds(self):
        """
        Test that out-of-bounds magnitudes return 0.0 selection.
        """

        # Draw random coordinates, both above and below dec = -30 degree line
        n_pix = 1000
        ra = -180. + 360.*np.random.random(n_pix)
        dec = -75. + 90.*np.random.random(n_pix)    # 45 degrees above/below
        G = -25+100*np.random.random(n_pix)
        c = Source(ra, dec, photometry={'gaia_g':G}, frame='icrs', unit='deg')

        sf_calc = self._cog_ii_dr2_sf(c)

        zero_below = sf_calc[G < 0.0]<1e-8
        zero_above = sf_calc[G > 25.0]<1e-8

        # print r'{:s}: {:.5f}% nan above dec=-25 deg.'.format(mode, 100.*pct_nan_above)

        self.assertTrue(np.all(zero_below))
        self.assertTrue(np.all(zero_above))

    def test_shape(self):
        """
        Test that the output shapes are as expected with input coordinate arrays
        of different shapes.
        """

        for reps in range(5):
            # Draw random coordinates, with different shapes
            n_dim = np.random.randint(1,4)
            shape = np.random.randint(1,7, size=(n_dim,))

            ra = -180. + 360.*np.random.random(shape)
            dec = -90. + 180. * np.random.random(shape)
            G = 1.7+19.8*np.random.random(shape)
            c = Source(ra, dec, photometry={'gaia_g':G}, frame='icrs', unit='deg')

            sf_calc = self._cog_ii_dr2_sf(c)

            np.testing.assert_equal(sf_calc.shape[:n_dim], shape)

            self.assertEqual(len(sf_calc.shape), n_dim) 


if __name__ == '__main__':
    unittest.main()
