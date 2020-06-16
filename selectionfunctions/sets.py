#!/usr/bin/env python
#
# sets.py
# Generates selection functions for combined samples (e.g. APOGEE + GaiaDR2)
#
# Copyright (C) 2020  Douglas Boubert & Andrew Everall
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

import os
import h5py
import numpy as np

import astropy.coordinates as coordinates
import astropy.units as units
import h5py
import healpy as hp
from scipy import interpolate, special

from .std_paths import *
from .map import SelectionFunction, ensure_flat_icrs, coord2healpix
from . import fetch_utils

from time import time


class intersect_sf(SelectionFunction):
    """
    Queries the Gaia DR2 selection function (Boubert & Everall, 2019).
    """

    def __init__(self, instances):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename of the BoubertEverall2019 selection function. Defaults to
                :obj:`None`, meaning that the default location is used.
            version (Optional[:obj:`str`]): The selection function version to download. Valid versions
                are :obj:`'modelT'` and :obj:`'modelAB'`
                Defaults to :obj:`'modelT'`.
            crowding (Optional[:obj:`bool`]): Whether or not the selection function includes crowding.
                Defaults to :obj:`'False'`.
            bounds (Optional[:obj:`bool`]): Whether or not the selection function is bounded to 0.0 < G < 25.0.
                Defaults to :obj:`'True'`.
        """

        t_start = time()

        self.sf_inst = instances

        t_sf = time()
        t_finish = time()

        print('t = {:.3f} s'.format(t_finish - t_start))
        print('          sf: {: >7.3f} s'.format(t_sf-t_start))

    def _selection_function(self,sources, sources_shape):


        _result=np.ones(sources_shape)

        for ii in range(len(self.sf_inst)):
            _result *= self.sf_inst[ii](sources)

        return _result


    @ensure_flat_icrs
    def query(self, sources):
        """
        Returns the selection function at the requested coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            Selection function at the specified coordinates, as a fraction.

        """

        # Extract sources array shape (if float, get's converted to 1D array)
        sources_shape = np.array(next(iter(sources.photometry.measurement.values()))).shape

        # Evaluate selection function
        selection_function = self._selection_function(sources, sources_shape)

        return selection_function


class union_sf(SelectionFunction):
    """
    Queries the Gaia DR2 selection function (Boubert & Everall, 2019).
    """

    def __init__(self, instances):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename of the BoubertEverall2019 selection function. Defaults to
                :obj:`None`, meaning that the default location is used.
            version (Optional[:obj:`str`]): The selection function version to download. Valid versions
                are :obj:`'modelT'` and :obj:`'modelAB'`
                Defaults to :obj:`'modelT'`.
            crowding (Optional[:obj:`bool`]): Whether or not the selection function includes crowding.
                Defaults to :obj:`'False'`.
            bounds (Optional[:obj:`bool`]): Whether or not the selection function is bounded to 0.0 < G < 25.0.
                Defaults to :obj:`'True'`.
        """

        t_start = time()

        self.sf_inst = instances

        t_sf = time()
        t_finish = time()

        print('t = {:.3f} s'.format(t_finish - t_start))
        print('          sf: {: >7.3f} s'.format(t_sf-t_start))

    def _selection_function(self,sources, sources_shape):


        _not_result=np.ones(sources_shape)

        for ii in range(len(self.sf_inst)):
            _not_result *= (1-self.sf_inst[ii](sources))

        _result = 1-_not_result

        return _result


    @ensure_flat_icrs
    def query(self, sources):
        """
        Returns the selection function at the requested coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            Selection function at the specified coordinates, as a fraction.

        """

        # Extract sources array shape (if float, get's converted to 1D array)
        sources_shape = np.array(next(iter(sources.photometry.measurement.values()))).shape

        # Evaluate selection function
        selection_function = self._selection_function(sources, sources_shape)

        return selection_function
