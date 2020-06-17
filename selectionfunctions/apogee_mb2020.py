#!/usr/bin/env python
#
# cog_ii.py
# Reads the Gaia DR2 selection function from Completeness
# of the Gaia-verse Paper II, Boubert & Everall (2020).
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
import h5py, tqdm
import numpy as np

import astropy.coordinates as coordinates
import astropy.units as units
import h5py
import healpy as hp
from scipy import interpolate, special, spatial

from .std_paths import *
from .map import SelectionFunction, ensure_flat_icrs, coord2healpix
from .source import ensure_tmass_hjk
from . import fetch_utils

from time import time

class apogee_sf(SelectionFunction):
    """
    Queries the Gaia DR2 selection function (Boubert & Everall, 2019).
    """

    def __init__(self, map_fname=None, multi_radius=True, bounds=True):
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

        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'apogee', 'mackerethbovy_2020.h5')

        t_start = time()

        with h5py.File(map_fname, 'r') as f:
            # Load auxilliary data
            print('Loading auxilliary data ...')
            self._nside = 128
            self._bounds = bounds

            self._h_grid = f['h_bins'][...]
            self._jk_grid = f['jk_bins'][...]

            _apogee_count = f['apogee'][...].astype(float)
            _tmass_count = f['tmass'][...].astype(float)

            _ra_field = f['racen'][...]
            _dec_field = f['deccen'][...]
            self._radius_field = f['radius'][...]
            self._unique_radii = np.unique(f['radius'][...])
            print(self._unique_radii)

        t_auxilliary = time()

        self.multi_radius = multi_radius

        # Calculate apogee selection function from counts
        self._sf_field = np.where(_tmass_count>_apogee_count, _apogee_count.astype(float)/_tmass_count.astype(float),
                                    np.where(_apogee_count==0, 0., 1.))

        # Create KDTree for field centers
        xyz_field = np.stack([np.cos(np.deg2rad(_ra_field))*np.cos(np.deg2rad(_dec_field)),
                               np.sin(np.deg2rad(_ra_field))*np.cos(np.deg2rad(_dec_field)),
                               np.sin(np.deg2rad(_dec_field))]).T

        if self.multi_radius: self.tree_field_radii = [spatial.cKDTree(xyz_field[self._radius_field==rad]) for rad in self._unique_radii]
        else: self.tree_field = spatial.cKDTree(xyz_field)

        if bounds == True:
            self._h_min = self._h_grid[0]; self._h_max = self._h_grid[-1];
            self._jk_min = self._jk_grid[0]; self._jk_max = self._jk_grid[-1];
        else:
            self._h_min = -np.inf; self._h_max = np.inf;
            self._jk_min = -np.inf; self._jk_max = np.inf;

        t_sf = time()

        t_finish = time()

        print('t = {:.3f} s'.format(t_finish - t_start))
        print('  auxilliary: {: >7.3f} s'.format(t_auxilliary-t_start))
        print('          sf: {: >7.3f} s'.format(t_sf-t_auxilliary))

    def _selection_function(self,_ra, _dec, _h, _jk):


        # Get magnitude ids
        Hid = np.zeros(_h.shape).astype(int)-1
        for ii in range(len(self._h_grid)): Hid += (_h>self._h_grid[ii]).astype(int)
        JKid = np.zeros(_jk.shape).astype(int)-1
        for ii in range(len(self._jk_grid)): JKid += (_jk>self._jk_grid[ii]).astype(int)

        xyz_source = np.stack([np.cos(np.deg2rad(_ra))*np.cos(np.deg2rad(_dec)),
                               np.sin(np.deg2rad(_ra))*np.cos(np.deg2rad(_dec)),
                               np.sin(np.deg2rad(_dec))]).T

        # Resultant Selection Function is union of overlapping fields (1 - product of non-selection)
        if self.multi_radius:
            _result = np.ones(len(_ra))
            for ii in range(len(self._unique_radii)):
                crossmatch = self.tree_field_radii[ii].query_ball_point(xyz_source, 2*np.sin(np.deg2rad(self._unique_radii[ii])/2))
                _result *= np.array([np.product(1 - self._sf_field[crossmatch[ii],Hid[ii],JKid[ii]]) for ii in range(len(crossmatch))])
            _result = 1-_result
        else:
            crossmatch = self.tree_field.query_ball_point(xyz_source, 2*np.sin(np.deg2rad(self._radius_field)/2))
            _result = 1 - np.array([np.product(1 - self._sf_field[crossmatch[ii],Hid[ii],JKid[ii]]) for ii in range(len(crossmatch))])

        return _result


    @ensure_flat_icrs
    @ensure_tmass_hjk
    def query(self, sources):
        """
        Returns the selection function at the requested coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            Selection function at the specified coordinates, as a fraction.

        """

        # Convert coordinates to healpix indices
        _ra = sources.coord.ra.deg; _dec = sources.coord.dec.deg

        # Extract 2MASS H magnitude and J-K colour
        H = sources.photometry.measurement['tmass_h']

        print(sources.photometry.measurement.keys())

        if 'tmass_j_tmass_k' in sources.photometry.measurement.keys():
            JK = sources.photometry.measurement['tmass_j_tmass_k']
        else:
            JK = sources.photometry.measurement['tmass_j'] - sources.photometry.measurement['tmass_k']

        shape = _ra.shape
        # Evaluate selection function
        selection_function = self._selection_function(_ra.flatten(), _dec.flatten(), H.flatten(), JK.flatten()).reshape(shape)

        if self._bounds == True:
            _outside_bounds = np.where( (H<self._h_min) | (H>self._h_max) | (JK<self._jk_min) | (JK>self._jk_max) )
            selection_function[_outside_bounds] = 0.0

        return selection_function


def fetch():
    """
    Downloads the specified version of the Bayestar dust map.

    Args:
        version (Optional[:obj:`str`]): The map version to download. Valid versions are
            :obj:`'bayestar2019'` (Green, Schlafly, Finkbeiner et al. 2019),
            :obj:`'bayestar2017'` (Green, Schlafly, Finkbeiner et al. 2018) and
            :obj:`'bayestar2015'` (Green, Schlafly, Finkbeiner et al. 2015). Defaults
            to :obj:`'bayestar2019'`.

    Raises:
        :obj:`ValueError`: The requested version of the map does not exist.

        :obj:`DownloadError`: Either no matching file was found under the given DOI, or
            the MD5 sum of the file was not as expected.

        :obj:`requests.exceptions.HTTPError`: The given DOI does not exist, or there
            was a problem connecting to the Dataverse.
    """

    doi = '10.7910/DVN/PDFOVC'

    requirements = {'filename': 'mackerethbovy_2020.h5'}

    local_fname = os.path.join(data_dir(), 'apogee', 'mackerethbovy_2020.h5')

    # Download the data
    fetch_utils.dataverse_download_doi(
        doi,
        local_fname,
        file_requirements=requirements)
