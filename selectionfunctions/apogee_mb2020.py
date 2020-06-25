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

    def __init__(self, map_fname=None, bounds=True, apogee=1, hemisphere='north'):
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
            if apogee==2:
                if hemisphere.lower()=='north': map_fname = os.path.join(data_dir(), 'apogee', 'mackerethbovy_2020_apo2north.h')
                elif hemisphere.lower()=='south': map_fname = os.path.join(data_dir(), 'apogee', 'mackerethbovy_2020_apo2south.h')
                else: raise ValueError('For apogee 2, hemisphere must be "north" or "south"')
            elif apogee==1:
                if hemisphere.lower()=='south': print('apogee 1 is only in northern hemisphere. If you want apogee 2, include apogee=2 in kwargs.')
                map_fname = os.path.join(data_dir(), 'apogee', 'mackerethbovy_2020_apo1.h')
            else: raise ValueError('apogee kwarg must be 1 or 2.')
        t_start = time()

        with h5py.File(map_fname, 'r') as f:
            # Load auxilliary data
            print('Loading auxilliary data ...')

            _frac4complete = f['_frac4complete'][...]
            minnspec = f['_minnspec'][...]

            _empty_field = np.isnan(f['radius'][:,0])

            self._locations = f['_locations'][...][~_empty_field]

            self._short_completion = f['_short_completion'][...][~_empty_field]
            self._short_completion[np.isnan(self._short_completion)] = 0.
            self._medium_completion = f['_medium_completion'][...][~_empty_field]
            self._medium_completion[np.isnan(self._medium_completion)] = 0.
            self._long_completion = f['_long_completion'][...][~_empty_field]
            self._long_completion[np.isnan(self._long_completion)] = 0.

            self._nspec_short = f['_nspec_short'][...][~_empty_field]
            self._nspec_medium = f['_nspec_medium'][...][~_empty_field]
            self._nspec_long = f['_nspec_long'][...][~_empty_field]
            self._nphot_short = f['_nphot_short'][...][~_empty_field]
            self._nphot_medium = f['_nphot_medium'][...][~_empty_field]
            self._nphot_long = f['_nphot_long'][...][~_empty_field]

            _ra_field = f['racen'][...][~_empty_field]
            _dec_field = f['deccen'][...][~_empty_field]
            self._radius_field = f['radius'][:,0][~_empty_field]
            self._unique_radii = np.unique(self._radius_field)

            self._h_bin = np.zeros((len(self._locations), 4))
            self._h_bin[:,0] = f['_short_hmin'][~_empty_field]
            self._h_bin[:,1] = f['_short_hmax'][~_empty_field]
            self._h_bin[:,2] = f['_medium_hmax'][~_empty_field]
            self._h_bin[:,3] = f['_long_hmax'][~_empty_field]
            self._h_bin[np.isnan(self._h_bin)] = np.inf

            self._jk_bin = np.zeros((len(self._locations), 6))
            self._jk_bin[:,:5] = f['_color_bins_jkmin'][...][~_empty_field]
            self._jk_bin[:,5] = np.inf

        t_auxilliary = time()

        # Selection function: Nfield x Nhbin x Njkbin
        self._sf_field = np.zeros((len(self._locations), 3, 5))
        for jj in range(5):
            self._sf_field[:,0,jj] = np.where((np.nanmax(self._short_completion, axis=1) >= _frac4complete)\
                                    & (np.nansum(self._nspec_short, axis=1) >= minnspec),
                                        self._nspec_short[:,jj]/self._nphot_short[:,jj], np.nan)
            self._sf_field[:,1,jj] = np.where((np.nanmax(self._medium_completion, axis=1) >= _frac4complete)\
                                    & (np.nansum(self._nspec_medium, axis=1) >= minnspec),
                                        self._nspec_medium[:,jj]/self._nphot_medium[:,jj], np.nan)
            self._sf_field[:,2,jj] = np.where((np.nanmax(self._long_completion, axis=1) >= _frac4complete)\
                                    & (np.nansum(self._nspec_long, axis=1) >= minnspec),
                                        self._nspec_long[:,jj]/self._nphot_long[:,jj], np.nan)

        # Create KDTree for field centers
        self.xyz_field = np.stack([np.cos(np.deg2rad(_ra_field))*np.cos(np.deg2rad(_dec_field)),
                               np.sin(np.deg2rad(_ra_field))*np.cos(np.deg2rad(_dec_field)),
                               np.sin(np.deg2rad(_dec_field))]).T

        self._bounds = bounds
        if bounds == True:
            self._h_min = -np.inf; self._h_max = np.inf;
            self._jk_min = -np.inf; self._jk_max = np.inf;

        t_sf = time()

        t_finish = time()

        print('t = {:.3f} s'.format(t_finish - t_start))
        print('  auxilliary: {: >7.3f} s'.format(t_auxilliary-t_start))
        print('          sf: {: >7.3f} s'.format(t_sf-t_auxilliary))

    def _selection_function(self,_ra, _dec, _h, _jk):

        # Build KDTree
        xyz_source = np.stack([np.cos(np.deg2rad(_ra))*np.cos(np.deg2rad(_dec)),
                               np.sin(np.deg2rad(_ra))*np.cos(np.deg2rad(_dec)),
                               np.sin(np.deg2rad(_dec))]).T
        tree_source = spatial.cKDTree(xyz_source)

        # Resultant Selection Function is union of overlapping fields (1 - product of non-selection)
        _result = np.ones(len(_ra))
        for radius in self._unique_radii:
            # location_ids with this radius
            loc_ids=np.argwhere(self._radius_field==radius)
            # Map fields to sources
            crossmatch=tree_source.query_ball_point(self.xyz_field[self._radius_field==radius], 2*np.sin(np.deg2rad(radius)/2))
            for ii in range(len(crossmatch)):

                # ID of H bin
                Hid = np.zeros(len(crossmatch[ii])).astype(int)-1
                for jj in range(self._h_bin.shape[1]): Hid += (_h[crossmatch[ii]]>self._h_bin[loc_ids[ii],jj]).astype(int)
                # ID of JK bin
                JKid = np.zeros(len(crossmatch[ii])).astype(int)-1
                for jj in range(self._jk_bin.shape[1]): JKid += (_jk[crossmatch[ii]]>self._jk_bin[loc_ids[ii],jj]).astype(int)

                # Prod(P(not selected))
                # Replace nan values with 0.
                _result[crossmatch[ii]] *= np.where(np.isnan(self._sf_field[loc_ids[ii],Hid,JKid]), \
                                                    1, 1 - self._sf_field[loc_ids[ii],Hid,JKid])

        # Union is 1-probability of not being selected on any field.
        _result = 1-_result

        #_result[np.isnan(_result)] = 0.

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


class apogeeCombinedSelect(SelectionFunction):

    def __init__(self):
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

        self.sf_inst = [apogee_sf(apogee=1),
                        apogee_sf(apogee=2, hemisphere='south'),
                        apogee_sf(apogee=2, hemisphere='north')]

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

    for version in ["apogee2north", "apogee2south", "apogee1"]:
        requirements = {'filename': 'mackerethbovy_2020_%s.h5' % version}

        local_fname = os.path.join(data_dir(), 'apogee', 'mackerethbovy_2020_%s.h5' % version)

        # Download the data
        fetch_utils.dataverse_download_doi(
            doi,
            local_fname,
            file_requirements=requirements)
