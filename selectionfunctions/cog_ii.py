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
import h5py
import numpy as np

import astropy.coordinates as coordinates
import astropy.units as units
import h5py
import healpy as hp
from scipy import interpolate, special

from .std_paths import *
from .map import SelectionFunction, ensure_flat_icrs, coord2healpix
from .source import ensure_gaia_g
from . import fetch_utils

from time import time


class dr2_sf(SelectionFunction):
    """
    Queries the Gaia DR2 selection function (Boubert & Everall, 2019).
    """

    def __init__(self, map_fname=None, version='modelAB', crowding=False, bounds=True):
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
            map_fname = os.path.join(data_dir(), 'cog_ii', 'cog_ii_dr2.h5')

        t_start = time()

        with h5py.File(map_fname, 'r') as f:
            # Load auxilliary data
            print('Loading auxilliary data ...')
            self._g_grid = f['g_grid'][...]
            self._n_field = f['n_field'][...]
            self._nside = hp.npix2nside(self._n_field.shape[0])
            self._crowding = crowding
            self._bounds = bounds
            if crowding == True:
                self._nside_crowding = 1024
                self._log10_rho_grid = f['log10_rho_grid'][...]
                self._log10_rho_field = np.log10(np.maximum(1.0,f['neighbour_field'][...])/hp.nside2pixarea(self._nside_crowding,degrees=True))

            t_auxilliary = time()

            # Load selection function
            print('Loading selection function ...')
            if version == 'modelT':
                if crowding == True:
                    self._theta = f['t_theta_percentiles'][:,:,2]
                else:
                    self._theta = f['t_theta_percentiles'][0,:,2]
            elif version == 'modelAB':
                if crowding == True:
                    self._alpha = f['ab_alpha_percentiles'][:,:,2]
                    self._beta = f['ab_beta_percentiles'][:,:,2]
                else:
                    self._alpha = f['ab_alpha_percentiles'][0,:,2]
                    self._beta = f['ab_beta_percentiles'][0,:,2]

            if bounds == True:
                self._g_min = 0.0
                self._g_max = 25.0
            else:
                self._g_min = -np.inf
                self._g_max = np.inf

            t_sf = time()

        # Create interpolator
        print('Creating selection function interpolator...')
        if version == 'modelT':
            if crowding == True:
                self._theta_interpolator = interpolate.RectBivariateSpline(self._log10_rho_grid,self._g_grid,self._theta)
                self._interpolator = lambda _log10_rho, _g : (self._theta_interpolator(_log10_rho, _g, grid = False),)
            else:
                self._theta_interpolator = interpolate.interp1d(self._g_grid,self._theta,kind='linear',fill_value='extrapolate')
                self._interpolator = lambda _g : (self._theta_interpolator(_g),)
        elif version == 'modelAB':
            if crowding == True:
                self._alpha_interpolator = interpolate.RectBivariateSpline(self._log10_rho_grid,self._g_grid,self._alpha)
                self._beta_interpolator = interpolate.RectBivariateSpline(self._log10_rho_grid,self._g_grid,self._beta)
                self._interpolator = lambda _log10_rho, _g : (self._alpha_interpolator(_log10_rho, _g, grid = False),self._beta_interpolator(_log10_rho, _g, grid = False))
            else:
                self._alpha_interpolator = interpolate.interp1d(self._g_grid,self._alpha,kind='linear',fill_value='extrapolate')
                self._beta_interpolator = interpolate.interp1d(self._g_grid,self._beta,kind='linear',fill_value='extrapolate')
                self._interpolator = lambda _g : (self._alpha_interpolator(_g),self._beta_interpolator(_g))

        t_interpolator = time()

        t_finish = time()

        print('t = {:.3f} s'.format(t_finish - t_start))
        print('  auxilliary: {: >7.3f} s'.format(t_auxilliary-t_start))
        print('          sf: {: >7.3f} s'.format(t_sf-t_auxilliary))
        print('interpolator: {: >7.3f} s'.format(t_interpolator-t_sf))

    def _selection_function(self,_n,_parameters):
        if len(_parameters) == 1:

            # This must be Model T, _parameters = (theta)
            _t = _parameters[0]

            # 0 < theta < 1, make it so!
            _t[_t<0.0] = 1e-6
            _t[_t>1.0] = 1-1e-6

            _result = special.betainc(5,_n-4,_t)

        elif len(_parameters) == 2:

            # This must be Model AB, _parameters = (alpha,beta)
            _a, _b = _parameters

            # 0.1 < alpha,beta < 10000, make it so!
            _a[_a<1e-1] = 1e-1
            _a[_a>1e+4] = 1e+4
            _b[_b<1e-1] = 1e-1
            _b[_b>1e+4] = 1e+4

            _result = np.ones(_a.shape)
            for _m in range(1,5)[::-1]:
                _result = 1.0 + _result*((_n-_m+1)/_m)*(_a+_m-1)/(_b+_n-_m)
            _result = 1.0 - _result*special.beta(_a,_b+_n)/special.beta(_a,_b)

        return _result


    @ensure_flat_icrs
    @ensure_gaia_g
    def query(self, sources):
        """
        Returns the selection function at the requested coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            Selection function at the specified coordinates, as a fraction.

        """

        # Convert coordinates to healpix indices
        hpxidx = coord2healpix(sources.coord, 'icrs', self._nside, nest=True)

        # Calculate the number of observations of each source
        n = self._n_field[hpxidx]

        # Extract Gaia G magnitude
        G = sources.photometry.measurement['gaia_g']

        if self._crowding == True:

            # Work out HEALPix index in crowding nside
            hpxidx_crowding = np.floor(hpxidx * hp.nside2npix(self._nside_crowding) / hp.nside2npix(self._nside)).astype(np.int)

            # Calculate the local density field at each source
            log10_rho = self._log10_rho_field[hpxidx_crowding]

            # Calculate parameters
            sf_parameters = self._interpolator(log10_rho,G)

        else:

            # Calculate parameters
            sf_parameters = self._interpolator(G)

        # Evaluate selection function
        selection_function = self._selection_function(n,sf_parameters)

        if self._bounds == True:
            _outside_bounds = np.where( (G<self._g_min) | (G>self._g_max) )
            selection_function[_outside_bounds] = 0.0

        return selection_function


class dr3_sf(dr2_sf):
    def __init__(self, map_fname_dr3=None,map_fname_dr2=None, version='modelAB', crowding=False, bounds=True):

        if map_fname_dr3 is None:
            map_fname_dr3 = os.path.join(data_dir(), 'cog_ii', 'n_field_dr3.h5')

        dr2_sf.__init__(self, map_fname_dr2, version, crowding, bounds)

        with h5py.File(map_fname_dr3, 'r') as f:
            # Load auxilliary data
            print('Loading auxilliary data ...')
            self._n_field = f['n_field'][...]
            self._nside = hp.npix2nside(self._n_field.shape[0])


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

    doi = {'dr2': '10.7910/DVN/PDFOVC',
           'edr3_nfield': '10.7910/DVN/I3TGTS'}

    requirements = {'dr2': {'filename': 'cog_ii_dr2.h5'},
                    'edr3_nfield':{'filename': 'n_field_dr3.h5'}}

    local_fname = os.path.join(data_dir(), 'cog_ii', requirements['dr2']['filename'])
    # Download the data
    fetch_utils.dataverse_download_doi(
        doi['dr2'],
        local_fname,
        file_requirements=requirements['dr2'])

    local_fname = os.path.join(data_dir(), 'cog_ii', requirements['edr3_nfield']['filename'])
    # Download the data
    fetch_utils.dataverse_download_doi(
        doi['edr3_nfield'],
        local_fname,
        file_requirements=requirements['edr3_nfield'])
