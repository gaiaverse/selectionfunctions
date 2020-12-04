#!/usr/bin/env python
#
# rrl_mateu_2020.py
# Implements querying of the completeness maps from
# Mateu, Holl, De Ridder & Rimoldini (2020).
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
import pandas as pd

from .std_paths import *
from .map import SelectionFunction, ensure_flat_galactic, coord2healpix
from .source import ensure_gaia_g, ensure_distance
from . import fetch_utils

from time import time


class rrl_sf(SelectionFunction):
    """
    Queries the Gaia DR2 selection function (Boubert & Everall, 2019).
    """

    def __init__(self, survey='gaiadr2', map_dimension='2d_bright', rrl_type='rrab', bounds_value=0.0, null_value=np.nan, map_fname=None):
        """
        Args:
            survey (Optional[:obj:`str`]): Which survey to return the completeness of. Valid surveys
                are :obj:`'gaiadr2'`, :obj:`'asas'` and :obj:`'ps1'`.
                Defaults to :obj:`'gaiadr2'`.
            map_dimension (Optional[:obj:`str`]): Whether to use the :obj:`'2d_bright'`, :obj:`'2d_faint'` or :obj:`'3d'` map.
                Defaults to :obj:`'2d_bright'`.
            rrl_type (Optional[:obj:`str`]): Whether to return the completeness of :obj:`'rrab'` or :obj:`'rrc'` RRL stars.
                Defaults to :obj:`'ab'`.
            bounds_value (Optional[:obj:`float`]): What value to use for the completeness of RRL that fall outside the map grid.
                Defaults to :obj:`'0.0'`.
            null_value (Optional[:obj:`float`]): What value to use for the completeness of RRL that fall outside the map grid.
                Defaults to :obj:`'0.0'`.
            map_fname (Optional[:obj:`str`]): Filename of the completeness map to use overriding the other keywords. Defaults to
                :obj:`None`.
        """

        # Identify which map we are loading
        if map_fname is None:

            # Validate input
            if survey == 'gaiadr2' or survey == 'asas' or survey == 'ps1':
                self._survey = survey
            else:
                raise ValueError('survey must be either "gaiadr2", "asas" or "ps1".')

            if map_dimension[:2] == '2d' and (map_dimension[3:] == 'bright' or map_dimension[3:] == 'faint'):
                self._dimension = '2d'
                self._bright_or_faint = map_dimension[3:]
                self._fname = 'completeness2d'+'.'+self._bright_or_faint
            elif map_dimension == '3d':
                self._dimension = '3d'
                self._fname = 'completeness3d'+'.'+self._survey
                if self._survey == 'gaiadr2':
                    self._fname = self._fname + '.vcsos'

            else:
                raise ValueError('map_dimension must be either "2d_bright", "2d_faint", or "3d".')


            if rrl_type == 'rrab' or rrl_type == 'rrc':
                self._rrl_type = rrl_type
                self._fname = self._fname + '.' + self._rrl_type + '.tab'
            else:
                raise ValueError('rrl_type must be either "rrab" or "rrc".')

            # Verify that our combination is valid
            if self._dimension == '2d' and self._bright_or_faint == 'bright' and self._survey == 'ps1':
                raise ValueError('Only "gaia" and "asas" have a "2d_bright" completeness map.')

            self._null_value = null_value
            map_fname = os.path.join(data_dir(), 'rrl_mateu_2020', self._fname)


        t_start = time()
        
        print('Loading data ...')
        map_data = pd.read_csv(map_fname)

        if self._dimension == '2d':

            if self._bright_or_faint == 'bright':
                self._nside = 4
                self._g_bins = np.array([11.0,13.0,15.0,17.0])

                if self._survey == 'gaiadr2':
                    self._completeness = np.stack([map_data['Gaia[11,13]'],map_data['Gaia[13,15]'],map_data['Gaia[15,17]']])
                elif self._survey == 'asas':
                    self._completeness = np.stack([map_data['ASAS[11,13]'],map_data['ASAS[13,15]'],map_data['ASAS[15,17]']])

            else:
                self._nside = 16
                self._g_bins = np.array([13.0,16.0,18.0,22.0])
                if self._survey == 'gaiadr2':
                    self._completeness = np.stack([map_data['Gaia[13,16]'],map_data['Gaia[16,18]'],map_data['Gaia[18,22]']])
                elif self._survey == 'asas':
                    self._completeness = np.stack([map_data['ASAS[13,16]'],map_data['ASAS[16,18]'],map_data['ASAS[18,22]']])
                elif self._survey == 'ps1':
                    self._completeness = np.stack([map_data['PS1[13,16]'],map_data['PS1[16,18]'],map_data['PS1[18,22]']])

            self._completeness_dict = {_n:{'bins':self._g_bins,'completeness':self._completeness[:,_n]} for _n in range(hp.nside2npix(self._nside))}

        elif self._dimension == '3d':
            self._nside = 4
            self._hpx = map_data['hpix'].values
            self._distance_lower = map_data['D_o'].values
            self._distance_upper = map_data['D_f'].values
            self._completeness = map_data['C'].values

            # Create the dictionary
            self._completeness_dict = {}
            counts = np.unique(self._hpx,return_counts=True)[1]
            _idx = 0
            for _n in range(hp.nside2npix(self._nside)):
                d_bins = np.concatenate([self._distance_lower[_idx:_idx+counts[_n]],[self._distance_upper[_idx+counts[_n]-1]]])
                self._completeness_dict[_n] = {'bins':d_bins,'completeness':self._completeness[_idx:_idx+counts[_n]]}
                _idx += counts[_n]



        t_data = time()
        
        t_finish = time()
        
        print('t = {:.3f} s'.format(t_finish - t_start))
        print('        data: {: >7.3f} s'.format(t_data-t_start))


    @ensure_flat_galactic
    def query(self, sources):
        """
        Returns the selection function at the requested coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            Selection function at the specified coordinates, as a fraction.

        """
        if self._dimension == '2d':
            @ensure_gaia_g
            def _query(self,sources):

                # Convert coordinates to healpix indices
                hpxidx = coord2healpix(sources.coord, 'galactic', self._nside, nest=False)

                # Extract Gaia G magnitude
                x = sources.photometry.measurement['gaia_g']

                n_coords = hpxidx.size
                selection_function = np.zeros(n_coords)
                for _idx in range(n_coords):
                    box = self._completeness_dict[hpxidx[_idx]]
                    _cidx = np.searchsorted(box['bins'],x[_idx])
                    if _cidx == 0 or _cidx == box['bins'].size:
                        selection_function[_idx] = self._null_value
                    else:
                        selection_function[_idx] = box['completeness'][_cidx-1]
                return selection_function

        elif self._dimension == '3d':
            @ensure_distance
            def _query(self,sources):

                # Convert coordinates to healpix indices
                hpxidx = coord2healpix(sources.coord, 'galactic', self._nside, nest=False)

                # Extract distance
                x = sources.coord.distance.kpc

                n_coords = hpxidx.size
                selection_function = np.zeros(n_coords)
                for _idx in range(n_coords):
                    box = self._completeness_dict[hpxidx[_idx]]
                    _cidx = np.searchsorted(box['bins'],x[_idx])
                    if _cidx == 0 or _cidx == box['bins'].size:
                        selection_function[_idx] = self._null_value
                    else:
                        selection_function[_idx] = box['completeness'][_cidx-1]
                return selection_function

        return _query(self,sources)


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

    doi = '10.7910/DVN/HX0RNB'
    for _fname in ['completeness2d.bright.rrab.tab',
                    'completeness2d.bright.rrc.tab',
                    'completeness2d.faint.rrab.tab',
                    'completeness2d.faint.rrc.tab',
                    'completeness3d.asas.rrab.tab',
                    'completeness3d.asas.rrc.tab',
                    'completeness3d.gaiadr2.vcsos.rrab.tab',
                    'completeness3d.gaiadr2.vcsos.rrc.tab',
                    'completeness3d.ps1.rrab.tab',
                    'completeness3d.ps1.rrc.tab']:

        requirements = {'filename': _fname}

        local_fname = os.path.join(data_dir(), 'rrl_mateu_2020', _fname)

        # Download the data
        fetch_utils.dataverse_download_doi(
            doi,
            local_fname,
            original = True,
            file_requirements=requirements)
