#!/usr/bin/env python
#
# healpix.py
# Implements a selection function that only varies with HEALPix index.
#
# Copyright (C) 2020  Berry Holl, Douglas Boubert & Andrew Everall
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
import healpy as hp

from .std_paths import *
from .map import SelectionFunction, ensure_flat_icrs, coord2healpix
from . import fetch_utils

from time import time


class healpix_sf(SelectionFunction):
    """
    Implements a general selection function based on a HEALPix grid.
    """

    def __init__(self, nside, indices=None, values=None, nest=True, coord='C', interp=False, default_value=1.0):
        """
        Args:
            nside ([:obj:`int`]): The nside of the HEALPix map that :obj:`'indices'` and/or :obj:`'values'` refer to.
            indices ([:obj:`int`]): Indices at which the :obj:`'values'` are true
            values (Optional[:obj:`float`]): The value of the selection function at the specified indices.
            nest (Optional[:obj:`bool`]): Whether the HEALPix map is nested.
                Defaults to :obj:`'True'`.
            coord (Optional[:obj:`str`]): The coordinate frame of the HEALPix map.
                Defaults to :obj:`'C'`.
            interp (Optional[:obj:`bool`]): Whether or not the HEALPix map will be bilinearly interpolated when it is queried.
                Defaults to :obj:`'False'`.
            default_value (Optional[:obj:`float`]): The value of the HEALPix map at indices that you have not specified a value for.
                Defaults to :obj:`'1.0'`.
        """

        # Verify if nside is power of two
        if ((nside & (nside-1) == 0) and nside != 0) == False:
            raise ValueError("The value of nside must be a power of two.")
        else:
            npix = hp.nside2npix(nside)
            self.nside = nside

        # If no values are specific, default to one everywhere.
        if values is None:
            values = default_value*np.ones(npix)
        else:
            values = np.asarray(values)

        # If no indices, then set the same value everywhere.
        if indices is None:
            if values.size == npix:
                self.values = values
            elif values.size == 1:
                self.values = values[0]*np.ones(npix)
            else:
                raise ValueError("You have not specified any indices and so the size of the parameter values must be equal to either one or the number of pixels in the HEALPix map.")
        else:
            indices = np.asarray(indices)

            if indices.size == values.size or values.size == 1:
                self.values = default_value*np.ones(npix)
                self.values[indices] = values
            else:
                raise ValueError("You have specified indices and values, but they are not of compatible lengths. The parameter values must either of the same size as indices or of length one.")

        self.nest = nest
        self.coord = coord
        self.interp = interp
        if self.coord == 'C':
            self.frame = 'icrs'
        elif self.coord == 'G':
            self.frame = 'galactic'
        elif self.coord == 'E':
            self.frame = 'barycentricmeanecliptic'
        else:
            raise ValueError("The parameter coord must be either C(elestial), G(alactic) or E(cliptic).")


    @ensure_flat_icrs
    def query(self, sources):
        """
        Returns the selection function at the requested coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            Selection function at the specified coordinates, as a fraction.

        """

        # Evaluate selection function
        if self.interp == False:
            # Convert coordinates to healpix indices
            hpxidx = coord2healpix(sources.coord, self.frame, self.nside, nest=self.nest)
            selection_function = self.values[hpxidx]
        else:
            if self.coord == 'C':
                lon,lat = sources.coord.ra.deg,sources.coord.dec.deg
            elif self.coord == 'G':
                lon,lat = sources.coord.galactic.l.deg,sources.coord.galactic.b.deg
            elif self.coord == 'E':
                lon,lat = coords.coord.barycentricmeanecliptic.lon.deg,coords.coord.barycentricmeanecliptic.lat.deg
            else:
                raise ValueError("The parameter coord must be either C(elestial), G(alactic) or E(cliptic).")

            selection_function = hp.get_interp_val(self.values, lon, lat, nest=self.nest, lonlat=True)

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

    requirements = {'filename': 'cog_ii_dr2.h5'}

    local_fname = os.path.join(data_dir(), 'cog_ii', 'cog_ii_dr2.h5')

    # Download the data
    fetch_utils.dataverse_download_doi(
        doi,
        local_fname,
        file_requirements=requirements)
