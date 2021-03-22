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
import tqdm
import healpy as hp
from scipy import interpolate, special

from .std_paths import *
from .map import SelectionFunction, ensure_flat_icrs, coord2healpix
import selectionfunctions.SelectionFunctionUtils as SelectionFunctionUtils
from .source import ensure_gaia_g
from . import fetch_utils

from time import time
from numba import njit
from scipy import sparse


class chisel(SelectionFunction):
    """
    Queries the Gaia DR2 selection function (Boubert & Everall, 2019).
    """
    basis_keyword='wavelet'

    def __init__(self, map_fname=None, bounds=True, basis_options={},
                       nside=128, lmax=100, M=17, C=1, lengthscale=0.3,
                       spherical_basis_directory='./SphericalBasis'):
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

        # Utilities
        self.order_to_nside = lambda order: 2**order
        self.nside_to_npix = lambda nside: 12*nside**2
        self.order_to_npix = lambda order: self.nside_to_npix(self.order_to_nside(order))


        map_fname = os.path.join(data_dir(), map_fname)

        t_start = time()

        self.nside=nside
        self.M=M
        self.C=C
        self.lmax=lmax
        self._bounds = bounds
        self.lengthscale = lengthscale

        self.L, self.H, self.R = 2 * self.lmax + 1, (self.lmax + 1) ** 2, 4 * self.nside - 1

        with h5py.File(map_fname, 'r') as f:
            # Load auxilliary data
            print('Loading auxilliary data ...')
            self.x = f['x'][...]
            self.b = f['b'][...]
            self.Mlim = f['Mlim'][...]
            self.Clim = f['Clim'][...]
        self.H = self.b.shape[0]
        t_auxilliary = time()

        if bounds == True:
            self._g_min = 5.0
            self._g_max = 22.0
        else:
            self._g_min = -np.inf
            self._g_max = np.inf

        t_sf = time()

        # Process basis-specific options
        self._process_basis_options(**basis_options)

        # Load spherical basis
        self.spherical_basis_directory=spherical_basis_directory
        self._load_spherical_basis()
        self.spherical_basis_file = f"{self.basis_keyword}_{self.needlet}_nside{self.nside}_B{self.B}_"+ (f"p{self.p}_" if self.needlet == 'chisquare' else '') + f"tol{self.wavelet_tol}_j[{','.join([str(_i) for _i in self.j])}].h5"

        # Initialise covariance kernel
        Mbins = np.linspace(self.Mlim[0],self.Mlim[1], M+1)
        self.Mcenters = (Mbins[1:]+Mbins[:-1])/2
        self._inv_KMM = np.linalg.inv(self.covariance_kernel(self.Mcenters, self.Mcenters, lengthscale=lengthscale) + 1e-15*np.eye(M))
        Cbins = np.linspace(self.Clim[0],self.Clim[1], C+1)
        self.Ccenters = (Cbins[1:]+Cbins[:-1])/2
        self._inv_KCC = np.linalg.inv(self.covariance_kernel(self.Ccenters, self.Ccenters, lengthscale=lengthscale) + 1e-15*np.eye(C))

        t_interpolator = time()

        t_finish = time()

        print('t = {:.3f} s'.format(t_finish - t_start))
        print('  auxilliary: {: >7.3f} s'.format(t_auxilliary-t_start))
        print('          sf: {: >7.3f} s'.format(t_sf-t_auxilliary))
        print('interpolator: {: >7.3f} s'.format(t_interpolator-t_sf))

    @ensure_flat_icrs
    @ensure_gaia_g
    def query(self, sources, chunksize=1000):
        """
        Returns the selection function at the requested coordinates.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.

        Returns:
            Selection function at the specified coordinates, as a fraction.

        """

        # Convert coordinates to healpix indices
        hpxidx = coord2healpix(sources.coord, 'icrs', self.nside, nest=True)

        # Extract Gaia G magnitude
        mag = sources.photometry.measurement['gaia_g']
        try: color = sources.photometry.measurement['gaia_bp_gaia_rp']
        except KeyError: color = np.zeros(len(mag))

        # Switch positions to ring ordering
        pix = hp.nest2ring(self.nside, hpxidx)

        # Evaluate selection function
        selection_function = self._selection_function(mag, color, pix)

        if self._bounds == True:
            _outside_bounds = np.where( (mag<self._g_min) | (mag>self._g_max) )
            selection_function[_outside_bounds] = 0.0

        return selection_function

    def _selection_function_pixel(self, mag, color, hpx):

        # Get mag and color IDs
        M_idx = ( self.M*(mag  -self.Mlim[0])/(self.Mlim[1]-self.Mlim[0]) ).astype(int)
        C_idx = ( self.C*(color-self.Clim[0])/(self.Clim[1]-self.Clim[0]) ).astype(int)

        # Load logit probability in pixels
        x = self.x[M_idx, C_idx, hpx]

        # Take expit
        p = self.expit(x)

        return p

    def _selection_function_bulky(self, mag, color, pix, chunksize=1000):

        # Batch up iterations:
        x = np.zeros(mag.shape)

        # Load in sparse matrix
        nmodes = 0
        for j in self.j:
            if j==-1: nmodes += 1
            else: nmodes += 12*(2**j)**2
        npix = self.nside_to_npix(self.nside)

        Y = sparse.csr_matrix((self.basis['wavelet_w'],self.basis['wavelet_v'],self.basis['wavelet_u']), shape=(npix, nmodes)).toarray()
        #print(Y.shape)
        Y = sparse.csr_matrix((self.basis['wavelet_w'],self.basis['wavelet_v'],self.basis['wavelet_u']), shape=(npix, nmodes)).toarray()[pix]

        for ii in tqdm.tqdm_notebook(range(x.shape[0]//chunksize + 1)):

            # Contstruct covariance kernel for new positions.
            KmM = self.covariance_kernel(mag[ii*chunksize:(ii+1)*chunksize], self.Mcenters, lengthscale=self.lengthscale)
            KcC = self.covariance_kernel(color[ii*chunksize:(ii+1)*chunksize], self.Ccenters, lengthscale=self.lengthscale)

            #print(KmM.shape, self._inv_KMM.shape, self.b.shape, KcC.shape, self._inv_KCC.shape)

            # Estimate alm using Gaussian Process
            _b_m = np.sum ( ((KmM @ self._inv_KMM) @ self.b) * (KcC @ self._inv_KCC)[None, :,:] , axis=2).T

            #print(_b_m.shape, Y.shape)

            #x[ii*chunksize:(ii+1)*chunksize] = Y.dot(_b_m.T)
            x[ii*chunksize:(ii+1)*chunksize] = np.sum(Y[ii*chunksize:(ii+1)*chunksize] * _b_m, axis=1)
            # x[m,c] = csr_matrix_times_vector(P, S, wavelet_w, wavelet_v, wavelet_u, _b_m);

        # Take expit
        p = self.expit(x)

        return p

    def _selection_function(self, mag, color, pix):

        # Load in sparse matrix
        nmodes = 0
        for j in self.j:
            if j==-1: nmodes += 1
            else: nmodes += 12*(2**j)**2
        npix = self.nside_to_npix(self.nside)
        n = len(mag)
        x = np.zeros(n)

        @njit
        def matrix_multiply(x, b, KM, KC, n, H, M, C, wavelet_w, wavelet_v, wavelet_u, pix):

            # Iterate over pixels
            for i, ipix in enumerate(pix):
                # Iterate over modes which are not sparsified in Y
                for iY, iH in enumerate(wavelet_v[wavelet_u[ipix]:wavelet_u[ipix+1]]):
                    x[i] += np.dot(np.dot(KM[i], b[iH]), KC[i]) * wavelet_w[wavelet_u[ipix]+iY]

        # Contstruct covariance kernel for new positions.
        KmM = self.covariance_kernel(mag, self.Mcenters, lengthscale=self.lengthscale)
        KcC = self.covariance_kernel(color, self.Ccenters, lengthscale=self.lengthscale)

        matrix_multiply(x, self.b, (KmM @ self._inv_KMM), (KcC @ self._inv_KCC), n, self.H, self.M, self.C,\
                  self.basis['wavelet_w'], self.basis['wavelet_v'], self.basis['wavelet_u'], pix)

        # Take expit
        p = self.expit(x)

        return p

    def covariance_kernel(self, x1, x2, lengthscale=None):

        # Initialise covariance kernel
        C = np.exp(-np.square(x1[:,None]-x2[None,:])/(2.0*lengthscale*lengthscale))

        return C

    def expit(self, x):
        return 1/(1+np.exp(-x))

    def logit(self, p):
        return np.log(p/(1-p))

    def _process_basis_options(self, needlet = 'littlewoodpaley', j=[0], B = 2.0, p = 1.0, wavelet_tol = 1e-10):

        if type(j) in [list,tuple,np.ndarray]:
            self.j = sorted([int(_j) for _j in j])

        else:
            self.j = [_j for _j in range(-1,j+1)]

        self.needlet, self.B, self.p, self.wavelet_tol = needlet, B, p, wavelet_tol

        self.spherical_basis_file = f"{self.basis_keyword}_{self.needlet}_nside{self.nside}_B{self.B}_"+ (f"p{self.p}_" if self.needlet == 'chisquare' else '') + f"tol{self.wavelet_tol}_j[{','.join([str(_i) for _i in self.j])}].h5"

        self.S = sum([self.order_to_npix(_j) if _j >= 0 else 1 for _j in self.j])

        assert self.B > 1.0
        assert self.wavelet_tol >= 0.0
        assert self.needlet in ['littlewoodpaley','chisquare']
        if self.needlet == 'chisquare':
            from selectionfunctions.SelectionFunctionUtils import chisquare
            self.weighting = chisquare(self.j, p = self.p, B = self.B)
        else:
            from selectionfunctions.SelectionFunctionUtils import littlewoodpaley
            self.weighting = littlewoodpaley(B = self.B)



    def _load_spherical_basis(self):
        """ Loads in the spherical basis file. If they don't exist, then generate them. The generator must be implemented in each child class. """

        _spherical_basis_files = self.spherical_basis_file

        if not os.path.isfile(self.spherical_basis_directory + self.spherical_basis_file):
            print('Spherical basis file does not exist, generating... (this may take some time!)')
            self._generate_spherical_basis(self.spherical_basis_directory + self.spherical_basis_file)

        # Load spherical wavelets
        with h5py.File(self.spherical_basis_directory + self.spherical_basis_file, 'r') as sbf:
            self.basis = {k:v[()] for k,v in sbf.items()}

        print('Spherical basis file loaded')

    def _generate_spherical_basis(self,gsb_file):

        # Import dependencies
        from numba import njit
        from math import sin, cos
        import sys

        nside = self.nside
        B = self.B
        needle_sparse_tol = self.wavelet_tol

        # Function to compute needlet across sky
        @njit
        def pixel_space (Y, cos_gamma, window, start, end, legendre):
            '''Return the value of a needlet at gamma radians from the needlet centre.'''

            legendre[0] = 1.0
            legendre[1] = cos_gamma
            for cur_l in range(2, end + 1):
                legendre[cur_l] = ((cos_gamma * (2 * cur_l - 1) * legendre[cur_l - 1] - (cur_l - 1) * legendre[cur_l - 2])) / cur_l

            Y[:] = np.dot(window,legendre[start:end+1])

        # Compute locations of pixels
        npix = self.nside_to_npix(nside)
        colat, lon = np.array(hp.pix2ang(nside=nside,ipix=np.arange(npix),lonlat=False))
        cos_colat, sin_colat = np.cos(colat), np.sin(colat)
        cos_lon, sin_lon = np.cos(lon), np.sin(lon)

        # Initialise variables
        running_index = 0
        needlet_w, needlet_v, needlet_u, needlet_j = [], [], [], []
        Y = np.zeros(npix)
        legendre = np.zeros((1+self.weighting.end(max(self.j)),npix))

        for j in self.j:

            print(f'Working on order {j}.')

            if j == -1:
                needlet_w.append(np.ones(npix))
                needlet_v.append(np.arange(npix))
                needlet_u.append(0)
                needlet_j.append(np.zeros(1))
                running_index += npix
                continue

            nside_needle = self.order_to_nside(j)
            npix_needle = self.nside_to_npix(nside_needle)

            start = self.weighting.start(j)
            end = self.weighting.end(j)
            modes = np.arange(start, end + 1, dtype = 'float')
            window = self.weighting.window_function(modes,j)*(2.0*modes+1.0)/np.sqrt(4.0*np.pi*npix_needle)

            for ipix_needle in tqdm.tqdm(range(npix_needle),file=sys.stdout):

                colat_needle, lon_needle = hp.pix2ang(nside=nside_needle,ipix=ipix_needle,lonlat=False)

                cos_gamma = cos(colat_needle) * cos_colat + sin(colat_needle) * sin_colat * (cos(lon_needle) * cos_lon + sin(lon_needle) * sin_lon)

                pixel_space(Y, cos_gamma = cos_gamma, window = window, start = start, end = end, legendre = legendre)

                _significant = np.where(np.abs(Y) > Y.max()*needle_sparse_tol)[0]
                needlet_w.append(Y[_significant])
                needlet_v.append(_significant)
                needlet_u.append(running_index)
                needlet_j.append(j*np.ones(self.order_to_npix(j)))
                running_index += _significant.size

        # Add the ending index to u
        needlet_u.append(running_index)

        # Concatenate the lists
        needlet_w = np.concatenate(needlet_w)
        needlet_v = np.concatenate(needlet_v)
        needlet_u = np.array(needlet_u)

        # Flip them round
        from scipy import sparse
        Y = sparse.csr_matrix((needlet_w,needlet_v,needlet_u)).transpose().tocsr()
        wavelet_w, wavelet_v, wavelet_u = Y.data, Y.indices, Y.indptr
        wavelet_j = np.concatenate(needlet_j).astype(int)
        wavelet_n = wavelet_w.size

        # Save file
        save_kwargs = {'compression':"lzf", 'chunks':True, 'fletcher32':False, 'shuffle':True}
        with h5py.File(gsb_file, 'w') as f:
            f.create_dataset('wavelet_w', data = wavelet_w, dtype = np.float64, **save_kwargs)
            f.create_dataset('wavelet_v', data = wavelet_v, dtype = np.uint64, scaleoffset=0, **save_kwargs)
            f.create_dataset('wavelet_u', data = wavelet_u, dtype = np.uint64, scaleoffset=0, **save_kwargs)
            f.create_dataset('wavelet_n', data = wavelet_n)
            f.create_dataset('modes', data = wavelet_j, dtype = np.uint64, scaleoffset=0, **save_kwargs)


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
