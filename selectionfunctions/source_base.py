#!/usr/bin/env python
#
# source_base.py
# Provides a new class that extends astropy SkyCoord.
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

from __future__ import print_function, division

import astropy.coordinates as coordinates
import astropy.units as units

from functools import wraps

class SourceCoord(coordinates.SkyCoord):
    def __init__(self,*args,photometry=None,photometry_error=None,**kwargs):
        coordinates.SkyCoord.__init__(self,*args,**kwargs)
        self.photometry = photometry
        self.photometry_error = photometry_error

    def transform_to(self,frame):
        _return_sourcecoord = coordinates.SkyCoord.transform_to(self,frame)
        _return_sourcecoord.photometry = self.photometry
        _return_sourcecoord.photometry_error = self.photometry_error
        return _return_sourcecoord

def ensure_gaia_g(f):
    """
    A decorator for class methods of the form

    .. code-block:: python

        Class.method(self, coords, **kwargs)

    where ``coords`` is an :obj:`astropy.coordinates.SkyCoord` object.

    The decorator ensures that the ``coords`` that gets passed to
    ``Class.method`` is a flat array of Equatorial coordinates. It also reshapes
    the output of ``Class.method`` to have the same shape (possibly scalar) as
    the input ``coords``. If the output of ``Class.method`` is a tuple or list
    (instead of an array), each element in the output is reshaped instead.

    Args:
        f (class method): A function with the signature
            ``(self, coords, **kwargs)``, where ``coords`` is a :obj:`SkyCoord`
            object containing an array.

    Returns:
        A function that takes :obj:`SkyCoord` input with any shape (including
        scalar).
    """

    @wraps(f)
    def _wrapper_func(self, coords, **kwargs):
        # t0 = time.time()

        has_photometry = hasattr(coords, 'photometry')
        has_photometry_error = hasattr(coords, 'photometry_error')
        if has_photometry:
            print(coords.photometry)
            has_gaia_g = 'gaia_g' in coords.photometry.keys()
            if has_gaia_g:
                print('Gaia G magnitude was passed.')
            else:
                print('No Gaia G passed, but transformation is not yet implemented.')
                raise ValueError('You need to pass in Gaia G-band photometric magnitudes to use this selection function.')
        else:
            raise ValueError('You need to pass in Gaia G-band photometric magnitudes to use this selection function.')

        out = f(self, coords, **kwargs)

        return out

    return _wrapper_func