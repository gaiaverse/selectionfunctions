#!/usr/bin/env python
#
# source.py
# Provides a new class that extends astropy SkyCoord.
#
# Copyright (C) 2020  Douglas Boubert
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

class Photometry():
    def __init__(self,photometry,photometry_error):
        self.measurement = {k:v for k,v in photometry.items()}
        self.error = {k:v for k,v in photometry_error.items()} if photometry_error is not None else None

class Source():
    def __init__(self,*args,photometry=None,photometry_error=None,**kwargs):
        self.coord = coordinates.SkyCoord(*args,**kwargs)
        self.photometry = Photometry(photometry,photometry_error) if photometry is not None else None

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
    def _wrapper_func(self, sources, **kwargs):
        # t0 = time.time()

        has_photometry = hasattr(sources, 'photometry')
        if has_photometry:
            has_gaia_g = 'gaia_g' in sources.photometry.measurement.keys()
            if has_gaia_g:
                pass
            else:
                print('No Gaia G passed, but transformation is not yet implemented.')
                raise ValueError('You need to pass in Gaia G-band photometric magnitudes to use this selection function.')
        else:
            raise ValueError('You need to pass in Gaia G-band photometric magnitudes to use this selection function.')

        out = f(self, sources, **kwargs)

        return out

    return _wrapper_func


def ensure_tmass_hjk(f):
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
    def _wrapper_func(self, sources, **kwargs):
        # t0 = time.time()

        has_photometry = hasattr(sources, 'photometry')
        if has_photometry:
            has_tmass_h = 'tmass_h' in sources.photometry.measurement.keys()
            has_tmass_jk = 'tmass_jk' in sources.photometry.measurement.keys()
            if has_tmass_h: pass
            else:
                print('No 2MASS H magnitude passed (tmass_h), but transformation is not yet implemented.')
                raise ValueError('You need to pass in 2MASS H photometric magnitudes to use this selection function.')
            if has_tmass_jk: pass
            else:
                print('No 2MASS J-K colour passed (tmass_jk), but transformation is not yet implemented.')
                raise ValueError('You need to pass in 2MASS J-K colour to use this selection function.')
        else:
            raise ValueError('You need to pass in 2MASS H photometric magnitudes and J-K colours to use this selection function.')

        out = f(self, sources, **kwargs)

        return out

    return _wrapper_func
