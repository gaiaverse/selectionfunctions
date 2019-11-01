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

class SourceCoord(coordinates.SkyCoord):
    def __init__(self,*args,photometry=None,photometry_errors=None,**kwargs):
        coordinates.SkyCoord.__init__(self,*args,**kwargs)
        self.photometry = photometry
        self.photometry_errors = photometry_errors