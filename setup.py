#!/usr/bin/env python
#
# setup.py
# Package "selectionfunctions" for pip.
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

from setuptools import setup, Extension
from setuptools.command.install import install
import distutils.cmd

import os
import json
import io


class InstallCommand(install):
    description = install.description
    user_options = install.user_options + [
        ('large-data-dir=', None, 'Directory to store large data files in.')
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.large_data_dir = None

    def finalize_options(self):
        if not self.large_data_dir is None:
            self.large_data_dir = os.path.abspath(os.path.expanduser(self.large_data_dir))

        install.finalize_options(self)

    def run(self):
        if not self.large_data_dir is None:
            print('Large data directory is set to: {}'.format(self.large_data_dir))
            with open(os.path.expanduser('~/.selectionfunctionsrc'), 'w') as f:
                json.dump({'data_dir': self.large_data_dir}, f, indent=2)

        # install.do_egg_install(self) # Due to bug in setuptools that causes old-style install
        install.run(self)


def fetch_cog_ii():
    import selectionfunctions.cog_ii
    selectionfunctions.cog_ii.fetch()


class FetchCommand(distutils.cmd.Command):
    description = ('Fetch selection functions from the web, and store them in the data '
                   'directory.')
    user_options = [
        ('map-name=', None, 'Which selection functions to load.')]

    map_funcs = {
        'cog_ii': fetch_cog_ii,
    }

    def initialize_options(self):
        self.map_name = None

    def finalize_options(self):
        try:
            import selectionfunctions
        except ImportError:
            print('You must install the package selectionfunctions before running the '
                  'fetch command.')
        if not self.map_name in self.map_funcs:
            print('Valid map names are: {}'.format(self.map_funcs.keys()))

    def run(self):
        print('Fetching map: {}'.format(self.map_name))
        self.map_funcs[self.map_name]()


def readme():
    with io.open('README.md', mode='r', encoding='utf-8') as f:
        return f.read()


setup(
    name='selectionfunctions',
    version='1.0.0',
    description='Uniform interface for the selection functions of astronomical surveys.',
    long_description=readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/gaiaverse/selectionfunctions',
    download_url='https://github.com/gaiaverse/selectionfunctions/archive/v1.0.0.tar.gz',
    author='Douglas Boubert',
    author_email='douglasboubert@gmail.com',
    license='GPLv2',
    packages=['selectionfunctions'],
    install_requires=[
        'numpy',
        'scipy',
        'astropy',
        'h5py',
        'healpy',
        'requests',
        'progressbar2',
        'six'
    ],
    include_package_data=True,
    test_suite='nose.collector',
    tests_require=['nose'],
    zip_safe=False,
    cmdclass = {
        'install': InstallCommand,
        'fetch': FetchCommand,
    },
)
