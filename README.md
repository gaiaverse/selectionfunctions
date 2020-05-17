[![DOI](http://joss.theoj.org/papers/10.21105/joss.00695/status.svg)](https://doi.org/10.21105/joss.00695)

selectionfunctions
==================

The ``selectionfunctions`` package aspires to provide a uniform interface to the selection functions of the major astronomical surveys.
This package is entirely derivative of the truly excellent ``dustmaps`` package created by Gregory M. Green.
The ``selectionfunctions`` package is a product of the [Completeness of the *Gaia*-verse (CoG)](https://www.gaiaverse.space/) collaboration.

Supported Selection Functions
-----------------------------

The currently supported selection functions are:

1. Gaia DR2 source catalogue (cog_ii.dr2_sf, Boubert & Everall 2020, submitted)

To request addition of another selection function in this package, [file an issue on
GitHub](https://github.com/gaiaverse/selectionfunctions/issues), or submit a pull request.


Installation
------------

Download the repository from [GitHub](https://github.com/gaiaverse/selectionfunctions) and
then run:

    python setup.py install --large-data-dir=/path/where/you/want/large/data/files/stored

Alternatively, you can use the Python package manager `pip`:

    pip install selectionfunctions


Getting the Data
----------------

To fetch the data for the GaiaDR2 selectionfunction, run:

    python setup.py fetch --map-name=cog_ii

You can download the other selection functions by changing "cog_ii" to (other selection functions will be added in future).

Alternatively, if you have used `pip` to install `selectionfunctions`, then you can
configure the data directory and download the data by opening up a python
interpreter and running:

    >>> from selectionfunctions.config import config
    >>> config['data_dir'] = '/path/where/you/want/large/data/files/stored'
    >>>
    >>> import selectionfunctions.cog_ii
    >>> selectionfunctions.cog_ii.fetch()


Querying the selection functions
-----------------

Selection functions are queried using Source objects, which are a variant on the 
[`astropy.coordinates.SkyCoord`](http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord)
object. This means that any coordinate system supported by `astropy` can be
used as input. For example, we can query the Gaia DR2 selection function as follows:

    >>> import selectionfunctions.cog_ii as CoGII
    >>> from selectionfunctions.source import Source
    >>>
    >>> dr2_sf = CoGII.dr2_sf()
    >>>
    >>> c = Source(
            '22h54m51.68s',
            '-51d11m44.19s',
            photometry={'gaia_g':16.02},
            frame='icrs')
    >>> print(dr2_sf(c))


Above, we have used the ICRS coordinate system (the inputs are RA and Dec). We
can use other coordinate systems, such as Galactic coordinates, and we can
provide coordinate arrays. The following example uses both features:

    >>> c = Source(
            [75.00000000, 130.00000000],
            [-89.00000000, 10.00000000],
            photometry={'gaia_g':[2.3,17.8]},
            frame='galactic',
            unit='deg')
    >>> print(dr2_sf(c))



Documentation
-------------

Read the full documentation at http://selectionfunctions.readthedocs.io/en/latest/.


Citation
--------

If you make use of this software in a publication, please always cite
[Green (2018) in The Journal of Open Source Software](https://doi.org/10.21105/joss.00695).

You should also cite the papers behind the selection functions you use.

1. cog_ii.dr2_sf - Please cite Completeness of the Gaia-verse [Paper I](https://ui.adsabs.harvard.edu/abs/2020arXiv200414433B/abstract) and Paper II.

Development
-----------

Development of `selectionfunctions` takes place on GitHub, at
https://github.com/gaiaverse/selectionfunctions. Any bugs, feature requests, pull requests,
or other issues can be filed there. Contributions to the software are welcome.
