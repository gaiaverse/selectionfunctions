[![DOI](http://joss.theoj.org/papers/10.21105/joss.00695/status.svg)](https://doi.org/10.21105/joss.00695)

selectionfunctions
========

The ``selectionfunctions`` package provides a uniform interface to the selection functions of the major astronomical surveys. This package is entirely derivative of the truly excellent ``dustmaps`` package created by Gregory M. Green. Please cite the ``dustmaps`` JOSS paper when using this package.

Supported Selection Functions
-------------------

The currently supported selection functions are:

1. Boubert & Everall (2019; GaiaDR2, in prep.)

To request addition of another selection function in this package, [file an issue on
GitHub](https://github.com/DouglasBoubert/selectionfunctions/issues), or submit a pull request.


Installation
------------

Download the repository from [GitHub](https://github.com/DouglasBoubert/selectionfunctions) and
then run:

    python setup.py install --large-data-dir=/path/where/you/want/large/data/files/stored

Alternatively, you can use the Python package manager `pip`:

    pip install dustmaps


Getting the Data
----------------

To fetch the data for the GaiaDR2 selectionfunction, run:

    python setup.py fetch --map-name=gaiadr2

You can download the other selection functions by changing "gaiadr2" to (other selection functions to be added in future).

Alternatively, if you have used `pip` to install `selectionfunctions` (not yet available), then you can
configure the data directory and download the data by opening up a python
interpreter and running:

    >>> from selectionfunctions.config import config
    >>> config['data_dir'] = '/path/where/you/want/large/data/files/stored'
    >>>
    >>> import selectionfunctions.gaiadr2
    >>> selectionfunctions.gaiadr2.fetch()


Querying the Maps
-----------------

Maps are queried using
[`astropy.coordinates.SkyCoord`](http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord)
objects. This means that any coordinate system supported by `astropy` can be
used as input. For example, we can query GaiaDR2 as follows:

    >>> from selectionfunctions.gaiadr2 import GaiaDR2Query
    >>> from astropy.coordinates import SkyCoord
    >>>
    >>> gaiadr2 = GaiaDR2Query()
    >>>
    >>> c = SkyCoord(
            '05h00m00.00000s',
            '+30d00m00.0000s',
            frame='icrs')
    >>> print gaiadr2(c)

Above, we have used the ICRS coordinate system (the inputs are RA and Dec). We
can use other coordinate systems, such as Galactic coordinates, and we can
provide coordinate arrays. The following example uses both features:

    >>> c = SkyCoord(
            [75.00000000, 130.00000000],
            [-89.00000000, 10.00000000],
            frame='galactic',
            unit='deg')
    >>> print gaiadr2(c)


Documentation
-------------

Read the full documentation at http://dustmaps.readthedocs.io/en/latest/.


Citation
--------

If you make use of this software in a publication, please cite
[Green (2018) in The Journal of Open Source Software](https://doi.org/10.21105/joss.00695):

    @ARTICLE{2018JOSS....3..695M,
           author = {{Green}, {Gregory M.}},
            title = "{dustmaps: A Python interface for maps of interstellar dust}",
          journal = {The Journal of Open Source Software},
             year = "2018",
            month = "Jun",
           volume = {3},
           number = {26},
            pages = {695},
              doi = {10.21105/joss.00695},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2018JOSS....3..695M},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }


Development
-----------

Development of `selectionfunctions` takes place on GitHub, at
https://github.com/DouglasBoubert/selectionfunctions. Any bugs, feature requests, pull requests,
or other issues can be filed there. Contributions to the software are welcome.
