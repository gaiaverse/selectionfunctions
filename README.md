[![DOI](http://joss.theoj.org/papers/10.21105/joss.00695/status.svg)](https://doi.org/10.21105/joss.00695)

selectionfunctions
==================

The ``selectionfunctions`` package provides a uniform interface to the selection functions of the major astronomical surveys. This package is entirely derivative of the truly excellent ``dustmaps`` package created by Gregory M. Green. Please cite the ``dustmaps`` JOSS paper when using this package.

Supported Selection Functions
-----------------------------

The currently supported selection functions are:

1. Boubert & Everall (2019; boubert_everall_2019, in prep.)

To request addition of another selection function in this package, [file an issue on
GitHub](https://github.com/DouglasBoubert/selectionfunctions/issues), or submit a pull request.


Installation
------------

Download the repository from [GitHub](https://github.com/DouglasBoubert/selectionfunctions) and
then run:

    python setup.py install --large-data-dir=/path/where/you/want/large/data/files/stored

Alternatively, you can use the Python package manager `pip`:

    pip install selectionfunctions


Getting the Data
----------------

To fetch the data for the GaiaDR2 selectionfunction, run:

    python setup.py fetch --map-name=boubert_everall_2019

You can download the other selection functions by changing "boubert_everall_2019" to (other selection functions to be added in future).

Alternatively, if you have used `pip` to install `selectionfunctions`, then you can
configure the data directory and download the data by opening up a python
interpreter and running:

    >>> from selectionfunctions.config import config
    >>> config['data_dir'] = '/path/where/you/want/large/data/files/stored'
    >>>
    >>> import selectionfunctions.boubert_everall_2019
    >>> selectionfunctions.boubert_everall_2019.fetch()


Querying the Maps
-----------------

Maps are queried using SourceCoord objects, which are a variant on the 
[`astropy.coordinates.SkyCoord`](http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord)
object. This means that any coordinate system supported by `astropy` can be
used as input. For example, we can query BoubertEverall2019 as follows:

    >>> from selectionfunctions.boubert_everall_2019 import BoubertEverall2019Query
    >>> from selectionfunctions.source_base import SourceCoord
    >>>
    >>> boubert_everall_2019 = BoubertEverall2019Query()
    >>>
    >>> c = SourceCoord(
            '05h00m00.00000s',
            '+30d00m00.0000s',
            photometry={'gaia_g':21.2},
            frame='icrs')
    >>> print(boubert_everall_2019(c))


Above, we have used the ICRS coordinate system (the inputs are RA and Dec). We
can use other coordinate systems, such as Galactic coordinates, and we can
provide coordinate arrays. The following example uses both features:

    >>> c = SourceCoord(
            [75.00000000, 130.00000000],
            [-89.00000000, 10.00000000],
            photometry={'gaia_g':[2.3,17.8]},
            frame='galactic',
            unit='deg')
    >>> print(boubert_everall_2019(c))



Documentation
-------------

Read the full documentation at http://selectionfunctions.readthedocs.io/en/latest/.


Citation
--------

If you make use of this software in a publication, please cite
[Green (2018) in The Journal of Open Source Software](https://doi.org/10.21105/joss.00695).

Development
-----------

Development of `selectionfunctions` takes place on GitHub, at
https://github.com/DouglasBoubert/selectionfunctions. Any bugs, feature requests, pull requests,
or other issues can be filed there. Contributions to the software are welcome.
