.. role:: python(code)
   :language: python

Examples
========

Getting Started
---------------

Here, we'll look up a selection function at a number of different locations on the sky and a number of different magnitudes.
The principal object in `selectionfunctions` the `Source` object, which has both a `SkyCoord` attribute giving the position and a `Photometry` attribute giving the photometric measurements.
We specify coordinates on the sky using
`astropy.coordinates.SkyCoord <http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html>`_
objects. This allows us a great deal of flexibility in how we specify sky
coordinates. We can use different coordinate frames (e.g.,
`Galactic <https://en.wikipedia.org/wiki/Galactic_coordinate_system>`_,
`equatorial <https://en.wikipedia.org/wiki/Equatorial_coordinate_system>`_,
`ecliptic <https://en.wikipedia.org/wiki/Ecliptic_coordinate_system>`_),
different units (e.g., degrees, radians,
`hour angles <https://en.wikipedia.org/wiki/Hour_angle>`_), and either
scalar or vector input.

For our first example, let's load the
`Boubert & Everall (2019, in prep.)`_
-- or "boubert_everall_2019" -- selection function for Gaia DR2, and then query the selection function of a G=21 source at one location
on the sky:

.. code-block :: python
    
    from __future__ import print_function
    from astropy.coordinates import SkyCoord
    from dustmaps.sfd import SFDQuery
    
    coords = SkyCoord('12h30m25.3s', '15d15m58.1s', frame='icrs')
    sfd = SFDQuery()
    ebv = sfd(coords)
    
    coords = SkyCoord('12h30m25.3s', '15d15m58.1s', frame='icrs')
    print('E(B-V) = {:.3f} mag'.format(ebv))
    
    >>> E(B-V) = 0.030 mag

A couple of things to note here:

1. In this example, we used :python:`from __future__ import print_function` in
   order to ensure compatibility with both Python 2 and 3.
2. Above, we used the
   `ICRS coordinate system <https://en.wikipedia.org/wiki/International_Celestial_Reference_System>`_,
   by specifying :python:`frame='icrs'`.
3. :python:`SFDQuery` returns reddening in a unit that is similar to magnitudes
   of **E(B-V)**. However, care should be taken: a unit of SFD reddening is not
   quite equivalent to a magnitude of **E(B-V)**. The way to correctly convert
   SFD units to extinction in various broadband filters is to use the
   conversions in
   `Table 6 of Schlafly & Finkbeiner (2011) <http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6>`_.

We can query the other maps in the :code:`dustmaps` package with only minor
modification to the above code. For example, here's how we would query the
Planck Collaboration (2013) dust map:

.. code-block :: python
    
    from __future__ import print_function
    from astropy.coordinates import SkyCoord
    from dustmaps.planck import PlanckQuery
    
    coords = SkyCoord('12h30m25.3s', '15d15m58.1s', frame='icrs')
    planck = PlanckQuery()
    ebv = planck(coords)
    
    print('E(B-V) = {:.3f} mag'.format(ebv))
    
    >>> E(B-V) = 0.035 mag


Querying Reddening at an Array of Coordinates
---------------------------------------------

We can also query an array of coordinates, as follows:


.. code-block :: python
    
    from __future__ import print_function
    import numpy as np
    from astropy.coordinates import SkyCoord
    from dustmaps.planck import PlanckQuery
    from dustmaps.sfd import SFDQuery
    
    l = np.array([0., 90., 180.])
    b = np.array([15., 0., -15.])
    
    coords = SkyCoord(l, b, unit='deg', frame='galactic')
    
    planck = PlanckQuery()
    planck(coords)
    >>> array([ 0.50170666,  1.62469053,  0.29259142])
    
    sfd = SFDQuery()
    sfd(coords)
    >>> array([ 0.55669367,  2.60569382,  0.37351534], dtype=float32)

The input need not be a flat array. It can have any shape -- the shape of the
output will match the shape of the input:

.. code-block :: python
    
    from __future__ import print_function
    import numpy as np
    from astropy.coordinates import SkyCoord
    from dustmaps.planck import PlanckQuery
    
    l = np.linspace(0., 180., 12)
    b = np.zeros(12, dtype='f8')
    l.shape = (3, 4)
    b.shape = (3, 4)
    
    coords = SkyCoord(l, b, unit='deg', frame='galactic')
    
    planck = PlanckQuery()
    
    ebv = planck(coords)
    
    print(ebv)
    >>> [[ 315.52438354   28.11778831   23.53047562   20.72829247]
         [   2.20861101   15.68559361    1.46233201    1.70338535]
         [   0.94013882    1.11140835    0.38023439    0.81017196]]
    
    print(ebv.shape)
    >>> (3, 4)


Plotting the Dust Maps
----------------------

We'll finish by plotting a comparison of the SFD, Planck Collaboration and
Bayestar Dust maps. First, we'll import the necessary modules:

.. code-block :: python
    
    from __future__ import print_function
    
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    
    import astropy.units as units
    from astropy.coordinates import SkyCoord
    
    from dustmaps.sfd import SFDQuery
    from dustmaps.planck import PlanckQuery
    from dustmaps.bayestar import BayestarQuery

Next, we'll set up a grid of coordinates to plot, centered on the Aquila South
cloud:

.. code-block :: python
    
    l0, b0 = (37., -16.)
    l = np.arange(l0 - 5., l0 + 5., 0.05)
    b = np.arange(b0 - 5., b0 + 5., 0.05)
    l, b = np.meshgrid(l, b)
    coords = SkyCoord(l*units.deg, b*units.deg,
                      distance=1.*units.kpc, frame='galactic')

Then, we'll load up and query three different dust maps:

.. code-block :: python
    
    sfd = SFDQuery()
    Av_sfd = 2.742 * sfd(coords)
    
    planck = PlanckQuery()
    Av_planck = 3.1 * planck(coords)
    
    bayestar = BayestarQuery(max_samples=1)
    Av_bayestar = 2.742 * bayestar(coords)

We've assumed :math:`R_V = 3.1`, and used the coefficient from
`Table 6 of Schlafly & Finkbeiner (2011) <http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6>`_
to convert SFD and Bayestar reddenings to magnitudes of :math:`A_V`.

Finally, we create the figure using :code:`matplotlib`:

.. code-block :: python
    
    fig = plt.figure(figsize=(12,4), dpi=150)
    
    for k,(Av,title) in enumerate([(Av_sfd, 'SFD'),
                                   (Av_planck, 'Planck'),
                                   (Av_bayestar, 'Bayestar')]):
        ax = fig.add_subplot(1,3,k+1)
        ax.imshow(
            np.sqrt(Av)[::,::-1],
            vmin=0.,
            vmax=2.,
            origin='lower',
            interpolation='nearest',
            cmap='binary',
            aspect='equal'
        )
        ax.axis('off')
        ax.set_title(title)
    
    fig.subplots_adjust(wspace=0., hspace=0.)
    plt.savefig('comparison.png', dpi=150)

Here's the result:

.. image :: figs/boubert_everall_2019.png
