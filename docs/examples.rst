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
`Boubert & Everall (2020, submitted)`
-- or "cog_ii" -- selection function for Gaia DR2, and then query the selection function of a G=21 source at one location
on the sky:

.. code-block :: python
    
    from selectionfunctions.source import Source
    from selectionfunctions import cog_ii
    
    coords = Source('12h30m25.3s', '15d15m58.1s', frame='icrs', photometry={'gaia_g':21.0})
    dr2_sf = cog_ii.dr2_sf(version='modelAB',crowding=True)
    prob_selection = dr2_sf(coords)
    
    print('Probability of selection = {:.3f}%'.format(prob_selection*100.0))
    
    >>> Probability of selection = 69.877%

A couple of things to note here:

1. Above, we used the
   `ICRS coordinate system <https://en.wikipedia.org/wiki/International_Celestial_Reference_System>`_,
   by specifying :python:`frame='icrs'`.
2. We specified the apparent Gaia G magnitude of the star through a Python dictionary. Photometric transformations have not yet been implemented.
3. We used the keywords :python:`version='modelAB'` and :python:`crowding=True` when constructing the selection function. By default, :python:`crowding=False`.

In future, you will be able to query other selection funtions from the :code:`selectionfunctions` package with only minor
modification to the above code.


Querying Selection Function at an Array of Coordinates
---------------------------------------------

We can also query an array of coordinates, as follows:


.. code-block :: python
    
    import numpy as np
    from selectionfunctions.source import Source
    from selectionfunctions import cog_ii
    
    l = np.array([0., 90., 180.])
    b = np.array([15., 0., -15.])
    g = np.array([20.8,21.0,21.2])
    
    coords = Source(l, b, unit='deg', frame='galactic', photometry={'gaia_g':g})
    
    dr2_sf = cog_ii.dr2_sf()
    dr2_sf(coords)
    >>> array([0.99997069, 0.96233884, 0.58957493])

The input need not be a flat array. It can have any shape -- the shape of the
output will match the shape of the input:

.. code-block :: python
    
    import numpy as np
    from selectionfunctions.source import Source
    from selectionfunctions import cog_ii
    
    l = np.linspace(0., 180., 12)
    b = np.zeros(12)
    g = 21.0*np.ones(12)
    l.shape = (3, 4)
    b.shape = (3, 4)
    g.shape = (3, 4)
    
    coords = Source(l, b, unit='deg', frame='galactic', photometry={'gaia_g':g})
    
    dr2_sf = cog_ii.dr2_sf()
    
    prob_selection = dr2_sf(coords)
    
    print(prob_selection)
    >>> [[0.74045863 0.69877491 0.74045863 0.94768624]
         [0.98794938 0.93834743 0.95561436 0.96803869]
         [0.99962099 0.97286789 0.91445208 0.59940653]]
    
    print(prob_selection.shape)
    >>> (3, 4)


Plotting a Selection Function
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
