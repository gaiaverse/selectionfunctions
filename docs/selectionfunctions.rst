Available Selection Functions
=============================


Gaia selection functions
----------------------------


Gaia DR2 (cog_ii.dr2_sf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The selection function of the Gaia source catalogue is principally determined by the scanning law, which gives the number of times that Gaia observes every location on the sky.
Sources must be detected at least five times to have made it into Gaia DR2, but sources are not detected every time Gaia observes them.
This selection function models the probability that a source is detected as a function of G magnitude and, optionally, the local source density.

There are four variants of this selection function implemented, which can be used by changing the `version` and `crowding` parameters. We recommend setting :python:`version='modelAB'` and :python:`crowding=True` for most applications, but these are not set by default.

* **Reference**: `Boubert & Everall (2020, submitted)`
