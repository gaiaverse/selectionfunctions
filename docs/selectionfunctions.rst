Available Selection Functions
=============================


Gaia selection functions
----------------------------


Gaia DR2 (boubert_everall_2019)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The selection function of the Gaia source catalogues is principally determined by the scanning law.
Sources must be detected at least five times to have made it into Gaia DR2, but sources are not detected every time Gaia observes them.
This selection function models the probability that a source is detected as a function of G magnitude and, optionally, the local source density.

There are four variants of this selection function implemented, which can be used by changing the `version` and `crowding` parameters.

* **Reference**: `Boubert & Everall (2019, in prep.)`_
