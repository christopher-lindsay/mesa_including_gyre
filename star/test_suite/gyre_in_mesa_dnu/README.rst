.. _gyre_in_mesa_dnu:

******************
gyre_in_mesa_dnu
******************

This test suite checks whether the delta_nu computed by MESA is within 1 percent of the one computed by GYRE.

The first inlist evolves a star of 1Msun from the pre-MS to near-TAMS.


The second inlist loads in this near-TAMS model and continues to evolve the star to log(Teff)=3.7 and uses GYRE to compute delta_nu and the fundamental frequency.

Last-Updated: 12Aug2022 (MESA 22.05.1).
