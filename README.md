Quasi-1D Euler
==============

Quasi 1D Euler solver for steady and unsteady problems using 2nd order MUSCL scheme with various limiters and two Riemann fluxes - local Lax-Freidrichs and Van Leer flux vector splitting. Explicit time stepping.

Two variants of MUSCL have been implemented - MUSCLReconstructionG based on least-squares slopes, and MUSCLReconstruction based on Blazek's method. The latter works better at the moment. With the former, note that only Van Albada limiter works.
