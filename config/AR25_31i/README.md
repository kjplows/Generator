# SBN tune with CRPA 1p1h

<<<<<<< HEAD:config/AR25_22i/README.md
This configuration is the same as SBN25_20i_00_000 except that the CRPA model is
used for charged-current 1p1h instead of the Valencia model
=======
This configuration is based upon AR23_20i_00_000 but has modifications
requested in summer 2025 by the Short Baseline Neutrino program experiments.

Notable components of the physics model are the following:
 - Valencia model for 1p1h, using z-expansion, with RPA turned off
 - SuSAv2 model for 2p2h
 - Spectral function for select nuclei, including one for 40Ar based on JLab
   measurements (see https://doi.org/10.1103/PhysRevD.105.112002
   and https://doi.org/10.1103/PhysRevD.107.012005)
 - Other nuclei use the AR23_20i_00_000 "spectral-function-like approach" for LFG
 - The parameters related to pion production are taken from the G18_10a_02_11b
   tune in order to ensure a better starting point.
 - De-exctitation photons are enabled for 40Ar
 - Minerva form factors from https://arxiv.org/pdf/2512.14097
>>>>>>> b142dcc... Jan 7 2025: rename minerva & lqcd tunes:config/AR25_31i/README.md
