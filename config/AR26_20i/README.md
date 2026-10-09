# AR26_20i_00_000: AR25_20i + MARLEY

Copy of AR25_20i_00_000 with two additional, dedicated MARLEY channels
(`ProcessInfo`: `kScMarley` + `kIntWeakCC` / `kIntWeakNC`) simulated with the external MARLEY generator:
 - **MARLEY-CC**: low-energy nu_e + 40Ar charged-current scattering (HF-CRPA continuum +
   Bhattacharya 2009 discrete transitions; see `$GENIE/data/evgen/marley/ar26_cc.js`).
 - **MARLEY-NC**: coherent elastic (anti)neutrino-40Ar scattering, all four flavours
   (CEvNS; see `$GENIE/data/evgen/marley/ar26_nc.js`).

How it works:
 - The cross section of the MARLEY-CC channel is the *total* MARLEY cross section.
   MARLEY chooses the reaction and the kinematics and de-excites the nucleus when it is called.
 - MARLEY is blended with the native GENIE CC channels with an energy-dependent probability
   (`ChannelBlend-Emin`/`ChannelBlend-Emax` in ModelConfiguration.xml): MARLEY only below
   60 MeV, native GENIE above 100 MeV, linear in probability in between. The weights multiply
   the channel cross sections (see `src/Framework/EventGen/ChannelBlend.h`), so the neutrino
   acceptance in GMCJDriver sees the blended total cross section.
 - Native GENIE channels are only suppressed for initial states for which a MARLEY channel exists.
 - Splines: the native GENIE splines of AR25_20i are unchanged. Only the MARLEY-CC splines are new
   (`gmkspl --tune AR26_20i_00_000 --event-generator-list MARLEY-CC ...`, and likewise for MARLEY-NC).

---

# SBN tune

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
