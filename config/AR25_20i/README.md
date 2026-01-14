# SBN tune

This configuration is based upon AR23_20i_00_000 but has modifications
requested in summer 2025 by the Short Baseline Neutrino program experiments.

Notable components of the physics model are the following:
 - Valencia model for 1p1h, using z-expansion, with RPA turned off
 - SuSAv2 model for 2p2h
 - Spectral function for select elements (C, O, Ar), including one for 40Ar based on JLab
   measurements (see https://doi.org/10.1103/PhysRevD.105.112002
   and https://doi.org/10.1103/PhysRevD.107.012005)
 - Isotopes of these elements use the SF of 12C, 16O, and 40 Ar, respectively.
 - Other nuclei use the AR23_20i_00_000 "spectral-function-like approach" for LFG
 - The parameters related to pion production are taken from the G18_10a_02_11b
   tune in order to ensure a better starting point.
 - De-exctitation photons are enabled for 40Ar

## Modifications

  Four new tunes have been added, taking different QEL z-expansion axial FF parameters following
  https://arxiv.org/abs/2512.14097 and T. Cai et al (MINERvA), Nature 614 (2023) 48-53.
  Note the arxiv paper uses the opposite sign convention to GENIE.

  AR25_20i_01_000: uses `ZExp_minervaFull`   parameters, with    the full sum rules, Eq. (33)
  AR25_20i_01_001: uses `ZExp_minerva`       parameters, without the full sum rules, Eq. (31)
  AR25_20i_01_002: uses `ZExp_minervaNature` parameters, from Nature paper Supp. Table 4
  AR25_20i_02_000: uses `ZExp_lqcdFull`      parameters, with    the full sum rules, Eq. (41)

  Note the `Q4limit` parameter, which applies the sum rules: 2512.14097 paper's kmax=6 is
  equivalent to GENIE kmax=2 with Q4limit set to true, which adds four extra constraints.
  This also means that Nature paper kmax=8 maps to kmax=4 in GENIE with Q4limit true

  Note also the Tcut values: The arxiv paper uses Tcut = 9 (m_pi0)^2 not 9 (m_pi+)^2

  Finally, note that AR25_20i_01_002 uses a different value for QEL-FA0 = -1.2723

  A summary table of the parameters and tunes is below:

  |-----------------------|------|---------|--------------|----------------|--------------------------|
  | Tune: AR25_20i_X0_000 | Kmax | Q4limit | T0 (GeV/c)^2 | Tcut (GeV/c)^2 |       A1,A2...           |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  |  Def: AR25_20i_00_000 |    4 | true    |        -0.28 |       0.1764   |      2.30,-0.6,-3.8,2.4  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  |  Mnv: AR25_20i_01_000 |    2 | true    |        -0.5  |       0.161604 |  1.64778080,-0.94181417  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  |  Mnv: AR25_20i_01_001 |    2 | false   |        -0.5  |       0.161604 |              1.65,-0.94  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  |  Mnv: AR25_20i_01_002 |    4 | true    |        -0.75 |       0.1764   |      1.50,-1.2,-0.1,0.2  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|
  | LQCD: AR25_20i_02_000 |    2 | true    |        -0.5  |       0.161604 |             1.721,-0.31  |
  |-----------------------|------|---------|--------------|----------------|--------------------------|