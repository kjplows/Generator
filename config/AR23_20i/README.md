# Liquid Argon Experiment tune

This configuration is the result of requests from SBN and DUNE collaborations.
It is designed to serve as baseline model for SBN analyses. 

The physics content is the following
 - Valencia model for 1p1h, using z-expansion
 - SuSAv2 model for 2p2h
 - Spectral function like approach for the Local Fermi Gas
    - Includes events with SRC-like missing momentum within the bulk of the $p_{miss}$ distribution, 
      but not extending pass the LFG cutoff in $p_{miss}$.
 - The parameters related to pion production are taken from the G18_10a_02_11b tune in order to ensure a better starting point. 
 - De-exitation photons are enabled for Argon

## Modifications

  Four new tunes have been added, taking different QEL z-expansion axial FF parameters following
  https://arxiv.org/abs/2512.14097

  AR23_20i_01_000: uses `ZExp_minerva` parameters, without the full sum rules, Eq. (31)
  AR23_20i_01_001: uses `ZExp_minerva` parameters, with    the full sum rules, Eq. (33)
  AR23_20i_02_000: uses `ZExp_lqcd`    parameters, without the full sum rules, Eq. (39)
  AR23_20i_02_001: uses `ZExp_lqcd`    parameters, with    the full sum rules, Eq. (41)

  A summary table of the parameters and tunes is below:

  |-----------------------|------|--------------|----------------|--------------------------|
  | Tune: AR23_20i_X0_000 | Kmax | T0 (GeV/c)^2 | Tcut (GeV/c)^2 |       A1,A2...           |
  |-----------------------|------|--------------|----------------|--------------------------|
  |  Def: AR23_20i_00_000 |    4 |        -0.28 |       0.1764   |      2.30,-0.6,-3.8,2.4  |
  |-----------------------|------|--------------|----------------|--------------------------| 
  |  Mnv: AR23_20i_01_000 |    2 |        -0.5  |       0.161604 |              1.65,-0.94  |
  |-----------------------|------|--------------|----------------|--------------------------|
  |  Mnv: AR23_20i_01_001 |    6 |        -0.5  |       0.161604 |  1.64778080,-0.94181417, |
  |                       |      |              |                | -0.41239729,-0.36611559, |
  |                       |      |              |                |  1.18722194,-0.49976799  |
  |-----------------------|------|--------------|----------------|--------------------------|
  | LQCD: AR23_20i_02_000 |    2 |        -0.5  |       0.161604 |             1.721,-0.31  |
  |-----------------------|------|--------------|----------------|--------------------------|
  | LQCD: AR23_20i_02_001 |    6 |        -0.5  |       0.161604 |  1.72089706,-0.30982708, |
  |                       |      |              |                | -1.62125837, 0.27506993, |
  |                       |      |              |                |  1.25297945,-0.60044079  |
  |-----------------------|------|--------------|----------------|--------------------------|