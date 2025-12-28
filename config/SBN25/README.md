# Testing tunes for SBN

The `SBN25` line of tunes is derived from AR23_20i, with certain changes.
The "baseline tune" `SBN25_00a` is another name for `AR23_20i` and forms the baseline.

## QEL-CC model changes
These live in the subdirectory `QELModelVariations`.
They are of the form `SBN25_xyz`.

| xyz | QEL model  | Form factor | RPA? |
|-----|------------|-------------|------|
| 00a | Nieves     | ZExp        | yes  |
| 00b | Nieves     | ZExp        | no   |
| 01a | Nieves     | ZExp        | yes  |
| 01b | Nieves     | Dipole      | no   |
| 02a | Nieves     | RunningMA   | yes  |
| 02b | Nieves     | RunningMA   | no   |
| 10a | LwlynSmith | ZExp        | n/a  |
| 11a | LwlynSmith | Dipole      | n/a  |
| 12a | LwlynSmith | RunningMA   | n/a  |
| 20a | Hybrid     | SuSAv2      | n/a  |
| 30a | Hybrid     | SuSAv2Blend | n/a  |
| 40a | Hybrid     | CRPA+SuSA   | n/a  |
| 50a | Hybrid     | CRPA        | n/a  |

## MEC model changes
These live in the subdirectory `MECModelVariations`.
They are of the form `SBN25_9ya`.

| 9ya | MEC model  | Fermi momentum table |
|-----|------------|----------------------|
| 00a | SuSAv2     | SuSA                 |
| 91a | SuSAv2     | Default              |
| 92a | NSValencia | SuSA                 |
| 93a | Empirical  | SuSA                 |

## Nuclear model changes
These live in the subdirectory `NuclearModelVariations`.
They are of the form `SBN25_99z`.

| 99z | Nuclear model   | Configuration |
|-----|-----------------|---------------|
| 00a | LocalFGM        | SFLike+Corr   |
| 99b | LocalFGM        | SFLike        |
| 99b | LocalFGM        | Correlated    |
| 99d | LocalFGM        | Default       |
| 99e | SpectralFunc    | Default       |
| 99f | SpectralFunc    | Coerced       |
| 99g | SpectralFunc1d  | Default       |
| 99h | EffectiveSF     | Default       |
| 99i | FGMBodekRitchie | Default       |