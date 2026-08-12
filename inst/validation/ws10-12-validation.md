# WS10-12 validation

Deterministic validation of governed SNIRF I/O. Reports contain no subject identifiers.

## Environment

- R: 4.3.3
- rhdf5: 2.46.1
- HDF5: 1.10.7
- SNIRF specification source: fNIRS/snirf tag v1.1
- Specification commit: 811c16ee730275d196b0ad910157ee632460bd57
- Official sample commit: e584d530a0903da250953df8a96affff547f039d
- Official sample SHA-256: 4673f295c2acba2f85beb80fcf9a1e5498a92c02b4c7d10b110c0331d30149db
- Official sample license: public domain
- pysnirf2: 0.7.3; NumPy 2 compatibility alias required
- MNE/MNE-NIRS: not installed; cross-check not run

## Result

- Generated cases: 100/100 pass
- Official sample: pass
- Maximum write/read data error: 0
- Maximum direct HDF5 data error: 0
- Maximum time error: 2.2204460492503131e-16
- Maximum distance error: 0
- Maximum event onset/duration error: 0 / 0
- Probe geometry and wavelength arrays: exact in every case
- Two-cycle semantic equality: pass in every case
- Package-owned fixture: two generations byte-identical
- pysnirf2 package-owned fixture: valid, zero errors, one coordinate-system warning
- pysnirf2 official sample: valid, zero errors, zero warnings

## Geometry classes

| Class | Cases | Pass |
|---|---:|---:|
| 2d | 46 | 46 |
| 3d | 27 | 27 |
| both | 27 | 27 |

## Time classes

| Class | Cases | Pass |
|---|---:|---:|
| compact | 35 | 35 |
| full | 35 | 35 |
| irregular | 30 | 30 |

## Measurement encodings

| Encoding | Cases | Pass |
|---|---:|---:|
| indexed | 50 | 50 |
| vectorized | 50 | 50 |

## Mutation gates

| Gate | Mutation | Detected |
|---:|---|:---:|
| 1 | data matrix transposed | yes |
| 2 | compact time interpreted as two samples | yes |
| 3 | milliseconds treated as seconds | yes |
| 4 | irregular time reported as uniform | yes |
| 5 | measurement list reordered | yes |
| 6 | source and detector indices swapped | yes |
| 7 | wavelength index treated as wavelength value | yes |
| 8 | indexed and vectorized lists combined | yes |
| 9 | missing measurement silently dropped | yes |
| 10 | probe coordinate axis dropped | yes |
| 11 | centimeters reported as meters | yes |
| 12 | 2-D geometry mislabeled as 3-D | yes |
| 13 | duplicate channel labels accepted | yes |
| 14 | stimulus duration/value columns swapped | yes |
| 15 | extra stimulus columns lost | yes |
| 16 | stale stored stimuli written after event edits | yes |
| 17 | unknown metadata tags lost | yes |
| 18 | subject identifier printed | yes |
| 19 | rank-1 singleton accepted as required scalar | yes |
| 20 | overwrite failure leaves a partial file | yes |
| 21 | official expectations inferred from R reader | yes |
| 22 | runtime Python/network dependency introduced | yes |

Detailed numeric results are in `ws10-12-validation.csv`; mutation results are in `ws10-12-mutation-gates.csv`.
