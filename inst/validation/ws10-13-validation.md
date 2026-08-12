# WS10-13 independent numeric validation

- Package: PhysioNIRS 0.2.0
- R: 4.3.3
- Oracle: independent base R equations; no PhysioNIRS helper is used
- Cases: 100/100
- Wavelength counts: 2, 3, 4
- Extinction table SHA-256: 261021720c0abece2b212cb35a8eb9cf4cc28d4fe0a00855143c1eb9d259d4d5
- Fixture SHA-256: 13842d1810bcdb669a33ce9daa4012ab47dd99a09484a584008543d4628e2ab6
- Fixture generator SHA-256: 49ea4b046c0360950d7eab6741e466313a5a9569ecd1f2163d836f44ce6706a8
- External comparator: unavailable (MNE/Homer3 not installed); pinned MNE data source is validated

## Tolerances

- OD max absolute error < 1e-12
- DPF max absolute error < 1e-12
- Extinction max relative error < 1e-12
- HbO/HbR max absolute error < 1e-8 uM
- HbT identity max absolute error < 1e-12 uM

## Results

- OD max absolute error: 9.5062846483529029e-16
- DPF max absolute error: 0e+00
- Extinction max relative error: 0e+00
- HbO max absolute error (uM): 1.1102230246251565e-15
- HbR max absolute error (uM): 1.3877787807814457e-15
- HbT identity max absolute error (uM): 1.1102230246251565e-16
- Sign cases: 100/100
- Mapping cases: 100/100
- Source immutability cases: 100/100
- Mutation gates: 30/30
- Gates: 10/10
- Final: PASS
