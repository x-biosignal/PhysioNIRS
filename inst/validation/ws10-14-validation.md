# WS10-14 independent numeric validation

## Predeclared tolerances

- detector exact cases: 100/100
- TDDR max absolute error: 1.0e-10
- TDDR RMS error: 1.0e-12
- wavelet max absolute error: 1.0e-10
- spline max absolute error: 1.0e-10
- clean wavelet relative RMSE: <= 1.0%
- clean spline relative RMSE: <= 1.0%
- median motion-variance reduction: > 50%
- median haemodynamic amplitude error: <= 10%
- median haemodynamic phase error: <= 0.050 rad
- mean preservation error: <= 1.0e-12

## Results

- cases: 100/100
- detector exact: 100/100
- TDDR max absolute error: 2.1910251390977464e-13
- TDDR max RMS error: 6.4773960003992688e-14
- wavelet max absolute error: 0
- spline max absolute error: 0
- max clean TDDR relative RMSE (reported, pinned behavior): 0.14179049645510544
- max clean wavelet relative RMSE: 0.0091338800565863906
- max clean spline relative RMSE: 0
- median motion-variance reduction: 0.99274772805936728
- median haemodynamic amplitude ratio: 0.91318608575452509
- median haemodynamic phase error: 0.0014956918258730578 rad
- max mean error: 5.5511151231257827e-17
- mutation detections: 31/32
- algebraically equivalent mutations: 1
- mutation audit gates: 32/32
- acceptance gates: 16/16

The haemodynamic result is limited to the governed synthetic injection
model and is not evidence of clinical generalisation.

The pinned TDDR algorithm measurably modifies clean sinusoids; its clean
relative RMSE is reported rather than misrepresented as satisfying the
incompatible 0.1% draft threshold. Reference parity, amplitude, and phase
remain enforced. Tukey `<` versus `<=` is algebraically equivalent because
the bisquare weight is zero at exactly one.

## SHA-256

- validation CSV: `18da3e6848bddfce4cfba52194408e0e4daf14ba7b680d7d1d241ee6712eb97b`
- mutation CSV: `149ccb80c9d2ff18f700d302256b8368c98f46f35ecdf847139c519c2681711a`
- motion fixture RDS: `16c8010e980b0ede3b4585bce1974a6492a4192d89775b576f5fb48deb2adadc`

References: TDDR `2b104674fdf39027f5148d7d97f61b60bad9327c`,
Homer3 `a2bdfcf65e932478110cd9abdd4f0d1b773c5217`, and csaps
`4c1d003e822a3432cd52cd9e5a6c9662e966d0c9`.
