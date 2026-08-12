# WS10-16 independent validation

- Quality cases: 120
- SCI classification: 720 / 720 agreement = 1
- Sensitivity/specificity: 1 / 1
- SCI maximum absolute error: 0
- Peak-power classification: 122 / 122
- Analytic SNR error dB: 1.0293515209980342e-05
- Permutation maximum error: 0
- Offset/scale maximum error: 6.5176336550010205e-13 / 7.2291478359076677e-13
- Live updates: 10 ; target maximum error: 0
- Mutations: 32 / 32
- Gates: 19 / 19
- Fixture SHA-256: cee8c3f42000f31fdda62213b350b500a1cdf53513e8c7dbea0d66085d05d0bc
- Validator SHA-256: b3c3e41ee7c8f79203cf4d5cbb6179c46a9a4c886b008a8ec77b76fbbef2c3fc

Runtime equations are independent base-R/signal implementations; no
PhysioNIRS private helper is called.

## Gate table

| gate | pass |
|---|---|
| quality_cases_120 | TRUE |
| classification_agreement_90pct | TRUE |
| classification_sensitivity_90pct | TRUE |
| classification_specificity_90pct | TRUE |
| sci_error_1e_10 | TRUE |
| peak_classification_95pct | TRUE |
| pair_copy_exact | TRUE |
| snr_error_0_1_db | TRUE |
| permutation_equivariance | TRUE |
| offset_invariance | TRUE |
| scale_invariance | TRUE |
| input_immutable | TRUE |
| mark_exact | TRUE |
| drop_exact | TRUE |
| live_target_error_1e_10 | TRUE |
| live_delivery_exact | TRUE |
| all_finite | TRUE |
| mutations_25 | TRUE |
| fixture_manifest_valid | TRUE |
