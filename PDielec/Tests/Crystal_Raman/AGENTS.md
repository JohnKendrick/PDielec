# Crystal Raman Test Guidelines

- Layered Raman amplitudes follow `eq-raman_at_zs` and depth integration follows `eq-raman_at_zs_integrated`.
- NAC direction tests should trace `eq-layer-qph`, `eq-layer-qhat`, `eq-layer-dnac`, and `eq-layer-nac-eigenproblem`.
- EO correction tests should cover zero `chi2`, zero Born charge, linear scaling, and tensor convention consistency.
- Focused command: `pytest PDielec/Tests/Crystal_Raman -v`.
