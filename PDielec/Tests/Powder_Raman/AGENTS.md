# Powder Raman Test Guidelines

- Test powder tensor algebra against `eq-invariants1` and `eq-Intensities`.
- Local-field tests should distinguish the stored bulk tensor, the physical `R_epsilon`, and the effective bulk tensor `R_eff`.
- Use deterministic orientation sampling or analytic invariants where possible.
- Focused command: `pytest PDielec/Tests/Powder_Raman -v`.
