### Fixed
- `_lc_lift!`: catastrophic cancellation in `r + x` when `x ≈ −r` (i.e. when the first LC pre-image coordinate is much smaller than the second). Direct subtraction `r + x` lost ~7 digits of precision; now uses `y²/(r − x)` (and symmetrically `y²/(r + x)` for the `r − x` branch), which is well-conditioned for all inputs. The bug surfaced as a flaky CI failure in the Levi-Civita round-trip identity test (`test_regularization.jl:215`, error ~3×10⁻⁷ vs tolerance 10⁻¹⁰).
- Levi-Civita round-trip test is now seeded (`Random.seed!(42)`) for deterministic reproduction.
