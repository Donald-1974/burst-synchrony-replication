# burst-synchrony-replication

Independent re-runs of simulation specifications from my research notes on
**burst-driven synchrony near the Kuramoto synchronization transition** — the idea that
bursts of partial synchrony can hold an oscillator network near its tipping point
("edge of synchronization"), producing neuronal-avalanche-like statistics.

The specifications were drafted in AI-assisted working sessions (Dec 2025). Several
of those sessions printed "expected outputs" that were never actually computed. This
repository runs the specs as written, with fixed seeds, so every number below can be
reproduced by anyone. **Negative results are reported as-is.**

## Contents

| File | What it tests |
|---|---|
| `rgg.py` | Shared helpers: 3D periodic random geometric graph (RGG), row-normalised Kuramoto velocity, order parameter. |
| `rgg_check.py` | Connectivity of 3D periodic RGGs with N = 10,000 at r0 = 0.03, 0.05, 0.07: mean degree, components, giant component, isolated nodes. |
| `exact_spec_rK.py` | The Dec 20, 2025 "Hyperplastic Criticality" r(K) spec exactly as written (N = 10,000, r0 = 0.03). The original session claimed r(K = 0.85) ≈ 0.378. |
| `kuramoto_avalanche_check.py` | The "Kuramoto avalanche layer" spec (N = 2,000, r0 = 0.1): order parameter, spike counts (forward phase-0 crossings) and avalanche counts at two bin widths. |
| `connected_sweep.py` | Corrected follow-up: fine K sweep (41 values, 0 to 8) on a connected RGG (N = 2,000, r0 = 0.1, giant component only), 5 seeds, Heun integration, order parameter, susceptibility, and avalanche statistics with Beggs-Plenz binning and maximum-likelihood tail fits. Writes `results/` (SVG/PNG plots, `sweep.csv`, and raw `.npz`). The committed figures are the SVG and CSV; running the script also writes PNG and `.npz`. |
| `dt_check_seed42_dt0.01.out` | Step-size check: `connected_sweep.py` for seed 42 rerun at dt = 0.01. |
| `*.out` | Exact stdout of each script with the default arguments (seed 42). |

## Model

dθᵢ/dt = ωᵢ + (K / kᵢ) Σⱼ Aᵢⱼ sin(θⱼ − θᵢ), ωᵢ ~ N(0, 1),

where A is the adjacency matrix of a random geometric graph in the unit cube with periodic
boundaries (nodes linked when closer than r0) and kᵢ is the degree of node i. Integration is
forward Euler (dt = 0.02 or 0.01) in the original checks and Heun's method (dt = 0.02) in
`connected_sweep.py`; the global order parameter is r = |⟨e^{iθ}⟩|.

## How to run

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python rgg_check.py                  # ~1 s
python exact_spec_rK.py              # ~30 s
python kuramoto_avalanche_check.py   # ~20 s
python connected_sweep.py            # ~3 min on 5+ cores (seeds run in parallel)
```

Each script accepts `--seed`, `--n` and the radius (`--r0`, or `--radii` for
`rgg_check.py`); run with `-h` for details. With the default arguments and the pinned
versions in `requirements.txt`, output matches the committed `.out` files exactly
(tested on Python 3.13.5, Debian Linux).

## Results

### 1. The Dec 2025 network is fragmented (`rgg_check.out`)

| r0 | mean degree (measured / expected N·4πr0³/3) | components | giant component | isolated nodes |
|---|---|---|---|---|
| 0.03 | 1.12 / 1.13 | 5,269 | 0.2% | 32.3% |
| 0.05 | 5.26 / 5.24 | 57 | 99.3% | 0.5% |
| 0.07 | 14.36 / 14.37 | 1 | 100.0% | 0.0% |

Continuum percolation in 3D requires a mean degree of roughly 2.7 (Lorenz & Ziff, 2001).
At r0 = 0.03 the mean degree is 1.12, so the "network" is thousands of tiny clusters and a
third of the nodes have no neighbours at all. Global synchronization is impossible.

### 2. The claimed r(K) curve does not reproduce (`exact_spec_rK.out`)

| K | 0.50 | 0.60 | 0.68 | 0.75 | 0.85 | 2.00 | 5.00 |
|---|---|---|---|---|---|---|---|
| ⟨r⟩ (t = 100–150) | 0.009 | 0.010 | 0.011 | 0.012 | 0.010 | 0.011 | 0.016 |

The order parameter stays at the finite-size noise floor (≈ 1/√N = 0.01) at every coupling,
including K = 5. The previously stated r(0.85) ≈ 0.378 is not supported; it follows from
result 1 that no synchronization transition exists in this spec.

### 3. Avalanche statistics are dominated by bin width (`kuramoto_avalanche_check.out`)

On a connected graph (r0 = 0.1, mean degree 8.5) the order parameter does rise with
coupling (r = 0.023 at K = 0.5 → 0.149 at K = 3 → 0.469 at K = 5), so this network does show
a synchronization transition. But the avalanche measure is fragile:

* With 0.10-time-unit bins, no bin is empty for K ≤ 2, so there are **zero** avalanches —
  the whole record is one unbroken event.
* With 0.01-time-unit bins there are ~2,100–2,500 avalanches for K ≤ 3. At K = 3 the two
  bin widths give 2,105 avalanches (max size 11) versus 35 (max size 555) from the same
  spike train, so the statistics depend strongly on the arbitrary bin width.

Avalanche analysis therefore needs a principled bin width (e.g. the mean inter-event
interval, as in Beggs & Plenz 2003), maximum-likelihood power-law fits tested against
alternatives, and a spike definition that does not ignore oscillators with negative ω.

### 4. Corrected, connected sweep (`connected_sweep.out`, `results/r_chi_vs_K.svg`, `results/sweep.csv`)

Setup: N = 2,000, r0 = 0.1 (mean degree 8.2 to 8.5), dynamics on the giant component only
(seed 43 drops 2 isolated nodes, the other four seeds are fully connected), ω ~ N(0, 1),
Heun integration with dt = 0.02, 100 time units of transient then 200 measured, 41 couplings
K = 0, 0.2, ..., 8, seeds 42 to 46.

* **There is a clear, continuous synchronization transition.** Seed-averaged ⟨r⟩ rises from
  0.020 (K = 0, the 1/√N noise floor) through 0.075 (K = 2), 0.164 (K = 3) and 0.317 (K = 4)
  to 0.584 (K = 5) and 0.833 (K = 8).
* **Critical coupling estimate: K_c ≈ 4.0.** The seed-averaged susceptibility χ = N·var(r)
  peaks at K = 4.0 (χ = 24 ± 17). Per-seed peaks are 4.2, 3.6, 4.6, 4.0 and 6.8; the 6.8 value
  is an isolated spike from one seed (likely a transient slip in a nearly locked state) and
  isolated spikes at K = 6.8 and 7.4 in the averaged curve come from single seeds. This is
  about 2.5 times the all-to-all mean-field value 2/(π g(0)) = 1.596, as expected for a sparse
  (mean degree ~8) spatially embedded network.
* **Step-size check.** Rerunning seed 42 with dt = 0.01 changes ⟨r⟩ by a median of 0.0004
  across all K, with a maximum difference of 0.048 at K = 3.8, right at the transition where
  fluctuations are largest. The χ peak for that seed moves from 4.2 to 3.8, so K_c should be
  read as 4.0 ± 0.4.
* **No avalanche criticality signature at K_c with this spike definition.** Using seed 42,
  spikes are phase-0 crossings in either direction (so negative-ω oscillators count), binned
  at the population mean inter-event interval. Away from the transition (K = 0, 2) the size
  tails are steeper than any power law (the fitted exponent hits the fit limit of 5) and are
  as well or better described by an exponential. At K = 4 the best power-law fit is
  α = 3.14 (xmin = 13, 122 tail events, KS = 0.054), but the likelihood-ratio test against an
  exponential tail is inconclusive (R = 0.5, Vuong p = 0.88), and α ≈ 3 is far from the
  ≈ 1.5 reported for neuronal avalanches (Beggs & Plenz 2003). Above the transition
  (K = 6, 8) the synchronized cluster rotates at the mean frequency (≈ 0), so it almost
  never crosses phase 0 (228 and 12 spikes in 200 time units) and avalanches cannot be
  measured with this definition.

Caveats: one network size only (no finite-size scaling), one seed for the avalanche fits,
avalanche maximum sizes are small (≤ 74), and the phase-crossing spike definition is itself
a modelling choice. The burst-driven variant (bursts that transiently raise K) was not run.
So this result says the plain connected Kuramoto model, measured this way, does **not**
show scale-free avalanches at its transition; it does not rule out the burst-driven
mechanism, which would need its own explicit model and a pre-stated falsification test.

## What this means for the hypothesis

These checks do **not** refute the burst-synchrony idea; they show that the specific
Dec 2025 specification could not test it. The idea itself sits close to published work in
which scale-free avalanches appear at the edge of a synchronization transition (di Santo
et al., 2018; Buendía et al., 2021). Next steps:

1. Done in result 4 (connected network, fine K sweep, five seeds, step-size check). Still
   open: finite-size scaling across several N.
2. Define bursts/avalanches following the methods in the references below and state, in
   advance, the result that would falsify the hypothesis.
3. Compare against the established Landau–Ginzburg / hybrid-transition results rather than
   claiming novelty for avalanches at the edge of synchronization itself.

## References

* di Santo, S., Villegas, P., Burioni, R., & Muñoz, M. A. (2018). Landau–Ginzburg theory of
  cortex dynamics: Scale-free avalanches emerge at the edge of synchronization. *PNAS*,
  115(7). https://doi.org/10.1073/pnas.1712989115
* Buendía, V., Villegas, P., Burioni, R., & Muñoz, M. A. (2021). Hybrid-type
  synchronization transitions: Where incipient oscillations, scale-free avalanches, and
  bistability live together. *Physical Review Research*, 3, 023224.
  https://doi.org/10.1103/PhysRevResearch.3.023224
* Peron, T. K. D., & Rodrigues, F. A. (2012). Explosive synchronization enhanced by
  time-delayed coupling. *Physical Review E*, 86, 016102.
  https://doi.org/10.1103/PhysRevE.86.016102
* Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-law distributions in
  empirical data. *SIAM Review*, 51(4), 661-703. https://doi.org/10.1137/070710111
* Lorenz, C. D., & Ziff, R. M. (2001). Precise determination of the critical percolation
  threshold for the three-dimensional "Swiss cheese" model using a growth algorithm.
  *Journal of Chemical Physics*, 114, 3659. https://doi.org/10.1063/1.1338506
* Beggs, J. M., & Plenz, D. (2003). Neuronal avalanches in neocortical circuits.
  *Journal of Neuroscience*, 23(35), 11167–11177.
  https://doi.org/10.1523/JNEUROSCI.23-35-11167.2003

## License

MIT — see [LICENSE](LICENSE).
