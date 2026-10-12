# Separate RadicalPy time-evolution simulation

The existing MARY files and `/simulate/mfe` route are retained unchanged.

## New files and changes

- `spin_models/time_evolution.py`: a separate Hilbert-space time-evolution calculation.
- `main.py`: adds `POST /simulate/time-evolution` while keeping the MARY endpoint.
- `web-demo/time_evolution.html`: new standalone Plotly page. The original `web-demo/index.html` is unchanged.

## Deployment

1. Copy **`spin_models/time_evolution.py`** into the `spin_models/` directory of your Render API repository.
2. Replace **`main.py`** with the provided version; redeploy Render.
3. Upload **`web-demo/time_evolution.html`** as e.g. `tutorials/time_evolution.html` to GitHub Pages.
4. Check that the page's API URL points to your Render service (`https://qubelab-spin-api.onrender.com`).
5. If desired, add a Simulations dropdown link to `tutorials/time_evolution.html`.

## Parameters

Retains all physically relevant controls from MARY: initial state (S, T, T+, T−, T0), observable, one-proton hyperfine coupling, fixed B field (rather than MARY maximum sweep field), J, D, exponential rate, duration and time step. MARY-only `field_points` is not applicable to a fixed-field time trajectory.

## Physical meaning

- Uses RadicalPy `HilbertSimulation.total_hamiltonian(B0, J, D)`, `time_evolution(initial_state, time, H)`, and `product_probability(state, density_matrices)`.
- The curves labelled *spin populations* show normalized, coherent spin-state populations, and **do not** depend on recombination rate.
- The curves labelled *surviving* are multiplied by `exp(-k*t)`, the probability a radical pair has not recombined, reproducing the **uniform exponential model** in the MARY example. This is **not spin-selective Haberkorn kinetics**.
- The total triplet population is the sum of T+, T0, T−. All four distinct populations sum to 1 in the unweighted calculation (subject to rounding).
- The selected rate of 10^8 s^-1 uses an adaptive sampling target; output is capped at 5001 samples. Long durations can therefore be undersampled, especially with strong magnetic fields; verify time resolution for scientific use.
- Model is educational and numerical/physical validation on your deployed RadicalPy environment is recommended.

## Test

`python -m compileall main.py spin_models/`

Run `uvicorn main:app --reload` and use `POST /simulate/time-evolution` or browse `/docs`.
