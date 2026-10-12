"""Small RadicalPy MARY calculation for the QubeLab tutorial API.

Model:
- two electron spins (one per radical)
- one 1H nucleus coupled isotropically to the first radical
- selectable initial singlet/triplet radical-pair state
- optional exchange J and dipolar D couplings (mT)
- exponential kinetic model
- magnetic-field sweep along the default z axis

This is a teaching model, not a fit to a particular molecule or experiment.
Check the installed RadicalPy version and validate results before research use.
"""

from __future__ import annotations

import numpy as np


def simulate_one_proton_mary(
    *,
    field_max_mT: float = 15.0,
    hyperfine_mT: float = 1.5,
    recombination_rate_s: float = 1.0e6,
    field_points: int = 31,
    duration_us: float = 2.0,
    time_step_ns: float = 10.0,
    observable: str = "singlet",
    initial_state: str = "singlet",
    exchange_mT: float = 0.0,
    dipolar_mT: float = 0.0,
) -> dict:
    try:
        import radicalpy as rp
        from radicalpy import kinetics
        from radicalpy.experiments import mary
        from radicalpy.simulation import HilbertSimulation, State
    except ImportError as exc:
        raise ImportError(
            "Install RadicalPy and its dependencies from requirements.txt."
        ) from exc

    # One isotropically coupled proton on radical A; radical B has no nuclei.
    radical_a = rp.simulation.Molecule.fromisotopes(
        isotopes=["1H"],
        hfcs=[float(hyperfine_mT)],
        name="Radical A",
    )
    radical_b = rp.simulation.Molecule("Radical B", [])
    sim = HilbertSimulation([radical_a, radical_b])

    # Time is in seconds; fields and hyperfine couplings are in mT.
    dt_s = float(time_step_ns) * 1e-9
    duration_s = float(duration_us) * 1e-6
    # Resolve the fastest decay even when k = 1e8 s^-1 (lifetime = 10 ns).
    # Avoid wasting thousands of samples beyond ~12 lifetimes at high rates.
    rate = float(recombination_rate_s)
    effective_dt_s = min(dt_s, 0.05 / rate)
    effective_duration_s = min(duration_s, 12.0 / rate)
    n_time = max(2, int(np.ceil(effective_duration_s / effective_dt_s)) + 1)
    time_s = np.linspace(0.0, effective_duration_s, n_time)
    fields_mT = np.linspace(0.0, float(field_max_mT), int(field_points))

    states = {
        "singlet": State.SINGLET,
        "triplet": State.TRIPLET,
        "triplet_plus": State.TRIPLET_PLUS,
        "triplet_minus": State.TRIPLET_MINUS,
        "triplet_zero": State.TRIPLET_ZERO,
    }
    if observable not in states:
        raise ValueError(f"Unknown observable: {observable}")
    if initial_state not in states:
        raise ValueError(f"Unknown initial state: {initial_state}")
    obs_state = states[observable]
    result = mary(
        sim,
        init_state=states[initial_state],
        obs_state=obs_state,
        time=time_s,
        B=fields_mT,
        D=float(dipolar_mT),
        J=float(exchange_mT),
        kinetics=[kinetics.Exponential(float(recombination_rate_s))],
    )

    # Recompute normalization from the kinetic-weighted populations rather than
    # relying on RadicalPy's (yield(B)-yield(0))/yield(0), which divides by zero
    # for uncoupled systems / unpopulated triplet substates.
    populations = np.asarray(result["time_evolutions"], dtype=float)
    if populations.shape != (len(fields_mT), len(time_s)):
        raise RuntimeError(f"Unexpected population array shape: {populations.shape}")
    # RadicalPy's Exponential kinetics weights the population by exp(-k*t).
    yields = rate * np.trapezoid(populations, x=time_s, axis=1) if hasattr(np, "trapezoid") else rate * np.trapz(populations, x=time_s, axis=1)
    if not np.all(np.isfinite(yields)):
        raise RuntimeError("Non-finite product yields")
    reference_yield = float(yields[0])
    if abs(reference_yield) > 1e-14:
        mary_values = 100.0 * (yields - reference_yield) / reference_yield
        response_unit = "percent relative to the zero-field yield"
    else:
        # A relative percentage is undefined if the B=0 yield is zero.
        # Show absolute changes (percentage points of yield) instead.
        mary_values = 100.0 * (yields - reference_yield)
        response_unit = "percentage points of product yield (zero-field yield is zero)"
    if len(mary_values) != len(fields_mT):
        raise RuntimeError(
            f"Unexpected MARY output length {len(mary_values)} "
            f"for {len(fields_mT)} fields."
        )

    return {
        "field_mT": fields_mT.tolist(),
        "mfe_percent": mary_values.tolist(),
        "observable": observable,
        "initial_state": initial_state,
        "model": f"RadicalPy one-proton radical pair; {initial_state.replace('_', ' ')}-born; exponential kinetics; J={exchange_mT:g} mT, D={dipolar_mT:g} mT",
        "parameters": {
            "hyperfine_mT": float(hyperfine_mT),
            "initial_state": initial_state,
            "exchange_mT": float(exchange_mT),
            "dipolar_mT": float(dipolar_mT),
            "recombination_rate_s": float(recombination_rate_s),
            "duration_us": float(duration_us),
            "time_step_ns": float(time_step_ns),
            "effective_time_step_ns": float((time_s[1] - time_s[0]) * 1e9),
            "effective_duration_us": float(time_s[-1] * 1e6),
            "field_points": int(field_points),
        },
        "units": {
            "field_mT": "mT",
            "mfe_percent": response_unit,
        },
        "note": (
            "Illustrative model. The observable, kinetic scheme, spin Hamiltonian, "
            "and parameter values should be checked against the intended physical system."
        ),
    }
