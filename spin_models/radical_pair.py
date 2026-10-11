"""Small RadicalPy MARY calculation for the QubeLab tutorial API.

Model:
- two electron spins (one per radical)
- one 1H nucleus coupled isotropically to the first radical
- singlet-born radical pair
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
    field_points: int = 61,
    duration_us: float = 2.0,
    time_step_ns: float = 10.0,
    observable: str = "singlet",
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
    n_time = max(2, int(round(duration_s / dt_s)) + 1)
    time_s = np.linspace(0.0, duration_s, n_time)
    fields_mT = np.linspace(0.0, float(field_max_mT), int(field_points))

    obs_state = State.SINGLET if observable == "singlet" else State.TRIPLET
    result = mary(
        sim,
        init_state=State.SINGLET,
        obs_state=obs_state,
        time=time_s,
        B=fields_mT,
        D=0.0,
        J=0.0,
        kinetics=[kinetics.Exponential(float(recombination_rate_s))],
    )

    # RadicalPy's MARY output is the normalised field-dependent response.
    # Convert arrays to native Python lists for JSON serialization.
    mary_values = np.asarray(result["MARY"], dtype=float).reshape(-1)
    if len(mary_values) != len(fields_mT):
        raise RuntimeError(
            f"Unexpected MARY output length {len(mary_values)} "
            f"for {len(fields_mT)} fields."
        )

    return {
        "field_mT": fields_mT.tolist(),
        "mfe_percent": mary_values.tolist(),
        "observable": observable,
        "model": "RadicalPy one-proton radical pair; singlet-born; exponential kinetics",
        "parameters": {
            "hyperfine_mT": float(hyperfine_mT),
            "recombination_rate_s": float(recombination_rate_s),
            "duration_us": float(duration_us),
            "time_step_ns": float(time_step_ns),
            "field_points": int(field_points),
        },
        "units": {
            "field_mT": "mT",
            "mfe_percent": "percent relative to the zero-field reference, as defined by RadicalPy MARY output",
        },
        "note": (
            "Illustrative model. The observable, kinetic scheme, spin Hamiltonian, "
            "and parameter values should be checked against the intended physical system."
        ),
    }
