"""Time-dependent radical-pair spin populations using RadicalPy Hilbert dynamics.

The exponential recombination factor is a uniform survival weight, not a
spin-selective Haberkorn reaction or an independently modelled product channel.
"""
from __future__ import annotations
import numpy as np

STATES = ("singlet", "triplet", "triplet_plus", "triplet_minus", "triplet_zero")


def simulate_time_evolution(*, field_mT: float = 5.0,
                            hyperfine_mT: float = 1.5,
                            recombination_rate_s: float = 1e6,
                            duration_us: float = 2.0,
                            time_step_ns: float = 10.0,
                            initial_state: str = "singlet",
                            observable: str = "singlet",
                            exchange_mT: float = 0.0,
                            dipolar_mT: float = 0.0) -> dict:
    import radicalpy as rp
    from radicalpy.simulation import HilbertSimulation, State

    mapping = {
        "singlet": State.SINGLET, "triplet": State.TRIPLET,
        "triplet_plus": State.TRIPLET_PLUS,
        "triplet_minus": State.TRIPLET_MINUS,
        "triplet_zero": State.TRIPLET_ZERO,
    }
    if initial_state not in mapping or observable not in mapping:
        raise ValueError("Unknown spin state")

    radical_a = rp.simulation.Molecule.fromisotopes(
        isotopes=["1H"], hfcs=[float(hyperfine_mT)], name="Radical A")
    radical_b = rp.simulation.Molecule("Radical B", [])
    sim = HilbertSimulation([radical_a, radical_b])

    # Unlike the MARY sweep, retain the full selected time interval on the x-axis.
    # The propagator performs exact fixed-step unitary evolution; no Euler
    # approximation is used. Limit returned samples to bound public API costs.
    duration_s = float(duration_us) * 1e-6
    requested_step = float(time_step_ns) * 1e-9
    kinetic_step = 0.05 / float(recombination_rate_s)
    target_step = min(requested_step, kinetic_step)
    samples = min(5001, max(2, int(np.ceil(duration_s / target_step)) + 1))
    time_s = np.linspace(0.0, duration_s, samples)

    H = sim.total_hamiltonian(B0=float(field_mT),
                              J=float(exchange_mT), D=float(dipolar_mT))
    rho = sim.time_evolution(mapping[initial_state], time_s, H)
    # Substates form a disjoint decomposition of the electron spin space.
    keys = ("singlet", "triplet_plus", "triplet_zero", "triplet_minus")
    coherent = {}
    for key in keys:
        arr = np.asarray(sim.product_probability(mapping[key], rho), dtype=float).reshape(-1)
        if len(arr) != samples or not np.all(np.isfinite(arr)):
            raise ValueError(f"Invalid population for {key}")
        coherent[key] = arr
    coherent["triplet"] = sum(coherent[k] for k in keys if k != "singlet")
    survival = np.exp(-float(recombination_rate_s) * time_s)
    # Observable convention: raw populations remain normalised within the
    # unrecombined ensemble; decaying populations account for uniform loss.
    all_keys = ("singlet", "triplet", "triplet_plus", "triplet_minus", "triplet_zero")
    return {
        "time_us": (time_s * 1e6).tolist(),
        "populations": {k: coherent[k].tolist() for k in all_keys},
        "surviving_populations": {k: (coherent[k] * survival).tolist() for k in all_keys},
        "survival": survival.tolist(),
        "observable": observable, "initial_state": initial_state,
        "selected_population": coherent[observable].tolist(),
        "selected_surviving_population": (coherent[observable] * survival).tolist(),
        "parameters": {
            "field_mT": field_mT, "hyperfine_mT": hyperfine_mT,
            "recombination_rate_s": recombination_rate_s,
            "duration_us": duration_us, "time_step_ns": time_step_ns,
            "effective_time_step_ns": (time_s[1] - time_s[0]) * 1e9,
            "initial_state": initial_state, "observable": observable,
            "exchange_mT": exchange_mT, "dipolar_mT": dipolar_mT,
        },
        "model": "One 1H, two electrons, coherent RadicalPy Hilbert dynamics; uniform exponential survival",
        "note": "Population is conditional on an unrecombined radical pair; surviving population is population × exp(-kt). Triplet is the sum of T+, T0, T−."
    }
