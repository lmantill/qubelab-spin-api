"""Fast 3D angular product yield using RadicalPy Hamiltonians and state projectors.

The exponential survival integral is evaluated analytically in the Hamiltonian
 eigenbasis, avoiding hundreds of slow time-step propagations on Render.
"""
from __future__ import annotations

import numpy as np
from .common import setup


def _yield_from_hamiltonian(H, rho0, Q, rate, duration):
    """k ∫_0^T exp(-kt) Tr[Q exp(-iHt) rho0 exp(iHt)] dt.

    H must be in angular-frequency (rad/s) units, as returned by RadicalPy.
    """
    energies, vectors = np.linalg.eigh(np.asarray(H, dtype=complex))
    rho = vectors.conj().T @ rho0 @ vectors
    obs = vectors.conj().T @ Q @ vectors
    delta = energies[:, None] - energies[None, :]
    z = rate + 1j * delta
    factor = -rate * np.expm1(-z * duration) / z
    value = np.sum(rho * obs.T * factor)
    return float(np.real(value))


def simulate_anisotropy(*, field_mT=5.0, hyperfine_mT=1.5,
                        recombination_rate_s=1e6, duration_us=2.0,
                        time_step_ns=10.0, initial_state='singlet',
                        observable='singlet', exchange_mT=0.0, dipolar_mT=0.0,
                        theta_points=13, phi_points=25):
    sim, states = setup(hyperfine_mT)
    theta = np.linspace(0, np.pi, int(theta_points))
    phi = np.linspace(0, 2*np.pi, int(phi_points))
    rate = float(recombination_rate_s)
    duration = float(duration_us) * 1e-6
    if rate <= 0 or duration <= 0:
        raise ValueError('Rate and duration must be positive')
    D = float(dipolar_mT)
    tensor = np.diag([-D/2, -D/2, D])
    Hbase = sim.total_hamiltonian(B0=0.0, J=float(exchange_mT),
                                   D=tensor, hfc_anisotropy=False)
    # Match RadicalPy's initialization/projection conventions, including
    # normalisation of the unpolarized triplet initial state.
    rho0 = np.asarray(sim.initial_density_matrix(states[initial_state], Hbase), dtype=complex)
    Q = np.asarray(sim.observable_projection_operator(states[observable]), dtype=complex)
    result = np.empty((len(theta), len(phi)), dtype=float)
    # The axial tensor is independent of azimuth. Singlet and unpolarized
    # triplet states are rotationally invariant, so phi is redundant for those.
    phi_invariant = initial_state in ('singlet', 'triplet') and observable in ('singlet', 'triplet')
    for i, th in enumerate(theta):
        for j, ph in enumerate(phi):
            if phi_invariant and j:
                result[i, j] = result[i, 0]
                continue
            H = Hbase + sim.convert(sim.zeeman_hamiltonian(float(field_mT), theta=float(th), phi=float(ph)))
            result[i, j] = _yield_from_hamiltonian(H, rho0, Q, rate, duration)
    if not np.isfinite(result).all():
        raise ValueError('The calculation generated non-finite yields')
    result = np.clip(result, 0., 1.)
    th, ph = np.meshgrid(theta, phi, indexing='ij')
    # Plot radius in [0.4,1] while retaining the actual yield in the color.
    top = float(np.max(result))
    radius = 0.4 + 0.6*result/top if top > 1e-15 else np.ones_like(result)
    return {
        'theta_rad': theta.tolist(), 'phi_rad': phi.tolist(),
        'x': (radius*np.sin(th)*np.cos(ph)).tolist(),
        'y': (radius*np.sin(th)*np.sin(ph)).tolist(),
        'z': (radius*np.cos(th)).tolist(),
        'yield': result.tolist(),
        'parameters': {
            'field_mT':field_mT,'hyperfine_mT':hyperfine_mT,
            'recombination_rate_s':recombination_rate_s,
            'duration_us':duration_us,'time_step_ns':time_step_ns,
            'initial_state':initial_state,'observable':observable,
            'exchange_mT':exchange_mT,'dipolar_mT':dipolar_mT,
        },
        'note': ('3D product-yield angular surface using RadicalPy Hamiltonians. '
                 'The dipolar coupling is modeled as an axial tensor '
                 'diag(-D/2,-D/2,D) in mT. Color is integrated yield; '
                 'radius is scaled for display. Exponential survival is '
                 'integrated exactly for the specified duration; requested '
                 'time step is retained for consistency but is not needed here.'),
    }
