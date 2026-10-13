"""RadicalPy 3D anisotropy using the official experiments.anisotropy API.

Uses one anisotropically coupled 14N nucleus (flavin N5 demonstration tensor)
for a fast web teaching model. It is NOT the two-nitrogen flavin_anion model.
The hyperfine slider scales the entire tensor; J is scalar, D is a molecular-
frame axial 3x3 tensor.  See RadicalPy's anisotropy_3d_polar_paper example.
"""
from __future__ import annotations

import numpy as np

# Flavin N5 tensor from RadicalPy's anisotropy_3d_polar_paper example (mT).
_N5_TENSOR_MT = np.array([
    [0.280, -0.138, 0.678],
    [-0.138, 0.043, -0.331],
    [0.678, -0.331, 1.412],
])
_STATES = {
    'singlet': 'SINGLET', 'triplet': 'TRIPLET',
    'triplet_plus': 'TRIPLET_PLUS', 'triplet_minus': 'TRIPLET_MINUS',
    'triplet_zero': 'TRIPLET_ZERO',
}


def simulate_anisotropy(*, field_mT=5.0, hyperfine_mT=1.5,
                        recombination_rate_s=1e6, duration_us=2.0,
                        time_step_ns=10.0, initial_state='singlet',
                        observable='singlet', exchange_mT=0.0,
                        dipolar_mT=0.0, theta_points=13, phi_points=25):
    import radicalpy as rp
    from radicalpy.experiments import anisotropy
    from radicalpy.simulation import State

    if initial_state not in _STATES or observable not in _STATES:
        raise ValueError('Unsupported spin state')
    rate = float(recombination_rate_s)
    if rate <= 0 or duration_us <= 0 or time_step_ns <= 0:
        raise ValueError('Rate, duration and time step must be positive')
    if not (5 <= int(theta_points) <= 25 and 9 <= int(phi_points) <= 49):
        raise ValueError('Angular resolution outside supported bounds')

    theta = np.linspace(0, np.pi, int(theta_points))
    phi = np.linspace(0, 2*np.pi, int(phi_points))

    # Use full angular tensor, with hyperfine_mT as its *isotropic component*.
    # Scaling the tensor preserves its anisotropy. A zero slider gives zero HFC.
    isotropic_n5 = np.trace(_N5_TENSOR_MT) / 3.0
    tensor = _N5_TENSOR_MT * (float(hyperfine_mT) / isotropic_n5)
    flavin = rp.simulation.Molecule.fromisotopes(
        isotopes=['14N'], hfcs=[tensor], name='Flavin N5 model'
    )
    partner = rp.simulation.Molecule('Partner', [])
    sim = rp.simulation.HilbertSimulation([flavin, partner])

    # The 3D dipolar tensor follows RadicalPy's tensor-in-mT convention.
    D = float(dipolar_mT)
    dipolar_tensor = D * np.diag([-0.5, -0.5, 1.0])

    duration = float(duration_us) * 1e-6
    # Once exp(-k*t) is below exp(-12), the yield contribution is negligible.
    # Truncate at 12 lifetimes to resolve 1e8/s without wasting points.
    integration_end = min(duration, 12.0/rate)
    step = min(float(time_step_ns) * 1e-9, 0.1/rate)
    n_full = max(2, int(np.ceil(integration_end/step)) + 1)
    n = min(n_full, 251)
    t = np.linspace(0.0, integration_end, n)

    result = anisotropy(
        sim,
        init_state=getattr(State, _STATES[initial_state]),
        obs_state=getattr(State, _STATES[observable]),
        time=t,
        theta=theta,
        phi=phi,
        B0=float(field_mT),
        D=dipolar_tensor,
        J=float(exchange_mT),
        kinetics=[rp.kinetics.Exponential(rate)],
    )
    Y = np.real(np.asarray(result['product_yield_sums']))
    expected = (len(theta), len(phi))
    if Y.shape != expected:
        raise ValueError(f'RadicalPy returned unexpected angular shape {Y.shape}; expected {expected}')
    if not np.isfinite(Y).all():
        raise ValueError('RadicalPy yielded non-finite values')

    average = float(np.real(rp.utils.spherical_average(Y, theta, phi)))
    delta, gamma = rp.utils.yield_anisotropy(Y, theta, phi)
    relative = Y - average
    th, ph = np.meshgrid(theta, phi, indexing='ij')
    span = float(np.max(np.abs(relative)))
    radius = (1.0 + 0.55 * relative/span) if span > 1e-12 else np.ones_like(Y)

    return {
        'theta_rad': theta.tolist(), 'phi_rad': phi.tolist(),
        'x': (radius*np.sin(th)*np.cos(ph)).tolist(),
        'y': (radius*np.sin(th)*np.sin(ph)).tolist(),
        'z': (radius*np.cos(th)).tolist(),
        'yield': Y.tolist(), 'yield_relative': relative.tolist(),
        'yield_average': average, 'delta_phi': float(delta),
        'gamma': float(gamma),
        'parameters': {
            'field_mT': field_mT, 'hyperfine_mT': hyperfine_mT,
            'recombination_rate_s': rate, 'duration_us': duration_us,
            'time_step_ns': time_step_ns, 'initial_state': initial_state,
            'observable': observable, 'exchange_mT': exchange_mT,
            'dipolar_mT': dipolar_mT, 'theta_points': len(theta),
            'phi_points': len(phi), 'effective_time_step_ns': float(t[1]-t[0])*1e9,
        },
        'note': ('RadicalPy experiments.anisotropy with one anisotropically coupled '
                 '14N (flavin N5 tensor). Hyperfine slider sets the tensor isotropic '
                 'component; J is scalar and D is the axial molecular-frame tensor '
                 'diag(-D/2,-D/2,D). Color is absolute product yield; the radius '
                 'represents yield minus its spherical average. '
                 f'Y_avg={average:.5g}, ΔΦ={float(delta):.5g}, γ={float(gamma):.5g}. '
                 'This is a teaching model, not the full two-nitrogen flavin model.'),
    }
