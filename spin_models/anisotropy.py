"""3D RadicalPy anisotropy using experiments.anisotropy, compatible with QUBE API.

Single anisotropic 14N nuclear HFC demonstration, not a complete flavin model.
Only the RadicalPy experiment computes spin dynamics; NumPy handles plotting data.
"""
from __future__ import annotations
import numpy as np

# Flavin N5 demonstration tensor (mT), taken from RadicalPy's example.
N5 = np.array([[0.280, -0.138, 0.678],
               [-0.138, 0.043, -0.331],
               [0.678, -0.331, 1.412]], dtype=float)
STATES = {
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

    if initial_state not in STATES or observable not in STATES:
        raise ValueError('Unknown initial state or observable')
    if not (5 <= int(theta_points) <= 25 and 9 <= int(phi_points) <= 49):
        raise ValueError('Angular grid is out of bounds')
    if recombination_rate_s <= 0 or duration_us <= 0 or time_step_ns <= 0:
        raise ValueError('Rate and integration times must be positive')

    # RadicalPy requires an even number of azimuthal samples.
    # Accept old frontend values (17/25/41), and enforce evenness server-side.
    phi_count = int(phi_points)
    if phi_count % 2:
        phi_count = min(phi_count + 1, 48)  # Keep within the API's 49-point limit.
    theta = np.linspace(0.0, np.pi, int(theta_points))
    phi = np.linspace(0.0, 2*np.pi, phi_count)

    # IMPORTANT: Hfc's documented dispatch accepts a nested Python list,
    # not necessarily a NumPy ndarray. Preserve anisotropic tensor symmetry.
    scale = float(hyperfine_mT) / (np.trace(N5)/3.0)
    hfc = (N5 * scale).tolist()
    flavin = rp.simulation.Molecule.fromisotopes(
        isotopes=['14N'], hfcs=[hfc], name='Flavin N5 demonstration')
    partner = rp.simulation.Molecule('Partner', [])
    sim = rp.simulation.HilbertSimulation([flavin, partner])

    # In this model, the angular anisotropy comes primarily from anisotropic HFC.
    # D is an axial molecular-frame tensor; J is an isotropic exchange parameter.
    D = float(dipolar_mT)
    dipolar = (D * np.diag([-0.5, -0.5, 1.0])).tolist()

    rate = float(recombination_rate_s)
    end = min(float(duration_us)*1e-6, 12.0/rate)
    dt = min(float(time_step_ns)*1e-9, 0.1/rate)
    # Cap work while maintaining sampling on the exponential decay scale.
    n_time = max(2, min(301, int(np.ceil(end/dt)) + 1))
    time = np.linspace(0, end, n_time)

    result = anisotropy(
        sim, init_state=getattr(State, STATES[initial_state]),
        obs_state=getattr(State, STATES[observable]), time=time,
        theta=theta, phi=phi, B0=float(field_mT),
        D=np.asarray(dipolar, dtype=float), J=float(exchange_mT),
        kinetics=[rp.kinetics.Exponential(rate)],
    )

    Y = np.asarray(result['product_yield_sums'], dtype=float)
    expected = (len(theta), len(phi))
    if Y.shape != expected:
        raise ValueError(f'Unexpected RadicalPy yield grid {Y.shape}; expected {expected}')
    if not np.isfinite(Y).all():
        raise ValueError('Non-finite yields returned from RadicalPy')

    # Compatible with RadicalPy's output and avoids a divide-by-zero from
    # yield_anisotropy() when the surface is exactly isotropic or zero.
    avg = float(rp.utils.spherical_average(Y, theta, phi))
    relative = Y - avg
    amplitude = float(np.max(Y) - np.min(Y))
    delta_phi = amplitude
    gamma = (amplitude / abs(avg)) if abs(avg) > 1e-14 else 0.0
    surface_amplitude = float(np.max(np.abs(relative)))
    radius = 1.0 + 0.55 * relative/surface_amplitude if surface_amplitude > 1e-12 else np.ones_like(Y)
    th, ph = np.meshgrid(theta, phi, indexing='ij')

    return {
        'theta_rad': theta.tolist(), 'phi_rad': phi.tolist(),
        'x': (radius*np.sin(th)*np.cos(ph)).tolist(),
        'y': (radius*np.sin(th)*np.sin(ph)).tolist(),
        'z': (radius*np.cos(th)).tolist(),
        'yield': Y.tolist(), 'yield_relative': relative.tolist(),
        'yield_average': avg, 'delta_phi': delta_phi, 'gamma': gamma,
        'parameters': {
            'field_mT': float(field_mT), 'hyperfine_mT': float(hyperfine_mT),
            'recombination_rate_s': rate, 'duration_us': float(duration_us),
            'time_step_ns': float(time_step_ns), 'initial_state': initial_state,
            'observable': observable, 'exchange_mT': float(exchange_mT),
            'dipolar_mT': D, 'theta_points': len(theta), 'phi_points': len(phi),
            'effective_time_step_ns': float(time[1]-time[0])*1e9,
        },
        'note': ('RadicalPy anisotropy experiment, simplified flavin N5 14N model. '
                 'Anisotropic hyperfine tensor is scaled by the hyperfine slider. '
                 'Colors show absolute yields, radial deformation shows deviations '
                 'from the spherical mean. ΔΦ is maximum minus minimum yield; '
                 'γ is ΔΦ divided by the absolute mean when nonzero.'),
    }
