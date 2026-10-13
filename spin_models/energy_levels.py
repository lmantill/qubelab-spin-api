"""RadicalPy Hamiltonian eigenvalues against field: eigenfrequencies in MHz."""
import numpy as np
from .common import setup

def simulate_energy_levels(*,field_max_mT=15.,field_mT=5.,hyperfine_mT=1.5,recombination_rate_s=1e6,duration_us=2.,time_step_ns=10.,initial_state='singlet',observable='singlet',exchange_mT=0.,dipolar_mT=0.,field_points=61):
    sim,_=setup(hyperfine_mT)
    field=np.linspace(0.,float(field_max_mT),int(field_points))
    spectra=[]
    for B in field:
        H=np.asarray(sim.total_hamiltonian(B0=float(B),J=float(exchange_mT),D=float(dipolar_mT)),dtype=complex)
        ev=np.linalg.eigvalsh(H)
        # RadicalPy Hamiltonians use angular frequencies in rad/s.
        spectra.append((ev/(2*np.pi*1e6)).tolist())
    return {'field_mT':field.tolist(),'levels_MHz':np.asarray(spectra).T.tolist(),'selected_field_mT':float(field_mT),'parameters':{'hyperfine_mT':hyperfine_mT,'exchange_mT':exchange_mT,'dipolar_mT':dipolar_mT,'initial_state':initial_state,'observable':observable,'recombination_rate_s':recombination_rate_s,'duration_us':duration_us,'time_step_ns':time_step_ns},'note':'Eigenvalues of the coherent Hamiltonian in MHz; initial state, observable, kinetics and timing do not change the energy spectrum. Eigenvalue branches are sorted by energy at each field, not tracked through crossings.'}
