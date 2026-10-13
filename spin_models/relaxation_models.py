"""Compare actual RadicalPy Liouville-space relaxation superoperators."""
import numpy as np
from .common import setup,times

def simulate_relaxation_models(*,field_mT=5.,hyperfine_mT=1.5,recombination_rate_s=1e6,duration_us=2.,time_step_ns=10.,initial_state='singlet',observable='singlet',exchange_mT=0.,dipolar_mT=0.,relaxation_rate_s=1e6,models=None):
    import radicalpy as rp
    from radicalpy.simulation import LiouvilleSimulation
    from radicalpy import relaxation
    available={'none':None,'singlet_triplet_dephasing':relaxation.SingletTripletDephasing,'random_fields':relaxation.RandomFields,'t1':relaxation.T1Relaxation,'t2':relaxation.T2Relaxation,'triplet_triplet_dephasing':relaxation.TripletTripletDephasing}
    if models is None: models=['none','singlet_triplet_dephasing','random_fields','t1','t2']
    if len(models)>6 or len(models)<1 or any(m not in available for m in models): raise ValueError('Invalid relaxation models')
    _,state=setup(hyperfine_mT)
    t=times(duration_us,time_step_ns,recombination_rate_s)
    series={}
    for key in models:
        a=rp.simulation.Molecule.fromisotopes(isotopes=['1H'],hfcs=[float(hyperfine_mT)],name='Radical A')
        b=rp.simulation.Molecule('Radical B',[])
        sim=LiouvilleSimulation([a,b])
        H=sim.total_hamiltonian(B0=float(field_mT),J=float(exchange_mT),D=float(dipolar_mT))
        if available[key] is not None:
            sim.apply_liouville_hamiltonian_modifiers(H,[available[key](rate_constant=float(relaxation_rate_s))])
        rho=sim.time_evolution(state[initial_state],t,H)
        population=np.real(np.asarray(sim.product_probability(state[observable],rho))).reshape(-1)
        if len(population)!=len(t) or not np.all(np.isfinite(population)): raise ValueError(f'Invalid series for {key}')
        series[key]=population.tolist()
    return {'time_us':(t*1e6).tolist(),'populations':series,'surviving_populations':{k:(np.asarray(v)*np.exp(-float(recombination_rate_s)*t)).tolist() for k,v in series.items()},'parameters':{'field_mT':field_mT,'hyperfine_mT':hyperfine_mT,'recombination_rate_s':recombination_rate_s,'relaxation_rate_s':relaxation_rate_s,'duration_us':duration_us,'time_step_ns':time_step_ns,'initial_state':initial_state,'observable':observable,'exchange_mT':exchange_mT,'dipolar_mT':dipolar_mT},'note':'Actual RadicalPy Liouville relaxation superoperators; recombination is separately shown as uniform exp(-kt) survival, not a spin-selective process.'}
