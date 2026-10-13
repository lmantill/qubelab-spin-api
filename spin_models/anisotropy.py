"""RadicalPy 3D orientation dependence from an axial dipolar tensor.

The original isotropic one-proton model is angularly isotropic at D=0.
A fixed molecular-frame axial dipolar tensor is introduced using the D control;
this is an explicit modelling choice rather than an anisotropic hyperfine tensor.
"""
import numpy as np
from .common import setup,times

def simulate_anisotropy(*,field_mT=5.,hyperfine_mT=1.5,recombination_rate_s=1e6,duration_us=2.,time_step_ns=10.,initial_state='singlet',observable='singlet',exchange_mT=0.,dipolar_mT=0.,theta_points=13,phi_points=25):
    from radicalpy.experiments import anisotropy_loop
    from radicalpy.kinetics import Exponential
    from radicalpy.simulation import HilbertSimulation
    # anisotropy() uses LiouvilleSimulation and orientational field rotation.
    import radicalpy as rp
    a=rp.simulation.Molecule.fromisotopes(isotopes=['1H'],hfcs=[float(hyperfine_mT)],name='Radical A')
    b=rp.simulation.Molecule('Radical B',[])
    sim=HilbertSimulation([a,b])
    _, states=setup(hyperfine_mT)
    theta=np.linspace(0,np.pi,int(theta_points))
    phi=np.linspace(0,2*np.pi,int(phi_points))
    D=float(dipolar_mT)
    tensor=np.diag([-D/2,-D/2,D])
    t=times(duration_us,time_step_ns,recombination_rate_s)
    if len(t) > 401:
        t = np.linspace(t[0], t[-1], 401)
    Hbase=sim.total_hamiltonian(B0=0.0,D=tensor,J=float(exchange_mT),hfc_anisotropy=False)
    # With an isotropic hyperfine interaction and no axial coupling, the
    # Hamiltonian is rotationally invariant. Avoid redundant angle sweeps.
    if abs(D) < 1e-15:
        H = sim.total_hamiltonian(B0=float(field_mT), D=0.0,
                                  J=float(exchange_mT), hfc_anisotropy=False)
        rho = sim.time_evolution(states[initial_state], t, H)
        curve = np.real(sim.product_probability(states[observable], rho))
        probabilities = np.broadcast_to(curve, (len(theta), len(phi), len(t)))
    else:
        probabilities=anisotropy_loop(sim,states[initial_state],states[observable],t,float(field_mT),Hbase,theta,phi)
    rate=float(recombination_rate_s)
    weighted=np.real(probabilities)*np.exp(-rate*t)[None,None,:]
    trapezoid=getattr(np,'trapezoid',np.trapz)
    yields=rate*trapezoid(weighted,x=t,axis=-1)
    expected=(len(theta),len(phi))
    if yields.shape != expected:
        raise ValueError(f'Unexpected RadicalPy anisotropy shape: {yields.shape} vs {expected}')
    if not np.all(np.isfinite(yields)):
        raise ValueError('Nonfinite anisotropy yield')
    # radius is nonnegative; color expresses actual product yield, including zeros.
    th,ph=np.meshgrid(theta,phi,indexing='ij')
    baseline=float(np.max(np.abs(yields)))
    r=0.4+0.6*np.abs(yields)/baseline if baseline>1e-15 else np.ones_like(yields)
    return {'theta_rad':theta.tolist(),'phi_rad':phi.tolist(),'x':(r*np.sin(th)*np.cos(ph)).tolist(),'y':(r*np.sin(th)*np.sin(ph)).tolist(),'z':(r*np.cos(th)).tolist(),'yield':yields.tolist(),'parameters':{'field_mT':field_mT,'hyperfine_mT':hyperfine_mT,'recombination_rate_s':recombination_rate_s,'duration_us':duration_us,'time_step_ns':time_step_ns,'initial_state':initial_state,'observable':observable,'exchange_mT':exchange_mT,'dipolar_mT':dipolar_mT},'note':'3D angular product-yield surface. Axial molecular-frame dipolar tensor diag(-D/2,-D/2,D) mT; isotropic hyperfine coupling. If D=0 the model is orientation independent. Surface radius is scaled for display; color represents yield.'}
