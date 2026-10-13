"""Shared two-electron/one-proton RadicalPy teaching model."""
import numpy as np
STATES = ('singlet','triplet','triplet_plus','triplet_minus','triplet_zero')
def setup(hyperfine_mT):
    import radicalpy as rp
    from radicalpy.simulation import HilbertSimulation, State
    a=rp.simulation.Molecule.fromisotopes(isotopes=['1H'],hfcs=[float(hyperfine_mT)],name='Radical A')
    b=rp.simulation.Molecule('Radical B',[])
    return HilbertSimulation([a,b]), {k:getattr(State,v) for k,v in [('singlet','SINGLET'),('triplet','TRIPLET'),('triplet_plus','TRIPLET_PLUS'),('triplet_minus','TRIPLET_MINUS'),('triplet_zero','TRIPLET_ZERO')]}
def times(duration_us,time_step_ns,rate):
    duration=float(duration_us)*1e-6
    end=min(duration,12/float(rate))
    step=min(float(time_step_ns)*1e-9,0.05/float(rate))
    n=min(2001,max(2,int(np.ceil(end/step))+1))
    return np.linspace(0,end,n)
