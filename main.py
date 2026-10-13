from __future__ import annotations

import os
import time
from typing import Literal

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from spin_models.radical_pair import simulate_one_proton_mary
from spin_models.time_evolution import simulate_time_evolution

app = FastAPI(
    title="QubeLab Spin Chemistry API",
    version="0.1.0",
    description="Remote RadicalPy simulations for QubeLab tutorials.",
)

# Add your production website origin(s) as a comma-separated environment variable.
# Example: https://qubelab.github.io,https://www.example.com
allowed_origins = [
    origin.strip()
    for origin in os.getenv(
        "ALLOWED_ORIGINS", "https://qubelab.github.io,http://localhost:8000"
    ).split(",")
    if origin.strip()
]
app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=False,
    allow_methods=["GET", "POST"],
    allow_headers=["Content-Type"],
)


class MARYRequest(BaseModel):
    field_max_mT: float = Field(default=15.0, ge=0.1, le=50.0)
    hyperfine_mT: float = Field(default=1.5, ge=0.0, le=5.0)
    recombination_rate_s: float = Field(default=1.0e6, ge=1.0e4, le=1.0e8)
    field_points: int = Field(default=31, ge=11, le=151)
    duration_us: float = Field(default=2.0, ge=0.1, le=5.0)
    time_step_ns: float = Field(default=10.0, ge=2.0, le=50.0)
    observable: Literal[
        "singlet", "triplet", "triplet_plus", "triplet_minus", "triplet_zero"
    ] = "singlet"
    initial_state: Literal[
        "singlet", "triplet", "triplet_plus", "triplet_minus", "triplet_zero"
    ] = "singlet"
    exchange_mT: float = Field(default=0.0, ge=-5.0, le=5.0)
    dipolar_mT: float = Field(default=0.0, ge=-5.0, le=5.0)


@app.get("/")
def root():
    return {
        "service": "QubeLab Spin Chemistry API",
        "status": "running",
        "docs": "/docs",
    }


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/simulate/mfe")
def simulate_mfe(request: MARYRequest):
    started = time.perf_counter()
    try:
        result = simulate_one_proton_mary(
            field_max_mT=request.field_max_mT,
            hyperfine_mT=request.hyperfine_mT,
            recombination_rate_s=request.recombination_rate_s,
            field_points=request.field_points,
            duration_us=request.duration_us,
            time_step_ns=request.time_step_ns,
            observable=request.observable,
            initial_state=request.initial_state,
            exchange_mT=request.exchange_mT,
            dipolar_mT=request.dipolar_mT,
        )
    except ImportError as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Simulation dependency unavailable: {exc}",
        ) from exc
    except Exception as exc:
        # Log full exception server-side; do not expose internals to public clients.
        print(f"Simulation failed: {type(exc).__name__}: {exc}")
        raise HTTPException(
            status_code=500,
            detail="The simulation failed. Check the server logs.",
        ) from exc

    result["elapsed_s"] = round(time.perf_counter() - started, 3)
    return result


class TimeEvolutionRequest(BaseModel):
    field_mT: float = Field(default=5.0, ge=0.0, le=50.0)
    hyperfine_mT: float = Field(default=1.5, ge=0.0, le=5.0)
    recombination_rate_s: float = Field(default=1e6, ge=1e4, le=1e8)
    duration_us: float = Field(default=2.0, ge=0.1, le=5.0)
    time_step_ns: float = Field(default=10.0, ge=2.0, le=50.0)
    observable: Literal["singlet", "triplet", "triplet_plus", "triplet_minus", "triplet_zero"] = "singlet"
    initial_state: Literal["singlet", "triplet", "triplet_plus", "triplet_minus", "triplet_zero"] = "singlet"
    exchange_mT: float = Field(default=0.0, ge=-5.0, le=5.0)
    dipolar_mT: float = Field(default=0.0, ge=-5.0, le=5.0)


@app.post("/simulate/time-evolution")
def simulate_time_evolution_endpoint(request: TimeEvolutionRequest):
    started = time.perf_counter()
    try:
        result = simulate_time_evolution(**request.model_dump())
    except ImportError as exc:
        raise HTTPException(status_code=503, detail=f"Simulation dependency unavailable: {exc}") from exc
    except Exception as exc:
        print(f"Time evolution failed: {type(exc).__name__}: {exc}")
        raise HTTPException(status_code=500, detail="Time evolution simulation failed. Check server logs.") from exc
    result["elapsed_s"] = round(time.perf_counter() - started, 3)
    return result


# Additional independent teaching simulations.
from spin_models.energy_levels import simulate_energy_levels
from spin_models.anisotropy import simulate_anisotropy
from spin_models.relaxation_models import simulate_relaxation_models

class EnergyRequest(TimeEvolutionRequest):
    field_max_mT: float = Field(default=15.0, ge=0.1, le=50.0)
    field_points: int = Field(default=61, ge=11, le=101)

class AnisotropyRequest(TimeEvolutionRequest):
    theta_points: int = Field(default=13, ge=5, le=25)
    phi_points: int = Field(default=25, ge=9, le=49)

class RelaxationRequest(TimeEvolutionRequest):
    relaxation_rate_s: float = Field(default=1e6, ge=1e3, le=1e8)
    models: list[Literal['none','singlet_triplet_dephasing','random_fields','t1','t2','triplet_triplet_dephasing']] = Field(default_factory=lambda:['none','singlet_triplet_dephasing','random_fields','t1','t2'],min_length=1,max_length=6)

def guarded(fn, values, label):
    start=time.perf_counter()
    try:
        result=fn(**values)
    except ImportError as e:
        raise HTTPException(status_code=503,detail=f'Simulation dependency unavailable: {e}') from e
    except Exception as e:
        print(f'{label} failed: {type(e).__name__}: {e}')
        raise HTTPException(status_code=500,detail=f'{label} failed. Check server logs.') from e
    result['elapsed_s']=round(time.perf_counter()-start,3)
    return result

@app.post('/simulate/energy-levels')
def energy_levels_endpoint(request:EnergyRequest):
    return guarded(simulate_energy_levels,request.model_dump(),'Energy levels')

@app.post('/simulate/anisotropy')
def anisotropy_endpoint(request:AnisotropyRequest):
    return guarded(simulate_anisotropy,request.model_dump(),'Anisotropy')

@app.post('/simulate/relaxation')
def relaxation_endpoint(request:RelaxationRequest):
    return guarded(simulate_relaxation_models,request.model_dump(),'Relaxation')
