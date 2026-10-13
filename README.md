# RadicalPy QUBE Lab: three additional simulations

Includes existing MARY and time-evolution files, plus independent energy-level, 3D anisotropy, and relaxation-model pages.

## Deployment

Copy `main.py` and the entire `spin_models/` folder into your Render backend repository; keep `requirements.txt` and `render.yaml`. Redeploy Render. Copy `web-demo/energy_levels.html`, `web-demo/anisotropy.html`, and `web-demo/relaxation_models.html` into your GitHub Pages site (such as `tutorials/`). Set the API URL on each page if your Render deployment URL differs.

## API endpoints

- `POST /simulate/energy-levels`: magnetic-field sweep of eigenvalues, MHz.
- `POST /simulate/anisotropy`: 3D angular yield via RadicalPy `experiments.anisotropy`.
- `POST /simulate/relaxation`: different RadicalPy Liouville relaxation models vs time.
- Original `POST /simulate/mfe` and `POST /simulate/time-evolution` unchanged.

## Physical caveats

Energy levels are coherent Hamiltonian eigenvalues; initial state and kinetics do not change them. This one-proton model uses isotropic hyperfine coupling. The anisotropy page explicitly models the scalar D control as an axial dipolar *tensor*, `diag(-D/2,-D/2,D)` in the molecular frame; this differs from using a scalar D in the previous simulations. If D=0, the model is isotropic and the 3D surface is spherical. The surface radius is normalized for visualization, while its colors show yields. Relaxation uses RadicalPy's Liouville-space operators; uniform recombination survival is applied separately in the displayed series. These are demonstration models, not validated fits.

## Verification

All Python modules were compiled for syntax, and the endpoints/page payloads were inspected. Numerical execution requires installing a compatible RadicalPy version; verify output after deployment, especially anisotropy output shape, and cross-check model conventions with package documentation.
