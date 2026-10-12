# QubeLab Spin Chemistry API — starter project

This project provides a small FastAPI service and an example webpage for a remote
RadicalPy magnetic-field-effect (MARY) tutorial.

## Important scientific note

The endpoint uses RadicalPy's `mary(...)` routine with a singlet-born pair, one
isotropically coupled `1H` nucleus on radical A, a second electron-only radical,
zero exchange/dipolar couplings, and an exponential kinetic process. This is a
teaching model, not a validated model for a specific molecule or experiment.
Validate package-version compatibility, numerical convergence, observable
normalisation, and kinetic assumptions before using results scientifically.

## Files

- `main.py`: API routes and request validation
- `spin_models/radical_pair.py`: RadicalPy simulation function
- `web-demo/index.html`: responsive page with sliders and Plotly graph
- `render.yaml`: optional Render deployment blueprint

## Run locally

Use Python 3.11 in a virtual environment:

```bash
python -m venv .venv
# macOS/Linux
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1

pip install -r requirements.txt
uvicorn main:app --reload
```

API docs: http://127.0.0.1:8000/docs

Health check: http://127.0.0.1:8000/health

## Test the endpoint

```bash
curl -X POST http://127.0.0.1:8000/simulate/mfe \
  -H "Content-Type: application/json" \
  -d '{"field_max_mT":15,"hyperfine_mT":1.5,"recombination_rate_s":1000000,"field_points":31,"duration_us":1,"time_step_ns":10,"observable":"singlet"}'
```

## Deploy on Render

1. Push this directory to a new GitHub repository.
2. In Render, create a Web Service from the repository.
3. Build command: `pip install -r requirements.txt`
4. Start command: `uvicorn main:app --host 0.0.0.0 --port $PORT`
5. Set `ALLOWED_ORIGINS` to your exact website origin, e.g.
   `https://qubelab.github.io`.
6. Deploy and test `/health` and `/docs`.

The `render.yaml` file is optional; review its plan setting against current
Render plans before using it.

## Connect the webpage

1. Open `web-demo/index.html`.
2. Replace `https://YOUR-SERVICE.onrender.com` with your deployed service URL.
3. Copy the file into your QubeLab website repository (for example,
   `tutorials/radical-pair-mfe/index.html`).
4. Make sure `ALLOWED_ORIGINS` includes `https://qubelab.github.io`.
5. Commit and push to the GitHub Pages publishing branch.

## Production hardening checklist

- Pin a tested RadicalPy version once you have confirmed the API works in the
  deployment environment.
- Set CPU/memory limits and rate limits; cap sweep points and duration.
- Add authentication or a quota if the API should not be publicly usable.
- Use a background job queue for long-running calculations.
- Add tests comparing backend output with a trusted local notebook.
- Avoid accepting arbitrary Python code or arbitrary serialized objects from users.

## Fixes for 10^8 s^-1 and triplet substate selection

- `triplet_plus`, `triplet_minus` and `triplet_zero` now map to RadicalPy's
  distinct `State.TRIPLET_PLUS`, `State.TRIPLET_MINUS` and `State.TRIPLET_ZERO`.
- The integration grid automatically shortens and refines for rapid exponential
  recombination (`dt <= 0.05/k`, duration <= 12/k), rather than sampling a 10 ns
  lifetime at only one point per lifetime.
- When the zero-field observed yield is exactly zero, a relative-percent MARY
  value is undefined. The API instead reports the **absolute yield change in
  percentage points** and marks the units explicitly; the UI labels it accordingly.
- The backend accepts the frontend's 50 mT maximum.

Deploy with `main.py` at repository root, `spin_models/radical_pair.py` in its
subfolder, and `web-demo/index.html` at your website's simulation page.
The included tests use an API stub and do **not** validate numerical outputs
against a real RadicalPy installation. Verify convergence and normalization
with the RadicalPy version deployed on Render.

## New controls

The simulation now accepts `initial_state` (`singlet`, `triplet`, `triplet_plus`, `triplet_minus`, or `triplet_zero`), `exchange_mT` (-5 to 5 mT), and `dipolar_mT` (-5 to 5 mT). Defaults retain the original singlet-born, J=0, D=0 model. `triplet` is RadicalPy's triplet manifold state, rather than an individual triplet sublevel. The scalar dipolar model and exchange sign/factor follow RadicalPy's conventions.

Example JSON fields: `{"initial_state":"triplet_zero","exchange_mT":0.5,"dipolar_mT":-0.2}`.

Deploy the backend changes before replacing the website HTML; otherwise the deployed API may silently ignore the new fields.
