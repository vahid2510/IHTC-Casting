# ihtccast
[![CI](https://img.shields.io/github/actions/workflow/status/yourname/ihtccast/ci.yml?branch=main)](../../actions)
[![PyPI](https://img.shields.io/pypi/v/ihtccast.svg)](https://pypi.org/project/ihtccast/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Interfacial Heat Transfer Coefficient (IHTC) estimator for sand-cast steel, built on a transparent inverse heat conduction workflow.

- 1D transient steel+mold model with interface Robin coupling
- Physics-aware IHTC parameterization (spike–decay–plateau)
- Regularized least squares with early-time weighting
- Optional bootstrap uncertainty bands
- Command-line tool `ihtc` for quick runs

## Install
```bash
# Clone the repository
git clone https://github.com/vahid2510/IHTC-Casting.git
# Navigate into the project directory
cd IHTC-Casting
# Install or upgrade the build tool
pip install -U build
# Build the package
python -m build
# Install the built wheel
pip install dist/ihtccast-*.whl
```
## Quickstart
```bash
ihtc --temps examples/temps.csv --L_c 0.04 --L_m 0.08 --t_end 200      --dx_c 0.0005 --dx_m 0.001 --bc_mold adiabatic      --init h_inf=900 h_peak=5200 tau=2.0 t0=0.0      --bounds h_inf=50:5000 h_peak=200:20000 tau=0.02:30 t0=0:1      --lambda 1e-3 --bootstrap 50 --out results
```

Outputs:
- `results/ihtc_curve.csv` — time vs h(t)
- `results/fig_h_curve.png` — IHTC plot (±95% CI if bootstrapped)
- `results/fig_temp_compare.png` — measured vs simulated temperatures

## Documentation
- [Methodology](docs/METHODOLOGY.md)
- [References (2022–2025)](docs/REFERENCES.md)
- [FAQ](docs/FAQ.md)

## Contributing
Pull requests are welcome. See [CONTRIBUTING](CONTRIBUTING.md) and the [Code of Conduct](CODE_OF_CONDUCT.md).

## Citation
If this helps your work, cite using the metadata in `CITATION.cff`.

## License
MIT. See [LICENSE](LICENSE).
