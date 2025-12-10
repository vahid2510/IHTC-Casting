# Contributing

Thanks for considering a contribution. Keep it simple and reproducible.

1. Fork and clone the repo
2. Create a branch: `git checkout -b feature/whatever`
3. Add tests for your change
4. Run tests locally
5. Open a pull request with a clear description

### Dev setup
```bash
python -m venv .venv && source .venv/bin/activate
pip install -U pip build pytest
pip install -e .
pytest -q
```
