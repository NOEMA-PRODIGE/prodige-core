# PRODIGE-core

![PyPI - License](https://img.shields.io/pypi/l/prodige_core?color=green)
![PyPI - Version](https://img.shields.io/pypi/v/prodige_core)
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2FNOEMA-PRODIGE%2Fprodige-core%2Fmaster%2Fpyproject.toml)
[![codecov](https://codecov.io/gh/NOEMA-PRODIGE/prodige-core/graph/badge.svg?token=OIQJYPIUO8)](https://codecov.io/gh/NOEMA-PRODIGE/prodige-core)


Collection of functions for a more common analysis/display of the [PRODIGE](NOEMA-PRODIGE.github.io) program.

This includes simple to use commands to make publication ready figures, as well as some convenience function to handle the data.


## Installation

Installation using pip:
```python
pip install prodige-core
```

From the repository root:
```bash
git clone https://github.com/NOEMA-PRODIGE/prodige-core.git
```

Install the package from the cloned repository with:
```bash
python -m pip install -e .
```

For development and testing, install the test dependencies as well:

```bash
python -m pip install -e ".[test]"
```


## Credits

---

This is developed by:

- Jaime E Pineda ([@jpinedaf](http://github.com/jpinedaf))
- Caroline Gieser ([@CarolineGieser](http://github.com/CarolineGieser))
- PRODIGE team

## Dependencies

---

- Python >=3.11
- scipy (>=1.10)
- astropy (>=5.0)
- matplotlib (>=3.6)
- radio-beam (>=0.3)
- spectral-cube (>=0.5.0)
- reproject

The plotting helpers expect FITS images and use the source metadata shipped with
the package. Data paths should include a trailing path separator because the
current API accepts directory strings.
