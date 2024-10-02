#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the pyproject.toml file.

import sys

# First provide helpful messages if contributors try and run legacy commands
# for tests or docs.

# Only import these if the above checks are okay
# to avoid masking the real problem with import error.
from setuptools import setup  # noqa: E402

setup()