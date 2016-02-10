#!/usr/bin/env python

import os
import sys
from setuptools import setup
import interpacf

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

setup(
    name="interpacf",
    version=interpacf.__version__,
    author="Brett Morris",
    author_email="bmmorris@uw.edu",
    url="https://github.com/bmorris3/interp-acf/",
    packages=["interpacf"],
    description="Take the autocorrelation function despite missing data.",
    long_description=open("README.md").read(),
    package_data={"": ["LICENSE"]},
    include_package_data=True,
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ])
