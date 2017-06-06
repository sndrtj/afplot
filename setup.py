"""
setup.py
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""
from os.path import abspath, dirname, join

from setuptools import setup

readme_file = join(abspath(dirname(__file__)), "README.rst")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()

setup(
    name="afplot",
    version="0.2",
    description="Plot allele frequencies in VCF files",
    long_description=long_desc,
    author="Sander Bollen",
    author_email="a.h.b.bollen@lumc.nl",
    url="https://github.com/sndrtj/afplot",
    license="MIT",
    packages=["afplot"],
    install_requires=[
        "click",
        "numpy",
        "matplotlib",
        "pandas",
        "seaborn",
        "progressbar2",
        "pysam",
        "pyvcf"
    ],
    entry_points={
        "console_scripts": [
            "afplot = afplot.cli:main"
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]

)