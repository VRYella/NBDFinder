#!/usr/bin/env python3
"""
Setup script for NBDFinder package
==================================
"""

from setuptools import setup, find_packages
import os

# Read README for long description
with open("README_Publication.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="nbdfinder",
    version="2.0.0",
    author="Venkata Rajesh Yella",
    author_email="yvrajesh_bt@kluniversity.in",
    description="Comprehensive computational framework for non-B DNA structure detection",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/VRYella/NBDFinder",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: Free for non-commercial use",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme>=0.5",
            "sphinx-autodoc-typehints>=1.11",
        ],
    },
    entry_points={
        "console_scripts": [
            "nbdfinder=app:main",
            "nbdfinder-batch=batch_analysis:main",
        ],
    },
    include_package_data=True,
    package_data={
        "nbdfinder": [
            "data/*.fasta",
            "examples/*.txt",
            "figures/*.png",
            "figures/*.svg",
        ],
    },
    keywords="bioinformatics genomics DNA structure motifs G-quadruplex Z-DNA R-loops",
    project_urls={
        "Documentation": "https://github.com/VRYella/NBDFinder/wiki",
        "Source": "https://github.com/VRYella/NBDFinder",
        "Tracker": "https://github.com/VRYella/NBDFinder/issues",
        "Citation": "https://doi.org/10.1093/nar/gkxx",
    },
)
