#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Setup script for the Python bindings of the bioseq library.
This is a thin wrapper around maturin, which is used for building the Rust extension.
"""

from setuptools import setup
from setuptools_rust import RustExtension, Binding

setup(
    name="biopython-rust",
    version="0.1.0",
    description="High-performance bioinformatics sequence operations library",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Xuliang Xu",
    author_email="xlxmobile@outlook.com",
    url="https://github.com/xlxmobile/Biopython-Rust",
    packages=["biopython-rust"],
    rust_extensions=[
        RustExtension(
            "biopython_rust._rust_bindings",
            "Cargo.toml",
            binding=Binding.PyO3,
            debug=False,
        )
    ],
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=23.0.0",
            "isort>=5.12.0",
            "mypy>=1.0.0",
        ],
        "viz": [
            "matplotlib>=3.5.0",
            "seaborn>=0.12.0",
        ],
    },
    python_requires=">=3.7",
    zip_safe=False,
    platforms="any",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Rust",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)