from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="spatial-needleman-wunsch",
    version="1.0.0",
    author="Joshua D. Curry",
    author_email="jcurry3428@smail.pcd.edu",
    description="A deterministic dynamic programming framework for 3D molecular docking",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JDCurry/spatial-needleman-wunsch",
    project_urls={
        "Source Code": "https://github.com/JDCurry/spatial-needleman-wunsch",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.812",
            "sphinx>=4.0",
            "jupyter>=1.0",
        ],
        "gpu": [
            "cupy>=9.0",
            "numba>=0.55",
        ],
        "rdkit": [
            "rdkit-pypi>=2021.09.4",
        ],
    },
    entry_points={
        "console_scripts": [
            "spatial-dock=spatial_docking.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "spatial_docking": [
            "data/*.json",
            "data/*.sdf",
            "examples/*.ipynb",
        ],
    },
    keywords=[
        "molecular docking",
        "computational biology",
        "dynamic programming",
        "structural biology",
        "drug discovery",
        "bioinformatics",
        "deterministic algorithms",
    ],
    zip_safe=False,
)
