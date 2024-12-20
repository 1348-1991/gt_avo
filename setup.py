from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='gt_avo',  # Geophysical Tools: AVO
    version='0.1.0',  # Versão Inicial
    description='Conjunto de Funções para Análise de AVO e afins',
    long_description=long_description,  # Use previously defined variable
    url="https://github.com/1348-1991/gt-avo",
    author='Vitor Azevedo dos Santos',
    author_email='vitors@id.uff.br',
    packages=find_packages(),  # Automatically find packages
    install_requires=[
        "numpy>=1.23.0",  # Specify version mínima
        "matplotlib>=3.5.0",
        "pandas>=1.4.0",
        "scipy>=1.9.0",  # Additional dependency for numerical functions
    ],
    extras_require={
        "dev": [
            "pytest",
            "black",
            "flake8",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Geosciences :: Seismology",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)