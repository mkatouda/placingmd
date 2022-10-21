from setuptools import setup

setup(
    name="placingmd",
    version="0.0.1",
    install_requires=[
        "pyyaml", "numpy", "scipy", "pandas", "rdkit", "meeko",
    ],
    entry_points={
        'console_scripts': [
            'placingmd=placingmd.placingmd:main',
        ],
    },
    author="Michio Katouda",
    author_email="katouda@rist.or.jp",
    description="Protein-Ligand Attached Complex Input Generator for Molecular Dynamics",
    url="https://github.com/mkatouda/placingmd",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.7',
)
