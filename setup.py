#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut für Eisenforschung, Dusseldorf, Germany 
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL). 
calphy is distributed in the hope that it will be useful for non-commercial academic research, 
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the ASL for more details. 

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
“Automated Free Energy Calculation from Atomistic Simulations.” Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de/yury.lysogorskiy@icams.rub.de
"""

from setuptools import setup, find_packages, Extension

with open('README.md') as readme_file:
    readme = readme_file.read()

test_requirements = ['pytest>=3', ]

setup(
    author="Sarath Menon, Yury Lysogorskiy, Ralf Drautz",
    author_email='sarathmenon@mailbox.org',
    python_requires='>=3.8',
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    description="free energy calculation for python",
    install_requires=['matplotlib', 'pytest',
    'pyyaml', 'mendeleev', 
    'tqdm', 'scipy', 'pydantic', 'pyscal3'],
    license="GNU General Public License v3",
    long_description=readme,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='calphy',
    name='calphy',
    packages=find_packages(include=['calphy', 'calphy.*']),
    test_suite='tests',
    url='https://github.com/ICAMS/calphy',
    version='1.4.2',
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'calphy = calphy.kernel:main',
            'calphy_kernel = calphy.queuekernel:main',
            'calphy_run_averaging = calphy.clitools:run_averaging',
            'calphy_process_averaging = calphy.clitools:process_averaging',
            'calphy_run_integration = calphy.clitools:run_integration',
            'calphy_process_integration = calphy.clitools:process_integration',
            'calphy_convert_input = calphy.clitools:convert_legacy_inputfile',
            'calphy_phase_diagram = calphy.clitools:phase_diagram',
        ],
    }
)

