#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages, Extension

with open('README.md') as readme_file:
    readme = readme_file.read()

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Sarath Menon",
    author_email='sarathmenon@mailbox.org',
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="a thermodynamic integration helper for python",
    install_requires=['matplotlib'],
    license="GNU General Public License v3",
    long_description=readme,
    include_package_data=True,
    keywords='pytint',
    name='pytint',
    packages=find_packages(include=['pytint', 'pytint.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/srmnitc/pytint',
    version='0.3.8',
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'tint = pytint.kernel:main',
            'tint_kernel = pytint.queuekernel:main',
        ],
    }
)
