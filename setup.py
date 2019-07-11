#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

setup_requirements = ['pytest-runner']

with open('requirements.txt') as f:
    requirements = list(f.readlines())

test_requirements = ['pytest']

setup(
    author="Andre Kahles",
    author_email='andre.kahles@inf.ethz.ch',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="This project contains genomics scripts and utilities used in the Ratschlab",
    entry_points={
        'console_scripts': [
            'count_expression=ratschlab_omicstools.counting.count_expression:main',
            'count_expression_prep_anno=ratschlab_omicstools.counting.count_expression_prep_anno:main',
            'splice_burden_project=ratschlab_omicstools.splice_burden.splice_burden_project:main',
            'splice_burden_compute=ratschlab_omicstools.splice_burden.splice_burden_compute:main',
            'splice_burden_plot_tcga=ratschlab_omicstools.splice_burden.splice_burden_plot_tcga:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='ratschlab_omicstools',
    name='ratschlab_omicstools',
    packages=find_packages(include=['ratschlab_omicstools']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/ratschlab/ratschlab_omicstools',
    version='0.1.0',
    zip_safe=False,
)
