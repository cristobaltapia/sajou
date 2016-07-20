# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='PyBar',
    version='0.0.1',
    description='Simple program to calculate bar structures',
    long_description=readme,
    author='Cristóbal Tapia Camú',
    author_email='crtapia@gmail.com',
    url='https://github.com/cristobaltapia/pybar',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

