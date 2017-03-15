# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='Sajou',
    version='0.1.0',
    description='Simple structural analysis library for Python',
    long_description=readme,
    author='Cristóbal Tapia Camú',
    author_email='crtapia@gmail.com',
    url='https://github.com/cristobaltapia/sajou',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

