# -*- coding: utf-8 -*-

# Learn more: https://github.com/dazdraperma/dezirteer/setup.py

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='Dezirteer',
    version='0.6.2020.04.28.01',
    description='Package for dezirteer.com',
    long_description=readme,
    author='Vladislav Powerman',
    author_email='vladislav.powerman@gmail.com',
    url='https://github.com/dazdraperma/dezirteer',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
	install_requires=["matplotlib","scipy"]
)