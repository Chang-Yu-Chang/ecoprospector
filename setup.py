#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='ecoprospector',
      version='0.0.1',
      description='Simulate community selection protocols',
      url='https://github.com/Chang-Yu-Chang/ecoprospector',
      author=['Chang-Yu Chang', 'Jean Villa'], 
      author_email=['chang-yu.chang@yale.edu'],
      license='MIT',
      packages = ['community_selection'],
      scripts = ['commandline_tool/ecoprospector'],
      include_package_data = True,
      package_data = {"": ["*.csv"]},
      install_requires=["community-simulator@ git+https://github.com/Emergent-Behaviors-in-Biology/community-simulator.git@master"]
      zip_safe=False)
