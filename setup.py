#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='Ecoprospector',
      version='0.0.9999',
      description='Simulate community selection protocols',
      url='https://github.com/Chang-Yu-Chang/Ecoprospector',
      author= ['Chang-Yu Chang', 'Jean Villa'], 
      author_email=['chang-yu.chang@yale.edu'],
      license='MIT',
      packages = ['community_selection'],
      scripts = ['commandline_tool/Ecoprospector'],
      zip_safe=False)

