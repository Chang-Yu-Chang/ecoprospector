#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from setuptools import setup


setup(name='community_selection',
      version='0.0.9999',
      description='Simulate community selection algorithms',
      url='https://github.com/Chang-Yu-Chang/community-selection',
      author='Chang-Yu Chang', 'Jean Vila',
      author_email='chang-yu.chang@yale.edu',
      license='MIT',
      packages=['community_selection'],
      py_modules = ['community_simulator', 'community_selection'],
      zip_safe=False)


