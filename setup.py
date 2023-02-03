# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 13:57:24 2023

@author: thuds
"""

from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

#requirements = ["glob2>=0.7", "pandas>=1.5.2", "openpyxl>=3", "FPDF>=1.7.2"]

setup(
    name="SCD_NGS_Analysis",
    version="0.0.1",
    author="Taylor Hudson, Marena Trinidad",
    author_email="thudson@berkeley.edu",
    description="A package for ",
    long_description=readme,
    long_description_content_type="Python, Bash, Text",
    url="https://github.com/ACTG802/SCD_NGS_Analysis",
    packages=find_packages(),
    #install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.9"
    ],
)