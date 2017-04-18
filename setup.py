#!/usr/bin/env python

from setuptools import setup

setup(name="alnview",
      version='0.0.1',
      description="To display alignments in commanline",
      author="Anmol M. Kiran\nJennifer E. Cornick",
      author_email="anmol@liv.ac.uk\nj.cornick@liv.ac.uk",
      url="https://github.com/codemeleon/AlnView",
      install_requires=['click',
                        'pandas',
                        'biopython'],
      license='GPLv3',
      scripts=['bin/alnview'],
      packages=["alnview"],
      zip_safe=False
      )
