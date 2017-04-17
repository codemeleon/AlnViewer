#!/usr/bin/env python

from setuptools import setup

setup(name="clalnview",
      version='0.0.1',
      description="To display alignments in commanline",
      author="Anmol M. Kiran\nJennifer E. Cornick",
      author_email="anmol@liv.ac.uk\nj.cornick@liv.ac.uk",
      url="https://github.com/codemeleon/ClAlnView",
      install_requires=['click',
                        'pandas',
                        'biopython'],
      license='GPLv3',
      scripts=['bin/clalnview'],
      packages=["clalnview"],
      zip_safe=False
      )
