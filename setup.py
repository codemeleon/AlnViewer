#!/usr/bin/env python

from setuptools import setup

setup(name="alnviewer",
      version='0.0.1',
      description="To display alignments in commanline",
      author="Anmol M. Kiran\nJennifer E. Cornick",
      author_email="anmol@liv.ac.uk\nj.cornick@liv.ac.uk",
      url="https://github.com/codemeleon/AlnViewer",
      install_requires=['click',
                        'pandas',
                        'biopython'],
      license='GPLv3',
      scripts=['bin/alnviewer'],
      packages=["alnviewer"],
      zip_safe=False
      )
