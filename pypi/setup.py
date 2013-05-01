#from distutils.core import setup
from setuptools import setup, find_packages
setup(name='Colonyzer2',
      version='1.0.52',
      packages=find_packages(),
      description='Image analysis for microbial cultures growing on solid agar surfaces',
      long_description=open('README.txt').read(),
      entry_points={"console_scripts":["independent = bin.independent","timecourse = bin.timecourse"]},
      author='Conor Lawless',
      author_email='conor.lawless@ncl.ac.uk',
      url='http://research.ncl.ac.uk/colonyzer/',
      py_modules=['colonyzer2'],
      scripts=['bin/independent.py','bin/timecourse.py'],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Image Recognition',
        'Intended Audience :: Science/Research'
        ]      
      )
