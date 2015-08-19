from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy

# Nice tutorial on packaging python code:
# http://www.scotttorborg.com/python-packaging/

setup(name='Colonyzer2',
      version='1.0.88',
      packages=['colonyzer2','scripts'],
      description='Image analysis for microbial cultures growing on solid agar surfaces',
      long_description=open('README.txt').read(),
      entry_points={"console_scripts":["independent = scripts.independent:main",
                                       "timecourse = scripts.timecourse:main",
                                       "fixedpos = scripts.timecourse_fixedpos:main",
                                       "parametryzer = scripts.parameteryzer_script:main"]},
      author='Conor Lawless',
      author_email='conor.lawless@ncl.ac.uk',
      url='http://research.ncl.ac.uk/colonyzer/',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Image Recognition',
        'Intended Audience :: Science/Research'
        ],
      install_requires=['numpy','scipy','pandas','matplotlib','pillow'],
      ext_modules=cythonize("colonyzer2/maskfill.pyx"),
      include_dirs=[numpy.get_include()]
      )
