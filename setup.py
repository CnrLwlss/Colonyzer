from distutils.core import setup
setup(name='Colonyzer2',
      version='1.0.2',
      description='Image analysis for microbial cultures growing on solid agar surfaces',
      long_description='A set of functions for spot location, image segmentation and spot size estimation, suitable for quantifying culture size/cell density for a library of microbial cultures growing on a solid agar surface.',
      author='Conor Lawless',
      author_email='conor.lawless@ncl.ac.uk',
      url='http://research.ncl.ac.uk/colonyzer/',
      py_modules=['Colonyzer2'],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Image Recognition',
        'Intended Audience :: Science/Research'
        ]
      )
