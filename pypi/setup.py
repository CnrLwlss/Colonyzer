from setuptools import setup, find_packages
setup(name='Colonyzer2',
      version='1.0.53',
      packages=find_packages(),
      description='Image analysis for microbial cultures growing on solid agar surfaces',
      long_description=open('README.txt').read(),
      entry_points={"console_scripts":["independent = colonyzer2.independent:main","timecourse = colonyzer2.timecourse:main",
                                       "parametryzer = colonyzer2.parameteryzer_script:main"]},
      author='Conor Lawless',
      author_email='conor.lawless@ncl.ac.uk',
      url='http://research.ncl.ac.uk/colonyzer/',
      py_modules=['colonyzer2'],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Image Recognition',
        'Intended Audience :: Science/Research'
        ]      
      )
